/*
 * BCA.cpp
 *
 *  Created on: Nov 24, 2016
 *      Author: cochez
 */

#include "BCA.h"
#include "doublePriorityQueue.h"

using namespace std;

void BCV::fixPaint(int ID, double amount) {
	double startAmount = this->GetDatWithDefault(ID, 0.0);
	double newAmount = startAmount + amount;
	this->AddDat(ID, newAmount);
}

string BCV::toString(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network) {
	string s = "{";
	string separator = "";

	for (THash<TInt, TFlt>::TIter iter = this->BegI(); iter < this->EndI(); iter++) {
		s += separator;
		const TInt k = iter.GetKey();
		TStr entity = network->GetNDat(k);
		s += entity.CStr();
		s += " = ";
		const TFlt v = iter.GetDat();
		s += v.GetStr().CStr();
		separator = ", ";
	}

	s.append("}");
	return s;
}

void BCV::removeEntry(int ID){
	this->DelKey(ID);
}

void BCV::normalizeInPlace() {
	double totalSum = 0.0;
	for (THash<TInt, TFlt>::TIter iter = this->BegI(); iter < this->EndI(); iter++) {
		totalSum += iter.GetDat();
	}
	for (THash<TInt, TFlt>::TIter iter = this->BegI(); iter < this->EndI(); iter++) {
		TFlt & value = iter.GetDat();
		double scaled = value / totalSum;
		//sets the value in the vector
		value = scaled;
	}
}

class BCAQueue: doublePriorityQueue<TInt> {
public:
	void addPaintTo(int toID, double paint) {
		double current = this->GetPriority(toID);
		this->SetPriority(toID, current + paint);
	}

	void addPaintToIfMoreAsEps(int toID, double paint, double eps) {
		double current = this->GetPriority(toID);
		double sum = current + paint;
		if (sum > eps) {
			this->SetPriority(toID, sum);
		}
	}

	bool empty() {
		return this->IsEmpty();
	}

	TPair<TInt, TFlt> pop() {
		double paint = this->GetMaxPriority();
		TInt ID = this->PopMax();
		return TPair<TInt, TFlt>(ID, paint);
	}

};

/**
 * Compute the bookmarking coloring algorithm (≃ personalized page rank) between node b_ID and all other nodes in the graph. Using teleportation parameter alpha and cut-off value eps.
 */
BCV computeBCA(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps) {
	BCAQueue Q;
	BCV p;
	Q.addPaintTo(b_ID, 1.0);
	while (!Q.empty()) {
		TPair<TInt, TFlt> element = Q.pop();
		int i = element.Val1;
		double w = element.Val2;
		p.fixPaint(i, alpha * w);
		if (w < eps) {
			continue;
		}
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI node_i = network->GetNI(i);
		int node_i_outdeg = node_i.GetOutDeg();

		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			int j = node_i.GetOutNId(outEdge);
			WeightedPredicate edgeData = node_i.GetOutEDat(outEdge);
			double edgeWeight = edgeData.W();
			double paintToJ = (1.0 - alpha) * w * edgeWeight;
			Q.addPaintTo(j, paintToJ);
		}
	}
	return p;

}

/**
 * Compute the bookmarking coloring algorithm (≃ personalized page rank) between node b_ID and all other nodes in the graph. Using teleportation parameter alpha and cut-off value eps.
 *
 * This version also adds pagerank for the edges, each time when BCA traverses them.
 *
 * The first element of the pair is the normal BCV, the second one the pagerank for the predicates
 *
 */
TPair<BCV, BCV> computeBCAIncludingEdges(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps) {
	BCAQueue Q;
	BCV p;
	BCV edgePagerank;
	Q.addPaintTo(b_ID, 1.0);

	while (!Q.empty()) {
		TPair<TInt, TFlt> element = Q.pop();
		int i = element.Val1;
		double w = element.Val2;
		p.fixPaint(i, alpha * w);
		if (w < eps) {
			continue;
		}
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI node_i = network->GetNI(i);
		int node_i_outdeg = node_i.GetOutDeg();

		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			int j = node_i.GetOutNId(outEdge);
			WeightedPredicate edgeData = node_i.GetOutEDat(outEdge);
			double edgeWeight = edgeData.W();
			double paintToJ = (1.0 - alpha) * w * edgeWeight;
			int outEdgeID = node_i.GetOutEId(outEdge);
			edgePagerank.fixPaint(outEdgeID, paintToJ);
			Q.addPaintTo(j, paintToJ);
		}
	}
	return TPair<BCV, BCV>(p, edgePagerank);

}

//From here on implementation of pushed Bookmark Coloring Algorithm

void PBCV::fixPaint(TPair<TInt, TInt> pred_obj_pair, double amount) {
	double startAmount = this->GetDatWithDefault(pred_obj_pair, 0.0);
	double newAmount = startAmount + amount;
	this->AddDat(pred_obj_pair, newAmount);
}

string PBCV::toString(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network) {
	string s = "{";
	string separator = "";

	for (THash<TPair<TInt, TInt>, TFlt>::TIter iter = this->BegI(); iter < this->EndI(); iter++) {
		s += separator;
		const TPair<TInt, TInt> k = iter.GetKey();
		WeightedPredicate edgeData = network->GetEDat(k.Val1);
		TStr predicate = edgeData.P();
		TStr entity = network->GetNDat(k.Val2);
		s += "(";
		s += predicate.CStr();
		s += ",";
		s += entity.CStr();
		s += ")";
		s += " = ";
		const TFlt v = iter.GetDat();
		s += v.GetStr().CStr();
		separator = ", ";
	}

	s.append("}");
	return s;
}

PBCV computePBCA(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps) {

	//Note!! This queue contains the amount of paint which still has to be moved OUT of the nodes.
	//Nothing of it should stay on the the node itself.
	//This is the main confusing difference between the ordering of the algorithms.
	PBCV p;
	BCAQueue Q;
	Q.addPaintTo(b_ID, 1.0);
	while (!Q.empty()) {
		TPair<TInt, TFlt> element = Q.pop();
		int i = element.Val1;
		double w = element.Val2;

//For all links i -> j
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI node_i = network->GetNI(i);
		int node_i_outdeg = node_i.GetOutDeg();
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			int j = node_i.GetOutNId(outEdge);

			WeightedPredicate edgeData = node_i.GetOutEDat(outEdge);
			double edgeWeight = edgeData.W();
			double paintToJ = w * edgeWeight;

			TPair<TInt, TInt> combinedNode = TPair<TInt, TInt>(node_i.GetOutEId(outEdge), j);
			//stand-in for p_i = p_i + alpha*w
			p.fixPaint(combinedNode, alpha * paintToJ);

			double extrapaintToBeMovedOutFromJ = (1 - alpha) * paintToJ;
			Q.addPaintToIfMoreAsEps(j, extrapaintToBeMovedOutFromJ, eps);
		}
	}
	return p;
}
