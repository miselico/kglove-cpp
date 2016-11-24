//============================================================================
// Name        : RDFConverter.cpp
// Author      : Michael Cochez
// Version     :
// Copyright   : (c) Michael Cochez
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "Snap.h"

#include "MurmurHashAdditions.h"

#include "doublePriorityQueue.h"

#include "nTripleParser.h"

using namespace std;

class WeightedPredicate: public TPair<TStr, TLFlt> {

public:
	WeightedPredicate() :
			TPair() {
	}

	WeightedPredicate(TStr predicate, double weight = 1.0) :
			TPair(predicate, weight) {
	}

	TStr P() {
		return this->Val1;
	}

	long double W() {
		return this->Val2.Val;
	}

};

TPt<TNodeEdgeNet<TStr, TStr> > buildRDFGraph(TStr filename) {

//The graph
	TPt<TNodeEdgeNet<TStr, TStr> > Net = TNodeEdgeNet<TStr, TStr>::New();

//temporary keep track of all the nodes
	THash<TStr, int> addedNodes = THash<TStr, int>();

	PSIn FInPt = TFIn::New(filename);
	TStr line;
	int count = 0;
//To reuse the strings in the node/edge labels
	THashSet<TStr> stringpool = THashSet<TStr>();

	while (FInPt->GetNextLn(line)) {
		Triple values = parsetripleLine(line);
		//This saves about 10% memory on a small test. Perhaps more on a larger one.
		values = Triple(stringpool.GetKey(stringpool.AddKey(values.S())), stringpool.GetKey(stringpool.AddKey(values.P())), stringpool.GetKey(stringpool.AddKey(values.O())));

		int subjectIndex = 0;
		if (addedNodes.IsKeyGetDat(values.S(), subjectIndex)) {
			//nothing, subjectIndex now contains the node ID
		} else {
			//add new Node, save index
			subjectIndex = Net->AddNode(-1, values.S());
			addedNodes.AddDat(values.S(), subjectIndex);
		}

		int objectIndex = 0;
		if (addedNodes.IsKeyGetDat(values.O(), objectIndex)) {
			//nothing, objectIndex now contains the node ID
		} else {
			//add new Node, save index
			objectIndex = Net->AddNode(-1, values.O());
			addedNodes.AddDat(values.O(), objectIndex);
		}
		//add edge
		Net->AddEdge(subjectIndex, objectIndex, -1, values.O());

		count++;
		if (count % 100000 == 0) {
			cout << "Processed " << count << " lines" << endl;
		}
	}

	stringpool.Clr(true, -1);

	return Net;
}

// weighing strategies
TPt<TNodeEdgeNet<TStr, WeightedPredicate> > uniformWeigher(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}

	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int outDegree = baseGraph->GetNI(EI.GetSrcNId()).GetOutDeg();
		double weight = 1.0 / (double) outDegree;
		TStr label = EI.GetDat();
		int ID = EI.GetId();
		int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		newNet->AddEdge(src, dst, ID, WeightedPredicate(label, weight));
	}

	return newNet;
}
//sparse vector representing the approx pagerank
class BCV: THash<TInt, TFlt> {
public:
	void fixPaint(int ID, double amount) {
		double startAmount = this->GetDatWithDefault(ID, 0.0);
		double newAmount = startAmount + amount;
		this->AddDat(ID, newAmount);
	}

	string toString(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network) {
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
};

class BCAQueue: doublePriorityQueue<TInt> {
public:
	void addPaintTo(int toID, double paint) {
		double current = this->GetPriority(toID);
		this->SetPriority(toID, current + paint);
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

BCV computeBCA(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps) {
	BCAQueue Q = BCAQueue();
	BCV p = BCV();
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
			double edgeWeight = node_i.GetOutEDat(outEdge).W();
			double paintToJ = (1.0 - alpha) * w * edgeWeight;
			Q.addPaintTo(j, paintToJ);
		}
	}
	return p;

}

void computeFrequencies(TStr filename, TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weighingStrategy(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph)) {
	TPt<TNodeEdgeNet<TStr, TStr> > graph = buildRDFGraph(filename);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy(graph);
	BCV bcv = computeBCA(weightedGraph, 2, 0.50, 0.0000001);

	cout << bcv.toString(weightedGraph);
	return;
}

int main() {
	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
	TStr file = "sample-wikidata-terms-fragment.nt";
	//TStr file = "sample-wikidata-terms.nt";

	computeFrequencies(file, uniformWeigher);

	return 0;
}

