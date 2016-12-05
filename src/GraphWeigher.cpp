/*
 * GraphWeigher.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWeigher.h"

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > UniformWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}

	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int outDegree = baseGraph->GetNI(EI.GetSrcNId()).GetOutDeg();
		double weight = 1.0 / (double) outDegree;
		TStr label = EI.GetDat();
		const int ID = EI.GetId();
		const int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		const WeightedPredicate pred(label, weight);
		newNet->AddEdge(src, dst, ID, pred);
	}

	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > InverseFrequencyWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}

	//count freq of each property:
	THash<TStr, TInt> absolute_freq = THash<TStr, TInt>();
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr prop = EI.GetDat();
		TInt start = 0;
		absolute_freq.IsKeyGetDat(prop, start);
		absolute_freq.AddDat(prop, start + 1);
	}
	//inverse all freq
	THash<TStr, TFlt> inverse_freq;
	for (THash<TStr, TInt>::TIter iter = absolute_freq.BegI(); iter < absolute_freq.EndI(); iter++) {
		inverse_freq.AddDat(iter.GetKey(), 1.0 / iter.GetDat());
	}

	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {

		int node_i_outdeg = NI.GetOutDeg();

		double totalWeight = 0;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			TStr edgeData = NI.GetOutEDat(outEdge);
			totalWeight += inverse_freq.GetDat(edgeData);
		}
		double totalWeightInverse = 1.0 / totalWeight;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			int j = NI.GetOutNId(outEdge);
			TStr label = NI.GetOutEDat(outEdge);
			double normalized_weight = inverse_freq.GetDat(label) * totalWeightInverse;
			const WeightedPredicate pred(label, normalized_weight);
			newNet->AddEdge(NI.GetId(), j, NI.GetOutEId(outEdge), pred);
		}

	}

	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > PushDownWeigher::weigh(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) const {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	//add all nodes
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}
	//add all edges with assigned weight
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		TFlt weight = this->defaultWeight;
		//overwrite default if a value is set
		nodeWeights.IsKeyGetDat(NI.GetDat(), weight);

		//add all *in* edges with the weight
		int node_i_indeg = NI.GetInDeg();
		for (int inEdge = 0; inEdge < node_i_indeg; ++inEdge) {
			int j = NI.GetInNId(inEdge);
			TStr label = NI.GetInEDat(inEdge);
			const WeightedPredicate pred(label, weight);
			newNet->AddEdge(j, NI.GetId(), NI.GetInEId(inEdge), pred);
		}
	}

	//normalize all weights in newNet (sum weight on outedges == 1)

	for (TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI NI = newNet->BegNI(); NI < newNet->EndNI(); NI++) {

		int node_i_outdeg = NI.GetOutDeg();

		double totalWeight = 0;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			WeightedPredicate edgeData = NI.GetOutEDat(outEdge);
			totalWeight += edgeData.W();
		}
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			WeightedPredicate& edgeData = NI.GetOutEDat(outEdge);
			double normalizedWeight = edgeData.W() / totalWeight;
			//overwrite
			edgeData.Val2 = normalizedWeight;
		}
	}

	return newNet;

}

