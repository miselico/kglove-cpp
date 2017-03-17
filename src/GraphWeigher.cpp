/*
 * GraphWeigher.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWeigher.h"

#include <iostream>

#include <string>
#include <unordered_map>
#include "graph/LabeledGraph.h"

using namespace std;

typedef boost::flyweight<std::string> flyString;

namespace { //helpers

//count freq of each property:
unordered_map<flyString, double> absolute_predicate_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	unordered_map<flyString, double> absolute_freq;
	std::vector<QuickGraph::Node> & nodes = baseGraph->nodes;
	for (auto nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
		std::vector<QuickGraph::Edge> & edges = nodeI->edges;
		for (auto edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
			flyString predictae = edgeI->label;
			auto oldValI = absolute_freq.find(predictae);
			if (oldValI == absolute_freq.end()) {
				absolute_freq[predictae] = 1;
			} else {
				absolute_freq[predictae] = oldValI->second + 1;
			}
		}
	}
	return absolute_freq;
}

//counts the absolute frequence of a the objects
unordered_map<string, double> absolute_object_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	unordered_map<string, double> absolute_freq;
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string object = EI.GetDstNDat();
		double start = 0.0;
		absolute_freq.IsKeyGetDat(object, start);
		absolute_freq.AddDat(object, start + 1);
	}
	return absolute_freq;
}

unordered_map<stringPr, double> absolute_predicate_object_freq(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) {
	unordered_map<stringPr, double> absolute_freq;
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		string obj = EI.GetDstNDat();
		stringPr PO(pred, obj);
		double start = 0.0;
		absolute_freq.IsKeyGetDat(PO, start);
		absolute_freq.AddDat(PO, start + 1);
	}
	return absolute_freq;
}

//inverse all freq
template<class type> unordered_map<type, double> inverse_the_frequency(unordered_map<type, double> & absolute_freq) {
	unordered_map<type, double> inverse_freq;
	for (typename unordered_map<type, double>::TIter iter = absolute_freq.BegI(); iter < absolute_freq.EndI(); iter++) {
		inverse_freq.AddDat(iter.GetKey(), 1.0 / iter.GetDat());
	}
	return inverse_freq;
}

//normalize all weights in unbalanced (sum weight on outedges == 1)
void normalize(TPt<TNodeEdgeNet<string, WeightedPredicate> > unbalanced) {
	for (TNodeEdgeNet<string, WeightedPredicate>::TNodeI NI = unbalanced->BegNI(); NI < unbalanced->EndNI(); NI++) {
		int node_i_outdeg = NI.GetOutDeg();
		double totalWeight = 0;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			totalWeight += NI.GetOutEDat(outEdge).W();
		}
		double totalWeightInverse = 1.0 / totalWeight;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			WeightedPredicate & wp = NI.GetOutEDat(outEdge);
			double normalized_weight = wp.W() * totalWeightInverse;
			wp.Val2 = normalized_weight;
		}
	}
}

}

void UniformWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {

	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		WeightedPredicate wpred(pred, 1);
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize (newNet);
	return newNet;
}

void PredicateFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	unordered_map<string, double> absolute_freq = absolute_predicate_freq(baseGraph);
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		WeightedPredicate wpred(pred, absolute_freq.GetDat(pred));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize (newNet);
	return newNet;
}

void InversePredicateFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	//count freq of each property:
	unordered_map<string, double> absolute_freq = absolute_predicate_freq(baseGraph);
	unordered_map<string, double> inverse_freq = inverse_the_frequency(absolute_freq);
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		WeightedPredicate wpred(pred, inverse_freq.GetDat(pred));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize (newNet);
	return newNet;
}

void ObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	unordered_map<string, double> absolute_freq = absolute_object_freq(baseGraph);
	PushDownWeigher subWeigher(absolute_freq, -1);
	return subWeigher.weigh(baseGraph);
}

void InverseObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	unordered_map<string, double> absolute_freq = absolute_object_freq(baseGraph);
	unordered_map<string, double> inverse_freq = inverse_the_frequency(absolute_freq);
	PushDownWeigher subWeigher(inverse_freq, -1);
	return subWeigher.weigh(baseGraph);
}

void InverseObjectFrequencyWeigherSplitDown::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	unordered_map<string, double> absolute_freq = absolute_object_freq(baseGraph);
	unordered_map<string, double> inverse_freq = inverse_the_frequency(absolute_freq);
	SplitDownWeigher subWeigher(inverse_freq, -1);
	return subWeigher.weigh(baseGraph);
}

void PredicateObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	unordered_map<stringPr, double> absolute_freq = absolute_predicate_object_freq(baseGraph);
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		string obj = EI.GetDstNDat();
		stringPr PO(pred, obj);
		WeightedPredicate wpred(pred, absolute_freq.GetDat(PO));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize (newNet);
	return newNet;
}

void InversePredicateObjectFrequencyWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	unordered_map<stringPr, double> absolute_freq = absolute_predicate_object_freq(baseGraph);
	unordered_map<stringPr, double> inverse_freq = inverse_the_frequency(absolute_freq);
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);
	for (TNodeEdgeNet<string, string>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		string pred = EI.GetDat();
		string obj = EI.GetDstNDat();
		stringPr PO(pred, obj);
		WeightedPredicate wpred(pred, inverse_freq.GetDat(PO));
		newNet->AddEdge(EI.GetSrcNId(), EI.GetDstNId(), EI.GetId(), wpred);
	}
	normalize (newNet);
	return newNet;
}

void PushDownWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);

	//add all edges with assigned weight
	for (TNodeEdgeNet<string, string>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		double weight = this->defaultWeight;
		//overwrite default if a value is set
		nodeWeights.IsKeyGetDat(NI.GetDat(), weight);
		int node_i_indeg = NI.GetInDeg();
		if (node_i_indeg > 0) {
			if (this->defaultWeight == -1.0 && weight == -1.0) {
				cerr << "Defaultweight was  -1 and node was not found in the given map -> ERROR, quitting";
				exit(1);
			}
			//add all *in* edges with the weight
			for (int inEdge = 0; inEdge < node_i_indeg; ++inEdge) {
				string label = NI.GetInEDat(inEdge);
				const WeightedPredicate pred(label, weight);
				newNet->AddEdge(NI.GetInNId(inEdge), NI.GetId(), NI.GetInEId(inEdge), pred);
			}
		}
	}
	normalize (newNet);

	return newNet;

}

void SplitDownWeigher::weigh(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph) const {
	TPt < TNodeEdgeNet<string, WeightedPredicate> > newNet = newNetCopyNodes(baseGraph);

	//add all edges with assigned weight
	for (TNodeEdgeNet<string, string>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		double weight = this->defaultWeight;
		//overwrite default if a value is set
		nodeWeights.IsKeyGetDat(NI.GetDat(), weight);
		int indeg = NI.GetInDeg();
		if (indeg > 0) {
			if (this->defaultWeight == -1.0 && weight == -1.0) {
				cerr << "Defaultweight was  -1 and node was not found in the given map -> ERROR, quitting";
				exit(1);
			}
			//add all *in* edges with the weight/indegree
			weight = weight / double(indeg);
			for (int inEdge = 0; inEdge < indeg; ++inEdge) {
				string label = NI.GetInEDat(inEdge);
				const WeightedPredicate pred(label, weight);
				newNet->AddEdge(NI.GetInNId(inEdge), NI.GetId(), NI.GetInEId(inEdge), pred);
			}
		}
	}
	normalize (newNet);

	return newNet;

}
