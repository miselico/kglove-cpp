/*
 * GraphWeigher.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef GRAPHWEIGHER_HA_
#define GRAPHWEIGHER_HA_

#include <memory>
#include "graph/LabeledGraph.h"
#include <unordered_map>

class GraphWeigher {
protected:
	GraphWeigher() {

	}

public:
	virtual ~GraphWeigher() {

	}
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const = 0;
};

class UniformWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class InversePredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class PredicateFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class ObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class InverseObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class PredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class InversePredicateObjectFrequencyWeigher: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

class InverseObjectFrequencyWeigherSplitDown: public GraphWeigher {
public:
	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

/**
 * Assigns to each inedge the weight assigned to the nodes.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class PushDownWeigher: public GraphWeigher {
	const std::unordered_map<std::string, double> nodeWeights;
	const double defaultWeight;

public:
	/**
	 * If defaultweight is set to -1, it indicates that all weights MUST be in the nodeWeights. If not, the program will be terminated.
	 */
	PushDownWeigher(const std::unordered_map<std::string, double> & nodeWeights, const double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

/**
 * Assigns to each inedge the weight assigned to the nodes divided by the number of in edges.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight divided by #inedge
 * Then all weights are normalized
 *
 *
 * First, each in edge gets the weight of the node / #inedge
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class SplitDownWeigher: public GraphWeigher {
	const std::unordered_map<std::string, double> nodeWeights;
	const double defaultWeight;

public:
	/**
	 * If defaultweight is set to -1, it indicates that all weights MUST be in the nodeWeights. If not, the program will be terminated.
	 */
	SplitDownWeigher(const std::unordered_map<std::string, double> & nodeWeights, const double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual void weigh(std::shared_ptr<QuickGraph::LabeledGraph>) const;
};

#endif /* GRAPHWEIGHER_HA_ */
