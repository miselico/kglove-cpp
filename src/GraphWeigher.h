/*
 * GraphWeigher.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef GRAPHWEIGHER_H_
#define GRAPHWEIGHER_H_

#include "Snap.h"
#include "WeightedPredicate.h"

class GraphWeigher {
protected:
	GraphWeigher() {

	}
	virtual ~GraphWeigher() {

	}
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const = 0;
};

class UniformWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

class InverseFrequencyWeigher: public GraphWeigher {
public:
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};

/**
 * Assigns a weight to the edges depending on the weight assigned to the nodes.
 * Nodes which are not in the nodeWeights provided get assigned the defaultWeight
 *
 *
 * First, each in edge gets the weight of the node
 * Then, each weight on the outedges of each node is normalized such that they sum to 1.
 *
 */
class PushDownWeigher: public GraphWeigher {
	const THash<TStr, TFlt>  nodeWeights;
	const double defaultWeight;

public:
	PushDownWeigher(const THash<TStr, TFlt>  nodeWeights, const double defaultWeight) :
			nodeWeights(nodeWeights), defaultWeight(defaultWeight) {
	}

	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) const;
};




/*
 Ideas for other weighers:

 - start by assigning weight/indegree to the edges
 - start by inverting all weights


 */




#endif /* GRAPHWEIGHER_H_ */
