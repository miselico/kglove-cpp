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
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >) = 0;
};

class UniformWeigher: public GraphWeigher {
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >);
};

class InverseFrequencyWeigher: public GraphWeigher {
	virtual TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weigh(TPt<TNodeEdgeNet<TStr, TStr> >);
};

#endif /* GRAPHWEIGHER_H_ */
