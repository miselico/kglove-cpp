/*
 * WeightedPredicate.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef WEIGHTEDPREDICATE_H_
#define WEIGHTEDPREDICATE_H_

#include "Snap.h"

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

#endif /* WEIGHTEDPREDICATE_H_ */
