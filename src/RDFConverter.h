/*
 * RDFConverter.h
 *
 *  Created on: Nov 24, 2016
 *      Author: cochez
 */

#ifndef RDFCONVERTER_H_
#define RDFCONVERTER_H_

#include <iostream>
#include "Snap.h"

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


#include "MurmurHashAdditions.h"
#include "doublePriorityQueue.h"
#include "nTripleParser.h"
#include "BCA.h"



#endif /* RDFCONVERTER_H_ */
