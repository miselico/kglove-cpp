/*
 * nTripleParser.h
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#ifndef NTRIPLEPARSER_H_
#define NTRIPLEPARSER_H_

#include <iostream>

#include "Snap.h" //this is a bit overkill, but includes everything needed for glib
class Triple: public TTriple<TStr, TStr, TStr> {
public:
	Triple(TStr S, TStr P, TStr O) :
			TTriple(S, P, O) {
	}

	TStr S() {
		return this->Val1;
	}
	TStr P() {
		return this->Val2;
	}
	TStr O() {
		return this->Val3;
	}

};

Triple parsetripleLine(TStr line);

#endif /* NTRIPLEPARSER_H_ */
