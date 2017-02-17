/*
 * BCA.h
 *
 *  Created on: Nov 24, 2016
 *      Author: cochez
 */

#ifndef BCA_H_
#define BCA_H_


#include "Snap.h"

#include "WeightedPredicate.h"
#include <string>

using namespace std;

//sparse vector representing the approx pagerank
class BCV: public THash<TInt, TFlt> {
public:
	string toString(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network);
	void fixPaint(int ID, double amount);
	void removeEntry(int ID);
	//This function normalizes the vector such that pageranks add up to 1 IN PLACE
	void normalizeInPlace();
};

BCV computeBCA(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps);


class PBCV: public THash<TPair<TInt, TInt>, TFlt> {
public:
	string toString(const TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network);
	void fixPaint(TPair<TInt, TInt> pred_obj_pair, double amount);

};


//PBCA = Pushed Bookmar cocloring algorithm
PBCV computePBCA(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps);

TPair<BCV, BCV> computeBCAIncludingEdges(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > network, int b_ID, double alpha, double eps);


#endif /* BCA_H_ */
