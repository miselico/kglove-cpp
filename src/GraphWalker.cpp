/*
 * GraphWalker.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWalker.h"

TVec<TStr> RandomProportionalWalker::walk(TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph) {
	TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI previousnode = graph->GetRndNI(this->random);
	TVec<TStr> path = TVec<TStr>();
	path.Add(previousnode.GetDat());
	for (int j = 0; j < this->walklength; ++j) {
		//select random neighbour taking into account the weights
		double d = this->random.GetUniDev();

		int outdegree = previousnode.GetOutDeg();
		if (outdegree < 1) {
			break;
		}
		for (int k = 0; k < outdegree; k++) {
			WeightedPredicate edgedata = previousnode.GetOutEDat(k);
			if (edgedata.W() > d) {
				//selected
				path.Add(edgedata.P());
				//make step in the walk  and get label of O
				previousnode = graph->GetNI(previousnode.GetOutNId(k));
				path.Add(previousnode.GetDat());
				break;
			} else {
				d = d - edgedata.W();
			}

		}
	}

	return path;
}

