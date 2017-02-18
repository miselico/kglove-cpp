/*
 * utils.h
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#ifndef UTILS_H_
#define UTILS_H_

#include "Snap.h"

//template<class NodeData, class EdgeData>
//TPt<TNodeEdgeNet<NodeData, EdgeData> > reverseGraph(TPt<TNodeEdgeNet<NodeData, EdgeData> > baseGraph);

TPt<TNodeEdgeNet<TStr, TStr> > reverseGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph);


#endif /* UTILS_H_ */
