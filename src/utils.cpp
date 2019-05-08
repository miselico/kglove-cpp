/*
 * utils.cpp
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#include "utils.h"
#include <memory>
#include <vector>
#include "graph/LabeledGraph.h"

using namespace std;

//template<class NodeData, class EdgeData>
//TPt<TNodeEdgeNet<NodeData, EdgeData> > reverseGraph(TPt<TNodeEdgeNet<NodeData, EdgeData> > baseGraph) {
//	TPt<TNodeEdgeNet<NodeData, EdgeData> > newNet = TNodeEdgeNet<NodeData, EdgeData>::New();
//	for (typename TNodeEdgeNet<NodeData, EdgeData>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
//		newNet->AddNode(NI.GetId(), NI.GetDat());
//	}
//	//add all edges reversed
//
//	for (typename TNodeEdgeNet<NodeData, EdgeData>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
//		newNet->AddEdge(EI.GetDstNId(), EI.GetSrcNId(), EI.GetId(), EI.GetDat());
//	}
//	return newNet;
//}


