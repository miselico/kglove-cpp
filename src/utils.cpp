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




std::shared_ptr<QuickGraph::LabeledGraph> reverseGraph(std::shared_ptr<const QuickGraph::LabeledGraph> baseGraph){
	std::shared_ptr<QuickGraph::LabeledGraph> newGraph(new QuickGraph::LabeledGraph);
	const std::vector<QuickGraph::Node> & originalNodes = baseGraph->nodes;
	std::vector<QuickGraph::Node> & newNodes = newGraph->nodes;

	//add all nodes:
	newNodes.reserve(baseGraph->nodes.size());
	for (vector<QuickGraph::Node>::const_iterator iter = originalNodes.begin(); iter != originalNodes.end(); ){
		newNodes.emplace_back(iter->label);
	}

	//add all edges reversed

	for (vector<QuickGraph::Node>::const_iterator oldSrcNode = originalNodes.begin(); oldSrcNode != originalNodes.end(); ){
		const int oldSourceIndex = oldSrcNode - originalNodes.begin();
		const int newTargetIndex = oldSourceIndex;
		const std::vector<QuickGraph::Node>::iterator newTargetNode = newNodes.begin() + newTargetIndex;

		std::vector<QuickGraph::Edge> oldEdges = oldSrcNode->edges;
		for(vector<QuickGraph::Edge>::const_iterator oldEdge = oldEdges.begin(); oldEdge != oldEdges.end(); oldEdge++){
			int oldTargetIndex = oldEdge->targetIndex;
			int newSourceIndex = oldTargetIndex;
			QuickGraph::Node & newSourceNode = newNodes[newSourceIndex];
			newSourceNode.edges.emplace_back(oldEdge->label, oldEdge->weight, newTargetIndex);
		}
	}
	return newGraph;
}
