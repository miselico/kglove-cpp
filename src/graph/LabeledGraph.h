/*
 * LabeledGraph.h
 *
 *  Created on: Mar 16, 2017
 *      Author: cochez
 */

#ifndef GRAPH_LABELEDGRAPH_H_
#define GRAPH_LABELEDGRAPH_H_

#include <vector>
#include <string>
#include <boost/flyweight.hpp>

namespace QuickGraph {

//fprward decl
class Node;

class Edge {
	//no copying of edges
	//Edge(const Edge & sedge);
public:
	boost::flyweight<std::string> label;
	long double weight;
	//the index of the target in the graph vector
	const unsigned int targetIndex;

	Edge(const std::string & label, long double weight, const unsigned int targetIndex) :
			label(label), weight(weight), targetIndex(targetIndex) {

	}

};

class Node {
	//no copying of nodes
	//Node(const Node & node);
public:
	Node(std::string label) : label(label){

	}
	Node(boost::flyweight<std::string> label) : label(label){

	}

	boost::flyweight<std::string> label;
	std::vector<Edge> edges;
};

class LabeledGraph {

public:
	LabeledGraph() {
	}
	virtual ~LabeledGraph() {
	}
	std::vector<Node> nodes;

	template<typename FuncType> void forEachEdge(FuncType operation) {
		std::vector<QuickGraph::Node> & nodes = this->nodes;
		for (std::vector<QuickGraph::Node>::iterator nodeI = nodes.begin(); nodeI != nodes.end(); nodeI++) {
			std::vector<QuickGraph::Edge> & edges = nodeI->edges;
			for (std::vector<QuickGraph::Edge>::iterator edgeI = edges.begin(); edgeI != edges.end(); edgeI++) {
				operation(*this, *edgeI);
			}
		}
	}

	//attempt to reduce the memory footprint of this graph
	void pack() {
		this->nodes.shrink_to_fit();
		for (std::vector<Node>::iterator node = this->nodes.begin(); node != this->nodes.end(); node++) {
			node->edges.shrink_to_fit();
		}
	}

};

//class ImmutableLabeledGraph {
//private:
//	vector<Node> nodes;
//public:
//	ImmutableLabeledGraph(const LabeledGraph &);
//	ImmutableLabeledGraph(const ImmutableLabeledGraph &);
//	virtual ~ImmutableLabeledGraph();
//
//};

}

#endif /* GRAPH_LABELEDGRAPH_H_ */

