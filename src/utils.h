/*
 * utils.h
 *
 *  Created on: Feb 18, 2017
 *      Author: cochez
 */

#ifndef UTILS_HA_
#define UTILS_HA_

#include <string>
#include "graph/LabeledGraph.h"

//template<class NodeData, class EdgeData>
//TPt<TNodeEdgeNet<NodeData, EdgeData> > reverseGraph(TPt<TNodeEdgeNet<NodeData, EdgeData> > baseGraph);

std::shared_ptr<QuickGraph::LabeledGraph> reverseGraph(std::shared_ptr<QuickGraph::LabeledGraph> baseGraph);

static std::string RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static std::string OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");

struct pairhash {
public:
	template<typename T, typename U>
	std::size_t operator()(const std::pair<T, U> &x) const {
		return std::hash<T>()(x.first) ^ std::hash<U>()(x.second);
	}
};

#endif /* UTILS_HA_ */
