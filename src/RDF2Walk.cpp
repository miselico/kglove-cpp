/*
 * RDF2Walk.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "RDF2Walk.h"
#include "nTripleParser.h"


TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraphIgnoringLiterals(TStr filename, GraphWeigher & weighingStrategy) {
	TPt<TNodeEdgeNet<TStr, TStr> > graph;
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	graph = graphAndNodeIndex.Val1;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return weightedGraph;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraph(TStr filename, GraphWeigher & weighingStrategy) {
	TPt<TNodeEdgeNet<TStr, TStr> > graph;
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraph(filename);
	graph = graphAndNodeIndex.Val1;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return weightedGraph;
}



void TextFileSink::consume(TVec<TStr> path) {
	int sepCount = 0;
	for (int partNr = 0; partNr < path.Len(); partNr++) {
		fwrite(" ", sizeof(char), sepCount, fout);
		sepCount = 1;
		TStr thePart = path[partNr];
		fwrite(thePart.CStr(), sizeof(char), thePart.Len(), fout);
	}
	fwrite("\n", sizeof(char), 1, fout);
}
