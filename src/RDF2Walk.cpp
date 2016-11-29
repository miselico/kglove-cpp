/*
 * RDF2Walk.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "RDF2Walk.h"
#include "nTripleParser.h"

namespace {
TPt<TNodeEdgeNet<TStr, WeightedPredicate> > buildWalkableGraph(TStr filename, GraphWeigher & weighingStrategy) {
	TPt<TNodeEdgeNet<TStr, TStr> > graph;
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = buildRDFGraph(filename);
	graph = graphAndNodeIndex.Val1;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);
	return weightedGraph;
}
}

void computeWalks(TStr filename, GraphWeigher & weighingStrategy, GraphWalker& walker, int amount, WalkSink & sink) {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraph(filename, weighingStrategy);

	TTmStopWatch w = TTmStopWatch(true);

	//these should come as parameters:
	for (int i = 0; i < amount; ++i) {
		TVec<TStr> path = walker.walk(graph);
		sink.consume(path);
	}
	w.Stop();
	printf("%d ms for %d random walks\n", w.GetMSecInt(), amount);

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
