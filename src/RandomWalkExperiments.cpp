/*
 * RandomWalkExperiments.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include "GraphWalker.h"
#include "GraphWeigher.h"
#include "RandomWalkExperiments.h"
#include "RDF2Walk.h"
#include <iostream>

//FIXME remove include
#include "nTripleParser.h"

#include "Snap.h"

namespace RandomWalkExperiments {

using namespace std;

void experiment0() {
	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
	//TStr file = "sample-wikidata-terms-fragment.nt";
	//TStr file = "sample-wikidata-terms.nt";
	//TStr file = "SmallTest3.nt";
	TStr file = "SmallTest5_duplicates.nt";

	//char* outfile = "frequencies_output.bin";

	//	computeFrequencies(file, inverseFrequencyWeigher, fopen(outfile, "w"));

	const char* walksoutfile = "walks";

	FILE* f = fopen(walksoutfile, "w");

	int length = 8;
	int seed = 45645;
	RandomProportionalWalker walker(length, seed);
	LengthEnforcingWalker walker2(walker, 2 * length + 1);
	InverseFrequencyWeigher weigher;
	TextFileSink sink = TextFileSink(f);

	int amount = 1E6;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraph(file, weigher);

	TRnd anRnd(24332);

	for (long int i = 0; i < amount; ++i) {
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI start = graph->GetRndNI(anRnd);
		TVec<TStr> path = walker.performWalk(graph, start);
		sink.consume(path);
	}

	fclose(f);
}

static TStr RDF_TYPE ("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static TStr OWL_THING ("<http://www.w3.org/2002/07/owl#Thing>");

void experiment1() {
//	* Let's start with 4 hops i.e., initial node + 4 predicates + 4 nodes.
//	* Shorter walks could be allowed, but I don't think that is going to happen in the DBpedia graph.
//	* All entities should be included (~5M). We can ignore the literals. Last time I used the following files from the DBpedia download server: instance_types, instance_types_transitive, categories, mappingbased_objects, Wikipedia_links and external_links. I think we can use the same subset again.
//	* let's try the uniform probability first. If that works, we can then use the inverse frequency of the predicate, and then the PageRank of the object.

	TStr file = "allData.nt";
//	TStr file = "SmallTest6_literals.nt";

	const char* walksoutfile = "walks.txt";
	FILE* f = fopen(walksoutfile, "w");
	int length = 4;
	int walksPerNode = 200;
	int seed = 45645;
	RandomProportionalWalker walker(length, seed);
	//LengthEnforcingWalker walker2 (walker, 2 * length + 1);
	UniformWeigher weigher;
	TextFileSink sink(f);

	cout << "Reading in data" << endl;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > graph = buildWalkableGraphIgnoringLiterals(file, weigher);


	long sumWalkLengths = 0;
	long generatedPaths = 0;
	int conceptCount = 0;
	long int nodeCount = graph->GetNodes();
	for(int i = 0; i < nodeCount; i++){
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI candidateNode = graph->GetNI(i);
		for (int outEdgeNr = 0; outEdgeNr < candidateNode.GetOutDeg(); ++outEdgeNr) {
			TStr predicate = candidateNode.GetOutEDat(outEdgeNr).P();
			if (predicate == RDF_TYPE){
				TStr object = candidateNode.GetOutNDat(outEdgeNr);
				if (object == OWL_THING){
					conceptCount++;
					//Generate paths
					for (int j = 0; j < walksPerNode; ++j) {
						TVec<TStr> path = walker.performWalk(graph, candidateNode);
						sumWalkLengths += path.Len();
						sink.consume(path);
						generatedPaths++;
					}
					//paths generated for this node
					break;
				}
			}
		}


	}
	fclose(f);

	cout << "Found " << conceptCount << " concepts out of " << nodeCount << " entities. Generated " << walksPerNode << " per concept, i.e., " << generatedPaths << " paths" << endl;

	cout << "All paths generated, Average path length (total length including start node and predicates) = " << double(sumWalkLengths) / double(generatedPaths) << endl;
}

void performExperiments() {
	experiment1();
}

}
