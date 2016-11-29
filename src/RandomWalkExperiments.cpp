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

#include "Snap.h"

namespace RandomWalkExperiments {
void performExperiments() {
	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
	//TStr file = "sample-wikidata-terms-fragment.nt";
	//TStr file = "sample-wikidata-terms.nt";
	//TStr file = "SmallTest3.nt";
	TStr file = "SmallTest4.nt";

	//char* outfile = "frequencies_output.bin";

	//	computeFrequencies(file, inverseFrequencyWeigher, fopen(outfile, "w"));

	const char* walksoutfile = "walks";

	FILE* f = fopen(walksoutfile, "w");

	int length = 8;
	int seed = 45645;
	RandomProportionalWalker walker = RandomProportionalWalker(length, seed);
	LengthEnforcingWalker walker2 = LengthEnforcingWalker(walker, 2 * length + 1);
	InverseFrequencyWeigher weigher = InverseFrequencyWeigher();
	TextFileSink sink = TextFileSink(f);

	int amount = 1E6;


	computeWalks(file, weigher, walker2, amount, sink);
	fclose(f);
}



}

