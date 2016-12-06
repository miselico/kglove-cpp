/*
 * RDF2Co_occurence.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef RDF2CO_OCCURENCE_H_
#define RDF2CO_OCCURENCE_H_

#include "GraphWeigher.h"

namespace RDF2CO {
void performExperiments() {
	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
	//TStr file = "sample-wikidata-terms-fragment.nt";
	//TStr file = "sample-wikidata-terms.nt";
	TStr file = "SmallTest4.nt";

	FILE* outfile = fopen("frequencies_output.bin", "w");

	InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

	computeFrequencies(file, weigher, outfile);

	fclose(outfile);
	fclose (f);
}
}

#endif /* RDF2CO_OCCURENCE_H_ */
