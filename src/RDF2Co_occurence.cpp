/*
 * RDF2Co_occurence.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include <iostream>
#include "Snap.h"
#include "WeightedPredicate.h"
#include "nTripleParser.h"
#include "BCA.h"
#include "GraphWeigher.h"

namespace {
/**
 * typedefs compatible with the input expected by glove
 */
typedef double real;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
} CREC;

/**
 * Get a table which assigns a unique index to each (predicate, object) pair.
 *
 * The index will be strictly greater as zero since that is what glove uses to skip an entry.
 */
THash<TPair<TStr, TStr>, TInt> createPairedWordIndexTable(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TPair<TStr, TStr>, TInt> table;
	int counter = 1;
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = graph->BegNI(); NI < graph->EndNI(); NI++) {

		int node_i_outdeg = NI.GetOutDeg();
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			TStr predicate = NI.GetOutEDat(outEdge);
			TStr object = NI.GetOutNDat(outEdge);
			TPair<TStr, TStr> pair = TPair<TStr, TStr>(predicate, object);
			if (!table.IsKey(pair)) {
				table.AddDat(pair, counter);
				counter++;
			}
		}
	}
	return table;
}


/**
 * Not finished!!
 *
 * Compute the BCA score for each pair in the graph under the given weiging strategy.
 *
 * Outputs the score as a sparse matrix which can be fed to glove.
 *
 *
 *
 */
void computeFrequencies(TStr filename, GraphWeigher& weighingStrategy, double bca_alpha, double bca_eps, FILE *fout) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	/*
	 * The wordIndexTable is indexed from 0 as it should be, but we need IDs which start from 1. So we create a new table with corrected gloveIDTable
	 */

	THash<TStr, int> gloveIDTable;
	{//scoping wordIndexTable to avoid mistakes
		THash<TStr, int> wordIndexTable = graphAndNodeIndex.Val2;
		for (THash<TStr, int>::TIter iter = wordIndexTable.BegI(); iter < wordIndexTable.EndI(); iter++) {
			gloveIDTable.AddDat(iter.GetKey(), iter.GetDat() + 1);
		}
	}

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);

	cout << "There are still some extra things done, remove if asserts are not triggered";

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {

		int focusWordGraphID = i;
		int focusWordGloveID = gloveIDTable.GetDat(weightedGraph->GetNDat(focusWordGraphID));

		BCV bcv = computeBCA(weightedGraph, focusWordGraphID, bca_alpha, bca_eps);

		for (THash<TInt, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
			int contextWordGraphID = iter.GetKey();
			int contextWordGloveID = gloveIDTable.GetDat(weightedGraph->GetNDat(contextWordGraphID));

			double freq = iter.GetDat();
			CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
			fwrite(&crec, sizeof(CREC), 1, fout);
		}
	}

	cerr << "TODO : write glove ID table to file";

}

void computeFrequenciesPushBCA(TStr filename, GraphWeigher& weighingStrategy, FILE *fout) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraph(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

	THash<TStr, int> wordIndexTable = graphAndNodeIndex.Val2;

	THash<TPair<TStr, TStr>, TInt> pairwordIndexTable = createPairedWordIndexTable(graph);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);

	graph.Clr();

	cout << "done weighing" << endl;

	cerr << "TODO : check - does the indexing for glove have to start from 1 or 0 ??";
	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {
		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
			continue;
		}
		PBCV bcv = computePBCA(weightedGraph, i, 0.10, 0.000001);
		int subjectIndex = wordIndexTable.GetDat(weightedGraph->GetNDat(i));
		for (THash<TPair<TInt, TInt>, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
			WeightedPredicate wpred = weightedGraph->GetEDat(iter.GetKey().Val1);
			TStr pred = wpred.P();
			TStr obj = weightedGraph->GetNDat(iter.GetKey().Val2);
			int offset = wordIndexTable.Len();
			int pairIndex = pairwordIndexTable.GetDat(TPair<TStr, TStr>(pred, obj)) + offset;
			double freq = iter.GetDat();

			CREC crec = CREC { word1:subjectIndex, word2:pairIndex, val: freq };
			fwrite(&crec, sizeof(CREC), 1, fout);

			cout << subjectIndex << " has " << pairIndex << " freq " << freq << endl;
		}
	}

//	TTmStopWatch w (true);
//	int needed = 10000;
//	int * selected = (int*) malloc(needed * sizeof(int));
//	int skipped = 0;
//	for (int i = 0; i < needed + skipped && i < weightedGraph->GetNodes(); ++i) {
//		if (weightedGraph->GetNDat(i).SearchCh('<') != 0) {
//			++skipped;
//			continue;
//		} else {
//			selected[i - skipped] = i;
//		}
//	}
//	for (int index = 0; index < needed && index < weightedGraph->GetNodes(); index++) {
//		int id = selected[index];
//		PBCV bcv = computePBCA(weightedGraph, id, 0.10, 0.000000000001);
//		//cout << bcv.toString(weightedGraph) << endl;
//		if (index % 1000 == 0) {
//			cout << "another 1000" << weightedGraph->GetNDat(id).CStr() << "->" << bcv.toString(weightedGraph) << endl;
//		}
//	}
//
//	w.Stop();
//	cout << w.GetMSecInt() << "ms" << endl;
	return;
}

} //end anonymous namespace.

namespace RDF2CO {
void performExperiments() {
//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
//TStr file = "sample-wikidata-terms-fragment.nt";
//TStr file = "sample-wikidata-terms.nt";
	TStr file = "SmallTest4.nt";

	FILE* outfile = fopen("frequencies_output.bin", "w");

	InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

	computeFrequenciesPushBCA(file, weigher, outfile);

	fclose(outfile);
//	fclose (f);
}

}

