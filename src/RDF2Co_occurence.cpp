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
#include "utils.h"

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

/*
 * The wordIndexTable is indexed from 0 as it should be, but we need IDs which start from 1. This function abstracts this away
 */

int graphIDToGloveID(int graphID) {
	return graphID + 1;
}

/**
 *
 * Compute the BCA score for each pair in the graph under the given weiging strategy.
 *
 * Outputs the score as a sparse matrix which can be fed to glove.
 *
 */
void computeFrequencies(TStr filename, GraphWeigher& weighingStrategy, double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {

		int focusWordGraphID = i;
		BCV bcv = computeBCA(weightedGraph, focusWordGraphID, bca_alpha, bca_eps);

		int focusWordGloveID = graphIDToGloveID(i);

		for (THash<TInt, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
			int contextWordGraphID = iter.GetKey();
			int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
			double freq = iter.GetDat();
			CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
			fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
		}

		TStr nodeLabel = weightedGraph->GetNDat(focusWordGraphID);

		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", nodeLabel.CStr());
		if (i % 1000 == 0) {
			cout << i / float(weightedGraph->GetNodes()) << endl;
		}
	}
}

/**
 *
 * Compute the BCA score for each pair in the graph under the given weighing strategy.
 * Additionally, adds a score for each edge as well.
 *
 * Outputs the score as a sparse matrix which can be fed to glove.
 *
 */
void computeFrequenciesIncludingEdges(TStr filename, GraphWeigher& weighingStrategy, double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy.weigh(graph);

	int predicateIDGloveOffset = weightedGraph->GetNodes();

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {

		const int focusWordGraphID = i;
		TPair<BCV, BCV> bcvs = computeBCAIncludingEdges(weightedGraph, focusWordGraphID, bca_alpha, bca_eps);
		const int focusWordGloveID = graphIDToGloveID(i);
		{ //scoping bcv
			BCV bcv = bcvs.Val1;
			if (normalize) {
				bcv.removeEntry(i);
				bcv.normalizeInPlace();
			}

			for (THash<TInt, TFlt>::TIter iter = bcv.BegI(); iter < bcv.EndI(); iter++) {
				int contextWordGraphID = iter.GetKey();
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}
		}
		{ //scoping bcv_predicates
			BCV bcv_predicates = bcvs.Val2;

			if (normalize) {
				//There is no need to remove anything, as the node under consideration is not in the predicate list.
				bcv_predicates.normalizeInPlace();
			}

			for (THash<TInt, TFlt>::TIter iter = bcv_predicates.BegI(); iter < bcv_predicates.EndI(); iter++) {
				int contextPredicateWordGraphID = iter.GetKey();
				int contextPredicateWordGloveID = graphIDToGloveID(contextPredicateWordGraphID + predicateIDGloveOffset);

				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextPredicateWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}

		}
		TStr nodeLabel = weightedGraph->GetNDat(focusWordGraphID);

		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", nodeLabel.CStr());
		if (i % 1000 == 0) {
			cout << i / float(weightedGraph->GetNodes()) << endl;
		}
	}

	//still need to write all predicates to the vocab file
	for (int i = 0; i < weightedGraph->GetEdges(); ++i) {
		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", weightedGraph->GetEDat(i).Val1.CStr());
	}
}

/**
 *
 * Compute the BCA score for each pair in the graph under the given weighing strategy.
 * Additionally, adds a score for each edge as well.
 *
 * Furthermore, it also performs a reverse walk and adds the result of that to the BCVs
 * The reverse walk can be performed according to a different weighing strategy
 *
 * Outputs the score as a sparse matrix which can be fed to glove.
 *
 */
void computeFrequenciesIncludingEdgesTheUltimate(TStr filename, GraphWeigher& weighingStrategy, GraphWeigher & reverseWeighingStrategy, double bca_alpha, double bca_eps, FILE * glove_input_file_out,
		FILE * glove_vocab_file_out, bool normalize) {

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph;
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedReverseGraph;

	{ //scoping to save on total memory
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
		{ //scoping to save on total memory
			TPt<TNodeEdgeNet<TStr, TStr> > reversed = reverseGraph(graph);
			weightedReverseGraph = reverseWeighingStrategy.weigh(reversed);
		}
		weightedGraph = weighingStrategy.weigh(graph);
	}

	int predicateIDGloveOffset = weightedGraph->GetNodes();

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {

		const int focusWordGraphID = i;
		TPair<BCV, BCV> bcvs = computeBCAIncludingEdges(weightedGraph, focusWordGraphID, bca_alpha, bca_eps);
		TPair<BCV, BCV> reversedBCVs = computeBCAIncludingEdges(weightedReverseGraph, focusWordGraphID, bca_alpha, bca_eps);

		const int focusWordGloveID = graphIDToGloveID(i);
		{ //scoping bcv to avoid programming err.
			BCV partialbcv = bcvs.Val1;
			BCV partialreversedBCV = reversedBCVs.Val1;
			partialbcv.add(partialreversedBCV);
			if (normalize) {
				partialbcv.removeEntry(i);
				partialbcv.normalizeInPlace();
			}
			BCV combined = partialbcv;

			for (THash<TInt, TFlt>::TIter iter = combined.BegI(); iter < combined.EndI(); iter++) {
				int contextWordGraphID = iter.GetKey();
				int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}
		}
		{ //scoping bcv_predicates to avoid programming err.
			BCV partialbcv_predicates = bcvs.Val2;
			BCV partial_reversed_bcv_predicates = reversedBCVs.Val2;
			partialbcv_predicates.add(partial_reversed_bcv_predicates);
			if (normalize) {
				//There is no need to remove anything, as the node under consideration is not in the predicate list.
				partialbcv_predicates.normalizeInPlace();
			}
			BCV combined_predictes = partialbcv_predicates;

			for (THash<TInt, TFlt>::TIter iter = combined_predictes.BegI(); iter < combined_predictes.EndI(); iter++) {
				int contextPredicateWordGraphID = iter.GetKey();
				int contextPredicateWordGloveID = graphIDToGloveID(contextPredicateWordGraphID + predicateIDGloveOffset);

				double freq = iter.GetDat();
				CREC crec = CREC { word1:focusWordGloveID, word2:contextPredicateWordGloveID, val: freq };
				fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
			}

		}
		TStr nodeLabel = weightedGraph->GetNDat(focusWordGraphID);

		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", nodeLabel.CStr());
		//if (i % 1000 == 0) {
			cout << i << "/" << weightedGraph->GetNodes() <<  " = " << i / float(weightedGraph->GetNodes()) << endl;
		//}
	}

	//still need to write all predicates to the vocab file
	for (int i = 0; i < weightedGraph->GetEdges(); ++i) {
		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", weightedGraph->GetEDat(i).Val1.CStr());
	}
}

/**
 * Implementation incomplete!!!
 */
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
	TStr file = "../../datasets/aifb_fixed_complete.nt";
	//	TStr file = "SmallTest4.nt";

	FILE* glove_input_file_out = fopen("glove_input_file_out.bin", "w");
	FILE* glove_vocab_file_out = fopen("glove_vocab_file_out", "w");

	//InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

	UniformWeigher weigher;

	//computeFrequenciesPushBCA(file, weigher, outfile);

	bool normalize = true;
	//computeFrequenciesIncludingEdges(file, weigher, 0.80, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize);
	//computeFrequencies(file, weigher, 0.50, 0.05, glove_input_file_out, glove_vocab_file_out);

	computeFrequenciesIncludingEdgesTheUltimate(file, weigher, weigher, 0.70, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize);

	fclose(glove_input_file_out);
	fclose(glove_vocab_file_out);
}

}

