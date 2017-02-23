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
 * For each unique predicate in the graph add a unique ID it to the returned hash.
 *
 * The returned vector contains the strings in the same order as their IDs
 *
 * The used IDs range from graph.Nodes() till graph.Nodes+(numberOfUniquePredicates=returnValue.Len()) exclusive
 */
TPair<TVec<TStr>, THash<TStr, int>> computePredicateIDs(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TStr, int> preds;
	TVec<TStr> labels;
	unsigned int currentID = graph->GetNodes();
	for (int i = 0; i < graph->GetEdges(); ++i) {
		TStr label = graph->GetEDat(i);
		if (!preds.IsKey(label)) {
			preds.AddDat(label, currentID);
			labels.Add(label);
			currentID++;
		}
	}
	return TPair<TVec<TStr>, THash<TStr, int>>(labels, preds);
}

namespace {
static TStr RDF_TYPE("<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>");
static TStr OWL_THING("<http://www.w3.org/2002/07/owl#Thing>");

/**
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 */
TVec<TInt> determineBCVcomputeOrder(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	TPt<TNEGraph> prunable = TNEGraph::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		prunable->AddEdge(EI.GetSrcNId(), EI.GetDstNId());
	}
	baseGraph = NULL; //making sure we do not accidentally write to the original
	TVec<TInt> result;

	while (prunable->GetNodes() > 0) {

		TVec<TInt> withZeroOutDegree;
		for (TNEGraph::TNodeI iter = prunable->BegNI(); iter < prunable->EndNI(); iter++) {
			if (iter.GetOutDeg() == 0) {
				withZeroOutDegree.Add(iter.GetId());
			}
		}
		if (withZeroOutDegree.Len() > 0) {
			//found some
			result.AddV(withZeroOutDegree);
			//remove from graph
			for (TVec<TInt>::TIter ID = withZeroOutDegree.BegI(); ID < withZeroOutDegree.EndI(); ID++) {
				prunable->DelNode(ID->Val);
			}
			continue;
		}
		//none with zero out degree. Take the one with highest indegree.
		int highestIndegree = -1;
		int withHighestInDegree = -1;
		for (TNEGraph::TNodeI iter = prunable->BegNI(); iter < prunable->EndNI(); iter++) {
			if (iter.GetInDeg() > highestIndegree) {
				highestIndegree = iter.GetInDeg();
				withHighestInDegree = iter.GetId();
			}
		}
		if (withHighestInDegree == -1) {
			cerr << "error, none found with in degree";
			exit(5);
		}
		result.Add(withHighestInDegree);
		prunable->DelNode(withHighestInDegree);
	}
	return result;
}

void testBCAComputeOrderSpeed(TStr filename) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	cout << "computing BCV compute order" << endl;
	TVec<TInt, int> order;
	order = determineBCVcomputeOrder(graph);
	cout << "end determining BCV compute order" << endl;
}

}

/**
 *
 * Compute the BCA score for each pair in the graph under the given weighing strategy.
 * Additionally, adds a score for each edge as well.
 *
 * Outputs the score as a sparse matrix which can be fed to glove.
 *
 *
 *
 */
void computeFrequenciesIncludingEdges(TStr filename, GraphWeigher& weighingStrategy, double bca_alpha, double bca_eps, FILE * glove_input_file_out, FILE * glove_vocab_file_out, bool normalize,
		bool onlyEntities) {

	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph;
	TPair<TVec<TStr>, THash<TStr, int>> predLabelsAndIDs;
	TVec<TInt> order;

	{ //scoping to save on total memory
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
		order = determineBCVcomputeOrder(graph);
		predLabelsAndIDs = computePredicateIDs(graph);
		weightedGraph = weighingStrategy.weigh(graph);
	}
	THash<TStr, int> predGraphIDs = predLabelsAndIDs.Val2;

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {
		TNodeEdgeNet<TStr, WeightedPredicate>::TNodeI candidateNode = weightedGraph->GetNI(i);
//		//only take specific one:
//		if (candidateNode.GetDat() != "<http://dbpedia.org/ontology/Province>"){
//			continue;
//		}

		if (onlyEntities) {
			bool accept = false;
			for (int outEdgeNr = 0; outEdgeNr < candidateNode.GetOutDeg(); ++outEdgeNr) {
				TStr predicate = candidateNode.GetOutEDat(outEdgeNr).P();
				if (predicate == RDF_TYPE) {
					TStr object = candidateNode.GetOutNDat(outEdgeNr);
					if (object == OWL_THING) {
						accept = true;
					}
				}
			}
			if (!accept) {
				continue;
			}
		}
		const int focusWordGraphID = i;
		BCV combinedbcv = computeBCAIncludingEdges(weightedGraph, focusWordGraphID, bca_alpha, bca_eps, predGraphIDs);
		const int focusWordGloveID = graphIDToGloveID(i);
		if (normalize) {
			combinedbcv.removeEntry(i);
			combinedbcv.normalizeInPlace();
		}
		for (THash<TInt, TFlt>::TIter iter = combinedbcv.BegI(); iter < combinedbcv.EndI(); iter++) {
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

//still need to write all predicates to the vocab file

	TVec<TStr> predicateLabels = predLabelsAndIDs.Val1;
	for (TVec<TStr>::TIter it = predicateLabels.BegI(); it < predicateLabels.EndI(); it++) {
		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", it->CStr());
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
	TPair<TVec<TStr>, THash<TStr, int>> predLabelsAndIDs;

	{ //scoping to save on total memory
		TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(filename);
		TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

		predLabelsAndIDs = computePredicateIDs(graph);

		{ //scoping to save on total memory
			TPt<TNodeEdgeNet<TStr, TStr> > reversed = reverseGraph(graph);
			weightedReverseGraph = reverseWeighingStrategy.weigh(reversed);
		}
		weightedGraph = weighingStrategy.weigh(graph);
	}

	THash<TStr, int> predGraphIDs = predLabelsAndIDs.Val2;

//int predicateIDGloveOffset = weightedGraph->GetNodes();

	for (int i = 0; i < weightedGraph->GetNodes(); ++i) {

		const int focusWordGraphID = i;
		BCV combinedbcv = computeBCAIncludingEdges(weightedGraph, focusWordGraphID, bca_alpha, bca_eps, predGraphIDs);
		BCV combinedreversedBCVs = computeBCAIncludingEdges(weightedReverseGraph, focusWordGraphID, bca_alpha, bca_eps, predGraphIDs);

		const int focusWordGloveID = graphIDToGloveID(i);
		combinedbcv.add(combinedreversedBCVs);

		if (normalize) {
			combinedbcv.removeEntry(i);
			combinedbcv.normalizeInPlace();
		}
		for (THash<TInt, TFlt>::TIter iter = combinedbcv.BegI(); iter < combinedbcv.EndI(); iter++) {
			int contextWordGraphID = iter.GetKey();
			int contextWordGloveID = graphIDToGloveID(contextWordGraphID);
			double freq = iter.GetDat();
			CREC crec = CREC { word1:focusWordGloveID, word2:contextWordGloveID, val: freq };
			fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
		}

		TStr nodeLabel = weightedGraph->GetNDat(focusWordGraphID);

		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", nodeLabel.CStr());
		if (i % 1000 == 0) {
			cout << i << "/" << weightedGraph->GetNodes() << " = " << i / float(weightedGraph->GetNodes()) << endl;
		}
	}
//still need to write all predicates to the vocab file

	TVec<TStr> predicateLabels = predLabelsAndIDs.Val1;
	for (TVec<TStr>::TIter it = predicateLabels.BegI(); it < predicateLabels.EndI(); it++) {
		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", it->CStr());
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


	testBCAComputeOrderSpeed(file);

//InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

//	UniformWeigher weigher;

//computeFrequenciesPushBCA(file, weigher, outfile);

//	bool normalize = true;
//	bool onlyEntities = false;
//	computeFrequenciesIncludingEdges(file, weigher, 0.80, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize, onlyEntities);
//computeFrequencies(file, weigher, 0.50, 0.05, glove_input_file_out, glove_vocab_file_out);

//computeFrequenciesIncludingEdgesTheUltimate(file, weigher, weigher, 0.70, 0.0039, glove_input_file_out, glove_vocab_file_out, normalize);

	fclose(glove_input_file_out);
	fclose(glove_vocab_file_out);
}

}

