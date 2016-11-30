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

THash<TPair<TStr, TStr>, TInt> createPairedWordIndexTable(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TPair<TStr, TStr>, TInt> table;
	int counter = 0;
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

typedef double real;

typedef struct cooccur_rec {
	int word1;
	int word2;
	real val;
} CREC;

void computeFrequencies(TStr filename, GraphWeigher& weighingStrategy, FILE *fout) {
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


