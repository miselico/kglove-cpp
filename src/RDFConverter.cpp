//============================================================================
// Name        : RDFConverter.cpp
// Author      : Michael Cochez
// Version     :
// Copyright   : (c) Michael Cochez
// Description : Hello World in C++, Ansi-style
//============================================================================

#include "RDFConverter.h"

using namespace std;

TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > buildRDFGraph(TStr filename) {

//The graph
	TPt<TNodeEdgeNet<TStr, TStr> > Net = TNodeEdgeNet<TStr, TStr>::New();

//temporary keep track of all the nodes
	THash<TStr, int> addedNodes = THash<TStr, int>();

	PSIn FInPt = TFIn::New(filename);
	TStr line;
	int count = 0;
//To reuse the strings in the node/edge labels
	THashSet<TStr> stringpool = THashSet<TStr>();

	while (FInPt->GetNextLn(line)) {
		Triple values = parsetripleLine(line);
		//This saves about 10% memory on a small test. Perhaps more on a larger one.
		values = Triple(stringpool.GetKey(stringpool.AddKey(values.S())), stringpool.GetKey(stringpool.AddKey(values.P())), stringpool.GetKey(stringpool.AddKey(values.O())));

		int subjectIndex = 0;
		if (addedNodes.IsKeyGetDat(values.S(), subjectIndex)) {
			//nothing, subjectIndex now contains the node ID
		} else {
			//add new Node, save index
			subjectIndex = Net->AddNode(-1, values.S());
			addedNodes.AddDat(values.S(), subjectIndex);
		}

		int objectIndex = 0;
		if (addedNodes.IsKeyGetDat(values.O(), objectIndex)) {
			//nothing, objectIndex now contains the node ID
		} else {
			//add new Node, save index
			objectIndex = Net->AddNode(-1, values.O());
			addedNodes.AddDat(values.O(), objectIndex);
		}
		//add edge
		Net->AddEdge(subjectIndex, objectIndex, -1, values.P());

		count++;
		if (count % 100000 == 0) {
			cout << "Processed " << count << " lines" << endl;
		}
	}

	stringpool.Clr(true, -1);

	return TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> >(Net, addedNodes);
}

// weighing strategies
TPt<TNodeEdgeNet<TStr, WeightedPredicate> > uniformWeigher(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}

	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int outDegree = baseGraph->GetNI(EI.GetSrcNId()).GetOutDeg();
		double weight = 1.0 / (double) outDegree;
		TStr label = EI.GetDat();
		const int ID = EI.GetId();
		const int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		const WeightedPredicate pred = WeightedPredicate(label, weight);
		newNet->AddEdge(src, dst, ID, pred);
	}

	return newNet;
}

TPt<TNodeEdgeNet<TStr, WeightedPredicate> > inverseFrequencyWeigher(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > newNet = TNodeEdgeNet<TStr, WeightedPredicate>::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		newNet->AddNode(NI.GetId(), NI.GetDat());
	}

	//count freq of each property:
	THash<TStr, TInt> absolute_freq = THash<TStr, TInt>();
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		TStr prop = EI.GetDat();
		TInt start = 0;
		absolute_freq.IsKeyGetDat(prop, start);
		absolute_freq.AddDat(prop, start + 1);
	}
	//inverse all freq
	THash<TStr, TFlt> inverse_freq = THash<TStr, TFlt>();
	for (THash<TStr, TInt>::TIter iter = absolute_freq.BegI(); iter < absolute_freq.EndI(); iter++) {
		inverse_freq.AddDat(iter.GetKey(), 1.0 / iter.GetDat());
	}

	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {

		int node_i_outdeg = NI.GetOutDeg();

		double totalWeight = 0;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			TStr edgeData = NI.GetOutEDat(outEdge);
			totalWeight += inverse_freq.GetDat(edgeData);
		}
		double totalWeightInverse = 1.0 / totalWeight;
		for (int outEdge = 0; outEdge < node_i_outdeg; ++outEdge) {
			int j = NI.GetOutNId(outEdge);
			TStr label = NI.GetOutEDat(outEdge);
			double normalized_weight = inverse_freq.GetDat(label) * totalWeightInverse;
			const WeightedPredicate pred = WeightedPredicate(label, normalized_weight);
			newNet->AddEdge(NI.GetId(), j, NI.GetOutEId(outEdge), pred);
		}

	}

	return newNet;
}

THash<TPair<TStr, TStr>, TInt> createPairedWordIndexTable(TPt<TNodeEdgeNet<TStr, TStr> > graph) {
	THash<TPair<TStr, TStr>, TInt> table = THash<TPair<TStr, TStr>, TInt>();
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

void computeFrequencies(TStr filename, TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weighingStrategy(TPt<TNodeEdgeNet<TStr, TStr> >), FILE *fout) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = buildRDFGraph(filename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;

	THash<TStr, int> wordIndexTable = graphAndNodeIndex.Val2;

	THash<TPair<TStr, TStr>, TInt> pairwordIndexTable = createPairedWordIndexTable(graph);
	TPt<TNodeEdgeNet<TStr, WeightedPredicate> > weightedGraph = weighingStrategy(graph);

	graph.Clr();

	cout << "done weighing" << endl;

	cerr << "TODO : check - does the indexing for glove have to start from 1 or 0 ??"
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

	fclose(fout);
//	TTmStopWatch w = TTmStopWatch(true);
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

int main() {

	//TStr file = "wikidata-simple-statements-1_000000-sample.nt";
//TStr file = "sample-wikidata-terms-fragment.nt";
//TStr file = "sample-wikidata-terms.nt";
	TStr file = "SmallTest3.nt";

	char* outfile = "frequencies_output.bin";

	computeFrequencies(file, inverseFrequencyWeigher, fopen(outfile,"w"));

	return 0;
}

