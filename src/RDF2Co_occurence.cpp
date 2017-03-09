/*
 * RDF2Co_occurence.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include <iostream>
#include <fstream>
#include "Snap.h"
#include "WeightedPredicate.h"
#include "nTripleParser.h"
#include "BCA.h"
#include "GraphWeigher.h"
#include "utils.h"
#include "boost/dynamic_bitset.hpp"
#include "boost/graph/graph_traits.hpp"
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
//#include <boost/graph/topological_sort.hpp>
#include <utility>                   // for std::pair
#include "PrintTime.h"

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
TVec<TInt> determineBCVcomputeOrderOLD(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

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

/**
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 */
TVec<TInt> determineBCVcomputeOrder(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
	//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
	//We also remove all self edges
	TPt<TNGraph> prunable = TNGraph::New();
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		if (src != dst) {
			prunable->AddEdge(src, dst);
		}
	}
	int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original
	TVec<TInt> finalOrder;
	THashSet<TInt> zeroOutDegNodes;
	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		if (node.GetOutDeg() == 0) {
			zeroOutDegNodes.AddKey(node.GetId());
		}
	}
	while (zeroOutDegNodes.Len() != 0) {
		const TInt n = *(zeroOutDegNodes.BegI());
		zeroOutDegNodes.DelKey(n);
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		TNGraph::TNodeI niter = prunable->GetNI(n);
		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
			int mid = niter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI m = prunable->GetNI(mid);
			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
			if (m.GetOutDeg() == 1) {
				zeroOutDegNodes.AddKey(m.GetId());
			}
		}
		//finalize
		prunable->DelNode(n);
		finalOrder.Add(n);
	}
	//now the more general case including loops is handled

	TMaxPriorityQueue<TInt> highestInDegree;
	//set-up the  highestInDegree PQ,
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
	}
	//algo start
	while (prunable->GetNodes() > 0) {
		while (zeroOutDegNodes.Len() != 0) {
			const TInt n = *(zeroOutDegNodes.BegI());
			zeroOutDegNodes.DelKey(n);
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			TNGraph::TNodeI niter = prunable->GetNI(n);
			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
				int mid = niter.GetInNId(inEdgeNumber);
				TNGraph::TNodeI m = prunable->GetNI(mid);
				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
				if (m.GetOutDeg() == 1) {
					zeroOutDegNodes.AddKey(m.GetId());
				}
			}
			//set indegree of n to 0 in PQueueu. PQueue does not support direct removal
			highestInDegree.SetPriority(n, 0.0);
			//finalize
			prunable->DelNode(n);
			finalOrder.Add(n);
		}
		if (prunable->GetNodes() == 0) {
			break;
		}
		IAssert(highestInDegree.Size() > 0);
		//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
		IAssert(highestInDegree.GetMaxPriority() > 0.0);
		cerr.flush();
		TInt k = highestInDegree.PopMax();
		//add all nodes which will get a zero out degree to set
		IAssert(prunable->IsNode(k));
		TNGraph::TNodeI kiter = prunable->GetNI(k);
		IAssert(kiter.GetOutDeg() > 0);
		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
			int lid = kiter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			if (l.GetOutDeg() == 1) {
				zeroOutDegNodes.AddKey(l.GetId());
			}
		}
		//update the priorities of all nodes k is pointing to
		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
			int lid = kiter.GetOutNId(outEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
		}
		//finalize
		prunable->DelNode(k);
		finalOrder.Add(k);
	}
	IAssert(startSize == finalOrder.Len());
	return finalOrder;
}

namespace {
class my_bitset: public boost::dynamic_bitset<> {

private:
	boost::dynamic_bitset<>::size_type lastFound = 0;
	boost::dynamic_bitset<>::size_type lastSetToTrue = 0;

public:
	boost::dynamic_bitset<>::size_type find_any() {
		if ((*this)[lastSetToTrue]) {
			return lastSetToTrue;
		}
		this->lastFound = this->find_next(lastFound);
		if (lastFound != this->npos) {
			return lastFound;
		} else {
			lastFound = find_first();
			return lastFound;
		}

	}

	dynamic_bitset& setTrueAndRecord(size_type n) {
		this->set(n, true);
		this->lastSetToTrue = n;
		return *this;
	}

};

//does the vector contain all numbers from 0 till all.Len()-1 ?
bool allNumbersIn(TVec<TInt> all) {
	TVec<TInt> copy(all);
	copy.Sort();
	for (int i = 0; i < all.Len(); i++) {
		if (copy[i] != i) {
			return false;
		}
	}
	return true;
}

/**
 * Attempts to determine the order in which most BCA computations can be reused.
 * This is using the heuristic that
 *
 * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
 * 2. if no such node exists, then the node with highest in degree is selected
 *
 * the input graph MUST have nodes indexed 0 till Nodes()-1
 *
 */
TVec<TInt> determineBCVcomputeOrderOptimized(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
	//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
	//We also remove all self edges
	cout << currentTime() << "copying graph" << endl;
	TPt<TNGraph> prunable = TNGraph::New();
	prunable->Reserve(baseGraph->GetNodes(), baseGraph->GetEdges());
	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = baseGraph->BegNI(); NI < baseGraph->EndNI(); NI++) {
		prunable->AddNode(NI.GetId());
	}
	cout << currentTime() << "copying graph - nodes copied" << endl;
	for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
		int src = EI.GetSrcNId();
		int dst = EI.GetDstNId();
		if (src != dst) {
			prunable->AddEdge(src, dst);
		}
	}
	const int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original
	cout << currentTime() << "graph copied" << endl;

	TVec<TInt> finalOrder;
	finalOrder.Reserve(startSize);

	my_bitset zeroOutDegrees;
	zeroOutDegrees.resize(startSize, false);

	boost::dynamic_bitset<> todo;
	todo.resize(startSize, true);

	//vector<bool> zeroOutDegrees;
	//zeroOutDegrees.

	//THashSet<TInt> zeroOutDegNodes;
	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		if (node.GetOutDeg() == 0) {
			zeroOutDegrees[node.GetId()] = true;
		}
	}
	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		//const int n = zeroOutDegrees.find_first();
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		TNGraph::TNodeI niter = prunable->GetNI(n);
		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
			int mid = niter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI m = prunable->GetNI(mid);
			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
			if (m.GetOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord(m.GetId());
				//zeroOutDegrees[m.GetId()] = true;
			}
		}
		//finalize
		prunable->DelNode(n);
		todo[n] = false;
		finalOrder.Add(n);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}
	}

	cout << currentTime() << "After first fast phase, " << finalOrder.Len() << "/" << startSize << " nodes are done, starting iterative phase" << endl;

	//now the more general case including loops is handled
	TMaxPriorityQueue<TInt> highestInDegree;
	//set-up the  highestInDegree PQ,
	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
	}
	//algo start
	while (prunable->GetNodes() > 0) {
		//while (zeroOutDegrees.any()) {
		//	int n = zeroOutDegrees.find_first();
		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			TNGraph::TNodeI niter = prunable->GetNI(n);
			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
				int mid = niter.GetInNId(inEdgeNumber);
				TNGraph::TNodeI m = prunable->GetNI(mid);
				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
				if (m.GetOutDeg() == 1) {
					zeroOutDegrees[m.GetId()] = true;
				}
			}
			//finalize
			prunable->DelNode(n);
			//We do not set indegree of n to 0 in PQueueu. PQueue does not support direct removal, but this is somewhat expensive. We use the to_do bitset instead.
			//			highestInDegree.SetPriority(n, 0.0);

			todo[n] = false;
			finalOrder.Add(n);
			if (finalOrder.Len() % infoFrequency == 0) {
				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
			}
		}
		if (prunable->GetNodes() == 0) {
			break;
		}
		TInt k = -1;
		do {
			IAssert(highestInDegree.Size() > 0);
			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
			k = highestInDegree.PopMax();
		} while (!todo[k]);
		//add all nodes which will get a zero out degree to set
		IAssert(prunable->IsNode(k));
		TNGraph::TNodeI kiter = prunable->GetNI(k);
		IAssert(kiter.GetOutDeg() > 0);
		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
			int lid = kiter.GetInNId(inEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			if (l.GetOutDeg() == 1) {
				zeroOutDegrees[l.GetId()] = true;
			}
		}
		//update the priorities of all nodes k is pointing to
		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
			int lid = kiter.GetOutNId(outEdgeNumber);
			TNGraph::TNodeI l = prunable->GetNI(lid);
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
		}
		//finalize
		prunable->DelNode(k);
		todo[k] = false;
		finalOrder.Add(k);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}

	}
	IAssert(startSize == finalOrder.Len());
	IAssert(todo.find_first() == todo.npos);
	IAssert(allNumbersIn(finalOrder));
	return finalOrder;
}

//int bla() {
//	using namespace boost;
//	// create a typedef for the Graph type
//	typedef adjacency_list<vecS, vecS, bidirectionalS> Graph;
//
//	// Make convenient labels for the vertices
//	int A = 0;
//	int B = 1;
//
//	const int num_vertices = 5;
//
//	// declare a graph object
//	Graph g(num_vertices);
//
//	// add the edges to the graph object
//	add_edge(A, B, g);
//	return 0;
//}
//
///**
// * Some potential future ideas for optimization:
// * do we need to know about both in AND out edges or only about one of them if we only need to know the count of the in/out edges, we can store that in an array!
// */
//
///**
// * Attempts to determine the order in which most BCA computations can be reused.
// * This is using the heuristic that
// *
// * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
// * 2. if no such node exists, then the node with highest in degree is selected
// *
// * the input graph MUST have nodes indexed 0 till Nodes()-1
// *
// * Uses boost graphs as the Snap ones seem to slow when they get large
// */
//TVec<TInt> determineBCVcomputeOrderOptimizedBOOST(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
//	bla();
//
//	using namespace boost;
//// TNGraph: directed graph (single directed edge between an ordered pair of nodes)
////We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
////We also remove all self edges
//
//	typedef adjacency_list<hash_setS, // Store out-edges of each vertex in a unordered set
//			listS,
//			// Store vertex set in a std::vector
//			bidirectionalS // The file dependency graph is directed
//	> Graph;
//
//	typedef graph_traits<Graph> GT;
//	typedef GT::vertex_descriptor Vertex;
//	typedef GT::edge_descriptor Edge;
//	typedef GT::vertex_iterator vertex_iter;
//
////Assumption: basegraph has nodes 0->nodes-1
//
//	cout << currentTime() << "copying graph" << endl;
//
//	Graph prunable_graph(baseGraph->GetNodes());
//
//	const unsigned int startSize = baseGraph->GetNodes();
//	{ //scoping allVertices
//	  //We keep iterators, they will even remain valid on removal, as we us a linked list!
//		std::vector<vertex_iter> allVertices(startSize);
//
//		int i = 0;
//		for (std::pair<vertex_iter, vertex_iter> vp = vertices(prunable_graph); vp.first != vp.second; ++vp.first) {
//			allVertices[i] = vp.first;
//			i++;
//		}
//		IAssert(allVertices.size() == startSize);
//		cout << currentTime() << "copying graph - nodes copied" << endl;
//		for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
//			int src = EI.GetSrcNId();
//			int dst = EI.GetDstNId();
//
//			Vertex srcV = *(allVertices[src]);
//			Vertex dstV = *(allVertices[dst]);
//
//			if (src != dst) {
//				add_edge(srcV, dstV, prunable_graph);
//				//prunable->AddEdge(src, dst);
//			}
//		}
//
//	} //end scoping allVertices
//	baseGraph = NULL; //making sure we do not accidentally write to the original
//	cout << currentTime() << "graph copied" << endl;
//
//	TVec<TInt> finalOrder;
//	finalOrder.Reserve(startSize);
//
//	my_bitset zeroOutDegrees;
//	zeroOutDegrees.resize(startSize, false);
//
//	boost::dynamic_bitset<> todo;
//	todo.resize(startSize, true);
//
//	//we do a first quicker check to eliminate all nodes which are not part of a cycle
//	typename property_map<Graph, vertex_index_t>::type index = get(vertex_index, prunable_graph);
//	for (std::pair<vertex_iter, vertex_iter> vp = vertices(prunable_graph); vp.first != vp.second; ++vp.first) {
//		unsigned int outDegree = out_degree(*vp.first, prunable_graph);
//
//		if (outDegree == 0) {
//			int a = index[*vp.first];
//
//			//zeroOutDegrees[]
//		}
//
//	}
//
////	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
////		if (node.GetOutDeg() == 0) {
////			zeroOutDegrees[node.GetId()] = true;
////		}
////	}
//
////	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;
////
////	unsigned long int n;
////	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
////		//const int n = zeroOutDegrees.find_first();
////		zeroOutDegrees[n] = false;
////		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
////		TNGraph::TNodeI niter = prunable->GetNI(n);
////		for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
////			int mid = niter.GetInNId(inEdgeNumber);
////			TNGraph::TNodeI m = prunable->GetNI(mid);
////			//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
////			if (m.GetOutDeg() == 1) {
////				zeroOutDegrees.setTrueAndRecord(m.GetId());
////				//zeroOutDegrees[m.GetId()] = true;
////			}
////		}
////		//finalize
////		prunable->DelNode(n);
////		todo[n] = false;
////		finalOrder.Add(n);
////		if (finalOrder.Len() % infoFrequency == 0) {
////			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////		}
////	}
////
////	cout << currentTime() << "After first fast phase, " << finalOrder.Len() << "/" << startSize << " nodes are done, starting iterative phase" << endl;
////
////	//now the more general case including loops is handled
////	TMaxPriorityQueue<TInt> highestInDegree;
////	//set-up the  highestInDegree PQ,
////	for (TNGraph::TNodeI node = prunable->BegNI(); node < prunable->EndNI(); node++) {
////		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
////		//We add one to indicate that even a node with 0 in degree is still a valid node.
////		highestInDegree.Insert(node.GetId(), node.GetInDeg() + 1);
////	}
////	//algo start
////	while (prunable->GetNodes() > 0) {
////		//while (zeroOutDegrees.any()) {
////		//	int n = zeroOutDegrees.find_first();
////		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
////			zeroOutDegrees[n] = false;
////			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
////			TNGraph::TNodeI niter = prunable->GetNI(n);
////			for (int inEdgeNumber = 0; inEdgeNumber < niter.GetInDeg(); inEdgeNumber++) {
////				int mid = niter.GetInNId(inEdgeNumber);
////				TNGraph::TNodeI m = prunable->GetNI(mid);
////				//it is enough to check the outdegree. The TNGraph type guarantees that there is only one directed edge between an ordered pair of nodes.
////				if (m.GetOutDeg() == 1) {
////					zeroOutDegrees[m.GetId()] = true;
////				}
////			}
////			//finalize
////			prunable->DelNode(n);
////			//We do not set indegree of n to 0 in PQueueu. PQueue does not support direct removal, but this is somewhat expensive. We use the to_do bitset instead.
////			//			highestInDegree.SetPriority(n, 0.0);
////
////			todo[n] = false;
////			finalOrder.Add(n);
////			if (finalOrder.Len() % infoFrequency == 0) {
////				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////			}
////		}
////		if (prunable->GetNodes() == 0) {
////			break;
////		}
////		TInt k = -1;
////		do {
////			IAssert(highestInDegree.Size() > 0);
////			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
////			k = highestInDegree.PopMax();
////		} while (!todo[k]);
////		//add all nodes which will get a zero out degree to set
////		IAssert(prunable->IsNode(k));
////		TNGraph::TNodeI kiter = prunable->GetNI(k);
////		IAssert(kiter.GetOutDeg() > 0);
////		for (int inEdgeNumber = 0; inEdgeNumber < kiter.GetInDeg(); inEdgeNumber++) {
////			int lid = kiter.GetInNId(inEdgeNumber);
////			TNGraph::TNodeI l = prunable->GetNI(lid);
////			if (l.GetOutDeg() == 1) {
////				zeroOutDegrees[l.GetId()] = true;
////			}
////		}
////		//update the priorities of all nodes k is pointing to
////		for (int outEdgeNumber = 0; outEdgeNumber < kiter.GetOutDeg(); outEdgeNumber++) {
////			int lid = kiter.GetOutNId(outEdgeNumber);
////			TNGraph::TNodeI l = prunable->GetNI(lid);
////			//We add one to indicate that even a node with 0 in degree is still a valid node.
////			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that lid != k
////			highestInDegree.SetPriority(lid, l.GetInDeg() - 1 + 1);
////		}
////		//finalize
////		prunable->DelNode(k);
////		todo[k] = false;
////		finalOrder.Add(k);
////		if (finalOrder.Len() % infoFrequency == 0) {
////			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
////		}
////
////	}
////	IAssert(startSize == finalOrder.Len());
////	IAssert(todo.find_first() == todo.npos);
////	IAssert(allNumbersIn(finalOrder));
//	return finalOrder;
//}

class Node {
public:
	boost::unordered_set<Node*> inedgeSources;
	boost::unordered_set<Node*> outedgeDestination;
	int ID;

	Node(int ID) :
			ID(ID) {			//Intentionally empty
	}

	int getInDeg() {
		return inedgeSources.size();
	}
	int getOutDeg() {
		return outedgeDestination.size();
	}

};

// directed graph (single directed edge between an ordered pair of nodes)
class MyGraph: public std::vector<Node> {
public:
	MyGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {
		cout << currentTime() << "copying graph" << endl;
		const int totalNodes = baseGraph->GetNodes();
		this->reserve(totalNodes);
		for (int i = 0; i < totalNodes; i++) {
			this->emplace_back(i);
		}
		cout << currentTime() << "copying graph - nodes copied" << endl;

		//We throw away the multiplicity of edges between two nodes. For the BCA caching it does not make a difference whether it is needed once or more times by the same other node
		//We also remove all self edges
		for (TNodeEdgeNet<TStr, TStr>::TEdgeI EI = baseGraph->BegEI(); EI < baseGraph->EndEI(); EI++) {
			int src = EI.GetSrcNId();
			int dst = EI.GetDstNId();
			if (src != dst) {
				//boost::unordered_set<Node*>::value_type  == Node*
				//eclipse does not find these defnitions, gcc does...
				boost::unordered_set<Node*>::value_type sourceNode = &(*this)[src];
				boost::unordered_set<Node*>::value_type destinationNode = &(*this)[dst];
				sourceNode->outedgeDestination.insert(destinationNode);
				destinationNode->inedgeSources.insert(sourceNode);
			}
		}
		cout << currentTime() << "graph copied" << endl;
	}
};

///**
// * Attempts to determine the order in which most BCA computations can be reused.
// * This is using the heuristic that
// *
// * 1. nodes with zero out degree (or one for which all outgoing edges points to nodes with precomputed BCAs) should always be computed first
// * 2. if no such node exists, then the node with highest in degree is selected
// *
// * the input graph MUST have nodes indexed 0 till Nodes()-1
// *
// */
TVec<TInt> determineBCVcomputeOrderOptimizedOwnGraph(TPt<TNodeEdgeNet<TStr, TStr> > baseGraph) {

	MyGraph prunable_graph(baseGraph);
	const int startSize = baseGraph->GetNodes();
	baseGraph = NULL; //making sure we do not accidentally write to the original

	TVec<TInt> finalOrder;
	finalOrder.Reserve(startSize);

	my_bitset zeroOutDegrees;
	zeroOutDegrees.resize(startSize, false);

	boost::dynamic_bitset<> todo;
	todo.resize(startSize, true);

	//we do a first quicker check to eliminate all nodes which are not part of a cycle
	for (std::vector<Node>::iterator node = prunable_graph.begin(); node < prunable_graph.end(); node++) {
		if (node->getOutDeg() == 0) {
			zeroOutDegrees[node->ID] = true;
		}
	}

	const int infoFrequency = startSize > 1000000 ? 100000 : 10000;

	unsigned long int n;
	while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
		zeroOutDegrees[n] = false;
		//n will be removed from the graph. Add all nodes which will get zero out degree to the set
		Node& toBeRemoved = prunable_graph[n];

		for (boost::unordered_set<Node*>::iterator dependant = toBeRemoved.inedgeSources.begin(); dependant != toBeRemoved.inedgeSources.end(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			boost::unordered_set<Node*>::value_type toBeRemovedAdress = &toBeRemoved;
			(*dependant)->outedgeDestination.erase(toBeRemovedAdress);
		}
		//finalize, delete node
		//the only things remaining are some bookkeeping:
#ifndef NDEBUG
		//probably not necessary anyhow
		toBeRemoved.inedgeSources.clear();
#endif
		todo[n] = false;
		finalOrder.Add(n);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}
	}

	cout << currentTime() << "After first fast phase, " << finalOrder.Len() << "/" << startSize << " nodes are done, starting iterative phase" << endl;

	//now the more general case including loops is handled
	TMaxPriorityQueue<TInt> highestInDegree;
	//set-up the  highestInDegree PQ,

	for (std::vector<Node>::iterator node = prunable_graph.begin(); node < prunable_graph.end(); node++) {
		//if the outdegree is 0, there is no need to get the node to the highestInDegree as it would be removed immediately again, but there are no zeroOutDegreeNodes at this point
		//We add one to indicate that even a node with 0 in degree is still a valid node.
		highestInDegree.Insert(node->ID, node->getInDeg() + 1);
	}

	//algo start
	while (finalOrder.Len() < startSize) {

		while ((n = zeroOutDegrees.find_any()) != zeroOutDegrees.npos) {
			zeroOutDegrees[n] = false;
			//n will be removed from the graph. Add all nodes which will get zero out degree to the set
			Node& toBeRemoved = prunable_graph[n];

			for (boost::unordered_set<Node*>::iterator dependant = toBeRemoved.inedgeSources.begin(); dependant != toBeRemoved.inedgeSources.end(); dependant++) {
				//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
				if ((*dependant)->getOutDeg() == 1) {
					zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
				}
				//start removing already
				(*dependant)->outedgeDestination.erase(&toBeRemoved);
			}
			//finalize, delete node
			//the only things remaining are some bookkeeping:
#ifndef NDEBUG
			//probably not necessary anyhow
			toBeRemoved.inedgeSources.clear();
#endif
			todo[n] = false;
			finalOrder.Add(n);
			if (finalOrder.Len() % infoFrequency == 0) {
				cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
			}
		}

		if (finalOrder.Len() == startSize) {
			break;
		}
		TInt k = -1;
		do {
			IAssert(highestInDegree.Size() > 0);
			//there are no nodes with zero out degree in the graph left. Attempt to break a cycle by removing the one with highest in degree
			k = highestInDegree.PopMax();
		} while (!todo[k]);
		//k will be removed add all nodes which will get a zero out degree to set
		Node& kNode = prunable_graph[k];

		for (boost::unordered_set<Node*>::iterator dependant = kNode.inedgeSources.begin(); dependant != kNode.inedgeSources.end(); dependant++) {
			//it is enough to check the outdegree being 1. The graph guarantees that there is only one directed edge between an ordered pair of nodes.
			if ((*dependant)->getOutDeg() == 1) {
				zeroOutDegrees.setTrueAndRecord((*dependant)->ID);
			}
			//start removing already
			(*dependant)->outedgeDestination.erase(&kNode);
		}

		//update the priorities of all nodes k is pointing to
		for (boost::unordered_set<Node*>::iterator dest = kNode.outedgeDestination.begin(); dest != kNode.outedgeDestination.end(); dest++) {
			//We add one to indicate that even a node with 0 in degree is still a valid node.
			//Here we use the fact that there are no self edges in the graph. Otherwise we have to make sure that dest.id != k
			highestInDegree.SetPriority((*dest)->ID, (*dest)->getInDeg() - 1 + 1);
			//start removing already
			(*dest)->inedgeSources.erase(&kNode);
		}

		//finalize, the removal is already done in the adjecancy lists
#ifndef NDEBUG
		//probably not necessary anyhow
		kNode.inedgeSources.clear();
		kNode.outedgeDestination.clear();
#endif
		todo[k] = false;
		finalOrder.Add(k);
		if (finalOrder.Len() % infoFrequency == 0) {
			cout << currentTime() << finalOrder.Len() << "/" << startSize << " done" << endl;
		}

	}
	IAssert(startSize == finalOrder.Len());
	IAssert(todo.find_first() == todo.npos);
	IAssert(allNumbersIn(finalOrder));
	return finalOrder;
}
} //end anonymous namspace

void testBCAComputeOrderSpeed(TStr inputFilename, TStr outputFileName) {
	TPair<TPt<TNodeEdgeNet<TStr, TStr> >, THash<TStr, int> > graphAndNodeIndex = n3parser::buildRDFGraphIgnoreLiterals(inputFilename);
	TPt<TNodeEdgeNet<TStr, TStr> > graph = graphAndNodeIndex.Val1;
	cout << currentTime() << " Now computing BCV compute order" << endl;
	TVec<TInt, int> order;
	//order = determineBCVcomputeOrderOptimized(graph);
	order = determineBCVcomputeOrderOptimizedOwnGraph(graph);
	cout << currentTime() << "end determining BCV compute order, writing to file" << endl;

	ofstream myfile(outputFileName.CStr());
	for (TVec<TInt>::TIter iter = order.BegI(); iter < order.EndI(); iter++) {
		myfile << int(*iter) << endl;
	}
	myfile.flush();
	if (!myfile.good()) {
		throw "WTF";
	}
	myfile.close();
	cout << currentTime() << "done writing" << endl;
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

	THash<TInt, BCV> bcvCache;

	int counter = 0;
	for (TVec<TInt>::TIter iter = order.BegI(); iter < order.EndI(); iter++) {
		int i = iter->Val;
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
		BCV combinedbcv = computeBCAIncludingEdgesCached(weightedGraph, focusWordGraphID, bca_alpha, bca_eps, predGraphIDs, bcvCache);
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
			//CREC crec = CREC { word1:contextWordGloveID, word2:focusWordGloveID, val: freq };
			fwrite(&crec, sizeof(CREC), 1, glove_input_file_out);
		}

		TStr nodeLabel = weightedGraph->GetNDat(focusWordGraphID);

		fprintf(glove_vocab_file_out, "%s fakeFrequency\n", nodeLabel.CStr());
		counter++;
		if ((counter % 10000) == 0) {
			cout << "Processed " << counter << "/" << order.Len() << " BCV computations" << endl;
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
//	TStr file = "wikidata-simple-statements-10_000000-sample.nt";
//TStr file = "sample-wikidata-terms-fragment.nt";
//TStr file = "sample-wikidata-terms.nt";
	TStr file = "../../datasets/dbPedia/allData27_30M.nt";
//TStr file = "SmallTest8_multiplePO.nt";

	FILE* glove_input_file_out = fopen("glove_input_file_out.bin", "w");
	FILE* glove_vocab_file_out = fopen("glove_vocab_file_out", "w");

	auto BCAOrderFile = "BCAComputeOrder.txt";
	testBCAComputeOrderSpeed(file, BCAOrderFile);

//InversePredicateFrequencyWeigher weigher = InversePredicateFrequencyWeigher();

	UniformWeigher weigher;

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

