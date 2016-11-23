//============================================================================
// Name        : RDFConverter.cpp
// Author      : Michael Cochez
// Version     :
// Copyright   : (c) Michael Cochez
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>

#include "Snap.h"

#include "MurmurHash3.h"

using namespace std;


namespace OLd{
class MyTNodeData {
private:
	TStr val;
public:
	MyTNodeData() {

	}

	MyTNodeData(TStr str) :
			val(str) {
	}
	MyTNodeData(TSIn& SIn) {
		SIn.GetNextLn(this->val);
	}
	void Save(TSOut& SOut) const {
		SOut.PutStr(this->val);
	}
};

//class Test {
//public:

//	void Save(TSOut& s) const {

//	}
///};

class PO {
public:
	TStr P;
	TStr O;
	PO() {

	}

	PO(TStr P, TStr O) :
			P(P), O(O) {
	}
	PO(TSIn& SIn) {
		SIn.GetNextLn(this->P);
		SIn.GetNextLn(this->O);
	}
	void Save(TSOut& SOut) const {
		SOut.PutStrLn(this->P);
		SOut.PutStrLn(this->O);
	}
};

bool DEBUGMODE = false;

// A blank node is at the start of this line, parses it of and returns the remaining of the line.
TPair<TStr, TStr> parseBlankNodeOf(TStr line) {
	int firstUnderScoreColon = line.SearchStr("_:", 0);
	if (firstUnderScoreColon == -1) {
		cerr << "_: expected not found " << line.CStr() << endl;
	}
	TStr atBNStart = line.GetSubStr(firstUnderScoreColon);
	//We assume it is ended by space.

	TStr blankNode = line.LeftOf(' ');
	TStr rest = line.RightOf(' ');
	return TPair<TStr, TStr>(blankNode, rest);
}

TPair<TStr, TStr> parseSubjectOf(TStr line) {
	int firstUnderScoreColon = line.SearchStr("_:", 0);
	int firstAngular = line.SearchCh('<', 0);

	if (firstUnderScoreColon != -1 && firstUnderScoreColon < firstAngular) {
		return parseBlankNodeOf(line);
	}
	//otherwise it is a resource
	line = line.RightOf('<');
	TStr S = line.LeftOf('>');
	return TPair<TStr, TStr>(S, line.RightOf('>'));
}

TStr parseObject(TStr line) {

	//the line now either contains one more resource or a literal
	//Note, there could be a '<' in the string and according to production of N3 also a " in resource

	while (line[0] == ' ') {
		line = line.GetSubStr(1);
	}

	if (line.GetCh(0) == '_') {
		return parseBlankNodeOf(line).Val1;
	} else if (line.GetCh(0) == '<') {
		line = line.RightOf('<');
		return line.LeftOf('>');
	} else if (line.GetCh(0) == '"') {
		//literal
		if (DEBUGMODE && ((line.SearchStr("^^", 0) != -1) || (line.SearchStr("@", 0) != -1))) {
			cerr << "TODO language tags and datatypes are not supported properly, assumed as part of the literal : " << line.CStr() << endl;
		}
		//cut of final dot and spaces
		while (line.GetCh(line.Len() - 1) == ' ' || line.GetCh(line.Len() - 1) == '.') {
			line = line.GetSubStr(0, line.Len() - 2);
		}
		return line;
	} else {
		cerr << " Parsing error, invalid object " << line.CStr() << endl;
		exit(1);
	}

}

TTriple<TStr, TStr, TStr> parsetripleLine(TStr line, THashSet<TStr> & Stringpool) {
	//line = <http://www.wikidata.org/ontology#gcPrecision> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#DatatypeProperty> .
	//line = _:node1ap4j3o3kx4 <http://example.com> <http://example.com#hasname> "Test"

	TPair<TStr, TStr> S_rest = parseSubjectOf(line);

	TStr S = S_rest.Val1;

	//predicate is easy
	line = S_rest.Val2.RightOf('<');
	TStr P = line.LeftOf('>');
	line = line.RightOf('>');

	TStr O = parseObject(line);

	//This saves about 10% memory on a small test. Perhaps more on a larger one.
	S = Stringpool.GetKey(Stringpool.AddKey(S));
	P = Stringpool.GetKey(Stringpool.AddKey(P));
	O = Stringpool.GetKey(Stringpool.AddKey(O));

	return TTriple<TStr, TStr, TStr>(S, P, O);
}

//TPt<TNodeEdgeNet<PO, TFlt> >
void buildGraph(TStr filename) {

//The graph
	TPt<TNodeEdgeNet<PO, TInt> > Net = TNodeEdgeNet<PO, TInt>::New();
//IDs of all nodes which have the given entity as O
	TStrHash<TVec<int> > nodesForObject = TStrHash<TVec<int> >();

//temporary keep track of all the nodes
	THash<TPair<TStr, TStr>, int> addedNodes = THash<TPair<TStr, TStr>, int>();

	PSIn FInPt = TFIn::New(filename);
	TStr line;
	int count = 0;
	TStr Xresource = TStr("<{}>");

	THashSet<TStr> stringpool = THashSet<TStr>();

	while (FInPt->GetNextLn(line)) {
		TTriple<TStr, TStr, TStr> values = parsetripleLine(line, stringpool);
		//each resource is at least added with an empty inlink
		{ //scoped to avoid mistakes
			TPair<TStr, TStr> XSstr = TPair<TStr, TStr>(Xresource, values.Val1);
			if (!addedNodes.IsKey(XSstr)) {
				int nodeID = Net->AddNode(-1, PO(Xresource, values.Val1));
				addedNodes.AddDat(XSstr, nodeID);
				//If understood correctly, this will instantiate the vector if it does not exist yet.
				nodesForObject.AddDat(values.Val1).Add(nodeID);
			}
		}
		{ //scoped to avoid mistakes
			TPair<TStr, TStr> POstr = TPair<TStr, TStr>(values.Val2, values.Val3);
			if (!addedNodes.IsKey(POstr)) {
				int nodeID = Net->AddNode(-1, PO(values.Val2, values.Val3));
				addedNodes.AddDat(POstr, nodeID);
				nodesForObject.AddDat(values.Val3).Add(nodeID);
			}
		}
		count++;
		if (count % 100000 == 0) {
			cout << "Processed " << count << " lines" << endl;
		}
	}

	//int a;
	//cin >> a;

	//for (TNodeEdgeNet<PO, TInt>::TNodeI NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
	//printf("node id %d P:%s O:%s with out-degree %d and in-degree %d\n", NI.GetId(), NI.GetDat().P.CStr(), NI.GetDat().O.CStr(), NI.GetOutDeg(), NI.GetInDeg());
	//}

//PNGraph Graph = TSnap::LoadEdgeList<PNGraph>(InFNm);

}

//TPt<TNodeEdgeNet<PO, TFlt> >
void buildSimpleGraph(TStr filename) {

//The graph
	TPt<TNodeEdgeNet<TStr, TStr> > Net = TNodeEdgeNet<TStr, TStr>::New();

//temporary keep track of all the nodes
	THash<TStr, int> addedNodes = THash<TStr, int>();

	PSIn FInPt = TFIn::New(filename);
	TStr line;
	int count = 0;

	THashSet<TStr> stringpool = THashSet<TStr>();

	while (FInPt->GetNextLn(line)) {
		TTriple<TStr, TStr, TStr> values = parsetripleLine(line, stringpool);
		int subjectIndex;
		if (addedNodes.IsKeyGetDat(values.Val1, subjectIndex)) {
			//nothing, subjectIndex now contains the node ID
		} else {
			//add new Node, save index
			subjectIndex = Net->AddNode(-1, values.Val1);
			addedNodes.AddDat(values.Val1, subjectIndex);
		}

		int objectIndex;
		if (addedNodes.IsKeyGetDat(values.Val3, objectIndex)) {
			//nothing, objectIndex now contains the node ID
		} else {
			//add new Node, save index
			objectIndex = Net->AddNode(-1, values.Val3);
			addedNodes.AddDat(values.Val3, objectIndex);
		}

		//add edge

		//Net->AddEdge(subjectIndex, objectIndex, -1, values.Val2);

		count++;
		if (count % 100000 == 0) {
			cout << "Processed " << count << " lines" << endl;
		}
	}

	stringpool.Clr(true, -1);

//	cout << Net->BegNI().GetDat().CStr() << endl;
//
//
//	for (TNodeEdgeNet<TStr, TStr>::TNodeI NI = Net->BegNI(); NI < Net->EndNI(); NI++) {
//		for(int i = 0; i < NI.GetOutDeg(); i++){
//			TStr outedge = NI.GetOutEDat(i);
//			cout << outedge.CStr();
//		}
//		//printf("node id %d P:%s O:%s with out-degree %d and in-degree %d\n", NI.GetId(), NI.GetDat().P.CStr(), NI.GetDat().O.CStr(), NI.GetOutDeg(), NI.GetInDeg());
//	}
//
//	int a;
//	cin >> a;

//PNGraph Graph = TSnap::LoadEdgeList<PNGraph>(InFNm);

}

const int MURMURSEED = 65765745;

const char base16[] = { '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'B', 'C', 'D', 'E' , 'F'};

TStr myhash(TStr in) {
	char* cstr = in.CStr();
	int len = in.Len();
	unsigned char* val = (unsigned char*) malloc(16);
	MurmurHash3_x86_128(cstr, len, MURMURSEED, val);
	char* stringVal = (char*) malloc(32 + 1);
	for (int i = 0; i < 16; ++i) {
		unsigned char c = val[i];
		unsigned char c1_index = c >> 4;
		unsigned char c2_index = c & 15;
		char c1 = base16[c1_index];
		char c2 = base16[c2_index];
		stringVal[2 * i] = c1;
		stringVal[2 * i + 1] = c2;
	}
	stringVal[32] = 0;
	return TStr(stringVal);
}

int mainOLD() {
	for (int var = 0; var < 100; ++var) {
		TStr v = myhash("test");

		cout << v.CStr() << endl;

	}

	return 0;

	buildSimpleGraph("wikidata-simple-statements-10_000000-sample.nt");

	//buildGraph("sample-wikidata-terms-fragment.nt");

	//buildGraph("sample-wikidata-terms.nt");

	return 0;

	cout << "!!!Hello World!!!" << endl; // prints !!!Hello World!!!
	cout << sizeof(TInt) << endl;
// create a graph
	PNGraph Graph = TNGraph::New();
	Graph->AddNode(1);
	Graph->AddNode(5);
	Graph->AddNode(32);
	Graph->AddEdge(1, 5);
	Graph->AddEdge(5, 1);
	Graph->AddEdge(5, 32);
	for (TNGraph::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
		printf("node id %d with out-degree %d and in-degree %d\n", NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
	}

	TPt<TNodeEdgeNet<TStr, TInt> > Net = TNodeEdgeNet<TStr, TInt>::New();

	TStr str = TStr("value");
//MyTNodeData data = MyTNodeData(str);
//int nodeID = Net->AddNode(-1, str);
//	Net->AddEdge(nodeID, nodeID, -1, data);

	return 0;
}
}
