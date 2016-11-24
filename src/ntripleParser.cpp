/*
 * ntripleParser.cpp
 *
 *  Created on: Nov 23, 2016
 *      Author: cochez
 */

#include "nTripleParser.h"

using namespace std;

namespace {

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
	TStr S = "<" + line.LeftOf('>') + ">";
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
		return "<" + line.LeftOf('>') + ">";
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

}

Triple parsetripleLine(TStr line) {
	//line = <http://www.wikidata.org/ontology#gcPrecision> <http://www.w3.org/1999/02/22-rdf-syntax-ns#type> <http://www.w3.org/2002/07/owl#DatatypeProperty> .
	//line = _:node1ap4j3o3kx4 <http://example.com> <http://example.com#hasname> "Test"

	TPair<TStr, TStr> S_rest = parseSubjectOf(line);

	TStr S = S_rest.Val1;

	//predicate is easy
	line = S_rest.Val2.RightOf('<');
	TStr P = line.LeftOf('>');
	P = "<" + P + ">";
	line = line.RightOf('>');

	TStr O = parseObject(line);



	return Triple(S, P, O);
}

