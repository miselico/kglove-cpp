/*
 * BCA.h
 *
 *  Created on: Nov 24, 2016
 *      Author: cochez
 */

#ifndef BCA_HA_
#define BCA_HA_

#include <string>
#include <unordered_map>
#include <memory>
#include "graph/LabeledGraph.h"
#include <utility>
#include "boost/flyweight.hpp"

typedef boost::flyweight<std::string> flyString;

//sparse vector representing the approx pagerank
class BCV: public std::unordered_map<unsigned int, double> {
//	BCV(BCV &other){
//
//	}
public:
	BCV(){

	}
	std::string toString(const std::shared_ptr<QuickGraph::LabeledGraph> network);
	void fixPaint(unsigned int ID, double amount);
	void removeEntry(unsigned int ID);
	//This function normalizes the vector such that pageranks add up to 1 IN PLACE
	void normalizeInPlace();
	void add(BCV & other);
};


BCV computeBCA(std::shared_ptr<QuickGraph::LabeledGraph> graph, int b_ID, double alpha, double eps);

BCV computeBCACached(std::shared_ptr<QuickGraph::LabeledGraph> network, int b_ID, double alpha, double eps, std::unordered_map<unsigned int, BCV> & bcvCache);

BCV computeBCAIncludingEdges(std::shared_ptr<QuickGraph::LabeledGraph> network, int b_ID, double alpha, double eps, const  std::unordered_map<flyString,unsigned int> & predIDs);

BCV computeBCAIncludingEdgesCached(std::shared_ptr<QuickGraph::LabeledGraph> network, int b_ID, double alpha, double eps, const std::unordered_map<flyString,unsigned int> & predIDs, std::unordered_map<unsigned int, BCV> & bcvCache);

//class PBCV : public std::unordered_map< std::pair<int, int> , double>{
//public:
//	std::string toString(const std::shared_ptr<QuickGraph::LabeledGraph> network);
//	void fixPaint(std::pair<int, int>  pred_obj_pair, double amount);
//};
//
//
////PBCA = Pushed Bookmar cocloring algorithm
//PBCV computePBCA(std::shared_ptr<QuickGraph::LabeledGraph> network, int b_ID, double alpha, double eps);

#endif /* BCA_HA_ */
