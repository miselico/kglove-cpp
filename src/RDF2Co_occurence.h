/*
 * RDF2Co_occurence.h
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#ifndef RDF2CO_OCCURENCE_HA_
#define RDF2CO_OCCURENCE_HA_

#include <string>
#include "GraphWeigher.h"
#include <tuple>

namespace RDF2CO {
void performExperiments();

namespace ParameterizedRun {

class Parameters {

public:
	/**
	 * The name of the graph, removeLiterals, should the vocabulary file be stored
	 */
	std::vector<std::tuple<std::string, bool, bool>> graphs;
	std::vector<std::pair<GraphWeigher&, GraphWeigher&>> weighers;
	std::vector<double> alphas;
	std::vector<double> epss;
	std::vector<bool> normalize;
	std::vector<bool> onlyEntities;

	void check() {
		if (!(graphs.size() > 0)) {
			throw "no graphs";
		}
		if (!(weighers.size() > 0)) {
			throw "no weighers";
		}
		if (!(alphas.size() > 0)) {
			throw "no alphas";
		}
		if (!(epss.size() > 0)) {
			throw "no epss";
		}

		if (!(normalize.size() == 1 || (normalize.size() == 2 && normalize[0] != normalize[1]))) {
			throw "normalize can only be [true], [false], [true, false], or [false,true]";
		}
		if (!(onlyEntities.size() == 1 || (onlyEntities.size() == 2 && onlyEntities[0] != onlyEntities[1]))) {
			throw "onlyEntities can only be [true], [false], [true, false], or [false,true]";
		}

	}

};

/**
 * Will perform an 'ultimate' run on each permutation of the parameters.
 *
 * Currently not implemented: use BCA compute order caching in a file. This speeds up if there are many runs and disk reading is faster then recomputing the order.
 *
 */
void parametrizedUltimateRun(Parameters & param);

}
}
#endif /* RDF2CO_OCCURENCE_HA_ */
