/*
 * Main.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

#include <utility>
#include <string>
#include <iostream>

#include "KGloVe.h"


int main(int argc, char *argv[]) {
    try {
        const std::string input_file = (argc > 1)? argv[1] : "./testInput/SmallTest2.nt";

        KGloVe::Parameters p;
        //p.graphs.push_back(std::tuple<string, bool, bool>("368303ALL_MergedMultiline_no-empty-lines_sort-uniq_error-boxer.nt", false, true));
        p.graphs.push_back(std::tuple<std::string, bool, bool>(input_file, true, true));
        weigher::UniformWeigher w;
        //PushDownWeigherMap w(readDBPediaPageRanks("pagerank_en_2016-04.tsv"), 0.2);
        p.weighers.push_back(&w);
        p.alphas.push_back(0.7);
        p.epss.push_back(0.00001);
        p.onlyEntities.push_back(false);
        p.includeEdges.push_back(true);
        p.outputPrefix = "output/";

        KGloVe::parametrizedRun(p);

    } catch (const std::exception& e) {
        std::cout << e.what() << std::endl;
        throw e;
    }
    return EXIT_SUCCESS;
}
