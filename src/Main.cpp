/*
 * Main.cpp
 *
 *  Created on: Nov 29, 2016
 *      Author: cochez
 */

//#include "RandomWalkExperiments.h"
#include "RDF2Co_occurence.h"
#include <iostream>

#include <utility>
#include <unordered_map>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

//int main(int argc, char **argv) {
//	//std::string filename = "wikidata-simple-statements-10_000000-sample.nt";
////	std::string filename = "../../datasets/dbPedia/allData27_30M.nt";
//	std::pair<std::shared_ptr<QuickGraph::LabeledGraph>, std::unordered_map<std::string, int>> g = n3parser::buildRDFGraph(filename);
//	std::cout << g.first->nodes.size() << std::endl;
//}

std::unordered_map<std::string, double> readDBPediaPageRanks(std::string tsvFile) {
	std::unordered_map<std::string, double> ranks;

	std::ifstream infile(tsvFile);

	if (!infile.is_open()) {
		// error! maybe the file doesn't exist.
		std::cerr << "Input file " << tsvFile << " not found, exiting!!" << std::endl;
		exit(7);
	}

	std::string line;

	while (std::getline(infile, line)) {
		if (line.find_first_not_of(' ') == std::string::npos) {
			continue;
		}
		if (line.find_first_of('#') == 0) {
			//comment
			continue;
		}

		std::vector<std::string> SplitVec; // #2: Search for tokens
		boost::split(SplitVec, line, boost::is_any_of("\t"), boost::token_compress_off);
		std::string resource = "<" + SplitVec[0] + ">";
		double rank = boost::lexical_cast<double>(SplitVec[1]);

		ranks[resource] = rank;
	}
	infile.close();
	return ranks;
}
GraphWeigher* GetWeigherPtr(std::string arg)
{
	if(arg == "UniformWeigher")
		return new UniformWeigher();
	else if(arg == "InversePredicateFrequencyWeigher")
		return new InversePredicateFrequencyWeigher();
	else if(arg == "PredicateFrequencyWeigher")
		return new PredicateFrequencyWeigher();
	else if(arg == "ObjectFrequencyWeigher")
		return new ObjectFrequencyWeigher();
	else if(arg == "InverseObjectFrequencyWeigher")
		return new InverseObjectFrequencyWeigher();
	else if(arg == "PredicateObjectFrequencyWeigher")
		return new PredicateObjectFrequencyWeigher();
	else if(arg == "InversePredicateObjectFrequencyWeigher")
		return new InversePredicateObjectFrequencyWeigher();
	else if(arg == "InverseObjectFrequencyWeigherSplitDown")
		return new InverseObjectFrequencyWeigherSplitDown();
//	else if(arg == "PushDownWeigher")
//		return new PushDownWeigher(readDBPediaPageRanks("pagerank_en_2016-04.tsv"));
	else if(arg == "PushDownWeigherMap")
		return new PushDownWeigherMap(readDBPediaPageRanks("../pagerank_en_2016-04.tsv"), 0.2);
//	else if(arg == "SplitDownWeigher")
//		return new SplitDownWeigher(readDBPediaPageRanks("pagerank_en_2016-04.tsv"));
	else
		throw std::runtime_error("Unknown Weigher");
}
int main(int argc, char **argv) {

	std::unique_ptr<GraphWeigher> weigher;

	// Declare the supported options.
	namespace po = boost::program_options;
	po::options_description desc("Options");
	desc.add_options()
		("help,h", "produce help message")
		("dataset", po::value<std::string>()->default_value("allData.nt"), "DataSet like file.nt")
		("remove-literals,l", po::bool_switch()->default_value(false), "Remove Literals")
		("store-vocab,s", po::bool_switch()->default_value(false), "Store vocabulary file")
		("normalize,n", po::bool_switch()->default_value(false), "Normalize the paint values")
		("only-entities", po::bool_switch()->default_value(false), "Include only the entities of the graph")
		("alpha,a", po::value<float>()->default_value(0.3), "Alpha of the BCA Algorithm")
		("epsilon,e", po::value<float>()->default_value(0.00001),"Epsilon of the BCA Algorithm")
		("weigher,w", po::value<std::string>()->default_value("PushDownWeigherMap"), "Choose Weigher")
		;
	po::positional_options_description p;
	p.add("dataset", -1);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).
			options(desc).positional(p).run(), vm);
	po::notify(vm);

	if (vm.count("help")) {
		std::cout << desc << "\n";
		return 1;
	}

	std::cout << "remove literals " << vm["remove-literals"].as<bool>() << std::endl;
	std::cout << "store vocab " << vm["store-vocab"].as<bool>() << std::endl;
	std::cout << "normalize " << vm["normalize"].as<bool>() << std::endl;
	std::cout << "only-entities " << vm["only-entities"].as<bool>() << std::endl;
	std::cout << "dataset: " << vm["dataset"].as<std::string>() << std::endl;
	std::cout << "alpha: " << vm["alpha"].as<float>() << std::endl;
	std::cout << "epsilon: " << vm["epsilon"].as<float>() << std::endl;
	std::cout << "weigher: " << vm["weigher"].as<std::string>() << std::endl;

	try {
		RDF2CO::ParameterizedRun::Parameters p;
		//p.graphs.push_back(std::tuple<std::string, bool, bool>("368303ALL_MergedMultiline_no-empty-lines_sort-uniq_error-boxer.nt", false, true));
		p.graphs.push_back(std::tuple<std::string, bool, bool>(vm["dataset"].as<std::string>(), vm["remove-literals"].as<bool>(), vm["store-vocab"].as<bool>()));
		//UniformWeigher w;
		std::unique_ptr<GraphWeigher> w(std::move(GetWeigherPtr(vm["weigher"].as<std::string>())));
		p.weighers.push_back(std::pair<GraphWeigher&, GraphWeigher&>(*w, *w));
		p.alphas.push_back(vm["alpha"].as<float>());
		p.epss.push_back(vm["epsilon"].as<float>());
		p.normalize.push_back(vm["normalize"].as<bool>());
		p.onlyEntities.push_back(vm["only-entities"].as<bool>());
		RDF2CO::ParameterizedRun::parametrizedUltimateRun(p);
	} catch (char const* str) {
		std::cout << str << std::endl;
		throw str;
	}
	return 0;
}

//int mainRandomWalk(int argc, char **argv) {
//
//
//
//
//
//	//RandomWalkExperiments::performExperiments(1);
//
//	//return 0;
//
//
//
//
//
//
//	if (argc < 2) {
//		std::cerr << "Usage: commandName [number], where number from" << std::endl;
//		std::cerr << "	1. Uniform" << std::endl << "	2. Predicate freq" << std::endl << "	3. Predicate inv. freq" << std::endl << "	4. Object freq" << std::endl << "	5. Object inv. freq" << std::endl << "	6. P-O freq" << std::endl
//				<< "	7. P-O inv. freq" << std::endl << "	8. Pagerank weight" << std::endl << "	9. inv Pagerank weight" << std::endl << "	10. Pagerank split weight" << std::endl << "	11. inv Pagerank split weight" << std::endl
//				<< "	12. Object inv. freq split" << std::endl << "	Object freq split == uniform"
//				<< std::endl;
//		return 1;
//	}
//
//	TStr strategy(argv[1]);
//	int strategyNumber = strategy.GetInt();
//
//	char * outFileName = NULL;
//	if (argc > 2){
//		outFileName = argv[2];
//	}
//
//	RandomWalkExperiments::performExperiments(strategyNumber, outFileName);
//
//
//
//	//RDF2CO::performExperiments();
//	return 0;
//}

