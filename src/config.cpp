/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <iostream>
#include <unordered_map>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>

#include "helpers.h"

//#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif

namespace cfg {
	std::unordered_map<char, size_t> AAmasses;
	std::map<char, std::vector<int> > allSitesMods;
	std::map<char, std::vector<int> > nTermMods;
	std::map<char, std::vector<int> > cTermMods;

	std::vector<size_t> AAmassDifferences;

	void loadConfig(std::string alphabetFilename, std::string modificationsFilename) {
		std::ifstream alphabetFile(alphabetFilename);
		std::ifstream modificationsFile(modificationsFilename);

		AAmasses.clear();
		AAmassDifferences.clear();
		allSitesMods.clear();
		nTermMods.clear();
		cTermMods.clear();

		if (!alphabetFile.good()) {
			std::cout << "ERROR: " << alphabetFilename << " not found. Switch to default amino acid masses." << std::endl;
			AAmasses = {
				{'X',     0},	// unknown amino acid
				{'*',     0},	// translation stop
				{'-',     0},	// gap of unknown length
				{'G',  5702},
				{'A',  7104},
				{'S',  8703},
				{'P',  9705},
				{'V',  9907},
				{'T', 10105},
				{'C', 10301},
				{'I', 11308},
				{'L', 11308},
				{'N', 11404},
				{'D', 11503},
				{'Q', 12806},
				{'K', 12809},
				{'E', 12904},
				{'M', 13104},
				{'H', 13706},
				{'F', 14707},
				{'U', 15004},
				{'R', 15610},
				{'Y', 16306},
				{'W', 18608}
			};
		} else {
			std::string line;
			char a;
			size_t mass;
			while (std::getline(alphabetFile, line)) {
				std::istringstream iss(line);
				if (line.front() == '#') continue;
				if (!(iss >> a >> mass)) {
					LOG("Load Alphabet Config File: not able to parse: " + line);
				} else
					AAmasses.insert(std::make_pair(a,mass));
			}
		}
		// compute AA mass differences;
		std::set<size_t> tmp;
		for ( auto a : AAmasses ) {
			if (a.second == 0) continue;
			for ( auto b : AAmasses ) {
				if (b.second == 0) continue;
				if (a.second > b.second)
					tmp.insert(a.second - b.second);
				else if (a.second < b.second)
					tmp.insert(b.second - a.second);
			}
		}
		for ( auto a : tmp )
			AAmassDifferences.push_back(a);
		std::sort(AAmassDifferences.begin(), AAmassDifferences.end());

		if (!modificationsFile.good()) {
			std::cout << "ERROR: " << modificationsFilename << " not found. No modifications used." << std::endl;
		} else {
			std::string line;
			char a;
			size_t delta;
			char site;
			while (std::getline(modificationsFile, line)) {
				std::istringstream iss(line);
				if (line.front() == '#') continue;
				if (	!(iss >> a >> delta >> site) || 
						AAmasses.find(a) == AAmasses.end() || 
						( site != 'N' && site != 'C' && site != '*' )
					) {
					LOG("Load Modifications Config File: not able to parse: " + line);
				} else {
					if (site == '*') {
						if (allSitesMods.find(a) == allSitesMods.end())
							allSitesMods.insert(std::make_pair(a,std::vector<int>()));
						allSitesMods.find(a)->second.push_back(-delta);
					} else if (site == 'N') {
						if (nTermMods.find(a) == nTermMods.end())
							nTermMods.insert(std::make_pair(a,std::vector<int>()));
						nTermMods.find(a)->second.push_back(-delta);
					} else if (site == 'C') {
						if (cTermMods.find(a) == cTermMods.end())
							cTermMods.insert(std::make_pair(a,std::vector<int>()));
						cTermMods.find(a)->second.push_back(-delta);
					} else
						LOG("Load Modifications Config File: not able to parse: " + line);
				}
			}
		}

#ifdef DEBUG
		LOG(std::to_string(AAmasses.size()) + " Amino Acids:");
		for ( auto a : AAmasses )
			LOG(std::string(1,a.first) + " with mass " + std::to_string(a.second));

		LOG(std::to_string(allSitesMods.size()) + " amino acid(s) with variable modifications (all sites):");
		for ( auto a : allSitesMods )
			for ( auto d : a.second )
				LOG("\t" + std::string(1,a.first) + " with mass difference " + std::to_string(d));
		LOG(std::to_string(nTermMods.size()) + " amino acid(s) with variable modifications (n-terminal):");
		for ( auto a : nTermMods )
			for ( auto d : a.second )
				LOG("\t" + std::string(1,a.first) + " with mass difference " + std::to_string(d));
		LOG(std::to_string(cTermMods.size()) + " amino acid(s) variable modifications (c-terminal):");
		for ( auto a : cTermMods )
			for ( auto d : a.second )
				LOG("\t" + std::string(1,a.first) + " with mass difference " + std::to_string(d));
#endif
	}
}
