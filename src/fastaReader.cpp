/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <vector>
#include <string>
#include <array>
#include <fstream>
#include <iostream>

#include "config.h"
#include "helpers.h"
#include "fastaReader.h"

//#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif


FastaReader::FastaReader(std::string f) : filename(f) {};
FastaReader::~FastaReader(){};
void split(const std::string& name, const std::string& seq, const size_t maxLength, std::vector<std::pair<std::array<std::string,2>,size_t>>& res) {
	for (auto s : seq) {
		if (cfg::AAmasses.find(s) == cfg::AAmasses.end()) {
			LOG("ERROR: Mass of amino acid " + std::string(1,s) +  " unknown. Ignoring protein " + name + ".");
			return;
		}
	}
	if (seq.size() < maxLength) {
		for (auto it = seq.begin(); it < seq.end(); it++) {
			std::array<std::string,2> id = { name, std::string(it,seq.end()) };
			res.push_back( std::make_pair( id, std::distance(seq.begin(), seq.end()) ) );
		}
		return;
	} else {
		auto start = seq.begin();
		auto end = seq.begin() + maxLength;
		while (end != seq.end()) {
			std::array<std::string,2> id = { name, std::string(start,end) };
			res.push_back( std::make_pair(id, std::distance(seq.begin(), end)) );
			start++;
			end++;
		}
		end--;
		while (start != end) {
			std::array<std::string,2> id = { name, std::string(start,end) };
			res.push_back( std::make_pair(id, std::distance(seq.begin(), end)) );
			start++;
		}
	}
}

size_t FastaReader::getPeptides(const size_t maxLength, std::vector<std::pair<std::array<std::string,2>,size_t>>& res) const {
	// res entry meaning: <protein ID, substring>, pos of last char
	size_t dbsize = 0;
	std::ifstream infile(filename);
	if (!infile.good()) {
		std::cout << "ERROR while reading fasta file " << filename << std::endl;
		return 0;
	}
	std::string line = "";
	std::string proteinName = "";
	std::string proteinSequence = "";
	while (std::getline(infile, line)) {
		if (line.size() == 0)	break;
		if (line[0] == '>') {
			// process last protein sequence
			if (proteinName.size() > 0)
				split(proteinName, proteinSequence, maxLength, res);
			// read next protein identifier
			proteinName = line.substr(1);
			dbsize += proteinSequence.size();
			proteinSequence.clear();
		} else
			proteinSequence += line;
	}
	if (proteinName.size() > 0) {
		split(proteinName, proteinSequence, maxLength, res);
		dbsize += proteinSequence.size();
	}
	return dbsize;
};
