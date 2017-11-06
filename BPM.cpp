#include <stdlib.h>
#include <chrono>
#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <stack>
#include <fstream>

#include "src/config.h"
#include "src/trie.h"
#include "src/helpers.h"
#include "src/minmaxpst.h"
#include "src/fastaReader.h"

#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif

int main(int argc, char *argv[]) {

	std::string bpFile;
	std::string dbFile;
	std::string modFile;
	std::string aaFile;
	
	if (argc < 2) {
		std::cout << "USAGE: " << argv[0] << " <DB file (fasta)> [<post-translational modifications file (cfg/modifications.cfg)> <AA masses file (cfg/aminoacids.cfg)>]" << std::endl;
		return 1;
	} else {
//		bpFile = argv[1];
		dbFile = argv[1];
	}
	if (argc > 2)
		modFile = argv[2];
	else
		modFile = "cfg/modifications.cfg";
	if (argc > 3)
		aaFile = argv[3];
	else
		aaFile = "cfg/aminoacids.cfg";
	LOG("Program: " + std::string(argv[0]));
	LOG("DB File: " + dbFile);
	LOG("PTM File: " + modFile);
	LOG("AAmasses File: " + aaFile);

	cfg::loadConfig(aaFile,modFile);

	// read DB File
	auto begin = std::chrono::high_resolution_clock::now();
	std::map< size_t, MinMaxPST > psts;
	MinMaxPST* leaves = nullptr;
	std::vector< std::pair<size_t,size_t> > trie;
	std::map< size_t, std::string > leafSeqs;

#ifdef MOD_TOLERANT
	std::map<size_t,std::map<char,size_t>> lastOcc;
	readDBFileMod( dbFile, psts, leaves, trie, leafSeqs, lastOcc );
#else
#ifdef MUT_TOLERANT
	std::unordered_map< size_t, std::vector<size_t> > links;
	readDBFileMut( dbFile, psts, leaves, trie, leafSeqs, links );
	if (links.size() == 0)
		std::cout << "Warnings: no links in DB index - probably index file has not been generated for mutation-tolerant BPM" << std::endl;
#else
	readDBFile( dbFile, psts, leaves, trie, leafSeqs );
#endif
#endif
	auto end = std::chrono::high_resolution_clock::now();
	auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	std::cout << "read DB file (in milliseconds): " << diff << std::endl;

#ifdef MOD_TOLERANT
	const bool modifications = ( cfg::allSitesMods.size() + cfg::nTermMods.size() + cfg::cTermMods.size() ) > 0;
	if (modifications) 
		std::cout << "Modification-tolerant blocked pattern matching" << std::endl;
	else
		std::cout << "No modifications -- standard blocked pattern matching" << std::endl;
#endif

#ifdef MUT_TOLERANT
	std::cout << "Mutation-tolerant blocked pattern matching" << std::endl;
#endif

	begin = std::chrono::high_resolution_clock::now();
	std::string line;
	std::vector<size_t> bp;
	size_t count = 0;
	while(getline(std::cin, line)) {
		if (line.size() == 0)
			break;
		if (line.at(0) == '#') {
			std::cout << line << std::endl;
			continue;
		}
		readBP(line,bp);
//		std::cout << "BP: ";
//		for ( auto a : bp )
//			std::cout << a << " ";
//		std::cout << std::endl;
		if (bp.size() == 0)
			continue;
		count++;
#ifdef MUT_TOLERANT
//		findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
		std::vector<std::string> results = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
//		if (results.size() > 0) {
//			std::cout << "Nr of hits: " << results.size() << std::endl;
//		}
		for ( auto r : results ) {
			std::cout << r << std::endl;
		}
#else
		std::vector<size_t> results;
#ifdef MOD_TOLERANT
		if (modifications)
			results = findBPMod(bp,psts,trie,lastOcc);
		else
			results = findBP(bp,psts);
#else
		results = findBP(bp,psts);
#endif
		for ( auto r : results ) {
			MinMaxPST_Node n = {r, trie.at(r).first};
			std::cout << getProteins( n, *leaves, trie, leafSeqs) << std::endl;
//			getProteins( n, *leaves, trie, leafSeqs); // only for benchmark
		}
#endif
	}
	end = std::chrono::high_resolution_clock::now();
	diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin).count();
	std::cout << count << " block pattern matchings (in milliseconds): " << diff << std::endl;
}

