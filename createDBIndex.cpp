#include <stdlib.h>
#include <ctime>
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

	std::string dbFile;
	std::string modFile;
	std::string aaFile;

	if (argc < 2) {
		std::cout << "USAGE: " << argv[0] << " <DB File (fasta)> [<post-translation    al modifications file (cfg/modifications.cfg)> <AA masses file (cfg/aminoacids.cfg)    >]" << std::endl;
		return 1;
	} else
		dbFile = argv[1];
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

	// read fasta file
	FastaReader f(dbFile);
	std::vector<std::pair<std::array<std::string,2>,size_t>> res;
	size_t dbsize = f.getPeptides(30,res);
	std::cout << "DB size: " << dbsize << std::endl;
	
	Trie t;
	for (auto pep : res)
		t.add(pep.first[1],pep.first[0]);
	res.clear();
	res.shrink_to_fit();
	std::cout << "Trie construction done" << std::endl;

	std::map< size_t, std::vector< MinMaxPST_Node > > points;
	t.getNodesByMass(points);
	std::cout << "points sorted by mass" << std::endl;
	
	std::ofstream outFile;
	outFile.open(dbFile + ".db");

	for ( auto m : points )
		outFile << m.first << ":" << MinMaxPST(m.second) << std::endl;
	points.clear();

	std::cout << "PSTs written" << std::endl;

	outFile << "LEAVES:" << MinMaxPST( t.getLeaves() ) << std::endl;
	std::cout << "Leaves PST written" << std::endl;

	outFile << "TRIE:" << t << std::endl;
	std::cout << "Trie structure written" << std::endl;
	std::cout << "done" << std::endl;
}
