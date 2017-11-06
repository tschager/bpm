/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <iostream>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <stack>
#include <set>
#include <cstdlib>
#include <random>
#include <algorithm>


#include <array>
#include <assert.h>

#include "../src/config.h"
#include "../src/trie.h"
#include "../src/minmaxpst.h"
#include "../src/helpers.h"
#include "../src/fastaReader.h"

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "../lib/doctest/doctest/doctest.h"
#include <random>

#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif

TEST_CASE("Trie tests") {
	cfg::loadConfig("cfg/aminoacids.cfg","tests/modifications.cfg");
	Trie t;
	CHECK(t.size() == 1);

	SUBCASE("add and find some words") {
		// add some words
		t.add("test","pep1");
		CHECK(t.find("test"));
		CHECK(!t.find("tes"));
		t.add("testa","pep2");
		CHECK(t.find("test"));
		CHECK(t.find("testa"));
		CHECK(!t.find("testb"));
		CHECK(!t.find("te"));
		CHECK(t.size() == 6);
	}
	
	SUBCASE("add and find some words") {
		// add some words
		t.add("test","pep1");
		CHECK(t.find("test"));
		CHECK(!t.find("tes"));
		t.add("testa","pep2");
		CHECK(t.find("test"));
		CHECK(t.find("testa"));
		CHECK(!t.find("testb"));
		CHECK(!t.find("te"));
		CHECK(t.size() == 6);
	}

	SUBCASE("add and find with empty word") {
		t.add("AGT","");
		CHECK(t.find("AGT"));
		CHECK(!t.find("AG"));
		CHECK(!t.find("A"));
		t.add("ADF","");
		CHECK(t.find("AGT"));
		CHECK(t.find("ADF"));
		CHECK(!t.find("ADFA"));
		CHECK(!t.find(""));
		t.add("","");
		CHECK(t.find(""));
		t.add("AGTAG","");
		CHECK(t.find("AGT"));
		CHECK(t.find("AGTAG"));
		CHECK(!t.find("AGTA"));
	}

	SUBCASE("preorder, postorder, parent") {
		t.add("AAA","");
		t.add("AG","");
		t.add("AAG","");
		t.add("GA","");
		t.add("GAG","");
		t.add("GG","");
		CHECK(t.size() == 10);

		TrieNode* cur = t.findByPreorder(2);
		CHECK(cur->getPreorder() == 2);
		CHECK(cur->getPostorder() == 2);
		CHECK(cur->getParentPreorder() == 1);

		cur = t.findByPreorder(3);
		CHECK(cur->getPreorder() == 3);
		CHECK(cur->getPostorder() == 0);
		CHECK(cur->getParentPreorder() == 2);

		cur = t.findByPreorder(0);
		CHECK(cur->getPreorder() == 0);
		CHECK(cur->getPostorder() == 9);
		CHECK(cur->getParentPreorder() == 0);

		cur = t.findByPreorder(7);
		CHECK(cur->getPreorder() == 7);
		CHECK(cur->getPostorder() == 6);
		CHECK(cur->getParentPreorder() == 6);
		
		cur = t.findByPreorder(9);
		CHECK(cur->getPreorder() == 9);
		CHECK(cur->getPostorder() == 7);
		CHECK(cur->getParentPreorder() == 6);
	}

	SUBCASE("getNodesByMass") {
		t.add("AAA","");
		t.add("AG","");
		t.add("AAG","");
		t.add("GA","");
		t.add("GAG","");
		t.add("GG","");
		CHECK(t.size() == 10);
		CHECK(!t.find(""));
		CHECK(!t.find("A"));
		CHECK(!t.find("AA"));
		CHECK(t.find("AAA"));
		CHECK(t.find("AG"));
		CHECK(t.find("AAG"));
		CHECK(t.find("GA"));
		CHECK(t.find("GAG"));
		CHECK(t.find("GG"));
		CHECK(!t.find("GAGA"));

		std::map< size_t, std::vector< MinMaxPST_Node > > points;
		t.getNodesByMass(points);

		CHECK(points.find(5702)->second.size() == 1);
		CHECK(points.find(7104)->second.size() == 1);
		// mass(AG) = 12806
		CHECK(points.find(12806)->second.size() == 2);
		// mass (AGG) = 18508
		CHECK(points.find(18508)->second.size() == 1);
	}
#ifdef MOD_TOLERANT
	SUBCASE("lastOccurence") {

		t.add("AAA","p1");
		t.add("AAC","p2");
		t.add("ACC","p3");
		t.add("CCA","p4");
		t.add("CCC","p5");
		std::map<size_t,std::map<char,size_t>> res = t.outputLastOcc();

		CHECK( res.find(0)->second.find('A')->second == 0 );
		CHECK( res.find(1)->second.find('A')->second == 1 );
		CHECK( res.find(2)->second.find('A')->second == 2 );
		CHECK( res.find(3)->second.find('A')->second == 3 );
		CHECK( res.find(4)->second.find('A')->second == 2 );
		CHECK( res.find(5)->second.find('A')->second == 1 );
		CHECK( res.find(6)->second.find('A')->second == 1 );
		CHECK( res.find(7)->second.find('A')->second == 0 );
		CHECK( res.find(8)->second.find('A')->second == 0 );
		CHECK( res.find(9)->second.find('A')->second == 9 );
		CHECK( res.find(10)->second.find('A')->second == 0 );
		CHECK( res.find(0)->second.find('C')->second == 0 );
		CHECK( res.find(1)->second.find('C')->second == 0 );
		CHECK( res.find(2)->second.find('C')->second == 0 );
		CHECK( res.find(3)->second.find('C')->second == 0 );
		CHECK( res.find(4)->second.find('C')->second == 4 );
		CHECK( res.find(5)->second.find('C')->second == 5 );
		CHECK( res.find(6)->second.find('C')->second == 6 );
		CHECK( res.find(7)->second.find('C')->second == 7 );
		CHECK( res.find(8)->second.find('C')->second == 8 );
		CHECK( res.find(9)->second.find('C')->second == 8 );
		CHECK( res.find(10)->second.find('C')->second == 10);
	}
#endif
#ifdef MUT_TOLERANT
	SUBCASE("Links information") {

		t.add("AAA","p1");
		t.add("AAC","p2");
		t.add("ACC","p3");
		t.add("CCA","p4");
		t.add("CCC","p5");

		std::unordered_map< size_t, std::vector<size_t> > links;
		t.computeLinks(links);

		CHECK(links.find(1)->second.size() == 3);
		CHECK(links.find(1)->second.at(0) == 1);
		CHECK(links.find(1)->second.at(1) == 2);
		CHECK(links.find(1)->second.at(2) == 8);
		CHECK(links.find(2)->second.size() == 1);
		CHECK(links.find(2)->second.at(0) == 1);
		CHECK(links.find(3) == links.end());
		CHECK(links.find(4) == links.end());
		CHECK(links.find(5)->second.size() == 1);
		CHECK(links.find(5)->second.at(0) == 1);
		CHECK(links.find(6) == links.end());
		CHECK(links.find(7)->second.size() == 5);
		CHECK(links.find(7)->second.at(0) == 2);
		CHECK(links.find(7)->second.at(1) == 1);
		CHECK(links.find(7)->second.at(2) == 5);
		CHECK(links.find(7)->second.at(3) == 7);
		CHECK(links.find(7)->second.at(4) == 8);
		CHECK(links.find(8)->second.size() == 2);
		CHECK(links.find(8)->second.at(0) == 1);
		CHECK(links.find(8)->second.at(1) == 7);
		CHECK(links.find(9) == links.end());
		CHECK(links.find(10) == links.end());
	}
#endif
}

TEST_CASE("Min-Max Priority Search Tree Tests") {
	std::vector< MinMaxPST_Node > points;

	SUBCASE("small example with 4 points") {
		points.push_back( std::make_pair(4,4) );
		points.push_back( std::make_pair(2,1) );
		points.push_back( std::make_pair(1,3) );
		points.push_back( std::make_pair(3,2) );
		MinMaxPST pst(points);
		std::vector< MinMaxPST_Node > res = pst.getArray();
		CHECK( res.size() == 4 );
		CHECK( res.at(0).first == 2 );
		CHECK( res.at(0).second == 1 );

		CHECK( res.at(1).first == 1 );
		CHECK( res.at(1).second == 3 );

		CHECK( res.at(2).first == 4 );
		CHECK( res.at(2).second == 4 );

		CHECK( res.at(3).first == 3 );
		CHECK( res.at(3).second == 2 );
	}

	SUBCASE("example with 8 points") {
		points.push_back( std::make_pair(2,6) );
		points.push_back( std::make_pair(0,0) );
		points.push_back( std::make_pair(1,2) );
		points.push_back( std::make_pair(5,7) );
		points.push_back( std::make_pair(3,4) );
		points.push_back( std::make_pair(6,5) );
		points.push_back( std::make_pair(7,3) );
		points.push_back( std::make_pair(8,1) );
		MinMaxPST pst(points);
		std::vector< MinMaxPST_Node > res = pst.getArray();
		CHECK( res.size() == 8 );
		CHECK( res.at(0).first == 0 );
		CHECK( res.at(0).second == 0 );

		CHECK( res.at(1).first == 5 );
		CHECK( res.at(1).second == 7 );

		CHECK( res.at(2).first == 6 );
		CHECK( res.at(2).second == 5 );

		CHECK( res.at(3).first == 1 );
		CHECK( res.at(3).second == 2 );

		CHECK( res.at(4).first == 3 );
		CHECK( res.at(4).second == 4 );

		CHECK( res.at(5).first == 7 );
		CHECK( res.at(5).second == 3 );

		CHECK( res.at(6).first == 8 );
		CHECK( res.at(6).second == 1 );

		CHECK( res.at(7).first == 2 );
		CHECK( res.at(7).second == 6 );
	}

	SUBCASE("leftmostNE") {
		points.push_back( std::make_pair(2,6) );
		points.push_back( std::make_pair(0,0) );
		points.push_back( std::make_pair(1,2) );
		points.push_back( std::make_pair(5,7) );
		points.push_back( std::make_pair(3,4) );
		points.push_back( std::make_pair(6,5) );
		points.push_back( std::make_pair(7,3) );
		points.push_back( std::make_pair(8,1) );
		MinMaxPST pst(points);

		MinMaxPST_Node p;
		p = pst.leftmostNE( std::make_pair(1,1) );
		CHECK(p.first == 1);
		CHECK(p.second == 2);
		p = pst.leftmostNE( std::make_pair(2,2) );
		CHECK(p.first == 2);
		CHECK(p.second == 6);
		p = pst.leftmostNE( std::make_pair(2,7) ); 
		CHECK(p.first == 5);
		CHECK(p.second == 7);
	}

	SUBCASE("leftmostSE") {
		points.push_back( std::make_pair(2,6) );
		points.push_back( std::make_pair(0,0) );
		points.push_back( std::make_pair(1,2) );
		points.push_back( std::make_pair(5,7) );
		points.push_back( std::make_pair(3,4) );
		points.push_back( std::make_pair(6,5) );
		points.push_back( std::make_pair(7,3) );
		points.push_back( std::make_pair(8,1) );
		MinMaxPST pst(points);

		MinMaxPST_Node p;
		p = pst.leftmostSE( std::make_pair(1,1) );
		CHECK(p.first == 8);
		CHECK(p.second == 1);
		p = pst.leftmostSE( std::make_pair(2,2) );
		CHECK(p.first == 8);
		CHECK(p.second == 1);
		p = pst.leftmostSE( std::make_pair(2,7) ); 
		CHECK(p.first == 2);
		CHECK(p.second == 6);
	}
}


TEST_CASE("Block pattern matching tests") {
	std::ifstream infile("tests/unittest.db");
	std::string line;
	Trie t;
	while (std::getline(infile, line)) {
		if (line.size() == 0)
			break;
		t.add(line,"");
	}
	std::map< size_t, std::vector< MinMaxPST_Node > > points;
	t.getNodesByMass(points);
	std::map< size_t, MinMaxPST > psts;
	for ( auto m : points )
		psts.insert( std::make_pair( m.first, MinMaxPST(m.second) ));

	const size_t nrPatterns = 10;
	const size_t patternLength = 10;
	size_t count = 0;
	std::vector<size_t> bp;
	std::ifstream infile2("tests/unittest.db");
	
	while (std::getline(infile2, line) && count < nrPatterns) {
		count++;
		if (line.size() == 0)
			break;
		bp.clear();
		for ( size_t j = 2; j < line.size() && bp.size() < patternLength; j+=2)
			bp.push_back( getMass(line.substr (j-2,2)) );
		
		CHECK( findBP(bp,psts).size() > 0);
	}
}


#ifdef MOD_TOLERANT
TEST_CASE("modification-tolerant BPM") {
	SUBCASE("one PTM start") {
		cfg::loadConfig("cfg/aminoacids.cfg","tests/mod2.cfg");
		std::map< size_t, MinMaxPST > psts;
		MinMaxPST* leaves = nullptr;
		std::vector< std::pair<size_t,size_t> > trie;
		std::map< size_t, std::string > leafSeqs;
		std::unordered_map< size_t, std::vector<size_t> > links;
		std::map<size_t,std::map<char,size_t>> lastOcc;
		readDBFileMod("tests/unittest2.fasta.db", psts, leaves, trie, leafSeqs, lastOcc);

		std::vector<size_t> bp = {
			35720+4321,			// M*LL
			28117				// APL
		};
		std::vector<size_t> res = findBPMod(bp,psts,trie,lastOcc);
		CHECK(res.size() == 1);
		MinMaxPST_Node n = {res.front(), trie.at(res.front()).first};
		CHECK(getProteins( n, *leaves, trie, leafSeqs) == "MLLAPL");
		
		bp = {
			54328+4321,			// WM*LL
			28117				// APL
		};
		res = findBPMod(bp,psts,trie,lastOcc);
		CHECK(res.size() == 0);
	}

	SUBCASE("one PTM mid") {
		cfg::loadConfig("cfg/aminoacids.cfg","tests/mod1.cfg");
		std::map< size_t, MinMaxPST > psts;
		MinMaxPST* leaves = nullptr;
		std::vector< std::pair<size_t,size_t> > trie;
		std::map< size_t, std::string > leafSeqs;
		std::unordered_map< size_t, std::vector<size_t> > links;
		std::map<size_t,std::map<char,size_t>> lastOcc;
		readDBFileMod("tests/unittest2.fasta.db", psts, leaves, trie, leafSeqs, lastOcc);

		std::vector<size_t> bp = {
			21312,					// GR
			34218,					// RW
			54328+1600,				// WM*LL
			28117					// APL
		};
		std::vector<size_t> res = findBPMod(bp,psts,trie,lastOcc);
		CHECK(res.size() == 1);
		MinMaxPST_Node n = {res.front(), trie.at(res.front()).first};
		CHECK(getProteins( n, *leaves, trie, leafSeqs) == "GRRWWMLLAPL");
	}

	SUBCASE("one PTM end") {
		cfg::loadConfig("cfg/aminoacids.cfg","tests/mod2.cfg");
		std::map< size_t, MinMaxPST > psts;
		MinMaxPST* leaves = nullptr;
		std::vector< std::pair<size_t,size_t> > trie;
		std::map< size_t, std::string > leafSeqs;
		std::unordered_map< size_t, std::vector<size_t> > links;
		std::map<size_t,std::map<char,size_t>> lastOcc;
		readDBFileMod("tests/unittest2.fasta.db", psts, leaves, trie, leafSeqs, lastOcc);

		std::vector<size_t> bp = {
			21312,				// GR
			34218,				// RW
			31712+1234,			// WM*
		};
		std::vector<size_t> res = findBPMod(bp,psts,trie,lastOcc);
		CHECK(res.size() == 1);
		MinMaxPST_Node n = {res.front(), trie.at(res.front()).first};
		CHECK(getProteins( n, *leaves, trie, leafSeqs) == "GRRWWM");
		
		bp = {
			21312,				// GR
			34218,				// RW
			40302+1234,			// WM*L
		};
		res = findBPMod(bp,psts,trie,lastOcc);
		CHECK(res.size() == 0);
	}

}
#endif

#ifdef MUT_TOLERANT
TEST_CASE("mutation-tolerant BPM") {
	cfg::loadConfig("cfg/aminoacids.cfg","cfg/modifications.cfg");
	std::map< size_t, MinMaxPST > psts;
	MinMaxPST* leaves = nullptr;
	std::vector< std::pair<size_t,size_t> > trie;
	std::map< size_t, std::string > leafSeqs;
	std::unordered_map< size_t, std::vector<size_t> > links;
	readDBFileMut("tests/unittest2.fasta.db", psts, leaves, trie, leafSeqs, links);

	// Insertion
	std::vector<size_t> bp = {
		32321,			// PLL
		48420,			// SPGWG
		19910 + 11404,	// AGA + (MUTATION +N)
		28416			// AGR
	};
	std::vector<std::string> res = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
	CHECK(res.size() == 1); 
	CHECK(res.front() == "PLLSPGWGAGAAGR");

	// Deletion
	bp = {
		32321,			// PLL
		48420 - 18608,	// SPGWG + (MUTATION -W)
		19910,			// AGA
		28416			// AGR
	};
	res = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
	CHECK(res.size() == 1); 
	CHECK(res.front() == "PLLSPGWGAGAAGR");

	// Mutation
	bp = {
		32321,					// PLL
		48420 - 18608 + 13706,	// SPGWG + (MUTATION W->H)
		19910,					// AGA
		28416					// AGR
	};
	res = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
	CHECK(res.size() == 1); 
	CHECK(res.front() == "PLLSPGWGAGAAGR");

	// Mutation in first block
	bp = {
		32321 - 11308 + 8703,	// PLL + Mutation L->S
		48420,					// SPGWG
		19910,					// AGA
		28416					// AGR
	};
	res = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
	CHECK(res.size() == 1); 
	CHECK(res.front() == "PLLSPGWGAGAAGR");

	// Mutation in last block
	bp = {
		32321,					// PLL
		48420,					// SPGWG
		19910,					// AGA
		28416 - 7104 + 9907		// AGR + Mutaton A->V
	};
	res = findBPMut(bp,psts,links,trie,*leaves,leafSeqs);
	CHECK(res.size() == 1); 
	CHECK(res.front() == "PLLSPGWGAGAAGR");
}
#endif

void checkAllinDB( std::string file ) {
		std::ifstream infile(file);
	std::string line;
	Trie t;
	while (std::getline(infile, line)) {
		if (line.size() == 0)
			break;
		t.add(line,"");
	}
	std::map< size_t, std::vector< MinMaxPST_Node > > points;
	t.getNodesByMass(points);
	std::map< size_t, MinMaxPST > psts;
	for ( auto m : points )
		psts.insert( std::make_pair( m.first, MinMaxPST(m.second) ));

	const size_t nrPatterns = 10;
	const size_t patternLength = 10;
	size_t c = 0;
	std::vector<size_t> bp;
	std::ifstream infile2(file);
	
	while (std::getline(infile2, line) && c < nrPatterns) {
		c++;
		if (line.size() == 0)
			break;
		bp.clear();
		for ( size_t j = 2; j < line.size() && bp.size() < patternLength; j+=2) {
			bp.push_back( getMass(line.substr (j-2,2)) );
		}
		CHECK( findBP(bp,psts).size() == 1);
	}
}
TEST_CASE("Block pattern matching test 2") {
	checkAllinDB("tests/unittestdb2.txt");
	checkAllinDB("tests/unittestdb3.txt");
}

TEST_CASE("Fasta Reader") {
	FastaReader f("tests/unittest.fasta");
	std::vector<std::pair<std::array<std::string,2>,size_t>> res;
	f.getPeptides(60,res);
	CHECK( res.front().first[0] == "pep1" );
	CHECK( res.front().first[1] == "MGAPLLSPGWGAGAAGRRWWMLLAPLLPALLL" );
	CHECK( res.front().second == 32 );
	CHECK( res.at(32).first[0] == "pep2" );
	CHECK( res.at(32).first[1] == "MSSHEGGKKKALKQPKKQAKEMDEEEKAFKQKQKEEQKKLEVLKAKVVGKGPLATGGIKK" );
	CHECK( res.at(32).second == 60 );
	CHECK( res.at(33).first[0] == "pep2" );
	CHECK( res.at(33).first[1] == "SSHEGGKKKALKQPKKQAKEMDEEEKAFKQKQKEEQKKLEVLKAKVVGKGPLATGGIKKS" );
	CHECK( res.at(33).second == 61 );
	CHECK( res.back().first[0] == "pep4" );
	CHECK( res.back().first[1] == "L" );
}
