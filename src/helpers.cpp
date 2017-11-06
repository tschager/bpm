/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <unordered_map>
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>
#include <assert.h>
#include <algorithm>

#include "config.h"
#include "minmaxpst.h"
#include "helpers.h"

//#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif


size_t getMass( const char s ) { return cfg::AAmasses.find(s)->second; }
size_t getMass( const std::string& seq) {
	int mass = 0;
	for ( auto s : seq ) {
		mass += cfg::AAmasses.find(s)->second;
	}
	return mass;
}

bool exploreDown( std::vector< MinMaxPST_Node >& curPath, const MinMaxPST& pst ) {
	MinMaxPST_Node nextPoint = pst.leftmostSE( curPath.back() );
	if ( 
		!(nextPoint.first == pst.getPosINF() && nextPoint.second == pst.getNegINF()) &&
		!(nextPoint.first == pst.getPosINF() && nextPoint.second == pst.getPosINF()) && 
		!(nextPoint.second > curPath.back().second)
	   ) {
		curPath.push_back(nextPoint);
		return true;
	} else
		return false;
}

bool exploreRight( std::vector< MinMaxPST_Node >& curPath, const MinMaxPST& pst ) {
	MinMaxPST_Node queryPoint = curPath.back();
	curPath.pop_back();
	queryPoint.second++;
	MinMaxPST_Node nextPoint = pst.leftmostNE( queryPoint );
	if (
			!(nextPoint.first == pst.getPosINF() && nextPoint.second == pst.getNegINF()) &&
			!(nextPoint.first == pst.getPosINF() && nextPoint.second == pst.getPosINF()) && 
			!(nextPoint.second >= curPath.back().second)
	   ) {
		curPath.push_back(nextPoint);
		return true;
	} else
		return false;
}

std::vector<size_t> findBP( std::vector< size_t >& bp, const std::map< size_t, MinMaxPST >& psts ) {
	// curPath stores the currently explored path. The mass of a vertex curPath[i] is masses[i];
	// we first explore the path downwards at curPath.back(); if no successor, explore to the right (siblings)
	std::vector<size_t> res;
	if (bp.size() == 0) return res;

	std::vector<size_t> masses;
	masses.reserve(bp.size());
	masses.push_back(bp.at(0));
	for (size_t i = 1; i < bp.size(); i++)
		masses.push_back(masses.back()+bp.at(i));

	std::vector< MinMaxPST_Node > curPath;
	curPath.push_back(std::make_pair(0,psts.begin()->second.getPosINF()));
	MinMaxPST_Node queryPoint, nextPoint;

	bool state = true; // state is true if we have to exploreDown from curPath.back() and otherwise false
	do {
		queryPoint = curPath.back();
		if (curPath.size() == masses.size()+1) {
			res.push_back(curPath.back().first);
			//curPath.pop_back();
			state = exploreRight( curPath, psts.find( masses[curPath.size()-2] )->second );

		} else {
			if ( psts.find( masses[curPath.size()-1] ) == psts.end() ) {
				curPath.pop_back();
				state = false;
			} else {
				if ( state ) {
					state = exploreDown( curPath, psts.find( masses[curPath.size()-1] )->second );
				}
				if ( !state ) {
					if (curPath.size() == 1) // no exploreRight at root
						curPath.pop_back();
					else 
						state = exploreRight( curPath, psts.find( masses[curPath.size()-2] )->second );
				}
			}
		}
	} while (curPath.size() > 0);
	return res;
}

std::vector<size_t> readIntArray(const std::string& s) {
	std::string::size_type pos = s.find(':');
	if (pos == std::string::npos)
		std::cout << "ERROR: parsing size_t array: " << std::endl << "> " << s << std::endl;
	std::stringstream sstream( s.substr(0,pos) );
	size_t resSize;
	sstream >> resSize;
	std::vector<size_t> res(resSize);
	for (size_t i = 0; i < resSize; i++) {
		std::string::size_type posnew = s.find(',', pos+1);
		std::stringstream ss( s.substr(pos+1,posnew-pos-1) );
		ss >> res.at(i);
		pos = posnew;
	}
	return res;
}

void readBP( std::string line, std::vector<size_t>& bp ) {
	std::string::size_type pos;
	size_t mass;
	bp.clear();
			
	if (line.at(0) != '[' && line.back() != ']') {
		std::cout << "ERROR parsing block pattern " << line << std::endl;
		return;
	}

	line.erase(0,1);
	line.pop_back();

	LOG("reading: " + line);
	while ((pos = line.find(", ")) != std::string::npos) {
		std::stringstream sstream( line.substr(0, pos) );
		sstream >> mass;
		bp.push_back(mass);
		LOG("block mass: " + std::to_string(mass));
		line.erase(0, pos + 2);
	}
	std::stringstream sstream( line );
	sstream >> mass;
	bp.push_back(mass);
	LOG("block mass: " + std::to_string(mass));
}

void readBPFile( std::string file, std::vector<std::vector<size_t> >& bps ) {
	std::ifstream db(file);

	if (!db.good()) {
		std::cout << "ERROR: " << file << " not found. Abort." << std::endl;
		return;
	}

	std::string line;
	std::vector<size_t> bp;
	while (std::getline(db, line)) {
		if (line.at(0) == '#')
			continue;
		readBP( line, bp );
		bps.push_back(bp);
	}
}

void readDBFile( std::string file, std::map<size_t, MinMaxPST>& psts, MinMaxPST*& leaves, std::vector< std::pair<size_t,size_t> >& trie, std::map<size_t,std::string>& leafSeqs ) {
	std::ifstream db(file);

	if (!db.good()) {
		std::cout << "ERROR: " << file << " not found. Abort." << std::endl;
		return;
	}

	std::string line;
	std::string::size_type pos = 0;
	std::vector<std::string> v;
	size_t mass;
	while (std::getline(db, line)) {
		pos = line.find(':');
		if (pos == std::string::npos) {
			std::cout << "ERROR: parsing line: " << std::endl << "> " << line << std::endl;
			continue;
		}
		if (line.substr(0,pos) == "LEAVES") // read of mass psts done
			break;
		std::stringstream sstream( line.substr(0,pos) );
		sstream >> mass;

		psts.insert(
				std::make_pair(mass, MinMaxPST( line.substr(pos+1,line.size()))
			));
	}

	leaves = new MinMaxPST(line.substr(pos+1,line.size()));

	std::getline(db, line);
	assert( line.substr(0,4) == "TRIE" );
	size_t size;
	std::stringstream sstream( line.substr(5,line.size()) );
	sstream >> size;
	trie.resize( size );
	size_t pre, post, parent_pre;
	std::string::size_type pos2,pos3,pos4;
	while(std::getline(db, line)) {
		if (line.substr(0,13) == "LASTOCC-ORDER")
			continue;
		if (line.size() == 0)
			break;
		pos = line.find(',');
		std::stringstream s1(line.substr(0,pos));
		s1 >> pre;
		
		pos2 = line.find(',',pos+1);
		std::stringstream s2(line.substr(pos+1,pos2-pos-1));
		s2 >> post;
		
		pos3 = line.find(',',pos2+1);
		std::stringstream s3(line.substr(pos2+1,pos3-pos2-1));
		s3 >> parent_pre;

		trie.at(pre) = std::make_pair(post,parent_pre);

		pos4 = line.find(',',pos3+1);
		if (pos4-pos3 > 1)
			leafSeqs.insert( std::make_pair(pre, line.substr(pos3+1,pos4-pos3-1)) );
	};
}


std::string getProteins(
		const MinMaxPST_Node& node,
		const MinMaxPST& leaves,
		const std::vector< std::pair<size_t,size_t> >& trie,
		const std::map<size_t,std::string>& leafSeqs ) {

	std::vector< MinMaxPST_Node > curPath = { node };
	
	bool status = exploreDown( curPath, leaves );
	assert( status );
	size_t pre = curPath.back().first;
	assert( leafSeqs.find(pre) != leafSeqs.end() );
	std::string seq = leafSeqs.find(pre)->second;
	LOG("FindProtein: start seq " + seq + " preorder(" + std::to_string(pre) + ")");
	while (pre != node.first) {
		assert( pre > node.first );
		pre = trie.at(pre).second;
		seq.pop_back();
		LOG("next " +  seq + " preorder(" + std::to_string(pre) + ")");
	}
	assert( pre == node.first );
	
	return seq;
}

/* 
 * MODIFICATION TOLERANT BLOCKED PATTERN MATCHING 
 */
#ifdef MOD_TOLERANT
// check if the last curPath-vertex is valid for curMod.back()
bool isValidVertex( 
		std::vector< MinMaxPST_Node >& curPath, 
		const std::array<std::vector<size_t>,3>& curMod,
		const std::vector< std::pair<size_t,size_t> >& trie,
		const std::map<size_t,std::map<char,size_t>>& lastOccAll ) {

	const size_t subtreeRoot = (curPath.size() > 1) ? curPath.at(curPath.size()-2).first : 0;
	const size_t lastNode = curPath.back().first;
	assert( lastOccAll.find(lastNode) != lastOccAll.end() );
	const std::map<char,size_t>& lastOcc = lastOccAll.find(lastNode)->second;

	auto mod = cfg::allSitesMods.begin();
	for ( auto a : curMod.at(0) ) {
			if ( a > 0 && 
				 (lastOcc.find( mod->first ) == lastOcc.end() ||
				  lastOcc.find( mod->first )->second <= subtreeRoot) )
			return false;
			mod++;
		}
		mod = cfg::nTermMods.begin();
		for ( auto a : curMod.at(1) ) {
			if ( a > 0 && 
				 (lastOcc.find( mod->first ) == lastOcc.end() ||
				  lastOcc.find( mod->first )->second <= subtreeRoot) ) {
				return false;
			} else if (a > 0 && lastOcc.find( mod->first ) != lastOcc.end() && lastOcc.find( mod->first )->second > subtreeRoot) {
				// check if lastOcc preorder p is really the first character
				// (i.e. parent is root)
				const size_t p = lastOcc.find(mod->first)->second;
				if (trie.at(p).second != 0)
					return false;
			}
			mod++;
		}
		mod = cfg::cTermMods.begin();
		for ( auto a : curMod.at(2) ) {
			if ( a > 0 && 
				 (lastOcc.find( mod->first ) == lastOcc.end() ||
				  lastOcc.find( mod->first )->second <= subtreeRoot) )
			return false;
			mod++;
		}

	return true;
}

bool exploreRightMod( 
	std::vector< MinMaxPST_Node >& curPath, 
		const std::array<std::vector<size_t>,3>& curMod, 
		const std::map<size_t,std::map<char,size_t>>& lastOccAll, 
		const std::vector< std::pair<size_t,size_t> >& trie,
		const MinMaxPST& pst ) {
	bool res = exploreRight(curPath,pst);
	if (res) {
		// check if new vertex is valid for current modification
		while( res && !isValidVertex(curPath, curMod, trie, lastOccAll) ) {
			res = exploreRightMod(curPath,curMod,lastOccAll,trie,pst);
		}
	}

	return res;
}

bool exploreDownMod(
		std::vector< MinMaxPST_Node >& curPath, 
		const std::array<std::vector<size_t>,3>& curMod, 
		const std::map<size_t,std::map<char,size_t>>& lastOccAll, 
		const std::vector< std::pair<size_t,size_t> >& trie,
		const MinMaxPST& pst ) {

	bool res = exploreDown(curPath,pst);
	if (res) {
		// check if new vertex is valid for current modification
		while( res && !isValidVertex(curPath, curMod, trie, lastOccAll) ) {
			res = exploreRightMod(curPath,curMod,lastOccAll,trie,pst);
		}
	}
	return res;
}


// return the modified mass for a set of modifications and a mass m 
// accumulate all modifications until modifications.at(until)
size_t getModMass( const std::vector<std::array<std::vector<size_t>,3> >& modifications, size_t until, size_t mass ) {

	for (size_t i = 0; i <= until; i++) {
		std::map< char, std::vector<int> >::iterator modA2 = cfg::allSitesMods.begin();
		//allSites modifications
		for (auto mod : modifications[i][0]) {
			if (mod == 0)
				modA2++;
			else if ( mod == std::numeric_limits<size_t>::max() )
				return std::numeric_limits<size_t>::max();
			else
				mass += (modA2++)->second.at(mod-1);
		}
		//n-term modification
		modA2 = cfg::nTermMods.begin();
		for (auto mod : modifications[i][1]) {
			if (mod == 0)
				modA2++;
			else if ( mod == std::numeric_limits<size_t>::max() )
				return std::numeric_limits<size_t>::max();
			else
				mass += (modA2++)->second.at(mod-1);
		}
		//c-term modification
		modA2 = cfg::cTermMods.begin();
		for (auto mod : modifications[i][2]) {
			if (mod == 0)
				modA2++;
			else if ( mod == std::numeric_limits<size_t>::max() )
				return std::numeric_limits<size_t>::max();
			else
				mass += (modA2++)->second.at(mod-1);
		}
	}
	return mass;
}
// returns the next Mass for the last vertex on curPath
size_t nextMass( std::vector< std::array<std::vector<size_t>,3> >& curMod, const std::vector< size_t >& masses ) {
	
	std::array<std::vector<size_t>,3> next =  curMod.back();
	bool incr = false;
	std::vector<size_t>::reverse_iterator it;
	std::map< char, std::vector<int> >::reverse_iterator modA;

	if (curMod.size() == 1 && next[1].size() > 0) {
		// next n-terminal modification
		it = next[1].rbegin();
		modA = cfg::nTermMods.rbegin();
		while(!incr && it != next[1].rend()) {
			if (*it == std::numeric_limits<size_t>::max())
				break;
			if (*it == modA->second.size()) {
				*it = 0;
				it++;
				modA++;
			} else {
				*it += 1;
				incr = true;
			}
		}
	} else if (curMod.size() == masses.size() && next[2].size() > 0) {
		// next c-terminal modification
		it = next[2].rbegin();
		modA = cfg::cTermMods.rbegin();
		while(!incr && it != next[2].rend()) {
			if (*it == std::numeric_limits<size_t>::max())
				break;
			if (*it == modA->second.size()) {
				*it = 0;
				it++;
				modA++;
			} else {
				*it += 1;
				incr = true;
			}
		}
	}
	if (!incr && next[0].size() > 0) {
		// next allSites modification
		it = next[0].rbegin();
		modA = cfg::allSitesMods.rbegin();
		while(!incr && it != next[0].rend()) {
			if (*it == std::numeric_limits<size_t>::max())
				break;
			if (*it == modA->second.size()) {
				*it = 0;
				it++;
				modA++;
			} else {
				*it += 1;
				incr = true;
			}
		}
	}
	if (!incr)
		return std::numeric_limits<size_t>::max();

	curMod.back() = next;

	return getModMass(curMod,curMod.size()-1,masses.at( curMod.size() - 1 ));
}


// return next mod. mass with a non-empty priority search tree
size_t nextMassWithPST(std::vector< std::array<std::vector<size_t>,3> >& curMod, const std::vector< size_t     >& masses, const std::map< size_t, MinMaxPST >& psts ) {

	size_t tmpMass = nextMass( curMod, masses );
	while ( psts.find(tmpMass) == psts.end() && tmpMass < std::numeric_limits<size_t>::max() )
		tmpMass = nextMass( curMod, masses );
	if (tmpMass == std::numeric_limits<size_t>::max()) {
		std::fill( curMod.back().at(0).begin(), curMod.back().at(0).end(), std::numeric_limits<size_t>::max() );
		std::fill( curMod.back().at(1).begin(), curMod.back().at(1).end(), std::numeric_limits<size_t>::max() );
		std::fill( curMod.back().at(2).begin(), curMod.back().at(2).end(), std::numeric_limits<size_t>::max() );
	}
		
	return tmpMass;
}


std::vector<size_t> findBPMod( std::vector< size_t >& bp,
							   const std::map< size_t, MinMaxPST >& psts,
							   const std::vector< std::pair<size_t,size_t> >& trie,
							   const std::map<size_t,std::map<char,size_t>>& lastOcc ) {
	// curPath stores the currently explored path. The mass of a vertex curPath[i] is masses[i];
	// we first explore the path downwards at curPath.back(); if no successor, explore to the right (siblings)
	std::vector< MinMaxPST_Node > curPath; 
	curPath.push_back(std::make_pair(0,psts.begin()->second.getPosINF()));

	// curMod stores the current modification for each mass of curPath
	std::vector< std::array<std::vector<size_t>,3> > curMod;
	const std::array< std::vector<size_t>, 3> firstMod = {
		std::vector<size_t>(cfg::allSitesMods.size(),0),
		std::vector<size_t>(cfg::nTermMods.size(),0),
		std::vector<size_t>(cfg::cTermMods.size(),0)
	};
	curMod.push_back( firstMod );

	std::vector<size_t> res;
	if (bp.size() == 0) return res;

	std::vector<size_t> masses;
	masses.reserve(bp.size());
	masses.push_back(bp.at(0));
	for (size_t i = 1; i < bp.size(); i++)
		masses.push_back(masses.back()+bp.at(i));

	MinMaxPST_Node queryPoint, nextPoint;
	bool state = true; // state is true if we have to exploreDown from curPath.back() and otherwise false
	size_t tmpMass;
	do {
		assert( curPath.size() == curMod.size() );
		queryPoint = curPath.back();
		if (curPath.size() == masses.size()+1) {
			res.push_back(curPath.back().first);
			curPath.pop_back();
			curMod.pop_back();
			state = false;
		} else {
			tmpMass = getModMass( curMod, curMod.size()-1, masses[curPath.size()-1] );
			if (tmpMass < std::numeric_limits<size_t>::max() && psts.find( tmpMass ) == psts.end())
				tmpMass = nextMassWithPST( curMod, masses, psts );

			if ( tmpMass ==  std::numeric_limits<size_t>::max() ) {
				if (curPath.size() == 1) break;
				curMod.pop_back();
				tmpMass = getModMass( curMod, curMod.size()-1, masses[curPath.size()-2] );
				assert( psts.find( tmpMass ) != psts.end() );
				state = exploreRightMod( curPath, curMod.back(), lastOcc, trie, psts.find(tmpMass)->second);
				if (state)
					curMod.push_back(firstMod);
				else {
					nextMassWithPST( curMod, masses, psts );
				}
				continue;
			}
			
			
			if (state) {
				state = exploreDownMod( curPath, curMod.back(), lastOcc, trie, psts.find( tmpMass )->second );
				if (state) {
					curMod.push_back(firstMod);
				}
			}
			if (!state) {
				if (curPath.size() == 1) {
					if ( nextMassWithPST( curMod, masses, psts ) == std::numeric_limits<size_t>::max() ) {
						// back at root, no explore right
						curMod.pop_back();
						curPath.pop_back();
					} else 
						state = true;
				} else {
					if ( nextMassWithPST( curMod, masses, psts ) == std::numeric_limits<size_t>::max() ) {
						tmpMass = getModMass( curMod, curMod.size()-2, masses[curPath.size()-2] );
						state = exploreRightMod( curPath, curMod.at(curMod.size()-2), lastOcc, trie, psts.find(tmpMass)->second );
						if (state)
							curMod.back() = firstMod;
						else
							curMod.pop_back();
					} else
						state = true;
				}
			}
		}
	} while (curPath.size() > 0);
	return res;
}

std::vector<char> readCharArray(const std::string& s) {
	std::string::size_type pos = s.find(':');
	if (pos == std::string::npos)
		std::cout << "ERROR: parsing char array: " << std::endl << "> " << s << std::endl;
	std::stringstream sstream( s.substr(0,pos) );
	pos++;
	size_t resSize;
	sstream >> resSize;
	std::vector<char> res(resSize);
	for (size_t i = 0; i < resSize; i++) {
		res.at(i) = *(s.begin() + pos);
		LOG("reading " + s.substr(pos,1));
		pos += 2;
	}
	return res;
}

std::map<char,size_t> readMap( std::string s, std::vector<char> order ) {
	std::map<char,size_t> m;
	std::string::size_type p1,p2;
	p1 = 0;
	for (size_t i = 0; i < order.size(); i++) {
		p2 = s.find(',',p1+1);
		std::stringstream sstream( s.substr(p1+1, p2-p1-1) );
		size_t p;
		sstream >> p;
		m.insert( std::make_pair( order.at(i), p ) );
		p1 = p2;
	}
	#ifdef DEBUG
	for (auto a : m)
		LOG( a.first + " : " + std::to_string(a.second) );
	#endif
	return m;
}

void readDBFileMod( std::string file,
					std::map<size_t, MinMaxPST>& psts,
					MinMaxPST*& leaves,
					std::vector< std::pair<size_t,size_t> >& trie,
					std::map<size_t,std::string>& leafSeqs,
					std::map< size_t, std::map<char,size_t> >& lastOcc ) {
	readDBFile(file, psts, leaves, trie, leafSeqs);
	
	std::ifstream db(file);
	if (!db.good()) return;
	std::string line;
	std::vector<char> lastOccOrder;
	while(std::getline(db, line)) {
		// read lastOcc order
		if( line.substr(0,13) != "LASTOCC-ORDER" )
			continue;
		lastOccOrder = readCharArray( line.substr(14,line.size()) );
		break;
	}
	std::string::size_type pos;
	size_t pre;
	while (std::getline(db, line)) {
		if( line.size() == 0 )
			break;
		pos = line.find(',');
		assert( pos != std::string::npos );
		std::stringstream sstream( line.substr(0,pos) );
		sstream >> pre;
		pos = line.find(',',pos+1); // skip postorder
		pos = line.find(',',pos+1); // skip parent_preorder
		pos = line.find(',',pos+1); // skip sequence
		std::map<char,size_t> m = readMap(line.substr(pos,line.size()), lastOccOrder);
		lastOcc.insert( std::make_pair( pre, m ) );
	}
}
#endif

/* 
 * MUTATION TOLERANT BLOCKED PATTERN MATCHING
 */
#ifdef MUT_TOLERANT
// check if mass can be explained by one mutation at seq
bool isPossibleModification( std::string seq, size_t mass ) {
	size_t diff;
	const size_t seqMass = getMass(seq);
	if (seqMass < mass) {
		// insertion possible
		diff = mass - seqMass;
		for (auto a : cfg::AAmasses) {
			if ( diff == a.second )
				return true;
		}
	} else {
		diff = seqMass - mass;
		for ( auto s : seq ) {
			// deletion possible
			if ( diff == cfg::AAmasses.find(s)->second )
				return true;
		}
	}
	// check modification possible in general (diff is equal to mass difference of pair of amino acids
	auto it = std::lower_bound(cfg::AAmassDifferences.begin(),cfg::AAmassDifferences.end(),diff);
	if (it == cfg::AAmassDifferences.end() || *it != diff) // modification not possible
		return false;

	// check modification for each character
	for ( auto s : seq ) {
		const size_t m = getMass(s);
		for (auto a : cfg::AAmasses) {
			if ( seqMass - m + a.second == mass )
				return true;
		}
	}

	return false;
}

// try to combine prefix match and suffix match (preorder number given) for a block mass,
// i.e. check if suffix has a link in prefix subtree and if seq can be explained by 
// one modification, and add to results if so;
// return true if combination is possible and false otherwise.
bool combine(
		size_t prefix,
		size_t suffix,
		size_t mass,
		const std::unordered_map< size_t, std::vector<size_t> >& links,
		const std::vector< std::pair<size_t,size_t> >& trie,
		const MinMaxPST& leaves,
		const std::map< size_t, std::string >& leafSeqs,
		std::vector<std::string>& results ) {

	if (links.find(suffix) == links.end())
		return false;
	bool found = false;
	const MinMaxPST_Node suf = {suffix, trie.at(suffix).first};
	const std::string sufseq = getProteins( suf, leaves, trie, leafSeqs);
	for ( auto a : links.find(suffix)->second ) {
		if ( a > prefix && trie.at(a).first < trie.at(prefix).first ) {
			// a is in subtree of prefix

			// extract sequence between prefix and a
			std::vector< MinMaxPST_Node > curPath = { std::make_pair(a,trie.at(a).first) };
			bool status = exploreDown( curPath, leaves );
			assert( status );
			std::string seq = leafSeqs.at(curPath.back().first);
			size_t cur = curPath.back().first;
			while( cur > a ) {
				seq.pop_back();
				cur = trie.at(cur).second;
			}
			assert( a == cur );
			size_t count = 0;
			while (cur > prefix) {
				count++;
				cur = trie.at(cur).second;
			}
			seq = seq.substr(seq.size()-count,seq.size());

			// check if mass can be explained by seq (w/ one mutation)
			if( isPossibleModification(seq, mass) ) {
				found = true;
				MinMaxPST_Node pre = {a, trie.at(a).first};
				results.push_back( getProteins( pre, leaves, trie, leafSeqs ) + sufseq );
			}
		}
	}
	return found;
}

std::vector<size_t> getLeavesInSubtree( const MinMaxPST_Node& n, const MinMaxPST& leaves ) {
	std::vector<size_t> res;
	std::vector< MinMaxPST_Node > curPath = { n };
	bool status = exploreDown( curPath, leaves );
	while( status ) {
		res.push_back( curPath.back().first );
		status = exploreRight( curPath, leaves );
	}
	return res;
}


bool combineWithLeaf(
		size_t prefix,
		size_t mass,
		const std::unordered_map< size_t, std::vector<size_t> >& links,
		const std::vector< std::pair<size_t,size_t> >& trie,
		const MinMaxPST& leaves,
		const std::map< size_t, std::string >& leafSeqs,
		std::vector<std::string>& results ) {

	bool found = false;
	for ( auto a : getLeavesInSubtree({prefix,trie.at(prefix).first}, leaves) ) {
		// extract sequence between prefix and a
		std::vector< MinMaxPST_Node > curPath = { std::make_pair(a,trie.at(a).first) };
		bool status = exploreDown( curPath, leaves );
		assert( status );
		std::string seq = leafSeqs.at(curPath.back().first);
		size_t cur = curPath.back().first;
		size_t count = 0;
		while (cur > prefix) {
			count++;
			cur = trie.at(cur).second;
		}

		for (size_t i = 1; i < seq.size()-count; i++) {
			// check if mass can be explained by seq (w/ one mutation)
			if( isPossibleModification(seq.substr(seq.size()-count,i), mass) ) {
				found = true;
				results.push_back( seq.substr(0,seq.size()-count+i) );
			}
		}
	}
	return found;
}

std::vector<std::string> findBPMut( std::vector< size_t >& masses,
							   const std::map< size_t,MinMaxPST >& psts,
							   const std::unordered_map< size_t, std::vector<size_t> >& links,
							   const std::vector< std::pair<size_t,size_t> >& trie,
							   const MinMaxPST& leaves,
							   const std::map< size_t, std::string >& leafSeqs
							   ) {
	std::vector<std::string> res;
	if (masses.size() == 0) return res;

	for ( size_t i = 0; i < masses.size(); i++ ) {
		// mutation in masses.at(i)
		std::vector<size_t> prefixMasses;
		std::vector<size_t> suffixMasses;
		if (i > 0)
			for ( size_t j = 0; j < i; j++ )
				prefixMasses.push_back(masses.at(j));
		if (i+1 < masses.size())
			for ( size_t j = i+1; j < masses.size(); j++ )
				suffixMasses.push_back(masses.at(j));

		if (i > 0 && i < masses.size() - 1) {
			// mutation neither in first nor in last block
			std::vector<size_t> prefixMatch = findBP(prefixMasses, psts);
			std::vector<size_t> suffixMatch = findBP(suffixMasses, psts);

			for (auto p : prefixMatch) {
				for (auto s : suffixMatch)
					combine(p,s,masses.at(i),links,trie,leaves,leafSeqs,res);
			}
		} else if ( i == 0 ) {
			// mutation in first block
			std::vector<size_t> suffixMatch = findBP(suffixMasses, psts);
			for (auto s : suffixMatch) 
				combine(0,s,masses.at(0),links,trie,leaves,leafSeqs,res);
		} else {
			// mutation in last block
			std::vector<size_t> prefixMatch = findBP(prefixMasses, psts);
			for (auto p : prefixMatch)
				combineWithLeaf(p,masses.back(),links,trie,leaves,leafSeqs,res);
		}
	}

	return res;
}

void readDBFileMut( std::string file,
					std::map<size_t, MinMaxPST>& psts,
					MinMaxPST*& leaves,
					std::vector< std::pair<size_t,size_t> >& trie,
					std::map<size_t,std::string>& leafSeqs,
					std::unordered_map< size_t, std::vector<size_t> >& links ) {
	readDBFile(file, psts, leaves, trie, leafSeqs);
	
	std::ifstream db(file);
	if (!db.good()) return;
	std::string line;
	std::vector<char> lastOccOrder;
	while(std::getline(db, line)) {
		// read lastOcc order
		if( line.substr(0,13) != "LASTOCC-ORDER" )
			continue;
		break;
	}
	std::string::size_type p1, p2;
	size_t pre, len;
	while (std::getline(db, line)) {
		if( line.size() == 0 )
			break;

		p1 = line.find(',');
		assert( p1 != std::string::npos );
		std::stringstream s1( line.substr(0,p1) );
		s1 >> pre;
		if (pre == 0) continue;
		p1 = line.find('|');
		if( p1 == std::string::npos ) continue; // no links (root or leaf)
		p2 = line.find(':', p1+1);
		assert( p2 != std::string::npos );
		std::stringstream sstream( line.substr(p1+1,p2-p1) );
		sstream >> len;
		links.insert( std::make_pair( pre, std::vector<size_t>() ) );
		std::vector<size_t>& linkvector = links.find(pre)->second;
		linkvector.reserve(len);

		size_t l;
		while( p2 != std::string::npos ) {
			p1 = p2;
			p2 = line.find(",",p1+1);
			std::stringstream s2( line.substr(p1+1,p2-p1) );
			s2 >> l;
			linkvector.push_back(l);
		}
		assert(linkvector.size() == len);
	}
}
#endif
