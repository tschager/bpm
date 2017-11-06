/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#ifndef HELPERS_H
#define HELPERS_H

#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include "minmaxpst.h"

size_t getMass(const std::string& seq);

std::vector<size_t> findBP( std::vector< size_t >& masses,
							const std::map< size_t,MinMaxPST >& psts);

void readBP( std::string line, std::vector<size_t>& bp );
void readBPFile( std::string file, std::vector<std::vector<size_t> >& bps );

void readDBFile( std::string file,
				 std::map<size_t, MinMaxPST>& psts,
				 MinMaxPST*& leaves,
				 std::vector< std::pair<size_t,size_t> >& trie,
				 std::map<size_t,std::string>& leafSeqs ); // read file and write psts, trie, leaves, and leafSeqs (for each leaf the corresponding sequence

std::string getProteins(
				const MinMaxPST_Node& node,
				const MinMaxPST& leaves,
				const std::vector< std::pair<size_t,size_t> >& trie,
				const std::map<size_t,std::string>& leafSeqs ); // return sequence and all proteins which contain this string

#ifdef MOD_TOLERANT
std::vector<size_t> findBPMod( std::vector< size_t >& masses,
							   const std::map< size_t, MinMaxPST >& psts,
							   const std::vector< std::pair<size_t,size_t> >& trie,
							   const std::map<size_t, std::map<char,size_t>>& lastOcc);

void readDBFileMod( std::string file,
					std::map<size_t, MinMaxPST>& psts,
					MinMaxPST*& leaves,
					std::vector< std::pair<size_t,size_t> >& trie,
					std::map< size_t, std::string >& leafSeqs,
					std::map< size_t, std::map<char,size_t> >& lastOcc );
#endif

#ifdef MUT_TOLERANT
std::vector<std::string> findBPMut( std::vector< size_t >& masses,
									const std::map< size_t, MinMaxPST >& psts,
									const std::unordered_map< size_t, std::vector<size_t> >& links,
									const std::vector< std::pair<size_t,size_t> >& trie,
									const MinMaxPST& leaves,
									const std::map< size_t, std::string >& leafSeqs
									);

void readDBFileMut( std::string file,
					std::map<size_t, MinMaxPST>& psts,
					MinMaxPST*& leaves,
					std::vector< std::pair<size_t,size_t> >& trie,
					std::map< size_t, std::string >& leafSeqs,
					std::unordered_map< size_t, std::vector<size_t> >& links );
#endif
#endif
