/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#ifndef FASTAREADER_H
#define FASTAREADER_H

#include <string>
#include <array>

class FastaReader { 
	private:
		const std::string filename;
		
	public:
		FastaReader(std::string f);
		~FastaReader();
		// push substrings of f to res (<protein ID, substring>, position of last char) and return size of db
		size_t getPeptides(const size_t maxLength, std::vector<std::pair<std::array<std::string,2>,size_t>>& res) const; 
};

#endif
