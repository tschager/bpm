/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#ifndef CONFIG_H
#define CONFIG_H

#include <string>
#include <vector>
#include <unordered_map>
#include <map>

namespace cfg {
	extern std::unordered_map<char, size_t> AAmasses;
	extern std::map<char, std::vector<int> > allSitesMods;
	extern std::map<char, std::vector<int> > nTermMods;
	extern std::map<char, std::vector<int> > cTermMods;
	
	extern std::vector<size_t> AAmassDifferences;

	void loadConfig(std::string alphabetFilename, std::string modificationsFilename);
}

#endif
