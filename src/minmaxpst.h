/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#ifndef MINMAXPST_H
#define MINMAXPST_H

#include <vector>
#include <functional>

typedef size_t coord_t;
typedef std::pair<coord_t,coord_t> MinMaxPST_Node;

class MinMaxPST {
	private:
		std::vector<MinMaxPST_Node> t;
		coord_t posINF;		// positive infinity
		coord_t negINF;		// negative infinity

		const MinMaxPST_Node& get(size_t index) const; // return entry for 1-base indexing
		size_t smallestYCoordIndex(size_t start, size_t end) const;
		size_t largestYCoordIndex(size_t start, size_t end) const;
		void swap(size_t i, size_t j);
		bool isLeaf(const size_t i) const;
		size_t leftChild(const size_t i) const;
		size_t rightChild(const size_t i) const;
		size_t parent(const size_t i) const;
		size_t nrChildren(const size_t i) const;

		static bool inNE(const MinMaxPST_Node& o, const MinMaxPST_Node& q);
		static bool inSE(const MinMaxPST_Node& o, const MinMaxPST_Node& q);
		static bool inNW(const MinMaxPST_Node& o, const MinMaxPST_Node& q);
		static bool inSW(const MinMaxPST_Node& o, const MinMaxPST_Node& q);
		static size_t getLevel(const size_t i);
		size_t getleftmost(const MinMaxPST_Node& o, const std::array<size_t,4>& next, std::function<bool(const MinMaxPST_Node&, const MinMaxPST_Node&)> inQuadrant, const size_t offset=0) const;

	public:
		MinMaxPST( const std::vector< MinMaxPST_Node >& points );
		~MinMaxPST();
		void printArray() const;
		std::vector< MinMaxPST_Node > getArray() const;
		std::vector< MinMaxPST_Node > enumerateUp(coord_t x0, coord_t x1, coord_t y) const;
		std::vector< MinMaxPST_Node > enumerateDown(coord_t x0, coord_t x1, coord_t y) const;
		MinMaxPST_Node leftmostNE(const MinMaxPST_Node& p) const;
		MinMaxPST_Node rightmostNW(const MinMaxPST_Node& p) const;
		MinMaxPST_Node highestNE(const MinMaxPST_Node& p) const;
		MinMaxPST_Node highestNW(const MinMaxPST_Node& p) const;
		MinMaxPST_Node leftmostSE(const MinMaxPST_Node& p) const;
		MinMaxPST_Node rightmostSW(const MinMaxPST_Node& p) const;
		MinMaxPST_Node lowestSE(const MinMaxPST_Node& p) const;
		MinMaxPST_Node lowestSW(const MinMaxPST_Node& p) const;
		MinMaxPST_Node lowest3SideUp(coord_t x0, coord_t x1, coord_t y) const;
		MinMaxPST_Node highest3SideDown(coord_t x0, coord_t x1, coord_t y) const;
		coord_t getNegINF() const;
		coord_t getPosINF() const;

		// serialization and deserialization
		MinMaxPST( const std::string );
		std::string serialize() const;
		friend std::ostream& operator<<(std::ostream& stream, const MinMaxPST& pst);
};

#endif

