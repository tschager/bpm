/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <vector>
#include <algorithm>
#include <iostream>
#include <math.h>
#include <functional>
#include <limits>
#include <assert.h>
#include <sstream>

#include "minmaxpst.h"

//#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif


typedef size_t coord_t;
typedef std::pair<coord_t,coord_t> MinMaxPST_Node;

constexpr size_t floorLog2(size_t n) { return ( (n<2) ? 0 : 1+floorLog2(n/2)); }
constexpr size_t pow2(size_t n) { return 1 << n; }


// 1-based
const MinMaxPST_Node& MinMaxPST::get(size_t i) const {
//	assert(i > 0 && i <= t.size());
	return t.at(i-1); 
}

// 1-based
size_t MinMaxPST::smallestYCoordIndex(size_t start, size_t end) const {
//	assert( start > 0 && end <= t.size());
	size_t res = start;
	for ( size_t i = start; i <= end; i++ )
		res = (get(res).second > get(i).second) ? i : res;
	return res;
}

// 1-based
size_t MinMaxPST::largestYCoordIndex(size_t start, size_t end) const {
//	assert( start > 0 && end <= t.size());
	size_t res = start;
	for ( size_t i = start; i <= end; i++ )
		res = (get(res).second > get(i).second) ? res : i;
	return res;
}

void MinMaxPST::swap(size_t i, size_t j) {
//	assert( i > 0 && j > 0 && i <= t.size() && j <= t.size());
	iter_swap(t.begin() + i-1, t.begin() + j-1);
}

MinMaxPST::MinMaxPST( const std::string s ) {
	posINF = std::numeric_limits<coord_t>::has_infinity ? std::numeric_limits<coord_t>::infinity() : std::numeric_limits<coord_t>::max();
	negINF = std::numeric_limits<coord_t>::has_infinity ? -std::numeric_limits<coord_t>::infinity() : std::numeric_limits<coord_t>::min();
	LOG("Deserialization of: " + s);
	std::string::size_type pos = s.find(':');
	if (pos == std::string::npos)
		std::cout << "ERROR in deserialization - size of PST not found." << std::endl;
	std::stringstream sstream( s.substr(0,pos) );
	size_t size;
	sstream >> size;
	LOG("Reserve: " + s.substr(0,pos) + " entries");
	t.reserve( size_t(size) );

	assert( s.at(pos+1) == '(' );
	std::string::size_type start = pos+1;
	std::string::size_type mid = s.find(',',pos+1);
	std::string::size_type end = s.find(')',pos+1);
	MinMaxPST_Node next;
	size_t first, second;
	while (end != std::string::npos) {
		LOG("next point: " + s.substr(start,mid-start) + " , " + s.substr(mid+1, end-mid));
		std::stringstream s1( s.substr(start+1,mid-start) );
		s1 >> first;
		std::stringstream s2( s.substr(mid+1,end-mid) );
		s2 >> second;

		next = std::make_pair( first, second );
		t.push_back(next);
		start = s.find('(',end+1);
		mid = s.find(',',end+1);
		end = s.find(')',end+1);
	}
}
MinMaxPST::MinMaxPST( const std::vector< MinMaxPST_Node >& points ) {

	posINF = std::numeric_limits<coord_t>::has_infinity ? std::numeric_limits<coord_t>::infinity() : std::numeric_limits<coord_t>::max();
	negINF = std::numeric_limits<coord_t>::has_infinity ? -std::numeric_limits<coord_t>::infinity() : std::numeric_limits<coord_t>::min();

	t.resize(points.size());
	std::copy ( points.begin(), points.end(), t.begin() );
//	t.reserve(points.size());
//	std::copy( points.begin(), points.end(), std::back_inserter(t) );
	const size_t h = floorLog2(points.size());			// tree height
	const size_t A = t.size() - (pow2(h) - 1);			// nr of leaf nodes
	
	// sort t by x coordinate
	std::sort(t.begin(),t.end());

	// build levels
	for ( size_t i = 0; i < h; i++ ) {
		// pow2(i) = nr of nodes in level i
		// k = nr of nodes in level i that are roots of subtrees of size k1
		// k1 = largest size of subtrees of nodes in level i
		// k2 = size of subtree of one node at level i
		// k3 = size of the subtrees of the remaining 2^i-1-k nodes in level i
		const size_t k = static_cast<size_t>(std::floor( A/pow2(h-i) )); 
		const size_t k1 = pow2(h+1-i) - 1;
		const size_t k2 = pow2(h-i) - 1 + A - k * pow2(h-i);
		const size_t k3 = pow2(h-i) - 1;

		// find the root (i.e. smallest/largest y coordinate) of the k subtrees of size k1 
		// and swap it to the correct position
		for ( size_t j = 1; j <= k; j++ ) {
			// find index with smallest/largest y coordinate
			size_t l;
			if ( i % 2 )
				l = largestYCoordIndex(pow2(i) + (j-1)*k1, pow2(i) + j*k1-1);
			else
				l = smallestYCoordIndex(pow2(i) + (j-1)*k1, pow2(i) + j*k1-1);

			// swap point with smallest/largest y coordinate
			swap(l,pow2(i) + j - 1);
		}

		// for the remaining subtrees: one has size k2, the remaining have size k3
		// again, find and swap the node with smallest/largest y coordinate
		if (k < pow2(i)) {
			size_t l;
			if ( i % 2 ) 
				l = largestYCoordIndex(pow2(i) + k*k1, pow2(i) + k*k1 + k2 - 1);
			else
				l = smallestYCoordIndex(pow2(i) + k*k1, pow2(i) + k*k1 + k2 - 1);
			swap(l,pow2(i)+k);
			const size_t m = pow2(i) + k*k1 + k2;
			for ( size_t j = 1; j < pow2(i)-k; j++ ) {
				if ( i % 2 ) 
					l = largestYCoordIndex(m+(j-1)*k3,m+j*k3-1);
				else
					l = smallestYCoordIndex(m+(j-1)*k3,m+j*k3-1);
				swap(l,pow2(i)+k+j);
			}
		}

		// finally, sort the remaining array with nodes of levels i+1...h
		std::sort(t.begin()+pow2(i+1)-1, t.end());
	}
}

MinMaxPST::~MinMaxPST() {}

void MinMaxPST::printArray() const {
	for (auto p : t) 
		std::cout << "(" << p.first << "," << p.second << ") ";
	std::cout << std::endl;
}

std::vector< MinMaxPST_Node > MinMaxPST::getArray() const { return t; }

bool MinMaxPST::isLeaf(const size_t i) const {
//	assert(i>0 && i <= t.size());
	return (2*i > t.size());
}
size_t MinMaxPST::getLevel(const size_t i) { return floorLog2(i); }
size_t MinMaxPST::parent(const size_t i) const { return i/2; }
size_t MinMaxPST::leftChild(const size_t i) const { return 2*i; }
size_t MinMaxPST::rightChild(const size_t i) const { return 2*i+1; }
size_t MinMaxPST::nrChildren(const size_t i) const {
//	assert(i>0 && i <= t.size());
	if (isLeaf(i))
		return 0;
	else if (leftChild(i) == t.size())
		return 1;
	else {
		//assert(rightChild(i) <= t.size());
		return 2;
	}
}
bool MinMaxPST::inNE(const MinMaxPST_Node& o, const MinMaxPST_Node& q) {
	return o.first <= q.first && o.second <= q.second;
}
bool MinMaxPST::inNW(const MinMaxPST_Node& o, const MinMaxPST_Node& q) {
	return o.first >= q.first && o.second <= q.second;
}
bool MinMaxPST::inSE(const MinMaxPST_Node& o, const MinMaxPST_Node& q) {
	return o.first <= q.first && o.second >= q.second;
}
bool MinMaxPST::inSW(const MinMaxPST_Node& o, const MinMaxPST_Node& q) {
	return o.first >= q.first && o.second >= q.second;
}
coord_t MinMaxPST::getPosINF() const { return posINF; }
coord_t MinMaxPST::getNegINF() const { return negINF; }

size_t MinMaxPST::getleftmost(const MinMaxPST_Node& o, const std::array<size_t,4>& next, std::function<bool(const MinMaxPST_Node&, const MinMaxPST_Node&)> inQuadrant, const size_t offset) const {
	auto it = next.begin() + offset;
	auto resPointer = next.end();
	while (it < next.end()) {
		if ( inQuadrant(o,get(*it)) ||
			((nrChildren(*it) >= 1) && inQuadrant(o,get(leftChild(*it))) ) ||
			((nrChildren(*it) == 2) && inQuadrant(o,get(rightChild(*it))) ) ) {
			if (resPointer != next.end())
				resPointer = (get(*resPointer).first > get(*it).first) ? it : resPointer;
			else
				resPointer = it;
		}
		it++;
	}
	if (resPointer == next.end()) 
		resPointer = next.begin() + offset;
	return *resPointer;
}

MinMaxPST_Node MinMaxPST::leftmostNE(const MinMaxPST_Node& o) const { 
	MinMaxPST_Node best = std::make_pair(posINF,posINF);
	size_t p = 1;
	size_t q = 1;

	while ( !isLeaf(p) ) {
		//assert( p <= q );

		// update best
		best = (inNE(o,get(p)) && best.first > get(p).first) ? get(p) : best;
		best = (inNE(o,get(q)) && best.first > get(q).first) ? get(q) : best;

		if (p == q) {
			q = (nrChildren(p) == 1) ? leftChild(p) : rightChild(p);
			p = leftChild(p);
		} else {
			if (isLeaf(q))
				q = p;
			else if (nrChildren(q) == 1) {
				if ( get(leftChild(q)).first < o.first ) {
					p = leftChild(q);
					q = leftChild(q);
				} else if ( get(rightChild(p)).first < o.first ) {
					p = rightChild(p);
					q = leftChild(q);
				} else {
					q = rightChild(p);
					p = leftChild(p);
				}
			} else {
				// next is sorted by x-coordinate
				const std::array<size_t,4> next = { 
					leftChild(p), rightChild(p), 
					leftChild(q), rightChild(q) 
				};

				if (o.first < get(leftChild(p)).first) {
					p = q = getleftmost(o,next,inNE);
				} else if (o.first >= get(rightChild(q)).first) {
					p = q = rightChild(q);
				} else {
					size_t i;
					if ( get(leftChild(p)).first <= o.first && o.first < get(rightChild(p)).first )
						i = 0;
					else if ( get(rightChild(p)).first <= o.first && o.first < get(leftChild(q)).first )
						i = 1;
					else if ( get(leftChild(q)).first <= o.first && o.first < get(rightChild(q)).first )
						i = 2;
					else
						i = 3;
					p = next.at(i);
					q = getleftmost(o,next,inNE,i);
				}
			}
		}
	}

	// update best
	best = (inNE(o,get(p)) && best.first > get(p).first) ? get(p) : best;
	best = (inNE(o,get(q)) && best.first > get(q).first) ? get(q) : best;
	return best;
};
MinMaxPST_Node MinMaxPST::leftmostSE(const MinMaxPST_Node& o) const { 
	MinMaxPST_Node best = std::make_pair(posINF,negINF);
	size_t p = 1;
	size_t q = 1;

	while ( !isLeaf(p) ) {
		//assert( p <= q );
		// update best
		best = (inSE(o,get(p)) && best.first > get(p).first) ? get(p) : best;
		best = (inSE(o,get(q)) && best.first > get(q).first) ? get(q) : best;

		if (p == q) {
			q = (nrChildren(p) == 1) ? leftChild(p) : rightChild(p);
			p = leftChild(p);
		} else {
			if (isLeaf(q))
				q = p;
			else if (nrChildren(q) == 1) {
				if ( get(leftChild(q)).first < o.first ) {
					p = leftChild(q);
					q = leftChild(q);
				} else if ( get(rightChild(p)).first < o.first ) {
					p = rightChild(p);
					q = leftChild(q);
				} else {
					q = rightChild(p);
					p = leftChild(p);
				}
			} else {
				// next is sorted by x-coordinate
				const std::array<size_t,4> next = { 
					leftChild(p), rightChild(p), 
					leftChild(q), rightChild(q) 
				};

				if (o.first < get(leftChild(p)).first) {
					p = q = getleftmost(o,next,inSE);
				} else if (o.first >= get(rightChild(q)).first) {
					p = q = rightChild(q);
				} else {
					size_t i;
					if ( get(leftChild(p)).first <= o.first && o.first < get(rightChild(p)).first )
						i = 0;
					else if ( get(rightChild(p)).first <= o.first && o.first < get(leftChild(q)).first )
						i = 1;
					else if ( get(leftChild(q)).first <= o.first && o.first < get(rightChild(q)).first )
						i = 2;
					else
						i = 3;
					p = next.at(i);
					q = getleftmost(o,next,inSE,i);
				}
			}
		}
	}

	// update best
	best = (inSE(o,get(p)) && best.first > get(p).first) ? get(p) : best;
	best = (inSE(o,get(q)) && best.first > get(q).first) ? get(q) : best;
	return best;
};

MinMaxPST_Node MinMaxPST::rightmostNW(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::highestNE(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::highestNW(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::rightmostSW(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::lowestSE(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::lowestSW(const MinMaxPST_Node& o) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::lowest3SideUp(coord_t x1, coord_t x2, coord_t y) const { return std::make_pair(0,0); };
MinMaxPST_Node MinMaxPST::highest3SideDown(coord_t x1, coord_t x2, coord_t y) const { return std::make_pair(0,0); };
std::vector< MinMaxPST_Node > MinMaxPST::enumerateUp(coord_t x1, coord_t x2, coord_t y) const { return std::vector< MinMaxPST_Node >(); };

std::string MinMaxPST::serialize() const {
	std::string res = std::to_string(t.size()) + ":";
	for (auto n : t)
		res += "(" + std::to_string(n.first) + "," + std::to_string(n.second) + ")";
	return res;
}
std::ostream& operator<<(std::ostream& stream, const MinMaxPST& pst)  {
	return stream << pst.serialize();
}
