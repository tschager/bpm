/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#include <vector>
#include <string>
#include <array>
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <stack>
#include <queue>
#include <map>
#include <set>

#include "config.h"
#include "trie.h"
#include "helpers.h"

//#define DEBUG
#ifdef DEBUG
#	define LOG(x) std::clog << "DEBUG: " << x << std::endl;
#else
#	define LOG(x) do {} while (0)
#endif


/*
 * TrieNode class implementation
 */

TrieNode::TrieNode(char c) : content(c) {
	end = false;
	preorder = -1;
	postorder = -1;
	mass = -1;
}

TrieNode::~TrieNode() {}
void TrieNode::setWordEnd(size_t protein) { end = true; }// proteins.push_back(protein); }
void TrieNode::setPreorder(size_t p) { preorder = p; }
void TrieNode::setPostorder(size_t p) { postorder = p; }
void TrieNode::setMass(size_t m) { mass = m; }
void TrieNode::setParentPreorder(size_t p) { parent_preorder = p; }
bool TrieNode::isWordEnd() const { return end; }
char TrieNode::getContent() const { return content; }
size_t TrieNode::getPreorder() const { return preorder; }
size_t TrieNode::getPostorder() const { return postorder; }
size_t TrieNode::getParentPreorder() const { return parent_preorder; }
size_t TrieNode::getMass() const { return mass; }
std::vector<TrieNode*> TrieNode::getChildren() const { return children; }
TrieNode* TrieNode::findChild(char c) const{
	for ( auto child : children ) {
		if (child->getContent() == c)
			return child;
	}
	return nullptr;
}

bool cmp( TrieNode const * l , TrieNode const * r ) {
	return l->getContent() < r->getContent();
}

bool TrieNode::addChild(TrieNode* p) {
	if ( findChild(p->getContent()) == nullptr ) {
		children.push_back(p);
		std::sort(children.begin(), children.end(), cmp);
		return true;
	} else
		return false;
}

// MODIFICATION-TOLERANT BPM
#ifdef MOD_TOLERANT
void TrieNode::setLastOcc( std::map<char,size_t> l ) { lastOcc = l; };
const std::map<char,size_t>& TrieNode::getLastOcc() const { return lastOcc; }
#endif

/*
 * Trie class implementation
 */

Trie::Trie() {
	root = new TrieNode( char(0) );
	nrNodes = 1;
	finalized = false;
	order = false;
	mass = false;
#ifdef MOD_TOLERANT
	lastOcc = false;
#endif
}

Trie::~Trie() {
	finalized = false;
	if (nrNodes > 0)
		deleteNode(root);
	assert(nrNodes == 0);
}

void Trie::add(std::string w, std::string protein) {
	assert( !finalized );
	TrieNode* cur = root;
	std::string::iterator it = w.begin();
	while ( it != w.end() ) {
		if ( cur->findChild(*it) == nullptr ) {
			TrieNode* next = new TrieNode(*it);
			cur->addChild(next);
			nrNodes++;
			cur = next;
		} else {
			cur = cur->findChild(*it);
		}
		it++;
	}
	cur->setWordEnd(proteinID.find(protein)->second);
}

bool Trie::find(std::string w) const {
	TrieNode* cur = root;
	std::string::iterator it = w.begin();
	while ( it != w.end() ) {
		if ( cur->findChild(*it) == nullptr )
			return false;
		else
			cur = cur->findChild(*it);
		LOG( "(" + std::to_string(cur->getPreorder()) + "," + std::to_string(cur->getPostorder()) + ") -- " + std::to_string(cur->getMass()) );
		it++;
	}
	return cur->isWordEnd();
}

size_t Trie::size() const { return nrNodes; }

std::map< size_t, std::vector< std::pair<size_t,size_t> > >& Trie::getNodesByMass( std::map< size_t, std::vector< std::pair<size_t,size_t> > >& res) {
	// [mass, preorder, postorder] for each node
	finalize();

	std::stack<TrieNode*> stack;
	stack.push(root);
	TrieNode* cur = nullptr;
	while (!stack.empty()) {
		cur = stack.top();
		stack.pop();
		if (res.find( cur->getMass() ) == res.end())
			res.insert( std::make_pair( 
						cur->getMass(), 
						std::vector< std::pair<size_t,size_t> >() ) );
		res.find(cur->getMass())->second.push_back( std::make_pair( cur->getPreorder(), cur->getPostorder()) );
		for ( auto c : cur->getChildren() )
			stack.push(c);
	}
	return res;
}

void Trie::deleteNode(TrieNode * cur) {
	assert( !finalized );
	for (auto c : cur->getChildren() )
		deleteNode(c);
	delete cur;
	nrNodes--;
}

void Trie::finalize() {
	finalized = true;
	if(!order) computeOrder();
	if(!mass) computeMass();
#ifdef MOD_TOLERANT
	if(!lastOcc) computeLastOcc();
#endif
}

void Trie::computeOrder() {
	assert( finalized );
	if (order) return;

	// preorder
	size_t i = 0;
	std::stack<TrieNode*> stack;
	stack.push(root);
	TrieNode* cur = nullptr;
	root->setParentPreorder(0);
	while(!stack.empty()) {
		TrieNode* cur = stack.top();
		stack.pop();
		cur->setPreorder(i);
		const std::vector<TrieNode*> children = cur->getChildren();
		for ( auto it = children.rbegin(); it != children.rend(); it++ ) {
			stack.push(*it);
			(*it)->setParentPreorder(i);
		}
		i++;
	}

	// postorder (using two stacks)
	i = 0;
	std::stack<TrieNode*> s1;
	std::stack<TrieNode*> s2;
	s1.push(root);
	cur = nullptr;
	while(!s1.empty()) {
		s2.push( s1.top() );
		cur = s1.top();
		s1.pop();
		for ( auto c : cur->getChildren() )
			s1.push(c);
	}
	while(!s2.empty()) {
		s2.top()->setPostorder(i);
		i++;
		s2.pop();
	}

	order = true;
}

void Trie::computeMass() {
	assert( finalized );
	if (mass) return;

	std::stack<TrieNode*> stack;
	stack.push(root);
	root->setMass(0);
	TrieNode* cur = nullptr;
	while(!stack.empty()) {
		cur = stack.top();
		stack.pop();
		assert( cur->getMass() >= 0 );
		for ( auto c : cur->getChildren() ) {
			assert(cfg::AAmasses.find(c->getContent()) != cfg::AAmasses.end());
			c->setMass( cur->getMass() + cfg::AAmasses.at( c->getContent() ) );
			stack.push(c);
		}
	}

	mass = true;
}

std::vector< std::pair<size_t,size_t> > Trie::getLeaves() const {
	assert( finalized && order );

	std::vector< std::pair<size_t,size_t> > res;

	std::stack<TrieNode*> stack;
	stack.push(root);
	TrieNode* cur = nullptr;
	while(!stack.empty()) {
		cur = stack.top();
		stack.pop();
		if (cur->isWordEnd())
			res.push_back( std::make_pair( cur->getPreorder(), cur->getPostorder() ) );
		const std::vector<TrieNode*> children = cur->getChildren();
		for ( auto it = children.rbegin(); it != children.rend(); it++ )
			stack.push(*it);
	}

	return res;
}

std::string Trie::getString( size_t preorder ) const {
	assert( finalized && order );
	assert( preorder < nrNodes );

	std::string s;
	TrieNode* cur = root;
	TrieNode* next = nullptr;
	while ( cur->getPreorder() != preorder ) {
		assert(  cur->getPreorder() < preorder );
		for (auto c : cur->getChildren()) {
			if (c->getPreorder() <= preorder)
				next = c;
		}
		assert( next != nullptr );
		s.push_back(next->getContent());
		cur = next;
		next = nullptr;
	}
	return s;
}

std::ostream& Trie::outputNodes (std::ostream& o) {
	assert( finalized && order );
	std::string res;

	LOG("Start outputNodes");
	o << std::to_string(nrNodes) + '\n';
	
	std::stack<TrieNode*> stack;
	stack.push(root);
	TrieNode* cur = nullptr;
#ifdef MOD_TOLERANT
	std::vector<char> AA;
	std::vector<char>::iterator it;
	for ( auto a : cfg::allSitesMods ) {
		it = std::lower_bound(AA.begin(), AA.end(), a.first);
		if (it == AA.end() || *it != a.first) {
			AA.push_back(a.first);
			std::sort(AA.begin(), AA.end());
		}
	}
	for ( auto a : cfg::cTermMods ) {
		it = std::lower_bound(AA.begin(), AA.end(), a.first);
		if (it == AA.end() || *it != a.first) {
			AA.push_back(a.first);
			std::sort(AA.begin(), AA.end());
		}
	}
	for ( auto a : cfg::nTermMods ) {
		it = std::lower_bound(AA.begin(), AA.end(), a.first);
		if (it == AA.end() || *it != a.first) {
			AA.push_back(a.first);
			std::sort(AA.begin(), AA.end());
		}
	}
	computeLastOcc();
	res += "LASTOCC-ORDER:" + std::to_string(AA.size()) + ":";
	for (auto a : AA)
		res += std::string(1,a) + ",";
	res.pop_back();
	res += '\n';
	o << res;
	res.clear();

	LOG("LastOcc annotation written");
#endif
#ifdef MUT_TOLERANT
	std::unordered_map< size_t, std::vector<size_t> > links;
	this->computeLinks(links);
	LOG("Link computation done");
#endif

	// output root information
	res += std::to_string(root->getPreorder()) + "," +
		   std::to_string(root->getPostorder()) + "," +
		   std::to_string(root->getParentPreorder()) + "," +
		   ",";
#ifdef MOD_TOLERANT
	for (auto a : AA)
		res += std::to_string(root->getLastOcc().find(a)->second) + ",";
	res.pop_back();
#endif
	res += ",";
	res += '\n';
	o << res;
	res.clear();

	while(!stack.empty()) {
		cur = stack.top();
		stack.pop();
		for ( auto it : cur->getChildren() ) {
			res += std::to_string(it->getPreorder()) + "," + 
				std::to_string(it->getPostorder()) + "," +
				std::to_string(it->getParentPreorder()) + ",";
			if (it->isWordEnd())
				res += getString(it->getPreorder()) + ",";
			else
				res += ",";
			// LastOcc
#ifdef MOD_TOLERANT
			for (auto a : AA)
				res += std::to_string(it->getLastOcc().find(a)->second) + ",";
			res.pop_back();
#endif
			res += ",";
			// Links
#ifdef MUT_TOLERANT
			if (links.find(it->getPreorder()) != links.end()) {
				res += "|" + std::to_string(links.find(it->getPreorder())->second.size()) + ":";
					for ( auto l : links.find(it->getPreorder())->second )
						res += std::to_string(l) + ",";
					res.pop_back();
			}
#endif
			res += '\n';
			o << res;
			res.clear();
			stack.push(it);
		}
	}
	return o;
}

std::ostream& operator<< (std::ostream& o, Trie& t) {
	return t.outputNodes(o);
}
/* 
 * MODIFICATION-TOLERANT BPM
 */
#ifdef MOD_TOLERANT
void Trie::computeLastOcc() {
	assert(finalized);
	if (lastOcc) return;

	std::map<char,size_t> lastOcc;
	for ( auto a : cfg::allSitesMods ) {
		if (lastOcc.find(a.first) == lastOcc.end())
			lastOcc.insert(std::make_pair(a.first,0));
	}
	for ( auto a : cfg::cTermMods ) {
		if (lastOcc.find(a.first) == lastOcc.end())
			lastOcc.insert(std::make_pair(a.first,0));
	}
	for ( auto a : cfg::nTermMods ) {
		if (lastOcc.find(a.first) == lastOcc.end())
			lastOcc.insert(std::make_pair(a.first,0));
	}

	std::queue<TrieNode*> queue;
	queue.push(root);
	root->setLastOcc(lastOcc);
	TrieNode* cur = nullptr;
	while(!queue.empty()) {
		cur = queue.front();
		queue.pop();
		lastOcc = cur->getLastOcc();
		if (lastOcc.find(cur->getContent()) != lastOcc.end())
			lastOcc.find(cur->getContent())->second = cur->getPreorder();
		cur->setLastOcc(lastOcc);
		const std::vector<TrieNode*> children = cur->getChildren();
		for ( auto it = children.rbegin(); it != children.rend(); it++ ) {
			(*it)->setLastOcc(lastOcc);
			queue.push(*it);
		}
	}

	this->lastOcc = true;
}

std::map<size_t,std::map<char,size_t>> Trie::outputLastOcc() {
	if(!finalized)
		finalize();
	assert(finalized);
	std::map<size_t,std::map<char,size_t>> res;
	computeOrder();
	computeLastOcc();

	std::stack<TrieNode*> stack;
	stack.push(root);
	TrieNode* cur = nullptr;
	while (!stack.empty()) {
		cur = stack.top();
		stack.pop();
		res.insert( std::make_pair(cur->getPreorder(), cur->getLastOcc()) );
		for ( auto c : cur->getChildren() )
			stack.push(c);
	}
	return res;
}
#endif


/*
 * MUTATION-TOLERANT BPM
 */
#ifdef MUT_TOLERANT
size_t Trie::findPreorder(std::string w) const {
	TrieNode* cur = root;
	std::string::iterator it = w.begin();
	while ( it != w.end() ) {
		if ( cur->findChild(*it) == nullptr )
			return false;
		else
			cur = cur->findChild(*it);
		it++;
	}
	return cur->getPreorder();
}

std::pair<std::vector<TrieNode*>,std::string> Trie::findPath( size_t p ) const {
	assert( finalized );
	assert( p < nrNodes );
	std::string s;
	std::vector<TrieNode*> stack;
	stack.push_back(root);
	TrieNode* cur = stack.back();
	TrieNode* next = nullptr;
	while ( cur->getPreorder() != p ) {
		assert(  cur->getPreorder() < p );
		for (auto c : cur->getChildren()) {
			if (c->getPreorder() <= p)
				next = c;
		}
		assert( next != nullptr );
		s.push_back(next->getContent());
		stack.push_back(next);
		cur = next;
		next = nullptr;
	}
	if (cur->getPreorder() == p)
		return std::make_pair(stack,s);
	else
		return std::make_pair(std::vector<TrieNode*>(),"");
}

std::unordered_map< size_t, std::vector<size_t> > Trie::computeLinks( std::unordered_map< size_t, std::vector<size_t> >& links ) {
	if(!finalized)
		finalize();
	assert( finalized && order);

	std::string s;
	std::vector<TrieNode*> stack;
	stack.push_back(root);
	TrieNode* cur = stack.back();
	TrieNode* next = nullptr;
	for (size_t i = 1; i < nrNodes; i++) {
		next = nullptr;
		for (auto c : cur->getChildren()) {
			if (c->getPreorder() == i) {
				next = c;
				break;
			}
		}
		if (next == nullptr) {
			std::pair<std::vector<TrieNode*>,std::string> tmp = findPath(i);
			stack = tmp.first;
			s = tmp.second;
			next = stack.back();
			stack.pop_back();
		} else 
			s.push_back(next->getContent());

		// compute links
		for (size_t j = 1; j < s.size(); j++) {
			size_t linkPartner = findPreorder( s.substr(s.size()-j,s.size()) );
			if (linkPartner == 0) continue;
			if (links.find(linkPartner) == links.end())
				links.insert(std::make_pair(linkPartner, std::vector<size_t>()));
			const size_t linkEnd = stack.at(stack.size()-j)->getPreorder();
			if (linkEnd == 0) continue;;
			links.find(linkPartner)->second.push_back(linkEnd);
		}
		stack.push_back(next);
		cur = next;
	}
	assert( links.size() < nrNodes );
	return links;
}
#endif

// ONLY FOR UNIT TESTS
TrieNode* Trie::findByPreorder(size_t p) {
	if (!finalized) finalize();
	assert(finalized);
	assert( p < nrNodes );
	TrieNode* cur = root;
	TrieNode* next = nullptr;
	while ( cur->getPreorder() != p ) {
		assert(  cur->getPreorder() < p );
		for (auto c : cur->getChildren()) {
			if (c->getPreorder() <= p)
				next = c;
		}
		assert( next != nullptr );
		cur = next;
		next = nullptr;
	}
	return (cur->getPreorder() == p) ? cur : nullptr;
}
