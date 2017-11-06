/*
 * This file is subject to the terms and conditions defined in
 * file 'LICENSE.txt', which is part of this source code package.
 */
#ifndef TRIE_H
#define TRIE_H

#include <vector>
#include <string>
#include <array>
#include <map>

#include "helpers.h"

class TrieNode {
	private:
		const char content;				// label of the edge from parent
		bool end;							// true if a word ends at this node
		std::vector<TrieNode*> children;	// list of children
		size_t preorder;
		size_t postorder;
		size_t parent_preorder;
		size_t mass;

		// MODIFICATION-TOLERANT BPM
#ifdef MOD_TOLERANT
		std::map<char,size_t> lastOcc;	// preorder label of last occurence of character i on path to root
#endif

	public:
		TrieNode(char c);
		~TrieNode();
		void setWordEnd(size_t p);
		void setPreorder(size_t p);
		void setPostorder(size_t p);
		void setMass(size_t m);
		void setParentPreorder(size_t p);
		char getContent() const;
		bool isWordEnd() const;
		size_t getPreorder() const;
		size_t getPostorder() const;
		size_t getMass() const;
		size_t getParentPreorder() const;
		TrieNode* findChild(char c) const;
		bool addChild(TrieNode* p);
		std::vector<TrieNode*> getChildren() const;

		// MODIFICATION-TOLERANT BPM
#ifdef MOD_TOLERANT
 		void setLastOcc( std::map<char,size_t> l );
		const std::map<char,size_t>& getLastOcc() const;
#endif
};

class Trie {
	private:
		TrieNode* root;
		size_t nrNodes;
		std::map< std::string, size_t > proteinID;
		std::map< size_t, std::string > idProtein;
	
		bool finalized;
		bool order;
		bool mass;

		void deleteNode(TrieNode* p);
		void computeOrder();
		void computeMass();
		void finalize();
		std::string getString( size_t preorder ) const;

		// MODIFICATION-TOLERANT BPM
#ifdef MOD_TOLERANT
		bool lastOcc;
		void computeLastOcc();
#endif

		// MUTATION-TOLERANT BPM
#ifdef MUT_TOLERANT
		std::pair<std::vector<TrieNode*>,std::string> findPath( size_t p ) const;
		size_t findPreorder(std::string w) const; // find word w and output preorder
#endif

	public:
		Trie();
		~Trie();
		void add(std::string w, std::string protein);					// add word w at line l and pos p
		bool find(std::string w) const;				// find word w
		size_t size() const;							// output number of nodes
		std::map< size_t, std::vector< std::pair<size_t,size_t> > >& getNodesByMass(std::map< size_t, std::vector< std::pair<size_t,size_t> > >& t);	// return <preorder, postorder> for each node grouped by mass, calls finalize()
		std::vector< std::pair<size_t,size_t> > getLeaves() const; // return <preorder,postorder> for each trie node with wordEnd == true; assumes that computeOrder has been called
		std::ostream& outputNodes(std::ostream&); // output <preorder,postorder,parent_preorder,isWordEnd() ? 0 : content,links,lastOcc-map> for each node
		friend std::ostream& operator<<(std::ostream& stream, Trie& t);

		// MODIFICATION-TOLERANT BPM
#ifdef MOD_TOLERANT
		std::map< size_t, std::map<char,size_t> > outputLastOcc(); // return map preorder->lastOcc
#endif
		// MUTATION-TOLERANT BPM
#ifdef MUT_TOLERANT
		std::unordered_map< size_t, std::vector<size_t> > computeLinks(std::unordered_map< size_t, std::vector<size_t> >& links);
#endif

		// ONLY FOR UNIT TESTS
		TrieNode* findByPreorder(size_t p);
};

#endif
