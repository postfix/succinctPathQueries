#ifndef INCLUDED_SUCCINCT_TREE
#define INCLUDED_SUCCINCT_TREE

#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "bp_support_algorithm.hpp"
#include "fast_cache.hpp"
#include <vector>
#include <stack>
#include <map>
#include <set>
#include <utility>
#include <stdexcept>
#ifndef NDEBUG
#include <algorithm>
#endif
#include <iostream>

class succinct_tree {
public:
	typedef bit_vector::size_type 	    size_type;
	typedef size_type				    node_type;
	typedef bit_vector::difference_type difference_type;

	// tree info
	virtual size_type 		  size() const ;

	// navigation
	virtual node_type 		  parent( node_type x ) const ;
	virtual node_type 		  ancestor( node_type x, size_type i ) const ;
	virtual vector<node_type> children( node_type x ) const ;
	virtual node_type 		  lca( node_type x, node_type y ) const ;

	// predicates
	virtual bool			  is_ancestor( node_type p, node_type x ) const ;
	virtual bool			  is_leaf( node_type x ) const ;
};

#endif
