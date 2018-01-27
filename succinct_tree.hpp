#ifndef INCLUDED_SUCCINCT_TREE
#define INCLUDED_SUCCINCT_TREE
#include "pq_types.hpp"

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
	typedef pq_types::size_type 	    size_type;
	typedef pq_types::node_type			node_type;

	// tree info
	virtual size_type 		  size() const ;

	// navigation
	virtual node_type 		  parent( const node_type x ) const ;
	virtual node_type 		  ancestor( const node_type x, const size_type i ) const ;
	virtual vector<node_type> children( const node_type x ) const ;
	virtual node_type 		  lca( const node_type x, const node_type y ) const ;
	virtual size_type		  depth( const node_type x ) const ;

	// predicates
	virtual bool			  is_ancestor( const node_type p, const node_type x ) const ;
	virtual bool			  is_leaf( const node_type x ) const ;
};

#endif
