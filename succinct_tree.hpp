#ifndef INCLUDED_SUCCINCT_TREE
#define INCLUDED_SUCCINCT_TREE
#include "pq_types.hpp"

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "sdsl/bp_support_algorithm.hpp"
#include "sdsl/fast_cache.hpp"
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
	virtual size_type 		  size() const = 0;
	virtual double			  size_in_bytes() const = 0;

	// navigation
	virtual node_type 		  parent( const node_type x ) const = 0;
	virtual node_type 		  ancestor( const node_type x, const size_type i ) const = 0;
	virtual std::vector<node_type> children( const node_type x ) const = 0;
	virtual node_type 		  lca( const node_type x, const node_type y ) const = 0;
	virtual size_type		  depth( const node_type x ) const = 0;

	// predicates
	virtual bool			  is_ancestor( const node_type p, const node_type x ) const = 0;
	virtual bool			  is_leaf( const node_type x ) const = 0;
};

#endif
