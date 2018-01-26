#ifndef INCLUDED_BP_TREE
#define INCLUDED_BP_TREE

#include "succinct_tree.hpp"
#include "bp_support_sada.hpp"

#include "int_vector.hpp"
#include "rank_support.hpp"
#include "select_support.hpp"
#include "bp_support_algorithm.hpp"
#include "fast_cache.hpp"
#include <cassert>
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

class bp_tree: public succinct_tree {

	std::unique_ptr<bp_support_sada> bp_sada = nullptr;

	/*! The position of the opening parenthesis of the node $x$
	 * \param x the pre-order number of the node
	 */
	size_type node2position( node_type x ) const {
		return bp_sada->select(x+1);
	}
	/*! The pre-order number of the node, whose opening parenthesis is at position i
	 * \param i position of the opening parenthesis
	 */
	node_type position2node( size_type i ) const {
		return bp_sada->rank(i)-1;
	}

	pair<size_type,size_type> interval( node_type x ) const {
		size_type i = node2position(x);
		return {i,bp_sada->find_close(i)};
	}

public:
	/*
	typedef bit_vector::size_type 	    size_type;
	typedef bit_vector::difference_type difference_type;
	typedef uint32_t					node_type;
	*/

	bp_tree( const bit_vector *bp ) {
		bp_sada = make_unique<bp_support_sada>(bp);
	}

	// tree info
	size_type size() const {
		return bp_sada->size()>>1;
	}

	// navigation
	node_type parent( node_type x ) const {
		return position2node( bp_sada->enclose( node2position(x) ) );
	}

	vector<node_type> children( node_type x ) const {
		auto ix = interval(x);
		vector<node_type> res;
		//TODO: assert res.size() == 0;
		for ( auto i = ix.first; i+1 < ix.second; i = bp_sada->find_close(i+1) ) {
			//TODO: assert( is_opening(i+1) );
			res.push_back(position2node(i+1));
		}
		return res;
	}

	node_type lca( node_type x, node_type y ) const {
		if ( is_ancestor(x,y) )
			return x;
		if ( is_ancestor(y,x) )
			return y;
		auto ix = interval(x),
			 iy = interval(y);
		assert( ix.second < iy.first || iy.second < ix.first );
		if ( ix.second < iy.first ) {
	label01:
			auto irmq = bp_sada->rmq(ix.second,iy.first);
			//TODO: assert that parenthesis as "irmq" is a closing one
			auto ipre_lca = bp_sada->find_open(irmq);
			return parent(position2node(ipre_lca));
		}
		swap(ix,iy);
		goto label01;
	}

	// predicates
	bool is_leaf( node_type x ) const {
		auto i = node2position(x);
		return bp_sada->find_close(i) == i+1;
	}
	bool is_ancestor( node_type p, node_type x ) const {
		auto ix = interval(x),
			 iy = interval(y);
		return ix.first <= iy.first && iy.second <= ix.second ||\
			   iy.first <= ix.first && ix.second <= iy.second;
	}
};

#endif
