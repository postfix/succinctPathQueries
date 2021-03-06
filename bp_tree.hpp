#ifndef INCLUDED_BP_TREE
#define INCLUDED_BP_TREE

#include "succinct_tree.hpp"
#include "sdsl/bp_support_sada.hpp"
#include "sdsl/bp_support_gg.hpp"

#include "sdsl/int_vector.hpp"
#include "sdsl/rank_support.hpp"
#include "sdsl/select_support.hpp"
#include "sdsl/bp_support_algorithm.hpp"
#include "sdsl/fast_cache.hpp"
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
#define SADA 0

class bp_tree: public succinct_tree {
public:
	typedef sdsl::bp_support_sada<> bp_support_sada;
	typedef sdsl::bp_support_gg<> bp_support_gg;
	typedef sdsl::bit_vector bit_vector;
	typedef pq_types::node_type node_type;
	typedef pq_types::size_type size_type;
private:
#if SADA
	bp_support_sada bp_sada;
#else
	bp_support_gg bp_sada;
#endif
	bit_vector *m_bv= NULL;
	size_type _sz= 0;

	/*! The position of the opening parenthesis of the node $x$
	 * \param x the pre-order number of the node
	 */
	size_type node2position( const node_type x ) const {
		return bp_sada.select(x+1);
	}
	/*! The pre-order number of the node, whose opening parenthesis is at position i
	 * \param i position of the opening parenthesis
	 */
	node_type position2node( const size_type i ) const {
		return bp_sada.rank(i)-1;
	}

	std::pair<size_type,size_type> interval( const node_type x ) const {
		size_type i = node2position(x);
		return {i,bp_sada.find_close(i)};
	}

public:
	
	// constructors
	bp_tree( const std::string &s ) {
		auto k= s.size();
		_sz= 0;
		assert( !(k&1) );
		m_bv= new bit_vector(k,0);
		assert( m_bv );
		for ( auto i = 0; i < k; ++i )
			if ( s[i] == '(' )
				(*m_bv)[i]= 1;
#if SADA
		bp_sada = bp_support_sada(m_bv);
#else
		bp_sada = bp_support_gg(m_bv);
#endif
		_sz+= sdsl::size_in_bytes(*m_bv);
	}
	
	bp_tree( bit_vector *bp ) {
		_sz= 0;
#if SADA
		bp_sada = bp_support_sada(bp);
#else
		bp_sada = bp_support_gg(bp);
#endif
		m_bv= bp;
		//_sz+= sdsl::size_in_bytes(*bp);
		//bp_sada.set_vector(bp);
	}

	// tree info
	size_type size() const {
		return bp_sada.size()>>1;
	}

	// size_in_bytes(bp_sada) does not take into consideration the backbone bitvector
	// now I fixed it
	double size_in_bytes() const {
		return sdsl::size_in_bytes(bp_sada)+sdsl::size_in_bytes(*m_bv);
	}

	// navigation
	node_type parent( const node_type x ) const {
		auto pos = node2position(x);
		if ( pos == 0 )
			return 0;
		//assert( bp_sada.is_opening(pos) );
		return position2node( bp_sada.enclose( node2position(x) ) );
	}

	node_type ancestor( const node_type x, const size_type i ) const {
		if ( x == 0 ) return 0;
#if SADA
		return position2node( bp_sada.enclose(node2position(x),i) );
#else
		return position2node( bp_sada.enclose(node2position(x)) );
#endif
	}

	std::vector<node_type> children( const node_type x ) const {
		std::vector<node_type> res;
		auto ix= interval(x);
		assert (res.size() == 0);
		for ( auto i= ix.first; i+1 < bp_sada.size() && bp_sada.is_opening(i+1); i= bp_sada.find_close(i+1) ) 
			res.push_back(position2node(i+1));
		assert( is_leaf(x) || res.size() > 0 );
		return res;
	}

	/*
	node_type lca( const node_type x, const node_type y ) const {
		if ( is_ancestor(x,y) )
			return x;
		if ( is_ancestor(y,x) )
			return y;
		auto ix = interval(x),
			 iy = interval(y);
		assert( ix.second < iy.first || iy.second < ix.first );
		if ( ix.second < iy.first ) {
	label01:
			auto irmq = bp_sada.rmq(ix.second,iy.first);
			assert( !bp_sada.is_opening(irmq) );
			auto ipre_lca = bp_sada.find_open(irmq);
			return parent(position2node(ipre_lca));
		}
		swap(ix,iy);
		goto label01;
		assert( false );
	}
	*/
	node_type lca( const node_type x, const node_type y ) const {
		if ( is_ancestor(x,y) )
			return x;
		if ( is_ancestor(y,x) )
			return y;
		auto ix = interval(x),
			 iy = interval(y);
		assert( ix.second < iy.first || iy.second < ix.first );
		if ( ix.second < iy.first ) {
	label01:
			return position2node(bp_sada.double_enclose(ix.first,iy.first));
		}
		swap(ix,iy);
		goto label01;
		assert( false );
	}

	size_type depth( const node_type x ) const {
		auto pos = node2position(x);
		if ( pos == 0 )
			return 0;
		return bp_sada.excess(pos)-1;
	}

	// predicates
	bool is_leaf( const node_type x ) const {
		auto i = node2position(x);
		return bp_sada.find_close(i) == i+1;
	}
	bool is_ancestor( const node_type x, const node_type y ) const {
		auto ix = interval(x),
			 iy = interval(y);
		return ix.first <= iy.first && iy.second <= ix.second;
	}
	~bp_tree() {
		if ( m_bv )
			delete m_bv;
	}
};

#endif
