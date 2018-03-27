#include <algorithm>
#include <cassert>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "path_query_processor.hpp"
#include "pq_types.hpp"
#include "succinct_tree.hpp"
#include "rs_bitvector01.hpp"
#include "sdsl/bp_support_gg.hpp"
#include "sdsl/wt_huff.hpp"
#include "sdsl/wt_int.hpp"
#include "bp_tree.hpp"
#include <cstring>
#include <sstream>
#include <cstdio>
#include "sdsl/io.hpp"
#define oo (0xfffffffful)
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;
typedef sdsl::bp_support_gg<> bp_support_gg;

template<  
		  class t_bp_support= bp_support_gg,
		  class t_wavelet_tree= wt_int>
class multiple_hierarchies: public path_query_processor {
		
	uint32_t num_par_types;
	size_type n;
	t_wavelet_tree *S;
	bit_vector *base_structure= nullptr, 
			   **sub_tree;
	t_bp_support *structure= nullptr,
				 **B;

	/*! The position of the opening parenthesis of the node $x$
	 * \param x the pre-order number of the node
	 */
	inline size_type node2position( const node_type x ) const {
		return structure->select(x+1);
	}
	inline size_type node2position( const t_bp_support *s, const node_type x ) const {
		return s->select(x+1);
	}
	/*! The pre-order number of the node, whose opening parenthesis is at position i
	 * \param i position of the opening parenthesis
	 */
	inline node_type position2node( const size_type i ) const {
		return structure->rank(i)-1;
	}
	inline node_type position2node( const t_bp_support *s, const size_type i ) const {
		return s->rank(i)-1;
	}

	inline std::pair<size_type,size_type> interval( const node_type x ) const {
		size_type i= node2position(x);
		return {i,structure->find_close(i)};
	}

	/*
	 * returns the depth including the dummy root
	 */
	inline size_type depth( node_type x ) const  { 
		size_type s= structure->excess(node2position(x));  
		assert( s >= 1 );
		return s-1;
	}
	inline size_type depth( const t_bp_support *st,  node_type x ) const  { 
		size_type s= st->excess(node2position(st,x));  
		assert( s >= 1 );
		return s-1;
	}

	inline node_type parent( node_type x ) const { 
		size_type pos= node2position(x);
		auto res= structure->enclose(pos);
		if ( res == structure->size() )
			return +oo;
		return position2node(res);
	}
	inline node_type parent( const t_bp_support *s, node_type x ) const { 
		size_type pos= node2position(s,x);
		auto res= s->enclose(pos);
		if ( res == s->size() )
			return +oo;
		return position2node(s,res);
	}

	node_type image_of( node_type x, value_type i ) const {
		if ( x == +oo )
			return x;

		auto pos= node2position(x);

		if ( (*S)[pos] == i ) {
			auto l= S->rank(pos,i);
			assert( B[i]->is_opening(l) );
			return position2node(B[i],l);
		}

		auto l= S->rank(pos,i);
		if ( !l ) return +oo;

		if ( B[i]->is_opening(l-1) )
			return position2node(B[i],l-1);
		l= B[i]->find_open(l-1);
		assert( B[i]->is_opening(l) );
		auto res= B[i]->enclose(l);
		if ( res == B[i]->size() )
			return +oo;
		return position2node(B[i],res);
	}

	bool is_ancestor( const node_type x, const node_type y ) const {
		auto ix= interval(x),
			 iy= interval(y);
		return ix.first <= iy.first && iy.second <= ix.second;
	}

	node_type _lca( const node_type x, const node_type y ) const {
		if ( is_ancestor(x,y) )
			return x;
		if ( is_ancestor(y,x) )
			return y;
		auto ix= interval(x),
			 iy= interval(y);
		assert( ix.second < iy.first || iy.second < ix.first );
		if ( ix.second < iy.first ) {
	label01:
			auto res= structure->double_enclose(ix.first,iy.first);
			if ( res == structure->size() )
				return +oo;
			return position2node(res);
		}
		swap(ix,iy);
		goto label01;
		assert( false );
	}

public:

	multiple_hierarchies( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2;
		value_type mi= 0;
		assert( s.length()/2 == wgt.size() );
		num_par_types= 0;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > num_par_types )
				num_par_types= wgt[i];
		size_type vc[num_par_types]= {0};

		sub_tree= new bit_vector *[num_par_types];
		B= new t_bp_support *[num_par_types];

		for ( auto i= 0; i < wgt.size(); ++vc[wgt[i++]-1] ) ;
		for ( auto i= 0; i < num_par_types; ++i )
			sub_tree[i]= new bit_vector(2*vc[i]), vc[i]= 0;

		base_structure= new bit_vector(2*size());
		int_vector weights= int_vector(size()*2);
		std::stack<node_type> st;
		node_type V= 0, kk= 0;
		for ( auto i= 0; i < 2*size(); ++i ) {
			(*base_structure)[i]= (s[i]=='('?1:0);
			if ( s[i] == '(' ) {
				st.push(V++);
				node_type ix= st.top();
				(*(sub_tree[wgt[ix]-1]))[vc[wgt[ix]-1]++]= 1;
				weights[kk++]= wgt[ix]-1;
			}
			else {
				node_type ix= st.top();
				(*(sub_tree[wgt[ix]-1]))[vc[wgt[ix]-1]++]= 0;
				weights[kk++]= wgt[ix]-1;
				st.pop();
			}
		}
		structure= new t_bp_support(base_structure);
		for ( auto i= 0; i < num_par_types; ++i )
			B[i]= new t_bp_support(sub_tree[i]);

		S= new t_wavelet_tree();
		construct_im(*S,weights);
	}

	size_type size() const { return n; }

	inline value_type weight( const node_type x ) const {
		return (*S)[node2position(x)]+1;
	}

	value_type query( const node_type x, const node_type y, const node_type z, const size_type k ) const {
		value_type wc= z<+oo ? weight_of(z):num_par_types+1;
		size_type acc= 0;
		for ( value_type i= 0; i < num_par_types; ++i ) {
			node_type ix= image_of(x,i),
					  iy= image_of(y,i),
					  iz= image_of(z,i);
			size_type dx= (ix<+oo?depth(B[i],ix)+1:0),
					  dy= (iy<+oo?depth(B[i],iy)+1:0),
					  dz= (iz<+oo?depth(B[i],iz)+1:0),
					  dw= dx+dy+(wc==i+1?1:0)-2*dz;
			if ( acc+dw > k )
				return (i+1);
			acc+= dw;
		}
		assert( false );
	}

	value_type query( const node_type x, const node_type y, const size_type k ) const {
		node_type z= _lca(x,y);
		value_type wc= z<+oo ? weight_of(z):num_par_types+1;
		size_type acc= 0;
		for ( value_type i= 0; i < num_par_types; ++i ) {
			node_type ix= image_of(x,i),
					  iy= image_of(y,i),
					  iz= image_of(z,i);
			size_type dx= (ix<+oo?depth(B[i],ix)+1:0),
					  dy= (iy<+oo?depth(B[i],iy)+1:0),
					  dz= (iz<+oo?depth(B[i],iz)+1:0),
					  dw= dx+dy+(wc==i+1?1:0)-2*dz;
			if ( acc+dw > k )
				return (i+1);
			acc+= dw;
		}
		assert( false );
	}

	value_type query( const node_type x, const node_type y ) const {
		node_type z= _lca(x,y);
		assert( z < +oo );
		value_type wc= weight_of(z);
		size_type acc= 0, k= (depth(x)+depth(y)+1-2*depth(z))>>1;
		for ( value_type i= 0; i < num_par_types; ++i ) {
			node_type ix= image_of(x,i),
					  iy= image_of(y,i),
					  iz= image_of(z,i);
			size_type dx= (ix<+oo?depth(B[i],ix)+1:0),
					  dy= (iy<+oo?depth(B[i],iy)+1:0),
					  dz= (iz<+oo?depth(B[i],iz)+1:0),
					  dw= dx+dy+(wc==i+1?1:0)-2*dz;
			if ( acc+dw > k )
				return (i+1);
			acc+= dw;
		}
		assert( false );
	}

	inline value_type weight_of( const node_type x ) const {
		return weight(x);
	}

	double bits_per_node() const {
		if ( !structure ) return 0.00;
		double ans= 8*(sdsl::size_in_bytes(*structure)+(sdsl::size_in_bytes(*S)));
		ans+= 8*(sizeof(base_structure) + sizeof(S) + sizeof(n));
		if ( base_structure )
			ans+= 8*sdsl::size_in_bytes(*base_structure);
		for ( auto i= 0; i < num_par_types; ++i ) {
			ans+= 8*sizeof(sub_tree[i])+8*sizeof(B[i]);
			if ( sub_tree[i] ) {
				ans+= 8*sdsl::size_in_bytes(*(sub_tree[i]));
				ans+= 8*sdsl::size_in_bytes(*(B[i]));
			}
		}
		return ans/size();
	}

	inline node_type lca( node_type x, node_type y ) const {
		return _lca(x,y);
	}

};

