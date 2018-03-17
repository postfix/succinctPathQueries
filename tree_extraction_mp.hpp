/*
 * multipar tree extraction without dummy node
 */
#include <algorithm>
#include <cassert>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "path_query_processor.hpp"
#include "pq_types.hpp"
#include "succinct_tree.hpp"
#include "rs_bitvector01.hpp"
#include "sdsl/bp_support_gg.hpp"
#include "bp_tree.hpp"
#include <cstring>
#include <sstream>
#include <cstdio>
#include "sdsl/io.hpp"
#define oo (0xfffffffful)
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;
typedef sdsl::bp_support_gg<> bp_support_gg;

template<class t_bp_support= bp_support_gg, uint8_t r=2>
class tree_extraction_mp: public path_query_processor {

private:

	value_type a,b;
	size_type n;
	bool is_homogeneous() const { return a==b; }

	bit_vector *base_structure= nullptr, *sub_tree[r];
	const t_bp_support *structure;
	rs_bitvector01 *S= nullptr;
	t_bp_support *B[r];

	tree_extraction_mp *t[r];

	/*! The position of the opening parenthesis of the node $x$
	 * \param x the pre-order number of the node
	 */
	size_type node2position( const node_type x ) const {
		return structure->select(x+1);
	}
	/*! The pre-order number of the node, whose opening parenthesis is at position i
	 * \param i position of the opening parenthesis
	 */
	node_type position2node( const size_type i ) const {
		return structure->rank(i)-1;
	}

	std::pair<size_type,size_type> interval( const node_type x ) const {
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

	inline node_type parent( node_type x ) const { 
		size_type pos= node2position(x);
		auto res= structure->enclose(pos);
		if ( res == structure->size() )
			return +oo;
		return position2node(res);
	}

	node_type image_of( node_type x, value_type i ) const {
		assert( x < +oo );

		auto pos= node2position(x);

		if ( (*S)[pos] == i ) {
			assert( t[i] );
			auto l= S->rank(pos,i);
			assert( B[i]->is_opening(l) );
			return t[i]->position2node(l);
		}

		auto l= S->rank(pos,i);
		if ( !l ) return +oo;

		if ( B[i]->is_opening(l-1) )
			return t[i]->position2node(l-1);
		l= B[i]->find_open(l-1);
		assert( B[i]->is_opening(l) );
		auto res= B[i]->enclose(l);
		if ( res == B[i]->size() )
			return +oo;
		return t[i]->position2node(res);
	}

	value_type
	_query( node_type x, node_type y, node_type z, const value_type wc, size_type k ) const {
		if ( this->is_homogeneous() ) {
			assert( k < this->size() );
			return this->a;
		}
		assert( this->a < this->b );
		value_type mid= (a+b)>>1;

		size_type accum= 0;
		for ( size_type i= 0; i < r; ++i ) {
			if ( !t[i] ) continue ;
			node_type ix,iy,iz;
			size_type dw,dx,dy,dz;
			if ( t[i]->size() >= 2 ) {
				ix= (x<+oo?image_of(x,i):x), iy= (y<+oo?image_of(y,i):y), iz= (z<+oo?image_of(z,i):z);
				dx= (ix<+oo?t[i]->depth(ix)+1:0), dy= (iy<+oo?t[i]->depth(iy)+1:0), dz= (iz<+oo?t[i]->depth(iz)+1:0);
				dw= dx+dy+(t[i]->a<=wc && wc<=t[i]->b?1:0)-2*dz;
			}
			else dw= (t[i]->a<=wc && wc<=t[i]->b?1:0);
			if ( accum+dw > k )
				return t[i]->_query(ix,iy,iz,wc,k-accum);
			accum+= dw;
		}
		std::cout << accum << " " << k << "\n";
		assert( false );
		return n;
	}

	void init( const std::string &s, const std::vector<value_type> &wgt ) {
		assert( size() == wgt.size() );

		for ( auto i= 0; i < r; t[i++]= NULL ) ;
		for ( auto i= 0; i < r; sub_tree[i++]= nullptr ) ;
		
		if ( is_homogeneous() ) return ;

		value_type mid= (a+b)>>1;

		size_type sizes[r], cur[r];
		memset(sizes,0,sizeof sizes), memset(cur,0,sizeof cur);
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] <= mid ) ++sizes[0];
			else ++sizes[1];

		for ( auto i= 0; i < r; ++i )
			sub_tree[i]= new bit_vector(2*sizes[i]);

		bit_vector *bv= new bit_vector(2*size());
		std::stack<size_type> st;
		size_type V= 0;
		for ( auto i= 0; i < 2*size(); ++i ) {
			if ( s[i] == '(' )
				st.push(V++);
			(*bv)[i]= wgt[st.top()]<=mid?0:1;
			if (wgt[st.top()] <= mid ) 
				(*(sub_tree[0]))[cur[0]++]= (s[i]=='('?1:0);
			else 
				(*(sub_tree[1]))[cur[1]++]= (s[i]=='('?1:0);
			if ( s[i] == ')' ) st.pop();
		}
		assert( st.empty() );
		S= new rs_bitvector01(*bv);
		for ( auto i= 0; i < r; ++i )
			B[i]= new t_bp_support(sub_tree[i]);

		std::stringstream str[r];
		std::vector<value_type> vec[r];
		for ( auto i= 0; i < r; vec[i++].clear() ) ;
		for ( auto i= 0; i < size(); ++i ) {
			size_type which_part= wgt[i]<=mid?0:1;
			vec[which_part].push_back(wgt[i]);
		}
		V= 0;
		for ( auto i= 0; i < 2*size(); ++i ) {
			if ( s[i] == '(' ) {
				st.push(V++);
				size_type which_part= wgt[st.top()]<=mid?0:1;
				str[which_part] << "(";
				continue ;
			}
			assert( !st.empty() );
			size_type which_part= wgt[st.top()]<=mid?0:1;
			st.pop();
			str[which_part] << ")";
		}
		assert( st.empty() );
		if ( !vec[0].empty() )
			t[0]= new tree_extraction_mp<>(str[0].str(),vec[0],a,mid,B[0]);
		else t[0]= NULL;
		if ( !vec[1].empty() )
			t[1]= new tree_extraction_mp<>(str[1].str(),vec[1],mid+1,b,B[1]);
		else t[1]= NULL;
	}

	tree_extraction_mp( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b, const t_bp_support *stru ) {
		this->a= a, this->b= b, n= s.length()/2;
		assert( !(s.length()&1) );
		assert( s.length()/2 == wgt.size() );
		structure= stru, init(s,wgt);
	}

	value_type original_weight( node_type x ) const {
		if ( is_homogeneous() ) 
			return a;
		size_type pos= node2position(x);
		value_type i= (value_type)(*S)[pos];
		assert( t[i] );
		node_type ix= image_of(x,i);
		assert( ix < +oo );
		return t[i]->original_weight(ix);
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

	tree_extraction_mp( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2;
		value_type mi= 0;
		assert( s.length()/2 == wgt.size() );
		value_type sigma= 0;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > sigma )
				sigma= wgt[i];
			else if ( wgt[i] < mi )
				mi= wgt[i];
		this->a= mi, this->b= sigma;

		base_structure= new bit_vector(2*size());
		for ( auto i= 0; i < 2*size(); ++i )
			(*base_structure)[i]= (s[i]=='('?1:0);
		structure= new t_bp_support(base_structure);

		init(s,wgt);
	}

	size_type size() const { return n; }

	value_type weight( const node_type x ) const {
		return original_weight(x);
	}

	value_type query( const node_type x, const node_type y ) const {
		node_type z= _lca(x,y);
		assert( z < +oo );
		size_type len= depth(x)+depth(y)+1-2*depth(z);
		return _query(x,y,z,original_weight(z),len>>1);
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x);
	}

	/*
	value_type a,b;
	size_type n;
	bool is_homogeneous() const { return a==b; }
	bit_vector *base_structure= nullptr, *sub_tree[r];
	const t_bp_support *structure;
	rs_bitvector01 *S= nullptr;
	t_bp_support *B[r];
	tree_extraction_mp *t[r];
	*/

	double bits_per_node() const {
		if ( !structure ) return 0.00;
		double ans= 8*(sdsl::size_in_bytes(*structure)+(is_homogeneous()?0:S->size_in_bytes())+sizeof(a)+sizeof(b));
		ans += 8*(sizeof(base_structure) + sizeof(S) + sizeof(n));
		if ( base_structure ) ans += 8*sdsl::size_in_bytes(*base_structure);
		for ( auto i= 0; i < r; ++i ) {
			ans += 8*sizeof(sub_tree[i]) + 8*sizeof(t[i]);
			if ( t[i] ) {
				ans+= t[i]->bits_per_node()*t[i]->size();
				if ( sub_tree[i] )
					ans += 8*sdsl::size_in_bytes(*(sub_tree[i]));
			}
		}
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		return _lca(x,y);
	}

};

