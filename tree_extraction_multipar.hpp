/*
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
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;
typedef sdsl::bp_support_gg<> bp_support_gg;

template<class t_bp_support= bp_support_gg, uint8_t r=2>
class tree_extraction_multipar: public path_query_processor {

private:

	const node_type dummy_root= 0;
	value_type a,b;
	size_type n;
	bool is_homogeneous() const { return a==b; }

	bit_vector *base_structure= nullptr, *sub_tree[r];
	const t_bp_support *structure;
	rs_bitvector01 *S= nullptr;
	t_bp_support *B[r];

	tree_extraction_multipar *t[r];

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
		return position2node(structure->enclose(pos));
	}

	node_type anc_nv( node_type x, value_type i ) const {
		if ( x == dummy_root )
			return x;
		size_type pos= node2position(x);
		if ( (*S)[pos] == i )
			return x;
		for ( ;x != dummy_root && (*S)[pos] != i; x= parent(x), pos= node2position(x) ) ;
		return x;
	}

	node_type image_of_nv( node_type x, value_type i ) const {
		if ( x == dummy_root ) return dummy_root;
		//node_type px= anc(x,i), ppx= anc_nv(x,i);
		//assert( px == dummy_root || (*B)[px] == i );
		//assert( px == ppx );
		//printf("%d %d\n",(int)px,(int)ppx);
		node_type px= anc_nv(x,i);
		if ( px == dummy_root )
			return dummy_root;
		auto l= S->rank(node2position(px),i)-S->rank(1,i);
		return B[i]->rank(l+1)-1;
	}

	node_type image_of( node_type x, value_type i ) const {
		/*
		 * since in bp_support_gg.hpp "select(i)" returnsthe index of the i-th
		 * opening parenthesis, we add +1 because the nodes are 0-based
		 */
		if ( x == dummy_root )
			return dummy_root;
		auto pos= node2position(x);

		assert( pos > 0 );

		if ( (*S)[pos] == i ) {
			assert( t[i] );
			auto aa= S->rank(pos,i), bb= S->rank(1,i);
			assert( aa >= bb );
			auto l= aa-bb;
			assert( B[i]->is_opening(l+1) );
			return t[i]->position2node(l+1);
		}

		auto aa= S->rank(pos,i), bb= S->rank(1,i);
		assert( aa >= bb );
		auto l= aa-bb;

		/*
		 * NOTE: in the official documentation:
		 *! Returns the number of opening parentheses up to and _including_ index i.
         * \pre{ \f$ 0\leq i < size() \f$ }
		 * NOTE: since "including", we need to subtract 1
		 */
		if ( B[i]->is_opening(l) )
			return t[i]->position2node(l);
		l= B[i]->find_open(l);
		assert( B[i]->is_opening(l) );
		return t[i]->position2node(B[i]->enclose(l));
	}

	value_type
	_query( node_type x, node_type y, node_type z, const value_type wc, size_type k ) const {
		if ( this->is_homogeneous() ) {
			assert( k < this->size() );
			return this->a;
		}
		assert( this->a < this->b );
		value_type mid= (a+b)>>1;

		//assert( T->depth(x)+T->depth(y)+1-2*T->depth(z) > k );

		size_type accum= 0;
		for ( size_type i= 0; i < r; ++i ) {
			if ( !t[i] ) continue ;
			size_type dw;
			node_type ix,iy,iz;
			if ( t[i]->size() >= 2 ) {
				ix= image_of(x,i), iy= image_of(y,i), iz= image_of(z,i);
				//ix= image_of_nv(x,i), iy= image_of_nv(y,i), iz= image_of_nv(z,i);
				dw= t[i]->depth(ix)+t[i]->depth(iy)+(t[i]->a<=wc && wc<=t[i]->b?1:0)-2*t[i]->depth(iz);
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
		assert( size() == wgt.size()+1 );

		for ( auto i= 0; i < r; t[i++]= NULL ) ;
		
		if ( is_homogeneous() ) return ;

		value_type mid= (a+b)>>1;

		size_type sizes[r], cur[r];
		memset(sizes,0,sizeof sizes), memset(cur,0,sizeof cur);
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] <= mid ) ++sizes[0];
			else ++sizes[1];

		for ( auto i= 0; i < r; ++i )
			sub_tree[i]= new bit_vector(2*(sizes[i]+1));

		bit_vector *bv= new bit_vector(2*size());
		std::stack<size_type> st;
		size_type V= 0;
		for ( auto i= 0; i < r; ++i )
			(*(sub_tree[i]))[cur[i]++]= 1;
		for ( auto i= 0; i < 2*(size()-1); ++i ) {
			if ( s[i] == '(' )
				st.push(V++);
			(*bv)[i+1]= wgt[st.top()]<=mid?0:1;
			if (wgt[st.top()] <= mid ) 
				(*(sub_tree[0]))[cur[0]++]= (s[i]=='('?1:0);
			else 
				(*(sub_tree[1]))[cur[1]++]= (s[i]=='('?1:0);
			if ( s[i] == ')' ) st.pop();
		}
		for ( auto i= 0; i < r; ++i ) {
			(*(sub_tree[i]))[cur[i]++]= 0;
			assert( cur[i] == 2*(sizes[i]+1) );
		}
		assert( st.empty() );
		S= new rs_bitvector01(*bv);
		for ( auto i= 0; i < r; ++i )
			B[i]= new t_bp_support(sub_tree[i]);

		std::stringstream str[r];
		std::vector<value_type> vec[r];
		for ( auto i= 0; i < r; vec[i++].clear() ) ;
		for ( auto i= 0; i < size()-1; ++i ) {
			size_type which_part= wgt[i]<=mid?0:1;
			vec[which_part].push_back(wgt[i]);
		}
		V= 0;
		for ( auto i= 0; i < 2*(size()-1); ++i ) {
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
			t[0]= new tree_extraction_multipar<>(str[0].str(),vec[0],a,mid,B[0]);
		else t[0]= NULL;
		if ( !vec[1].empty() )
			t[1]= new tree_extraction_multipar<>(str[1].str(),vec[1],mid+1,b,B[1]);
		else t[1]= NULL;
	}

	tree_extraction_multipar( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b, const t_bp_support *stru ) {
		this->a= a, this->b= b, n= s.length()/2+1;
		assert( !(s.length()&1) );
		assert( s.length()/2 == wgt.size() );
		structure= stru, init(s,wgt);
	}

	//TODO: image_of_nv --> image_of
	value_type original_weight( node_type x ) const {
		if ( is_homogeneous() ) 
			return a;
		assert( x != dummy_root );
		size_type pos= node2position(x);
		assert( pos > 0 );
		value_type i= (value_type)(*S)[pos];
		assert( t[i] );
		return t[i]->original_weight(image_of(x,i));
	}

	bool is_ancestor( const node_type x, const node_type y ) const {
		auto ix = interval(x),
			 iy = interval(y);
		return ix.first <= iy.first && iy.second <= ix.second;
	}

	node_type _lca( const node_type x, const node_type y ) const {
		if ( is_ancestor(x,y) )
			return x;
		if ( is_ancestor(y,x) )
			return y;
		auto ix = interval(x),
			 iy = interval(y);
		assert( ix.second < iy.first || iy.second < ix.first );
		if ( ix.second < iy.first ) {
	label01:
			return position2node(structure->double_enclose(ix.first,iy.first));
		}
		swap(ix,iy);
		goto label01;
		assert( false );
	}

public:

	tree_extraction_multipar( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2+1;
		assert( s.length()/2 == wgt.size() );
		value_type sigma= 0;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > sigma )
				sigma= wgt[i];
		this->a= 0, this->b= sigma;

		base_structure= new bit_vector(2*size());
		(*base_structure)[0]= 1, (*base_structure)[2*size()-1]= 0;
		for ( auto i= 0; i < 2*(size()-1); ++i )
			(*base_structure)[i+1]= (s[i]=='('?1:0);
		structure= new t_bp_support(base_structure);

		init(s,wgt);
	}

	size_type size() const { return n; }

	value_type weight( const node_type x ) const {
		return original_weight(x+1);
	}

	value_type query( const node_type x, const node_type y ) const {
		node_type z= _lca(x+1,y+1);
		assert( z != dummy_root );
		size_type len= depth(x+1)+depth(y+1)+1-2*depth(z);
		return _query(x+1,y+1,z,original_weight(z),len>>1);
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x+1);
	}

	double bits_per_node() const {
		if ( !structure ) return 0.00;
		double ans= 8*(sdsl::size_in_bytes(*structure)+(is_homogeneous()?0:S->size_in_bytes())+2*r*sizeof(int)+2*sizeof(a));
		for ( auto i= 0; i < r; ++i )
			if ( t[i] ) {
				ans+= t[i]->bits_per_node()*t[i]->size();
				//ans+= 8*sdsl::size_in_bytes(*B[i]);
			}
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		return _lca(x+1,y+1)-1;
	}

};

