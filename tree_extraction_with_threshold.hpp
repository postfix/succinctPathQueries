/*
 * the goal is to get rid of the dummy node
 */
#include <algorithm>
#include <cassert>
#include <sdsl/int_vector.hpp>
#include <sdsl/bit_vectors.hpp>
#include "path_query_processor.hpp"
#include "pq_types.hpp"
#include "succinct_tree.hpp"
#include "rs_bitvector01.hpp"
#include "bp_tree.hpp"
#include "sdsl/bp_support_gg.hpp"
#include <cstring>
#include <sstream>
#include <cstdio>
#define oo (0xfffffffful)
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;
typedef sdsl::bp_support_gg<> bp_support_gg;

extern uint32_t K;

template<class t_bp_support= bp_support_gg,uint8_t r=2>
class tree_extraction_with_threshold: public path_query_processor {
private:
	//static const uint32_t K= 0x400;
	bool is_final_tree;
	bit_vector *base_structure= nullptr;
	value_type a,b;
	rs_bitvector01 *B;
	bool is_homogeneous() const { return a==b; }
	size_type n;
	value_type *weights= nullptr;

	t_bp_support *structure= nullptr;
	tree_extraction_with_threshold *t[r];

	/*! The position of the opening parenthesis of the node $x$
	 * \param x the pre-order number of the node
	 */
	inline size_type node2position( const node_type x ) const {
		return structure->select(x+1);
	}
	/*! The pre-order number of the node, whose opening parenthesis is at position i
	 * \param i position of the opening parenthesis
	 */
	inline node_type position2node( const size_type i ) const {
		return structure->rank(i)-1;
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

	inline node_type parent( node_type x ) const { 
		size_type pos= node2position(x);
		auto res= structure->enclose(pos);
		if ( res == structure->size() )
			return +oo;
		return position2node(res);
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
			//position2node(structure->double_enclose(ix.first,iy.first));
			auto res= structure->double_enclose(ix.first,iy.first);
			if ( res == structure->size() )
				return +oo;
			return position2node(res);
		}
		swap(ix,iy);
		goto label01;
		assert( false );
	}

	node_type image_of( node_type x, value_type i ) const {
		if ( (*B)[x] == i )
			return B->rank(x,i);
		size_type howManyPrecede= B->rank(x,i);
		if ( !howManyPrecede )
			return +oo;
		node_type u,v,z;
		u= B->select(howManyPrecede,i);

		if ( is_ancestor(u,x) )
			return B->rank(u,i);
		v= _lca(u,x);
		if ( v < +oo && (*B)[v] == i || v == +oo )
			return v<+oo?B->rank(v,i):v;
		howManyPrecede= B->rank(v,i);
		z= B->select(howManyPrecede+1,i);
		node_type image_of_z= B->rank(z,i), pz;
		pz= image_of_z==0?+oo:t[i]->parent(image_of_z);
		return pz;
	}

	value_type
	_query( node_type x, node_type y, node_type z, const value_type wc, size_type k ) const {
		if ( this->is_homogeneous() ) {
			assert( k < this->size() );
			return this->a;
		}
		assert( this->a < this->b );

		if ( is_final_tree ) {
			std::vector<value_type> path;
			assert( weights );
			path.clear();
			for ( ;x < +oo && x != z; path.push_back(weights[x]), x= parent(x) ) ;
			for ( ;y < +oo && y != z; path.push_back(weights[y]), y= parent(y) ) ;
			if ( z < +oo && wc == weights[z] )
				path.push_back(wc);
			std::sort(path.begin(),path.end());
			return path[k];
		}

		value_type mid= (a+b)>>1;

		size_type accum= 0;
		for ( size_type i= 0; i < r; ++i ) {
			if ( !t[i] ) continue ;
			node_type ix,iy,iz;
			size_type dw,dx,dy,dz;
			if ( t[i]->size() >= 1 ) {
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
		B= NULL;

		for ( auto i= 0; i < r; t[i++]= NULL ) ;
		
		if ( is_homogeneous() ) { is_final_tree= true; return ; }

		is_final_tree= false;

		size_type dt= 0;
		for ( auto x= 0; x < size(); ++x ) {
			auto dx= depth(x);
			if ( dt < dx ) dt= x;
		}
		if ( dt <= K ) {
			is_final_tree= true ;
			weights= new value_type[size()];
			for ( auto x= 0; x < size(); ++x )
				weights[x]= wgt[x];
			return ;
		}

		value_type mid= (a+b)>>1;

		bit_vector *bv = new bit_vector(size());
		for ( auto l= 0; l < size(); ++l ) {
			assert( a <= wgt[l] && wgt[l] <= b );
			(*bv)[l]= (wgt[l]<=mid?0:1);
		}
		B= new rs_bitvector01(*bv);

		std::stringstream str[r];
		std::vector<value_type> vec[r];
		for ( auto i= 0; i < r; vec[i++].clear() ) ;
		for ( auto i= 0; i < size(); ++i ) {
			size_type which_part= wgt[i]<=mid?0:1;
			vec[which_part].push_back(wgt[i]);
		}
		std::stack<size_type> st;
		size_type V= 0;
		//st.push(V++);
		//for ( auto i= 0; i < r; str[i++] << "(" ) ;
		for ( size_type i= 0; i < 2*size(); ++i ) {
			if ( s[i] =='(' ) {
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
		//st.pop();
		//for ( auto i= 0; i < r; str[i++] << ")" ) ;
		//assert( V == size() );
		assert( st.empty() );
		if ( !vec[0].empty() )
			t[0]= new tree_extraction_with_threshold<>(str[0].str(),vec[0],a,mid);
		else t[0]= NULL;
		if ( !vec[1].empty() )
			t[1]= new tree_extraction_with_threshold<>(str[1].str(),vec[1],mid+1,b);
		else t[1]= NULL;
	}

	tree_extraction_with_threshold( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b ) {
		this->a= a, this->b= b, n= s.length()/2;
		assert( !(s.length()&1) );
		assert( s.length()/2 == wgt.size() );

		base_structure= new bit_vector(2*size());
		for ( auto i= 0; i < 2*size(); ++i )
			(*base_structure)[i]= (s[i]=='('?1:0);
		structure= new t_bp_support(base_structure);

		init(s,wgt);
	}

	value_type original_weight( node_type x ) const {
		if ( is_homogeneous() ) 
			return a;
		if ( is_final_tree ) {
			assert( weights );
			return weights[x];
		}
		assert( x < +oo );
		value_type i= (value_type)(*B)[x];
		assert( t[i] );
		node_type ix= image_of(x,i);
		assert( ix < +oo );
		return t[i]->original_weight(ix);
		/*
		for ( auto i= 0; i < r; ++i )
			if ( anc(x,i) == x ) {
				assert( t[i] );
				return t[i]->original_weight(image_of(x,i));
			}
		assert( false );
		return b+1;
		*/
	}

public:

	tree_extraction_with_threshold( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2;
		assert( s.length()/2 == wgt.size() );
		value_type sigma= 0;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > sigma )
				sigma= wgt[i];
		this->a= 0, this->b= sigma;

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
		auto res= _query(x,y,z,original_weight(z),len>>1);
		return res;
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x);
	}

	double bits_per_node() const {
		if ( !structure ) return 0.00;

		double ans= 8*( sdsl::size_in_bytes(*structure) + (is_final_tree?0:B->size_in_bytes()) );
		ans+= 8*(sizeof(a)+sizeof(b)+sizeof base_structure+sizeof structure+sizeof B);
		ans+= 8*(sizeof weights);
		ans+= 8*(sizeof is_final_tree);
		if ( is_final_tree ) 
			ans+= 8*size()*sizeof *weights;
		if ( base_structure ) ans+= 8*sdsl::size_in_bytes(*base_structure);
		for ( auto i= 0; i < r; ++i ) {
			ans+= 8*sizeof t[i];
			if ( t[i] )
				ans+= t[i]->bits_per_node()*t[i]->size();
		}
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		return _lca(x,y);
	}
};

