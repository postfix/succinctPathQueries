/*
 */
#include <algorithm>
#include <cassert>
#include <sdsl/wavelet_trees.hpp>
#include <sdsl/int_vector.hpp>
#include <sdsl/wt_int.hpp>
#include "path_query_processor.hpp"
#include "pq_types.hpp"
#include "succinct_tree.hpp"
#include "rs_bitvector.hpp"
#include "bp_tree.hpp"
#include <cstring>
#include <sstream>
#include <cstdio>
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;

template<size_type r= 2,class t_succinct_tree= bp_tree,class wavelet_tree= wt_int>
class tree_extraction: public path_query_processor {
private:
	const node_type dummy_root= 0;
	value_type a,b;
	wavelet_tree *B;
	bool is_homogeneous() const { return a==b; }
	size_type n;

	t_succinct_tree *T;
	tree_extraction *t[r];

	size_type depth( node_type x ) const  { return T->depth(x);  }
	node_type parent( node_type x ) const { return T->parent(x); }

	/* 
	 * find the lowest ancestor of "x" with weight from [a_i..b_i]
	 */

	node_type anc( node_type x, value_type i ) const {
		if ( x == dummy_root ) return x;
		if ( (*B)[x] == i )
			return x;
		size_type howManyPrecede= B->rank(x,i)+1-B->rank(1,i);
		assert( howManyPrecede );
		if ( !howManyPrecede )
			return size();
		node_type u,v,z;
		if ( i > 0 )
			u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede-1,i));
		else u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede,i));

		assert( u == dummy_root || (*B)[u] == i );

		if ( T->is_ancestor(u,x) )
			return u;
		v= T->lca(u,x);
		if ( v == dummy_root || (*B)[v] == i )
			return v;
		howManyPrecede= B->rank(v,i)+1-B->rank(1,i);
		if ( i > 0 )
			z= B->select(howManyPrecede,i);
		else z= B->select(howManyPrecede+1,i);
		assert( T->is_ancestor(v,z) );
		assert( z != dummy_root && (*B)[z] == i );
		node_type image_of_z= B->rank(z,i)+1-B->rank(1,i), pz;
		pz= image_of_z==dummy_root?dummy_root:t[i]->parent(image_of_z);
		if ( pz == dummy_root )
			return pz;
		if ( i > 0 )
			return B->select(pz,i);
		else return B->select(pz+1,i);
	}

	node_type anc_nv( node_type x, value_type i ) const {
		if ( x == dummy_root )
			return x;
		if ( (*B)[x] == i )
			return x;
		for ( ;x != dummy_root && (*B)[x] != i; x= T->parent(x) ) ;
		return x;
	}

	node_type image_of( node_type x, value_type i ) const {
		if ( x == dummy_root ) return dummy_root;
		node_type px= anc(x,i), ppx= anc_nv(x,i);
		assert( px == dummy_root || (*B)[px] == i );
		//printf("%d %d\n",(int)px,(int)ppx);
		assert( px == ppx );
		if ( px == dummy_root )
			return dummy_root;
		return B->rank(px,i)+1-B->rank(1,i);
	}

	value_type
	_query( node_type x, node_type y, node_type z, const value_type wc, size_type k ) const {
		if ( this->is_homogeneous() ) {
			assert( k < this->size() );
			return this->a;
		}
		assert( this->a < this->b );
		value_type quot= (b-a)/r+1;
		assert( quot <= (b-a)/(r-1) );
		assert( a+quot*r > b );
		assert( a+quot*(r-1)+quot-1 >= b );
		assert( a+quot*(r-1) <= b );

		assert( T->depth(x)+T->depth(y)+1-2*T->depth(z) > k );

		/*
		if ( size() <= 1000 ) {
			std::vector<value_type> path(depth(x)+depth(y)+1-2*depth(z));
			size_type l= 0;
			for ( auto cx= x; cx != z; path[l++]= original_weight(cx), cx= parent(cx) );
			for ( auto cy= y; cy != z; path[l++]= original_weight(cy), cy= parent(cy) );
			path[l++]= original_weight(z);
			std::sort(path.begin(),path.end());
			assert( l == path.size() );
			return path[k];
		}
		*/

		size_type accum= 0;
		for ( size_type i= 0; i < r; ++i ) {
			size_type dw;
			node_type ix= image_of(x,i), iy= image_of(y,i), iz= image_of(z,i);
			//printf("ix = %d, iy = %d, iz = %d\n",(int)ix,(int)iy,(int)iz);
			//printf("x = %d, y = %d, z = %d\n",(int)x,(int)y,(int)z);
			dw= t[i]->depth(ix)+t[i]->depth(iy)+(t[i]->a<=wc && wc<=t[i]->b?1:0)-2*t[i]->depth(iz);
			//printf("%d %d %d\n",(int)t[i]->depth(ix),(int)t[i]->depth(iy),(int)t[i]->depth(iz));
			if ( accum+dw > k )
				return t[i]->_query(ix,iy,iz,wc,k-accum);
			//std:: cout << "dw = " << dw << "\n";
			accum+= dw;
		}
		std::cout << accum << " " << k << "\n";
		assert( false );
		/*
		std::cout << x << " " << y << " " << z << "\n";
		printf("a= %d, b= %d\n",(int)a,(int)b);
		*/
		return n;
	}

	void init( const std::string &s, const std::vector<value_type> &wgt ) {
		//printf("[%d,%d]\n",(int)a,(int)b);
		assert( size() == wgt.size()+1 );
		T= new t_succinct_tree("("+s+")"); //holds the overall structure
		if ( is_homogeneous() ) return ;

		for ( auto i= 0; i < r; t[i++]= NULL ) ;

		value_type quot= (b-a)/r+1;
		assert( quot <= (b-a)/(r-1) );
		assert( a+quot*r > b );
		assert( a+quot*(r-1)+quot-1 >= b );
		assert( a+quot*(r-1) <= b );

		B= new wavelet_tree();
		int_vector w(size(),0,32);
		for ( size_type good, bad, mid, l= 0; l < size()-1; w[++l]= good ) {
			for ( good= 0, bad= r, mid; good+1 < bad; ) {
				mid= (good+bad)>>1;
				if ( a+mid*quot <= wgt[l] )
					good= mid;
				else bad= mid;
			}
			assert( a+quot*good <= wgt[l] );
		    assert( wgt[l] <= a+quot*good+quot-1 );
		}
		construct_im(*B,w);

		std::stringstream str[r];
		std::vector<value_type> vec[r];
		for ( auto i= 0; i < r; vec[i++].clear() ) ;
		/*
		for ( auto i= 0; i < r; ++i )
			vec[i].push_back(a+i*quot);
		*/
		for ( auto i= 0; i < size()-1; ++i ) {
			size_type which_part= (wgt[i]-a)/quot;
			if ( !(0 <= which_part && which_part < r) )
				std::cout << wgt[i] << " " << which_part << " " << r << " " << quot << " " << a << " " << b << "\n";
			assert( 0 <= which_part && which_part < r );
			vec[which_part].push_back(wgt[i]);
		}
		std::stack<size_type> st;
		size_type V= 0;
		//st.push(V++);
		//for ( auto i= 0; i < r; str[i++] << "(" ) ;
		for ( size_type i= 0; i < 2*(size()-1); ++i ) {
			if ( s[i] =='(' ) {
				st.push(V++);
				size_type which_part= (wgt[st.top()]-a)/quot;
				assert( 0 <= which_part && which_part < r );
				str[which_part] << "(";
				continue ;
			}
			assert( !st.empty() );
			size_type which_part= (wgt[st.top()]-a)/quot;
			assert( 0 <= which_part && which_part < r );
			st.pop();
			str[which_part] << ")";
		}
		//st.pop();
		//for ( auto i= 0; i < r; str[i++] << ")" ) ;
		//assert( V == size() );
		assert( st.empty() );
		for ( size_type i= 0; i < r; ++i ) {
			if ( !vec[i].empty() )
				t[i]= new tree_extraction<>(str[i].str(),vec[i],a+quot*i,i==r-1?b:a+quot*i+quot-1);
			else t[i]= NULL;
		}
	}

	tree_extraction( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b ) {
		this->a= a, this->b= b, n= s.length()/2+1;
		assert( !(s.length()&1) );
		assert( s.length()/2 == wgt.size() );
		init(s,wgt);
	}

	value_type original_weight( node_type x ) const {
		if ( is_homogeneous() ) 
			return a;
		for ( auto i= 0; i < r; ++i )
			if ( anc(x,i) == x ) {
				assert( t[i] );
				return t[i]->original_weight(image_of(x,i));
			}
		assert( false );
		return b+1;
	}

public:

	tree_extraction( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2+1;
		assert( s.length()/2 == wgt.size() );
		value_type sigma= 1;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > sigma )
				sigma= wgt[i];
		this->a= 1, this->b= sigma;
		init(s,wgt);
	}

	size_type size() const { return n; }

	value_type weight( const node_type x ) const {
		return original_weight(x+1);
	}

	value_type query( const node_type x, const node_type y ) const {
		node_type z= T->lca(x+1,y+1);
		size_type len= T->depth(x+1)+T->depth(y+1)+1-2*T->depth(z);
		return _query(x+1,y+1,z,original_weight(z),len>>1);
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x+1);
	}

	double bits_per_node() const {
		double ans= 8 * ( T->size_in_bytes() + (is_homogeneous()?0:sdsl::size_in_bytes(*B)) );
		for ( auto i= 0; i < r; ++i )
			if ( t[i] )
				ans += t[i]->bits_per_node()*t[i]->size();
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		return T->lca(x+1,y+1)-1;
	}
};

