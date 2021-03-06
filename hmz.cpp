/*
 * this tree uses succinctly represented tree during the
 * construction of the explicit 0/1 extraction;
 * O(1) LCA is constructed (Bender & Farach-Colton-type)
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
#include <cstring>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include "bender_farach_colton.hpp"
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;

template<class t_succinct_tree= bp_tree,uint8_t r=2>
class hmz: public path_query_processor {
private:
	const node_type dummy_root= 0, LEFT= 0, RIGHT= 0;
	node_type *explict_parent;
	value_type a,b;
	rs_bitvector01 *B;
	bool is_homogeneous() const { return a==b; }
	size_type n;
	mutable size_type *level;

	t_succinct_tree *T;
	hmz *t[r];
	node_type *ptr[2];

	inline size_type depth( node_type x ) const  { return level[x];  }
	node_type parent( node_type x ) const { assert( T ); return T->parent(x); }

	/* 
	 * find the lowest ancestor of "x" with weight from [a_i..b_i]
	 */
	node_type anc( node_type x, value_type i ) const {
		if ( x == dummy_root ) return x;
		if ( (*B)[x] == i )
			return x;
		//size_type howManyPrecede= B->rank(x,i)+1-B->rank(1,i);
		size_type howManyPrecede= B->rank(x,i)+1-(i==0?1:0);
		//assert( howManyPrecede );
		if ( !howManyPrecede )
			return size();
		node_type u,v,z;
		if ( i > 0 )
			u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede-1,i));
		else u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede,i));

		//assert( u == dummy_root || (*B)[u] == i );

		if ( T->is_ancestor(u,x) )
			return u;
		v= T->lca(u,x);
		if ( v == dummy_root || (*B)[v] == i )
			return v;
		//howManyPrecede= B->rank(v,i)+1-B->rank(1,i);
		howManyPrecede= B->rank(v,i)+1-(i==0?1:0);
		if ( i > 0 )
			z= B->select(howManyPrecede,i);
		else z= B->select(howManyPrecede+1,i);
		//assert( T->is_ancestor(v,z) );
		//assert( z != dummy_root && (*B)[z] == i );
		node_type image_of_z= B->rank(z,i)+1-B->rank(1,i), pz;
		pz= image_of_z==dummy_root?dummy_root:t[i]->explict_parent[image_of_z];
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
		/*node_type px= anc(x,i), ppx= anc_nv(x,i);
		assert( px == dummy_root || (*B)[px] == i );
		assert( px == ppx );
		*/
		node_type px= anc(x,i);
		if ( px == dummy_root )
			return dummy_root;
		return B->rank(px,i)+1-B->rank(1,i);
	}

	inline node_type get_ptr( node_type x, value_type i ) const {
		return ptr[i][x];
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
			if ( !t[i] ) continue ;
			size_type dw;
			node_type ix,iy,iz;
			if ( t[i]->size() >= 2 ) {
				ix= ptr[i][x], iy= ptr[i][y], iz= ptr[i][z];
				dw= t[i]->level[ix]+t[i]->level[iy]+(t[i]->a<=wc && wc<=t[i]->b?1:0)-2*t[i]->level[iz];
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

	void init( const std::string &s, const std::vector<value_type> &wgt, bool delete_T= true ) {
		assert( size() == wgt.size()+1 );
		B= NULL, T= new t_succinct_tree("("+s+")"); //holds the overall structure
		for ( auto i= 0; i < r; t[i++]= NULL ) ;
		

		value_type mid= (a+b)>>1;
		level= new size_type[size()];
		for ( auto l=0; l < size(); ++l )
			level[l]= T->depth((node_type)l);
		for ( auto l= 0; l < r; ++l )
			ptr[l]= new node_type[size()];

		explict_parent= new node_type[size()];
		for ( node_type x= 1; x < size(); ++x )
			explict_parent[x]= T->parent(x);

		if ( is_homogeneous() ) {
			if ( delete_T )
				delete T, T= NULL;
			return ;
		}

		bit_vector *bv = new bit_vector(size());
		for ( auto l= 0; l < size()-1; ++l ) {
			assert( a <= wgt[l] && wgt[l] <= b );
			(*bv)[l+1]= (wgt[l]<=mid?0:1);
		}
		B= new rs_bitvector01(*bv);

		std::stringstream str[r];
		std::vector<value_type> vec[r];
		for ( auto i= 0; i < r; vec[i++].clear() ) ;
		for ( auto i= 0; i < size()-1; ++i ) {
			size_type which_part= wgt[i]<=mid?0:1;
			vec[which_part].push_back(wgt[i]);
		}
		std::stack<size_type> st;
		size_type V= 0;
		for ( size_type i= 0; i < 2*(size()-1); ++i ) {
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
		assert( st.empty() );
		if ( !vec[0].empty() )
			t[0]= new hmz<>(str[0].str(),vec[0],a,mid);
		else t[0]= NULL;
		if ( !vec[1].empty() )
			t[1]= new hmz<>(str[1].str(),vec[1],mid+1,b);
		else t[1]= NULL;

		for ( auto v= 0; v < size(); ++v ) 
			for ( auto l= 0; l < r; ++l )
				ptr[l][v]= image_of(v,l);
		if ( delete_T )
			delete T, T= NULL;
		delete B;
	}

	hmz( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b ) {
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

	lca_processor *lca_proc= NULL;

	value_type *oweight= NULL;

public:

	hmz( const std::string &s, const std::vector<value_type> &wgt ) {
		n= s.length()/2+1;
		assert( s.length()/2 == wgt.size() );
		value_type sigma= 0;
		for ( auto i= 0; i < wgt.size(); ++i )
			if ( wgt[i] > sigma )
				sigma= wgt[i];
		this->a= 0, this->b= sigma;
		init(s,wgt,false);
		lca_proc= new lca_processor(T);
		oweight= new value_type[size()];
		for ( auto i= 0; i < wgt.size(); ++i )
			oweight[i+1]= wgt[i];
	}

	size_type size() const { return n; }

	value_type weight( const node_type x ) const {
		return original_weight(x+1);
	}

	value_type query( const node_type x, const node_type y ) const {
		//node_type z= T->lca(x+1,y+1);
		node_type z= (*lca_proc)(x+1,y+1);
		size_type len= depth(x+1)+depth(y+1)+1-2*depth(z);
		//return _query(x+1,y+1,z,original_weight(z),len>>1);
		return _query(x+1,y+1,z,oweight[z],len>>1);
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x+1);
	}

	double bits_per_node() const {
		if ( !T ) return 0.00;
		//double ans= 8*( T->size_in_bytes() + (is_homogeneous()?0:B->size_in_bytes()) );
		double ans= 0;
		for ( auto i= 0; i < r; ++i ) {
			ans+= 8*size()*sizeof *ptr[i];
			if ( t[i] )
				ans+= t[i]->bits_per_node()*t[i]->size();
		}
		ans+= 8*size()*sizeof *level;
		ans+= 8*sizeof a + 8*sizeof b;
		ans+= 8*size()*sizeof *explict_parent;
		if ( lca_proc )
			ans+= 8*T->size_in_bytes(), ans+= 8*lca_proc->size_in_bytes(), ans+= 8*size()*sizeof *oweight;
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		//return T->lca(x+1,y+1)-1;
		return (*lca_proc)(x+1,y+1)-1;
	}

	/*
	~hmz() {
		if ( is_homogeneous() ) return ;
		for ( auto i= 0; i < r; ++i )
			if ( t[i] )
				delete t[i];
		delete lca_proc;
		delete explict_parent;
	}
	*/
};
