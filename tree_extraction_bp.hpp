/*
 * The plan is to get rid of the T succinct_tree and replace it with bare bp_support instead
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
using namespace pq_types;
typedef sdsl::int_vector<> int_vector;
typedef sdsl::bp_support_gg<> bp_support_gg;

template<class t_bp_support= bp_support_gg,uint8_t r=2>
class tree_extraction_bp: public path_query_processor {
private:
	const node_type dummy_root= 0;
	bit_vector *base_structure= nullptr;
	value_type a,b;
	rs_bitvector01 *B;
	bool is_homogeneous() const { return a==b; }
	size_type n;

	t_bp_support *structure= nullptr;
	tree_extraction_bp *t[r];

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

		if ( is_ancestor(u,x) )
			return u;
		v= _lca(u,x);
		if ( v == dummy_root || (*B)[v] == i )
			return v;
		//howManyPrecede= B->rank(v,i)+1-B->rank(1,i);
		howManyPrecede= B->rank(v,i)+1-(i==0?1:0);
		if ( i > 0 )
			z= B->select(howManyPrecede,i);
		else z= B->select(howManyPrecede+1,i);
		assert( is_ancestor(v,z) );
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
		for ( ;x != dummy_root && (*B)[x] != i; x= parent(x) ) ;
		return x;
	}

	node_type image_of( node_type x, value_type i ) const {
		if ( x == dummy_root ) return x;
		if ( (*B)[x] == i )
			return B->rank(x,i)+1-B->rank(1,i);
		//size_type howManyPrecede= B->rank(x,i)+1-B->rank(1,i);
		size_type howManyPrecede= B->rank(x,i)+1-(i==0?1:0);
		assert( howManyPrecede );
		if ( !howManyPrecede )
			return size();
		node_type u,v,z;
		if ( i > 0 )
			u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede-1,i));
		else u= (howManyPrecede==1?dummy_root:B->select(howManyPrecede,i));

		//assert( u == dummy_root || (*B)[u] == i );

		if ( is_ancestor(u,x) )
			return u==dummy_root?dummy_root:B->rank(u,i)+1-B->rank(1,i);
		v= _lca(u,x);
		if ( v == dummy_root || (*B)[v] == i )
			return v==dummy_root?dummy_root:B->rank(v,i)+1-B->rank(1,i);
		//howManyPrecede= B->rank(v,i)+1-B->rank(1,i);
		howManyPrecede= B->rank(v,i)+1-(i==0?1:0);
		if ( i > 0 )
			z= B->select(howManyPrecede,i);
		else z= B->select(howManyPrecede+1,i);
		assert( is_ancestor(v,z) );
		assert( z != dummy_root && (*B)[z] == i );
		node_type image_of_z= B->rank(z,i)+1-B->rank(1,i), pz;
		pz= image_of_z==dummy_root?dummy_root:t[i]->parent(image_of_z);
		if ( pz == dummy_root )
			return pz;
		return pz;
	}

	node_type image_of_old( node_type x, value_type i ) const {
		if ( x == dummy_root ) return dummy_root;
		/*
		node_type px= anc(x,i), ppx= anc_nv(x,i);
		assert( px == dummy_root || (*B)[px] == i );
		assert( px == ppx );*/
		//printf("%d %d\n",(int)px,(int)ppx);
		node_type px= anc_nv(x,i);
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
				ix= image_of(x,i), iy= image_of(y,i), iz= image_of(z,i);
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
		B= NULL;

		for ( auto i= 0; i < r; t[i++]= NULL ) ;
		
		if ( is_homogeneous() ) return ;

		value_type mid= (a+b)>>1;

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
		//st.push(V++);
		//for ( auto i= 0; i < r; str[i++] << "(" ) ;
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
		//st.pop();
		//for ( auto i= 0; i < r; str[i++] << ")" ) ;
		//assert( V == size() );
		assert( st.empty() );
		if ( !vec[0].empty() )
			t[0]= new tree_extraction_bp<>(str[0].str(),vec[0],a,mid);
		else t[0]= NULL;
		if ( !vec[1].empty() )
			t[1]= new tree_extraction_bp<>(str[1].str(),vec[1],mid+1,b);
		else t[1]= NULL;
	}

	tree_extraction_bp( const std::string &s, const std::vector<value_type> &wgt, value_type a, value_type b ) {
		this->a= a, this->b= b, n= s.length()/2+1;
		assert( !(s.length()&1) );
		assert( s.length()/2 == wgt.size() );

		base_structure= new bit_vector(2*size());
		(*base_structure)[0]= 1, (*base_structure)[2*size()-1]= 0;
		for ( auto i= 0; i < 2*(size()-1); ++i )
			(*base_structure)[i+1]= (s[i]=='('?1:0);
		structure= new t_bp_support(base_structure);

		init(s,wgt);
	}

	value_type original_weight( node_type x ) const {
		if ( is_homogeneous() ) 
			return a;
		assert( x != dummy_root );
		value_type i= (value_type)(*B)[x];
		assert( t[i] );
		node_type ix= image_of(x,i);
		assert( ix != dummy_root );
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

	tree_extraction_bp( const std::string &s, const std::vector<value_type> &wgt ) {
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
		auto res= _query(x+1,y+1,z,original_weight(z),len>>1);
		return res;
	}

	value_type weight_of( const node_type x ) const {
		return original_weight(x+1);
	}

	double bits_per_node() const {
		if ( !structure ) return 0.00;
		double ans= 8 * ( sdsl::size_in_bytes(*structure) + (is_homogeneous()?0:B->size_in_bytes()) );
		ans+= 8*(sizeof(a) + sizeof(b));
		for ( auto i= 0; i < r; ++i )
			if ( t[i] )
				ans += t[i]->bits_per_node()*t[i]->size();
		return ans/size();
	}

	node_type lca( node_type x, node_type y ) const {
		return _lca(x+1,y+1)-1;
	}
};

