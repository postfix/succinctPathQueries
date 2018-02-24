#include <cstdio>
#include "succinct_tree.hpp"
#include "bp_tree.hpp"
#include "path_query_processor.hpp"
#include <cassert>
#include <cmath>
#include <vector>
#define MAXLOG (30)
#define MAXN (1<<MAXLOG)
#include <random>
using namespace pq_types;

class path_query_processor_impl: public path_query_processor {
private:
	succinct_tree *T;
	std::vector<value_type> weights;
	size_type m;
	value_type sigma;
	double log_sigma;

inline size_type next_int( size_type l, size_type r ) const {
	size_type res= rand()%(r-l+1)+l;
	assert( l <= res );
	assert( res <= r );
	return res;
}

size_type random_partition( std::vector<value_type> &arr, size_type l, size_type r ) const {
    srand(time(NULL));
  	size_type pivotIdx= next_int(l,r);
    size_type pivot= arr[pivotIdx];
    size_type i= l-1,j= l;
	std::swap(arr[pivotIdx],arr[r]); // move pivot element to the end
    pivotIdx= r;

    for ( ;j < r; ++j )
        if ( arr[j] <= pivot )
            ++i, std::swap(arr[i], arr[j]);
	std::swap(arr[i+1], arr[pivotIdx]);
	//check(arr,l,r,i+1,pivot);
    return i+1;
}

value_type randomized_select( std::vector<value_type> &A, size_type l, size_type r, size_type i ) const {
	assert( 0 <= i );
	assert( i < r-l+1 );
	if ( l == r )
		return A[l];
	auto p= random_partition(A,l,r);
	auto k= p-l+1;
	if ( k == i+1 )
		return A[p];
	if ( i+1 < k )
		return randomized_select(A,l,p-1,i);
	return randomized_select(A,p+1,r,i-k);
}

public:

	path_query_processor_impl( succinct_tree *s, std::vector<value_type> &w ) {
		this->T= s;
		this->m= s->size();
		weights= std::vector<value_type>(w);
		assert( weights.size() == s->size() );
		sigma= 0;
		for ( auto x: weights )
			sigma= std::max(sigma,x);
		log_sigma= log(sigma)/log(2.00);
	}

	size_type size() const { return m; }
	value_type weight( const node_type x ) const { return weights[x]; }
	value_type weight_of( const node_type x ) const { return weight(x); }

	value_type
	query( const node_type x, const node_type y ) const {
		auto z= T->lca(x,y);
		std::vector<value_type> path(T->depth(x) + T->depth(y) - 2*T->depth(z) + 1);
		auto cur= 0;
		for ( auto cx= x; cx != z; path[cur++]= weight_of(cx), cx= T->parent(cx) ) ;
		for ( auto cy= y; cy != z; path[cur++]= weight_of(cy), cy= T->parent(cy) ) ;
		path[cur++]= weight_of(z);
		assert( cur == path.size() );
		return randomized_select(path,0,path.size()-1,path.size()/2);
	}

	double bits_per_node() const {
		return 8*T->size_in_bytes()/(size()+.0) + log_sigma;
	}

};

