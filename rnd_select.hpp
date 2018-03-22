#ifndef RND_SELECT
#define RND_SELECT
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include "pq_types.hpp"

using namespace pq_types;

template<typename t_value_type= value_type,
		 typename t_size_type = size_type>
class rnd_select {
private:
	static t_size_type 
	next_int( t_size_type l, t_size_type r )  {
		assert( r >= l );
		t_size_type res= (t_size_type)(rand()%(r-l+1)+l);
		assert( l <= res );
		assert( res <= r );
		return res;
	}
	static t_size_type 
	partition( t_value_type *A, t_size_type p, t_size_type r )  {
		auto x= A[r];
		bool flag= false;
		long long i= ((long long)p)-1;
		assert( r > p );
		for ( auto j= p; j <= r-1; ++j ) {
			if ( A[j] <= x ) {
				if ( !flag ) flag= true, i= p;
				else ++i;
				std::swap(A[i],A[j]);
			}
		}
		std::swap(A[i+1],A[r]);
		return i+1;
	}
	static t_size_type
	randomized_partition( t_value_type *A, t_size_type p, t_size_type r )  {
		auto i= next_int(p,r);
		assert( p <= i );
		assert( i <= r );
		std::swap(A[r],A[i]);
		return partition(A,p,r);
	}
	static t_size_type 
	partition( std::vector<t_value_type> &A, t_size_type p, t_size_type r )  {
		auto x= A[r];
		bool flag= false;
		long long i= ((long long)p)-1;
		assert( r > p );
		for ( auto j= p; j <= r-1; ++j ) {
			if ( A[j] <= x ) {
				if ( !flag ) flag= true, i= p;
				else ++i;
				std::swap(A[i],A[j]);
			}
		}
		std::swap(A[i+1],A[r]);
		return i+1;
	}
	static t_size_type
	randomized_partition( std::vector<t_value_type> &A, t_size_type p, t_size_type r )  {
		auto i= next_int(p,r);
		assert( p <= i );
		assert( i <= r );
		std::swap(A[r],A[i]);
		return partition(A,p,r);
	}

public:
	static t_value_type select( t_value_type *A, t_size_type p, t_size_type r, t_size_type i )  {
		assert( p <= r );
		assert( i < r-p+1 );
		if ( p == r )
			return A[p];
		auto q= randomized_partition(A,p,r);
		assert( q >= p );
		assert( q <= r );
		auto k= q-p+1;
		if ( i+1 == k )
			return A[q];
		else if ( i+1 < k )
			return select(A,p,q-1,i);
		else return select(A,q+1,r,i-k);
	}
	static t_value_type select( std::vector<t_value_type> &A, t_size_type p, t_size_type r, t_size_type i )  {
		assert( p <= r );
		assert( i < r-p+1 );
		if ( p == r )
			return A[p];
		auto q= randomized_partition(A,p,r);
		assert( q >= p );
		assert( q <= r );
		auto k= q-p+1;
		if ( i+1 == k )
			return A[q];
		else if ( i+1 < k )
			return select(A,p,q-1,i);
		else return select(A,q+1,r,i-k);
	}
	static t_value_type naive_select( t_value_type *A, t_size_type p, t_size_type r, t_size_type i )  {
		std::sort(A+p,A+r+1);
		return A[p+i];
	}
};

#endif
