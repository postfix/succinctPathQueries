#ifndef MY_LCA
#define MY_LCA

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
#ifdef M
#undef M
#endif

using namespace pq_types;

class small_block_rmq {
	size_type BS, ***M= nullptr;
	void preprocess() {
		for ( auto u= 0; u < (1<<(BS-1)); ++u )
			for ( auto i= 0; i < BS; ++i )
				M[u][i][i]= i;
		for ( auto k= 2; k <= BS; ++k )
			for ( auto u= 0; u < (1<<(BS-1)); ++u ) 
				for ( auto i= 0, j=i+k-1; (j=i+k-1) < BS; ++i ) {
					M[u][i][j]= i;
					for ( auto curmin= 0, cursum= 0, t= i+1; t <= j; ++t ) {
						cursum+= ((u>>(t-1))&1)?1:-1;
						if ( cursum < curmin )
							curmin= cursum,  M[u][i][j]= t;
					}
				}
	}
public:
	small_block_rmq( size_type bs ) : BS{bs} {
		M= (size_type ***)malloc((1<<(BS-1))*sizeof *M);
		assert( M );
		for ( auto u= 0; u < (1<<(BS-1)); ++u ) {
			M[u]= (size_type **)malloc(BS*sizeof *M[u]);
			assert( M[u] );
			for ( auto i= 0; i < BS; ++i ) {
				M[u][i]= (size_type *)malloc(BS*sizeof *M[u][i]);
				assert( M[u][i] );
			}
		}
		preprocess();
	}
	~small_block_rmq() {
		if ( M ) {
			for ( auto u= 0; u < (1<<(BS-1)); free(M[u++]) ) 
				for ( auto i= 0; i < BS-1; free(M[u][i++]) ) ;
			free(M);
		}
	}
	size_type operator () ( const size_type u, size_type i, size_type j ) const {
		if ( i > j ) { auto t= i; i= j; j= t; } 
		assert( u < (1<<(BS-1)) );
		assert( j-i+1 <= BS );
		assert( i < BS && j < BS );
		return M[u][i][j];
	}
	double size_in_bytes() const {
		return (1<<(BS-1))*BS*BS*sizeof ***M;
	}
};

class sparse_table_rmq {
	const long long *A;
	size_type **M=NULL, n, K;
public:
	sparse_table_rmq( const long long *AA, size_type nn ) : A{AA}, n{nn} {
		for ( K= 0; (1<<K) <= n; ++K ) ;
		M= (size_type **)malloc(n*sizeof *M);
		assert( M );
		for ( size_type i= 0; i < n; ++i ) {
			M[i]= (size_type *)malloc(K*sizeof *M[i]);
			assert( M[i] );
		}
		for ( size_type i= 0; i < n; M[i][0]= i, ++i ) ;
		for ( size_type k= 1; k < K; ++k ) 
			for ( size_type i= 0; i+(1<<k) < n; ++i ) {
				assert( M[i][k-1] < n );
				assert( M[i+(1<<(k-1))][k-1] < n );
				M[i][k]= A[M[i][k-1]]<A[M[i+(1<<(k-1))][k-1]]?M[i][k-1]:M[i+(1<<(k-1))][k-1];
			}
	}
	~sparse_table_rmq() {
		if ( M ) {
			for ( size_type i= 0; i < n; free(M[i++]) ) ;
			free(M);
		}
	}
	size_type operator() ( const size_type i, const size_type j ) const {
		assert( i <= j );
		if ( i == j ) 
			return i;
		assert( i < j );
		size_type k,l,r;
		for ( k= 0; i+(1<<k) <= j; ++k ) ;
		assert( i+(1<<(k-1)) <= j );
		return A[l=M[i][k-1]] < A[r=M[j-(1<<(k-1))+1][k-1]] ? l:r;
	}
	double size_in_bytes() const {
		return n*K*sizeof **M;
	}
};

/*
 * ±1RMQ 
 */
class pm_one_rmq {
	const long long *L;
	long long *A;
	size_type n,K,BS,border;
	size_type *min_for_block= nullptr;
	size_type *B= nullptr;
	sparse_table_rmq *st= nullptr;
	small_block_rmq *smbl= nullptr;
	size_type *block_type;

	size_type _rmq( size_type i, size_type j ) const {
		if ( i > j ) { auto t= i; i= j; j= t; }
		size_type ans, bi= i/BS, bj= j/BS;

		if ( bi == bj ) 
			return bi*BS+(*smbl)(block_type[bi],i-bi*BS,j-bj*BS);

		auto idxi = bi*BS+(*smbl)(block_type[bi],i-bi*BS,BS-1),
			 idxj = bj*BS+(*smbl)(block_type[bj],0,j-bj*BS);
		ans= L[idxi]<L[idxj]?idxi:idxj;
		if ( bi+1 <= bj-1 ) {
			auto t= (*st)(bi+1,bj-1);
			if ( A[t] < L[ans] )
				ans= t*BS+B[t];
		}
		return ans;
	}

public:

	pm_one_rmq( const long long *a, size_type nn ) : L{a}, n{nn} {
		for ( K= 0; (1<<K) <= n; ++K ) ;
		BS= (K>>1);
		for ( border= n; n%BS; ++n ) ;
		A= new long long[n/BS+3];
		B= new size_type[n/BS+3];
		smbl= new small_block_rmq(BS);
		puts("Small block built");
		for ( auto t= 0; t < n/BS; ++t ) {
			A[t]= L[t*BS], B[t]= 0;
			for ( auto i= t*BS+1; i < t*BS+BS && i < border; ++i ) 
				if ( L[i] < A[t] )
					A[t]= L[i], B[t]= i-t*BS;
		}
		st= new sparse_table_rmq(A,n/BS);
		block_type= new size_type[n/BS+1];
		for ( auto t= 0; t < n/BS; ++t ) {
			size_type u= (1<<(BS-1))-1;
			for ( auto i= t*BS+1; i < t*BS+BS && i < border; ++i ) 
				if ( L[i] != L[i-1]+1 ) 
					u &= ~(1<<(i-t*BS-1));
			block_type[t]= u;
		}
	}
	size_type operator ()( const size_type i, const size_type j ) const {
		return _rmq(i,j);
	}
	~pm_one_rmq() {
		if ( min_for_block )
			delete min_for_block;
		if ( A ) delete A;
		if ( B ) delete B;
		if ( st ) delete st;
		if ( smbl ) delete smbl;
		delete block_type;
	}
	double size_in_bytes() const {
		return (n/BS+3)*(sizeof *A+sizeof *B+sizeof *block_type)+st->size_in_bytes()+smbl->size_in_bytes();
	}
};

class lca_processor {
	size_type *R;
	long long *L;
	node_type *E;
	const succinct_tree *T;
	size_type n;
	void euler_tour( node_type x, size_type d, size_type &cur ) {
		L[R[E[cur]= x]= cur]= d, ++cur;
		auto children= T->children(x);
		for ( auto y: children ) {
			euler_tour(y,1+d,cur);
			L[cur]= d, E[cur++]= x;
		}
	}
	const pm_one_rmq *rmq;
public:
	lca_processor( const succinct_tree *t ) : T(t), n(T->size()) {
		L= new long long[2*n];
		R= new size_type[2*n];
		E= new node_type[2*n];
		size_type cur= 0;
		euler_tour(0,0,cur);
		assert( cur == 2*n-1 );
		rmq= new pm_one_rmq(L,2*n-1);
	}
	~lca_processor() {
		delete L, delete R, delete E;
		delete rmq;
	}
	node_type operator ()( node_type x, node_type y ) const {
		return E[(*rmq)(R[x],R[y])];
	}
	double size_in_bytes() const {
		return 2*n*(sizeof *L+sizeof *R+sizeof *E)+rmq->size_in_bytes();
	}
};

#endif
