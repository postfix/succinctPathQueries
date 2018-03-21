#ifndef NAIVE_TREE_LCA_INCLUDED
#define NAIVE_TREE_LCA_INCLUDED

#include "succinct_tree.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stack>
#include <vector>
#include <random>
#include "bender_farach_colton.hpp"
#define M (1<<MLOG)
#define MLOG (30)


class naive_tree_lca: public succinct_tree {
public:
	typedef pq_types::size_type size_type;
	typedef pq_types::node_type node_type;
	typedef pq_types::value_type value_type;

private:

	enum { NONE= -1, ROOT= 0 };
	int n,p[M],d[M];
	std::vector<pq_types::node_type> adj[M];

	/*
	Randomized_Partition(A, l, u)
	i) Choose a random index i from [l..u]. Assign p ← A[i].
	ii) Rearrange the elements of A[l..u] so that there exists an integer j s.t.
	    elements in A[l..j-1] are <= p, 
			     A[j] = p (i.e. p is moved from A[i] to A[j]), and
						     elements in A[j+1..u] are >= p.
				 iii) return j
	*/
	/*
	 Randomized_Select(A, l, u, i) // find the ith smallest element from A[l..u]
	 if l = r then
	   return A[l]

	   p ← Randomized_Partition(A, l, u)
	   k ← p - l + 1 // note that a[p] is the k-th smallest element in A[l..u]

	   if k = i then
	     return A[p]
		 else if i < k then
		        return Randomized_Select(A, l, p-1, i)
				     else
					        return Randomized_Select(A, p+1, u, i-k)
	*/
inline size_type next_int( size_type l, size_type r ) {
	size_type res= rand()%(r-l+1)+l;
	assert( l <= res );
	assert( res <= r );
	return res;
}

size_type random_partition( std::vector<value_type> &arr, size_type l, size_type r )
{
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

value_type randomized_select( std::vector<value_type> &A, size_type l, size_type r, size_type i ) {
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
	lca_processor *prc;

public:

	naive_tree_lca( const std::string &s ) {
		//is >> n;
		int i,j,k,x,y,V= 0;
		//for ( i = 0; i < n; ++i ) is >> k; //read weights, no need to store them
		assert( !(s.size() & 1) );
		n= s.size()/2;
		for ( i= 0; i < n; adj[i++].clear() ) ;
		std::stack<int> st;
		assert( st.empty() );
		for ( i= 0; i < s.size(); ++i ) {
			if ( s[i] == '(' ) {
				if ( !st.empty() ) {
					adj[st.top()].push_back(V);
					p[V]= st.top();
				}
				d[V]= st.size();
				st.push(V++);
				continue ;
			}
			st.pop();
		}
		assert( V == n );
		prc= new lca_processor(this);
	}

	// tree info
	size_type size() const {
		return n;
	}

	// navigation
	node_type parent( const node_type x ) const {
		if ( x == ROOT ) return ROOT;
		return p[x];
	}

	// FIXME
	node_type ancestor( const node_type x, const size_type i ) const {
		return ROOT;
	}

	std::vector<node_type> children( const node_type x ) const {
		return adj[x];
	}

	node_type lca( const node_type cx, const node_type cy ) const {
		return (*prc)(cx,cy);
		/*
		node_type x= cx, y= cy;
		for ( ;d[x] > d[y]; x= parent(x) ) ;
		for ( ;d[x] < d[y]; y= parent(y) ) ;
		assert( d[x] == d[y] );
		for ( ;x != y; x= parent(x), y= parent(y) ) ;
		return x;
		*/
	}

	//FIXME
	size_type depth( const node_type x ) const {
		return d[x];
	}

	// predicates
	// FIXME
	bool is_ancestor( const node_type p, const node_type x ) const {
		//return in[p] <= in[x] && out[x] <= out[p];
		return false;
	}
	bool is_leaf( const node_type x ) const {
		return adj[x].size() == 0;
	}

	value_type query( node_type x, node_type y, std::vector<value_type> &w ) {
		std::vector<value_type> path;
		node_type z = lca(x,y);
		for ( auto cur = x; cur != z; path.push_back(w[cur]), cur = parent(cur) );
		for ( auto cur = y; cur != z; path.push_back(w[cur]), cur = parent(cur) );
		path.push_back(w[z]);
		//std::sort(path.begin(),path.end());
		/*
		for ( auto p: path )
			std::cout << p << " ";
		std::cout << std::endl;*/
		//return path[path.size()/2];
		return randomized_select(path,0,path.size()-1,path.size()/2);
	}

	//FIXME
	double size_in_bytes() const {
		return 0;
	}
};

#endif
