#ifndef RAW_TREE_INCLUDED
#define RAW_TREE_INCLUDED

#include "succinct_tree.hpp"
#include <iostream>
#include <algorithm>
#include <cassert>
#include <stack>
#include <vector>
#define M (1<<MLOG)
#define MLOG (25)

class raw_tree: public succinct_tree {

private:

	enum { NONE= -1, ROOT = 0 };
	int n,p[M],d[M],in[M],out[M],tick,K;
	bool seen[M];
	std::vector<pq_types::node_type> adj[M];
	int anc[M][MLOG+1];

	int up( int x, unsigned int u ) const {
		int k;
		if ( u > d[x] )
			return ROOT;
		for ( k = 0; u; ++k, u >>= 1 )
			if ( u & 1 )
				x = anc[x][k];
		return x;
	}

	int _lca( int x, int y ) const {
		if ( d[x] > d[y] )
			return _lca(up(x,d[x]-d[y]),y);
		if ( d[x] < d[y] )
			return _lca(y,x);
		assert( d[x] == d[y] );
		if ( x == y )
			return x;
		for ( int k = K-1; k; --k ) {
			assert( anc[x][k] == anc[y][k] );
			if ( anc[x][k-1] != anc[y][k-1] )
				x = anc[x][k-1], y = anc[y][k-1];
		}
		return anc[x][0];
	}

	void dfs( int x ) {
		int i,y,k;
		assert( !seen[x] );
		if ( x == ROOT ) d[x] = 0;
		else d[x] = d[p[x]]+1;
		for ( in[x]= ++tick, seen[x]= true, i= 0; i< (int)adj[x].size(); ++i ) {
			y= adj[x][i];
			assert( !seen[y] );
			for ( p[y]= anc[y][0]= x, k= 1; anc[y][k-1]!= NONE; anc[y][k]= anc[anc[y][k-1]][k-1], ++k ) ;
			dfs(y);
		}
		out[x]= ++tick;
	}

public:

	typedef pq_types::size_type size_type;
	typedef pq_types::node_type node_type;
	typedef pq_types::value_type value_type;

	raw_tree( const std::string &s) {
		//is >> n;
		int i,j,k,x,y,V = 0;
		//for ( i = 0; i < n; ++i ) is >> k; //read weights, no need to store them
		assert( !(s.size() & 1) );
		n = s.size()/2;
		for ( K = 0; (1<<K) <= n; ++K ) ;
		for ( tick = -1, i = 0; i < n; ++i ) 
			for ( adj[i].clear(), seen[i] = false, k = 0; k < K; anc[i][k++]= NONE) ;
		std::stack<int> st;
		assert( st.empty() );
		for ( i = 0; i < s.size(); ++i ) {
			if ( s[i] == '(' ) {
				if ( !st.empty() ) 
					adj[st.top()].push_back(V);
				st.push(V++);
				continue ;
			}
			st.pop();
		}
		assert( V == n );
		dfs(ROOT);
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

	node_type ancestor( const node_type x, const size_type i ) const {
		return up(x,i);
	}

	std::vector<node_type> children( const node_type x ) const {
		return adj[x];
	}

	node_type lca( const node_type x, const node_type y ) const {
		return _lca(x,y);
	}

	size_type depth( const node_type x ) const {
		return d[x];
	}

	// predicates
	bool is_ancestor( const node_type p, const node_type x ) const {
		return in[p] <= in[x] && out[x] <= out[p];
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
		std::sort(path.begin(),path.end());
		/*
		for ( auto p: path )
			std::cout << p << " ";
		std::cout << std::endl;*/
		return path[path.size()/2];
	}
};

#endif
