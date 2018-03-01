/*
 * given a BPS (i.e. a forest, but in our case it will be a tree),
 * produce another BPS, in which the preorder traversal corresponds to
 * first going along the heavy path
 */
#include <cassert>
#include <vector>
#include <unordered_map>
#include <sstream>
#include "pq_types.hpp"
#define MAXLOG (22)
#define MAXV (1<<MAXLOG)
using value_type = pq_types::value_type;
using size_type = pq_types::size_type;
using node_type = pq_types::node_type;

static int cnt= 0;

class hpd_remapper {
private:
	enum { ROOT = 0, oo = MAXV+1 };
	std::vector<node_type> adj[MAXV];
	size_type n,chain_id;
	node_type p[MAXV], best_son[MAXV];
	size_type card[MAXV], which_chain[MAXV];
	size_type head[MAXV],len[MAXV];
	node_type chain[MAXV];
	int cur;
	size_type dfs( node_type x ) {
		assert( 0 <= x && x < n );
		assert( card[x] == +oo );
		card[x]= 1;
		++cnt;
		auto &c= card[x];
		for ( size_type l= 0; l < adj[x].size(); ++l ) {
			node_type y= adj[x][l];
			c+= dfs(y);
			if ( best_son[x] == +oo || card[y] > card[best_son[x]] )
				best_son[x]= y;
		}
		return card[x];
	}
	void hld( node_type x, bool new_chain= false ) {
		if ( new_chain ) {
			if ( chain_id == n )
				chain_id = 0;
			else ++chain_id;
			len[chain_id] = 0;
		}
		chain[cur++]= x, which_chain[x]= chain_id, ++len[chain_id];
		if ( best_son[x] < +oo )
			hld(best_son[x]);
		for ( auto y: adj[x] )
			if ( best_son[x] == +oo || y != best_son[x] )
				hld(y,true);
	}

	void dfs( node_type x, char *ss, size_type &cur ) {
		ss[cur++]= '(';
		if ( best_son[x] < +oo )
			dfs(best_son[x],ss,cur);
		for ( auto y: adj[x] )
			if ( best_son[x] == +oo || y != best_son[x] )
				dfs(y,ss,cur);
		ss[cur++]= ')';
	}

public:

	hpd_remapper() {}

	std::tuple<std::string,std::vector<value_type>>
	convert( const std::string &s, const std::vector<value_type> &w ) {
		size_type i,j,k,V= 0;
		assert( !(s.size() & 1) );
		n= s.size()/2;
		for ( cur= 0, i= 0; i < n; ++i ) {
			adj[i].clear();
			card[i]= best_son[i]= +oo, len[i]= 0;
		}
		std::stack<int> st;
		assert( st.empty() );
		for ( i = 0; i < s.size(); ++i ) {
			if ( s[i] == '(' ) {
				if ( !st.empty() ) 
					adj[st.top()].push_back(V);
				st.push(V++);
				continue ;
			}
			assert( !st.empty() );
			st.pop();
		}
		assert( st.empty() );
		assert( V == n );
		dfs(ROOT);
		chain_id= n, hld(ROOT,true);
		std::vector<value_type> weights(w.size());
		for ( i = 0; i < w.size(); ++i )
			weights[i] = w[chain[i]];
		char *str= new char[2*n+1];
		str[2*n]= '\0';
		size_type cur= 0;
		dfs(ROOT,str,cur);
		assert( cur == 2*n );
		auto res = std::make_tuple<>(std::string(str),weights);
		delete str;
		return res;
	}
};

