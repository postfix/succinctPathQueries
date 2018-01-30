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
using value_type = pq_types::value_type;
using size_type = pq_types::size_type;
using node_type = pq_types::node_type;

class hpd_remapper {
private:
	enum { ROOT = 0 };
	std::unordered_map<node_type,std::vector<node_type>> adj;
	size_type n,chain_id;
	std::unordered_map<node_type,node_type> p, best_son;
	std::unordered_map<node_type,size_type> card, which_chain;
	std::unordered_map<size_type,size_type> head;
	std::vector<node_type> chain;
	std::vector<size_type> len;
	size_type dfs( node_type x ) {
		assert( !card.count(x) );
		card[x]= 1;
		auto &c = card[x];
		for ( auto y: adj[x] ) {
			c += dfs(y);
			if ( !best_son.count(x) || card[y] > card[best_son[x]] )
				best_son[x] = y;
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
		chain.push_back(x), which_chain[x]= chain_id, ++len[chain_id];
		if ( best_son.count(x) )
			hld(best_son[x]);
		for ( auto y: adj[x] )
			if ( !best_son.count(x) || y != best_son[x] )
				hld(y,true);
	}

	void dfs( node_type x, std::stringstream &ss ) {
		ss << "(";
		if ( best_son.count(x) )
			dfs(best_son[x],ss);
		for ( auto y: adj[x] )
			if ( !best_son.count(x) || y != best_son[x] )
				dfs(y,ss);
		ss << ")";
	}

public:

	hpd_remapper() {}

	std::tuple<std::string,std::vector<value_type>>
	convert( const std::string &s, const std::vector<value_type> &w ) {
		size_type i,j,k,V = 0;
		assert( !(s.size() & 1) );
		n= s.size()/2;
		adj.clear(), p.clear(), card.clear(), best_son.clear();
		chain.clear(), which_chain.clear(), len.reserve(n), head.clear();
		std::stack<int> st;
		assert( st.empty() );
		for ( i = 0; i < s.size(); ++i ) {
			if ( s[i] == '(' ) {
				if ( !st.empty() ) {
					if ( !adj.count(st.top()) )
						adj[st.top()].clear();
					adj[st.top()].push_back(V);
				}
				st.push(V++);
				continue ;
			}
			st.pop();
		}
		assert( V == n );
		dfs(ROOT);
		chain_id= n, hld(ROOT,true);
		std::vector<value_type> weights(w.size());
		for ( i = 0; i < w.size(); ++i )
			weights[i] = w[chain[i]];
		std::stringstream ss;
		dfs(ROOT,ss);
		return std::make_tuple<>(ss.str(),weights);
	}
};
