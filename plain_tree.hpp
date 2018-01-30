/*
 */
#include <cassert>
#include "sdsl/int_vector.hpp"
#include "raw_tree.hpp"
#include "succinct_tree.hpp"
#include "pq_types.hpp"
#include <vector>
#include <unordered_map>
#include <map>
#include <set>

class plain_tree {
public:
	typedef pq_types::size_type	size_type;
	typedef pq_types::node_type	node_type;
	typedef pq_types::value_type value_type;
	typedef sdsl::bit_vector bit_vector;
	typedef sdsl::int_vector<> int_vector;
private:
	std::map<node_type,std::set<node_type>> adj;
	bool str_repr= false;
	size_type num_of_edges= 0;
	sdsl::bit_vector bv;
	void dfs( node_type x, size_type &cur ) {
		bv[cur++]= 1;
		for ( auto y: adj[x] )
			dfs(y,cur);
		bv[cur++]= 0;
	}
public:
	plain_tree() { adj.clear(); }
	inline void add_arc( node_type x, node_type y ) {
		if ( !adj.count(x) )
			adj[x].clear();
		assert( !adj[x].count(y) );
		adj[x].insert(y), ++num_of_edges;
	}
	operator sdsl::bit_vector() {
		if ( !str_repr ) {
			bv= sdsl::bit_vector(2*(nedges()+1),0);
			size_type cur= 0;
			dfs(0,cur);
			assert( cur == 2*size() );
			str_repr= true ;
		}
		return bv;
	}
	inline size_type size() const {
		return adj.size();
	}
	inline size_type nedges() const {
		return num_of_edges;
	}
};

class hpd {
public:
	typedef pq_types::size_type	size_type;
	typedef pq_types::node_type	node_type;
	typedef pq_types::value_type value_type;
	typedef sdsl::bit_vector bit_vector;
	typedef sdsl::int_vector<> int_vector;
private:
	enum { NONE = (1<<22) };
	const succinct_tree *T;
	std::vector<node_type> chain;
	std::vector<size_type> len;
	mutable std::unordered_map<node_type,size_type> which_chain, pos_in_chain;
	std::unordered_map<node_type,size_type> card;
	std::unordered_map<node_type,node_type> best_son, parent;
	mutable std::unordered_map<size_type,size_type> head;
	size_type n, chain_id;

	size_type dfs( node_type x ) {
		if ( card.count(x) )
			return card[x];
		card[x]= 1;
		std::vector<node_type> children= T->children(x);
		for ( auto y: children ) {
			parent[y]= x, card[x]+= dfs(y);
			if ( !best_son.count(x) || card[y] > card[best_son[x]] )
				best_son[x]= y;
		}
		return card[x];
	}

	void hld( node_type x, bool new_chain= false ) {
		assert( card.count(x) );
		//TODO: this is ugly, needs FIX
		if ( new_chain ) {
			if ( chain_id == n )
				chain_id= 0;
			else ++chain_id;
		}
		chain.push_back(x), which_chain[x]= chain_id;
		if ( best_son.count(x) ) {
			assert( best_son[x] == x+1 );
			hld(best_son[x]);
		}
		std::vector<node_type> children= T->children(x);
		for ( auto y: children )
			if ( !best_son.count(x) || y != best_son[x] )
				hld(y,true);
	}

	node_type ref( node_type x ) const {
		return chain[head[which_chain[x]]];
	}

public:

	hpd( const succinct_tree *t ) { 
		n= (T= t)->size();
	};

	std::tuple<bit_vector, bit_vector, std::vector<node_type>> 
	operator()() {

		chain_id= n, chain.clear(), which_chain.clear(), pos_in_chain.clear(), len.reserve(n);
		best_son.clear(), dfs(0), hld(0,true);

		assert( chain.size() == n );
		assert( which_chain.size() == n );

		size_type i,j,k,ch,prev= chain_id+1;
		node_type x,y;
		for ( i= 0; i < n; ++i ) {
			if ( prev != (ch=which_chain[x=chain[i]]) ) 
				head[ch]= i, len[ch]= 0;
			pos_in_chain[x]= i, ++len[ch], prev= ch;
		}

		plain_tree pt{};
		for ( x= 1; x < n; ++x ) 
			pt.add_arc(ref(parent[x]),x);
		assert( pt.nedges()+1 == n );

		bit_vector B= bit_vector(2*n,0);
		for ( k= 0, x= 0; x < n; ++x ) {
			B[k++]= 1, ch= which_chain[x];
			if ( chain[head[ch]] == x ) 
				k+= len[ch];
		}

		assert( k == 2*n );

		return std::make_tuple(bit_vector(pt),B,chain);

	}
};

