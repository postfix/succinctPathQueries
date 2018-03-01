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
#define  MAXLG (22)
#define  MAXN   (1<<MAXLG)
#define  infty  (MAXN)

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
	const succinct_tree *T;
	mutable size_type which_chain[MAXN],pos_in_chain[MAXN];
	size_type card[MAXN],len[MAXN];
	node_type best_son[MAXN],parent[MAXN],chain[MAXN];
	mutable size_type head[MAXN];
	size_type n, chain_id;

	size_type dfs( node_type x ) {
		if ( card[x] )
			return card[x];
		card[x]= 1;
		auto &c= card[x];
		assert( T );
		assert( T->size() > x );
		std::vector<node_type> children= T->children(x);
		for ( auto y: children ) {
			parent[y]= x, c+= dfs(y);
			if ( best_son[x] == infty || card[y] > card[best_son[x]] )
				best_son[x]= y;
		}
		return card[x];
	}

	void hld( node_type x, size_type &cur, bool new_chain= false ) {
		assert( card[x] );
		//TODO: this is ugly, needs FIX
		if ( new_chain ) {
			if ( chain_id == n )
				chain_id= 0;
			else ++chain_id;
		}
		chain[cur++]= x, which_chain[x]= chain_id;
		if ( best_son[x] < infty ) {
			assert( best_son[x] == x+1 );
			hld(best_son[x],cur);
		}
		std::vector<node_type> children= T->children(x);
		for ( auto y: children )
			if ( best_son[x] == infty || y != best_son[x] )
				hld(y,cur,true);
	}

	node_type ref( node_type x ) const {
		return chain[head[which_chain[x]]];
	}

public:

	hpd( const succinct_tree *t ) { 
		assert( t );
		n= (T= t)->size();
	};

	std::tuple<bit_vector, bit_vector, std::vector<node_type>> 
	operator()() {

		chain_id= n;
		size_type cur= 0;
		for ( auto x= 0; x < n; best_son[x]= infty, card[x++]= 0 ) ;
		dfs(0), hld(0,cur,true);
		assert( cur == n );

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
		std::vector<node_type> rchain(n);
		for ( auto i= 0; i < n; rchain[i]= chain[i], ++i ) ;

		return std::make_tuple(bit_vector(pt),B,rchain);
	}
};

