/*
 */
#include <bits/stdc++>
#define BIT(k) (1ULL<<(x))
#define TST(a,x) (a[(x)>>3] & BIT((x)&7))
#define CLR(a,x) (a[(x)>>3] &= ~BIT((x)&7))
#define SET(a,x) (a[(x)>>3] |= BIT((x)&7))
#define N (1<<22)

#include "sdsl/bit_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "raw_tree.hpp"
#include "succinct_tree.hpp"
#include "pq_types.hpp"
#include <vector>
#include <unordered_map>

class hpd {
public:
	typedef pq_types::size_type	size_type;
	typedef pq_types::node_type	node_type;
	typedef pq_types::value_type value_type;
private:
	enum { NONE = N };
	const succinct_tree &T;
	std::vector<node_type> chain;
	std::vector<size_type> len;
	std::unordered_map<node_type,size_type> which_chain, pos_in_chain;
	std::unordered_map<node_type,size_type> card;
	std::unordered_map<node_type,node_type> best_son, parent;
	std::unordered_map<size_type,size_type> head;
	size_type n, chain_id;

	size_type dfs( node_type x ) {
		if ( card.count(x) )
			return card[x];
		card[x] = 1;
		for ( const_iterator<int> it = T->begin(x); it != T->end(x); ++it ) {
			auto y = *it;
			parent[y] = x, card[x] += dfs(y);
			if ( !best_son.count(x) || card[y] > card[best_son[x]] )
				best_son[x]= y;
		}
		return card[x];
	}
	void hld( node_type x, bool new_chain = false ) {
		if ( new_chain ) 
			++chain_id;
		chain.push_back(x), which_chain[x] = chain_id;
		if ( best_son.count(x) )
			hld(best_son[x]);
		for ( const_iterator<int> it = T->begin(x); it != T->end(); ++it )
			if ( !best_son.count(x) || *it != best_son[x] )
				hld(*it,true);
	}

	node_type ref( node_type x ) const {
		return chain[head[which_chain[x]]];
	}

public:

	hpd( const succinct_tree &t ) { 
		n= (t= T).size();
	};

	std::tuple<bit_vector,bit_vector,bit_vector,int_vector<value_type>> 
	operator()() {

		chain_id = -1, chain.clear(), which_chain.clear(), pos_in_chain.clear(), len.reserve(n);
		hld(0,true);

		int ch,x,y,i,j,k,prev = -1;
		for ( k = 0, i = 0; i < n; ++i ) {
			if ( prev != (ch=which_chain[x=chain[i]]) ) 
				head[ch] = i, len[ch] = 0;
			pos_in_chain[x] = i, ++len[ch], prev = ch;
		}

		// we need the plain_tree to support add_arcs operation
		plain_tree pt(n);
		for ( x = 1; x < n; ++x ) 
			pt.add_arc(ref(parent[x]),x);
		pt.normalize();

		bit_vector B = bit_vector(2*n,0);
		for ( k = 0, x = 0; x < n; ++x ) {
			B[k++]= 1, ch= which_chain[x];
			if ( chain[head[ch]] == x ) 
				k += len[ch];
		}
		assert( k == 2*n );
		int_vector<value_type> C(n);
		for ( i = 0; i < n; C[i] = static_cast<value_type>(T->weight(chain[i])), ++i );

		return std::make_tuple(*T,pt,B,C);

	}
};

