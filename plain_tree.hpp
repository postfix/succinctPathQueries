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
#define  MAXLG (30)
#define  MAXN   (1<<MAXLG)
#define  infty  (MAXN)

/*
 * this class is useless, and is here only for historic reasons
 */
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

/* now, this class is the real thing */
class hpd {
public:
	typedef pq_types::size_type	size_type;
	typedef pq_types::node_type	node_type;
	typedef pq_types::value_type value_type;
	typedef sdsl::bit_vector bit_vector;
	typedef sdsl::int_vector<> int_vector;
private:
	const succinct_tree *T;
	mutable size_type *which_chain;
	node_type *best_son,*chain,*len,*head_of_chain;
	size_type n,num_chains;

	size_type dfs( node_type x ) {
		size_type c= 1, cc, bc= 0;
		assert( T );
		assert( T->size() > x );
		std::vector<node_type> children= T->children(x);
		best_son[x]= n;
		for ( auto y: children ) {
			c+= (cc= dfs(y));
			if ( cc > bc )
				bc= cc, best_son[x]= y;
		}
		return c;
	}

	void hld( node_type x, bool is_start= false ) {
		which_chain[x]= num_chains-1;
		if ( is_start )
			head_of_chain[num_chains-1]= x;
		if ( best_son[x] < n ) 
			hld(best_son[x]);
		std::vector<node_type> children= T->children(x);
		for ( auto y: children )
			if ( y != best_son[x] )
				++num_chains, hld(y,true);
	}

	node_type ref( node_type x ) const {
		return head_of_chain[which_chain[x]];
	}

public:

	hpd( const succinct_tree *t ) { 
		assert( t );
		n= (T= t)->size();
		best_son= new node_type[n];
		head_of_chain= new node_type[n];
		which_chain= new size_type[n];
	};

	std::tuple<bit_vector, bit_vector, std::vector<node_type>> 
	operator()() {

		dfs(0), num_chains= 1, hld(0,true);

		size_type i,j,k,ch;
		node_type x,y;

		for ( len= best_son, ch= 0; ch < num_chains; len[ch++]= 0 ) ;
		for ( x= 0; x < n; ++len[which_chain[x++]] ) ;

		plain_tree pt{};
		for ( x= 1; x < n; ++x ) 
			pt.add_arc(ref(T->parent(x)),x);
		assert( pt.nedges()+1 == n );

		bit_vector B= bit_vector(2*n,0);
		for ( k= 0, x= 0; x < n; ++x ) {
			B[k++]= 1, ch= which_chain[x];
			if ( head_of_chain[ch] == x ) 
				k+= len[ch];
		}

		assert( k == 2*n );
		std::vector<node_type> rchain(n);
		node_type *counts= best_son;
		for ( x= 0; x < n; counts[x++]= 0 ) ;
		for ( x= 0; x < n; ++counts[ref(x++)] ) ;
		for ( x= 1; x < n; counts[x]+= counts[x-1], ++x ) ;
		for ( long long z= n-1; z >= 0; --z ) {
			assert( counts[ref((node_type)z)] );
			rchain[--counts[ref((node_type)z)]]= (node_type)z;
		}
		return std::make_tuple(bit_vector(pt),B,rchain);
	}

	~hpd() {
		delete best_son;
		delete head_of_chain;
		delete which_chain;
	}
};

