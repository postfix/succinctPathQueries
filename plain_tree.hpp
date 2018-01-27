/*
 * plain, pointer-based tree
 */
#include <bits/stdc++>
#define BIT(k) (1ULL<<(x))
#define TST(a,x) (a[(x)>>3] & BIT((x)&7))
#define CLR(a,x) (a[(x)>>3] &= ~BIT((x)&7))
#define SET(a,x) (a[(x)>>3] |= BIT((x)&7))
#define N (1<<21)

#include "bit_vector.hpp"
#include "int_vector.hpp"

class plain_tree {
	int n,w[N];
	std::vector<int> adj[N];
	char bp_seq[2*N+1];
	unsigned char a[(N>>3)+8];
	bool is_weighted;
	bit_vector<> b;

#define seen(x) TST(a,x)
#define set_seen(x) SET(a,x)

	void dfs( int x, int &cur ) {
		int i,y;
		assert( !seen(x) );
		for (set_seen(x),bp_seq[cur]= '1',b[cur++]= 1,i=0;i<(int)adj[x].size();++i) {
			y = adj[x][i];
			assert( !seen(y) );
			dfs(y,cur);
		}
		bp_seq[cur]= '0', b[cur++] = 0;
	}

public:

	int size() const { return n; }

	void add_arc( int from, int to ) {
		adj[from].push_back(to);
	}

	plain_tree( int nn, bool flag = false ) n{nn} {
		int i;
		for ( i = 0; i < n; CLR(a,i), adj[i++].clear() ) ;
		is_weighted = flag, bp_seq[0] = bp_seq[2*nn] = '\0';
		b = bit_vector(2*n,0);
	};

	void normalize() {
		for ( int x = 0; x < n; ++x )
			sort(adj[x].begin(),adj[x].end());
		int cur = 0; dfs(0,cur);
	}

	istream &operator >> ( istream &is ) {
		if ( is_weighted )
			for ( auto j = 0; j < n; is >> w[j++] ) ;
		for ( int x,y,i = 0; i < n-1; ++i ) {
			is >> x >> y;
			add_arc(x,y);
		}
		return is;
	}

	int weight( int x ) const { return w[x]; }

	ostream &operator << ( ostream &os ) const {
		os << n << "\n";
		if ( is_weighted )
			for ( auto i = 0; i < n; os << w[i++] << " " ) ;
		os << "\n" << bp_seq << "\n";
	}

	operator bit_vector<>() {
		return b;
	}

	const_iterator<int> begin( int x ) const {
		return adj[x].begin();
	}

	const_iterator<int> end( int x ) const {
		return adj[x].end();
	}

};

class hpd {
	int best_son[N],card[N],head[N],parent[N];
	unsigned char a[(N>>3)+8];
#define seen(x) TST(a,x)
#define set_seen(x) SET(a,x)
	const plain_tree *T;
	vector<int> chain, which_chain, pos_in_chain, len;
	int chain_id;
	int n;

	int dfs( int x ) {
		assert( !seen(x) );
		set_seen(x), best_son[x] = -1, card[x] = 1;
		for ( const_iterator<int> it = T->begin(x); it != T->end(x); ++it ) {
			auto y = *it;
			assert( !seen(y) );
			parent[y] = x, card[x] += dfs(y);
			if ( best_son[x]==-1 || card[y] > card[best_son[x]] )
				best_son[x]= y;
		}
		return card[x];
	}
	void hld( int x, bool on_track ) {
		if ( !on_track ) 
			++chain_id;
		chain.push_back(x), which_chain[x] = chain_id;
		if ( best_son[x] != -1 )
			hld(best_son[x],true);
		for ( const_iterator<int> it = T->begin(x); it != T->end(); ++it )
			hld(*it,false);
	}

	int ref( int x ) const {
		return chain[head[which_chain[x]]];
	}

public:

	hpd( const plain_tree *t ) { 
		T = t, n = T->size();
		for ( int i = 0; i < n; CLR(a,i), ++i ) ;
	};

	std::tuple<bit_vector,bit_vector,bit_vector,int_vector<value_type>> 
	operator()() {

		chain_id = -1, chain.clear(), which_chain.reserve(n), pos_in_chain.reserve(n), len.reserve(n);
		hld(0,false);

		int ch,x,y,i,j,k,prev = -1;
		for ( k = 0, i = 0; i < n; ++i ) {
			if ( prev != (ch=which_chain[x=chain[i]]) ) 
				head[ch] = i, len[ch] = 0;
			pos_in_chain[x] = i, ++len[ch], prev = ch;
		}

		plain_tree pt(n);
		for ( x = 1; x < n; ++x ) 
			pt.add_arc(ref(parent[x]),x);
		pt.normalize();

		bit_vector B = bit_vector(2*n,0);
		for ( k = 0, x = 0; x < n; ++x ) {
			B[k++] = 1, ch = which_chain[x];
			if ( chain[head[ch]] == x ) 
				for ( j = len[ch]; j--; B[k++] = 0 ) ;
		}
		assert( k == 2*n );
		int_vector<value_type> C(n);
		for ( i = 0; i < n; C[i] = static_cast<value_type>(T->weight(chain[i])), ++i );

		return std::make_tuple(*T,pt,B,C);

	}
};

