/*
 * Tree Extraction technique for path selection queries
 * should work in O(max(logn,log\sigma))
 */
#include <cassert>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <set>
#include <map>
#define N (1<<25)
#define MAXE (N<<1)
#define BIT(k) (1ULL<<(k))
#define MASK(k) (BIT(k)-1ULL)
#define L(k) ((k)&((~(k))+1ULL))
#define NONE (-1)
#include <chrono>
using namespace std;
typedef unsigned long long u64;
typedef long long i64;

/* the T_{a,b} storing
 * also left = T_{a,mid}, right = T_{mid+1,b} as fields */
class Tree {
#define ROOT 0
private:
	int n,a,b,mid,cardinality;
	int *ptr_left, *ptr_right, 
		*weight, *d,
		**anc,K;
	vector<int> *adj; // adj[x] lists all the children of "x"
	Tree *left, *right;
	int up( int x, unsigned int u ) const { // find u^th ancestor of "x"
		for ( int k = 0; u; u >>= 1, ++k )
			if ( u & 1 ) x = anc[x][k];
		return x;
	}
	int lca( int x, int y ) const {
		if ( d[x] > d[y] )
			return lca(up(x,d[x]-d[y]),y);
		if ( d[x] < d[y] )
			return lca(y,x);
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
	/* 
	 * "gpl" and "gpr" are "growth points" in left and right tree, resp.
	 * if, in the current tree T, the current node "x" goes to the left subtree,
	 * i.e. its weight is in [a..mid], then we attach it to "gpl"
	 * and make "ngpl" ("new" gpl) the node in "left" that corresponds to "x";
	 * this way, we conduct two extractions in parallel
	 */
	void dfs( int gpl, int gpr, int x, int depth ) {
		int ngpl = gpl, ngpr = gpr, y, k;
		assert( x == ROOT || a <= weight[x] && weight[x] <= b );
		if ( a < b && a <= weight[x] && weight[x] <= b ) 
			weight[x]<=mid?(ngpl=left->add_son(gpl,weight[x])):(ngpr=right->add_son(gpr,weight[x]));
		ptr_left[x] = ngpl, ptr_right[x] = ngpr, d[x] = depth;
		for ( int i = 0; i < (int)adj[x].size(); dfs(ngpl,ngpr,y,1+depth) ) 
			for ( anc[y = adj[x][i++]][0] = x, k = 1; anc[y][k-1] != ROOT; anc[y][k] = anc[anc[y][k-1]][k-1], ++k ) ;
	}
	void init( int card ) {
		assert( a <= b );
		assert( card >= 0 );
		cardinality = card;
		ptr_left = (int *)malloc((card+1) * sizeof *ptr_left);
		assert( ptr_left );
		ptr_right = (int *)malloc((card+1) * sizeof *ptr_right);
		assert( ptr_right );
		for ( int x = 0; x <= card; ++x )
			ptr_left[x] = ptr_right[x] = ROOT;
		weight = (int *)malloc((card+1) * sizeof *weight);
		assert( weight );
		d = (int *)calloc((card+1), sizeof *d);
		assert( d );
		adj = new vector<int>[card+1];
		assert( adj );
		for ( int x = 0; x <= card; adj[x++].clear() ) ;
		for ( n = 1, K = 0; (1<<K) <= card; ++K ) ;
		anc = (int **)malloc((card+1)*sizeof *anc);
		assert( anc );
		for ( int x = 0; x <= card; ++x ) {
			anc[x] = (int *)malloc(K*sizeof *anc[x]);
			assert( anc[x] );
		}
		for ( int x = 0; x <= card; ++x )
			for ( int k = 0; k < K; anc[x][k++] = ROOT ) ;
		d[ROOT] = 0, weight[ROOT] = 0;
	}
	int get_depth( const int x ) const { return d[x]; }
	int select( int x, int y, int z, const int c, const int wc, int k ) const {
		int xleft, xright, yleft, yright, zleft, zright;
		assert( d[x]+d[y]-2*d[z]+1 > k );
		if ( a == b ) {
			assert( k < cardinality );
			return a;
		}
		xleft  = ptr_left[x], yleft   = ptr_left[y], zleft   = ptr_left[z];
		xright = ptr_right[x], yright = ptr_right[y], zright = ptr_right[z];
		int cn = left?left->get_depth(xleft)+left->get_depth(yleft)-2*left->get_depth(zleft)+(a <= wc && wc <= mid?1:0):0;
		if ( k >= cn ) 
			assert( right );
		return k < cn ? left->select(xleft,yleft,zleft,c,wc,k) : right->select(xright,yright,zright,c,wc,k-cn);
	}
public:
	Tree( int a, int b, int card ) : a(a), b(b) { 
		assert( this->a == a );
		assert( this->b == b );
		mid = ((a+b)>>1), init(card);
		left = right = NULL;
	};
	/* add <x,y> arc, where "x" is the parent of "y";
	 * "x" can be NONE -- when "y" is root */
	void add_arc( int x, int y ) {
		if ( x >= n ) n = x+1;
		if ( y >= n ) n = y+1;
		assert( 0 <= y && y < n );
		if ( (anc[y][0] = x) != NONE ) 
			adj[x].push_back(y);
	}
	void set_weight( int x, int w ) { assert( a <= w && w <= b ); weight[x] = w; }
	/* "grow" the tree starting at "px": attach a new vertex with weight "w"
	 * "px" can be NONE */
	int add_son( int px, int w ) {
		assert( a <= w );
		assert( w <= b );
		int x = n++;
		assert( x <= cardinality );
		add_arc(px,x), set_weight(x,w);
		return x;
	}
	/* initialize "left" and "right" trees for ranges
	 * [a..mid] and [mid+1..b], resp.;
	 * launch a DFS with two NONEs as "growth points" 
	 * in the two trees */
	void preprocess() {
		if ( !cardinality ) return ;
		if ( a < b ) {
			int ll = 0, rr = 0;
			for ( int i = ROOT+1; i < n; ++i ) {
				assert( a <= weight[i] && weight[i] <= b );
				if ( weight[i] <= mid ) ++ll; else ++rr;
			}
			left = ll?new Tree(a,mid,ll):NULL;
		   	right = rr?new Tree(mid+1,b,rr):NULL;
		}
		else left = right = NULL;
		dfs(ROOT,ROOT,ROOT,d[ROOT]);
		if ( left ) left->preprocess();
		if ( right ) right->preprocess();
	}
	int kth_weight( int u, int v, int k ) const {
		int c = lca(u,v);
		return select(u,v,c,c,weight[c],k);
	}
	int median( int x, int y ) const {
		int len = get_depth(x) + get_depth(y) - 2*get_depth(lca(x,y)) + 1;
		return kth_weight(x,y,len>>1);
	}
	int get_size() const { return n; }
	~Tree() {
		if ( !n || n == 1 ) return ;
		free(ptr_left), free(ptr_right), free(d), free(weight);
		for ( int i = 0; i < n; free(anc[i++]) ) ;
		free(anc);
		if ( left ) delete left;
		if ( right ) delete right;
	}
} *T;

int n,sigma,weight[N],
	to[MAXE],_next[MAXE],last[N],E,
	seen[N],yes;

void add_arcs( int x, int y ) {
	int i = E++, j = E++;
	to[i] = y, _next[i] = last[x], last[x] = i;
	to[j] = x, _next[j] = last[y], last[y] = j;
}

/* orient the given undirected graph via a DFS
 * and insert the obtained oriented edges to the tree */
void dfs( int x ) {
	int i,y;
	assert( seen[x] != yes );
	for ( T->set_weight(x,weight[x]), seen[x] = yes, i = last[x]; i != NONE; i = _next[i] )
		if ( seen[y = to[i]] != yes ) 
			T->add_arc(x,y), T->set_weight(y,weight[y]), dfs(y);
}

int main( int argc, char **argv ) {
	int i,j,k,ts,cs = 0,qr,oqr;
	double ax = 0;
	FILE *fp = fopen(argv[1],"w");
	for ( ;2==scanf("%d %d",&n,&sigma) && ++yes; ) {
		for ( E = 0, i = 1; i <= n; last[i] = NONE, scanf("%d",&weight[i++]) ) ;
		for ( k = 0; k < n-1 && 2 == scanf("%d %d",&i,&j); add_arcs(i+1,j+1), ++k ) ;
		T = new Tree(1,sigma,n), T->add_arc(ROOT,1), T->set_weight(1,weight[1]), dfs(1);
		assert( T->get_size() == n+1 );
		T->preprocess();
		auto start = std::chrono::high_resolution_clock::now();
		for ( ++cs, scanf("%d",&qr), oqr = qr; qr-- && 2 == scanf("%d %d",&i,&j); )
			printf("%d\n",T->median(i+1,j+1));
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		ax += elapsed.count()/oqr;
		fprintf(fp,"%.6lf\n",elapsed.count()/oqr);
		delete T;
	}
	fprintf(fp,"%.6lf\n",ax/cs);
	fclose(fp);
	return 0;
}

