#include "succinct_tree.hpp"
#include "bp_tree.hpp"
#include "path_query_processor.hpp"
#include "rs_bitvector.hpp"
#include "sdsl/wavelet_trees.hpp"
#include <cassert>
#include <vector>

using namespace pq_types;
typedef pq_types::wt_int wt_int;
typedef sdsl::bit_vector bit_vector;


// wt_hpd: implements the path_query_processor interface
class wt_hpd: public path_query_processor {
private:
	size_type m;
	const succinct_tree *original= nullptr, 
				  		*condensed= nullptr;
	const wt_int *wavelet_tree;
	rs_bitvector B;

	size_type ref_count( node_type x ) const {
		return B.select(x+2)-B.select(x+1)-1;
	}

	node_type ref( node_type x ) const {
		if ( is_head_of_chain(x) )
			return x;
		return condensed->parent(x);
	}

	bool is_head_of_chain( node_type x ) const {
		return ref_count(x) > 0;
	}

	size_type position_in_chain( const node_type x ) const {
		node_type px= ref(x);
		return B.select(px+1) - px + original->depth(x) - original->depth(px);
	}

	// half-open segments exclusive of "p" itself
	void get_intervals( node_type p, node_type x, 
						std::vector<std::pair<size_type,size_type>> res, 
						bool inclusive = false ) const {
		// assert( original.is_ancestor(p,x) );
		auto pp = ref(p);
		for ( auto px = ref(x); px != pp; ) {
			res.push_back({position_in_chain(px),position_in_chain(x)+1});
			x = condensed->parent(px), px = ref(x);
		}
		res.push_back({position_in_chain(p)+(inclusive?0:1),position_in_chain(x)+1});
	}

	value_type query( std::vector<std::pair<size_type,size_type>> &vec, size_type k ) const {
		return wavelet_tree->range_quantile(vec,k);
	}

public:

	typedef pq_types::node_type		node_type;
	typedef pq_types::value_type 	value_type;
	using stree = succinct_tree;

	value_type weight( const node_type x ) const {
		return (*wavelet_tree)[position_in_chain(x)];
	}

	size_type size() const { return m; }

	wt_hpd( const stree *o, const stree *c, const wt_int *w, const bit_vector &b ) {
		original= o;
		condensed= c;
		wavelet_tree= w;
		B= rs_bitvector(b);
	}

	value_type 
	query( const node_type x, const node_type y ) const {
		auto z = original->lca(x,y);
		std::vector<std::pair<size_type,size_type>> segments;
		assert( segments.size() == 0 );
		get_intervals(z,x,segments);
		get_intervals(z,z,segments,true);
		get_intervals(z,y,segments);
		size_type k = original->depth(x)+original->depth(y)-2*original->depth(z);
		return query(segments,k>>1);
	}

};

