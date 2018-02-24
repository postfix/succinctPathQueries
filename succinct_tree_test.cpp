#include "sdsl/int_vector.hpp"
#include "sdsl/bits.hpp"
#include "sdsl/util.hpp"
#include "rs_bitvector.hpp"
#include "raw_tree.hpp"
#include "bp_tree.hpp"
#include "succinct_tree.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <random>
#include <cstring>

namespace
{

TEST(trees,allOperations)
{
	std::string s;
	std::cin >> s;
	s = "("+s+")";
	succinct_tree *raw = new raw_tree(s);	
	int n = s.size()/2,i,j,k,x,y;
	sdsl::bit_vector b = sdsl::bit_vector(2*n,0);
	for ( i = 0; i < s.size(); ++i )
		if ( s[i] == '(' )
			b[i] = 1;
	succinct_tree *bp = new bp_tree(&b);

	ASSERT_EQ(bp->size(),raw->size());
	std::cout << "sizes match"<<"\n";

	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(0,n-1);
	auto dice = std::bind(distribution,generator);

	std::cout << "children: "<<"\n";
	// children
	for ( x = 0; x < n; ++x ) {
		auto cx = bp->children(x);
		auto rx = raw->children(x);
		ASSERT_EQ(cx.size(),rx.size());
		for ( i = 0; i < (int)cx.size(); ++i )
			ASSERT_EQ(cx[i],rx[i]);
	}
	std::cout << "children match"<<"\n";

	std::cout << "checking: \"parent\" and \"depth\" and \"leaf\""<<"\n";
	// parent, depth, leaf
	for ( x = 0; x < n; ++x ) 
		ASSERT_EQ(bp->parent(x),raw->parent(x));
	std::cout << "\"parent\" OK"<<"\n";
	for( x = 0; x < n; ++x )
		ASSERT_EQ(bp->depth(x),raw->depth(x)) << "x = " << x << "\n";
	std::cout << "\"depth\" OK"<<"\n";
	for ( x = 0; x < n; ++x )
		ASSERT_EQ(bp->is_leaf(x),raw->is_leaf(x));
	std::cout << "\"is_leaf\" OK"<<"\n";


	std::cout << "Checking: \"lca\", \"is_ancestor\""<<"\n";
	// lca, is_ancestor
	for ( k = 0; k < n/4; ++k ) {
		x = dice(), y = dice();
		ASSERT_EQ(bp->is_ancestor(x,y),raw->is_ancestor(x,y)) << x << " " << y << std::endl;
		ASSERT_EQ(bp->lca(x,y),raw->lca(x,y)) << x << " " << y << " " << raw->is_ancestor(y,x) << std::endl;
	}
	std::cout << "\"lca\" and \"is_ancestor\" OK"<<"\n";

	// i-th ancestor
	for ( x = 0; x < n; ++x ) {
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0,bp->depth(x));
		auto dice = std::bind(distribution,generator);
		for ( k = 0; k < 32; ++k ) {
			i = dice();
			ASSERT_EQ(bp->ancestor(x,i),raw->ancestor(x,i));
		}
	}

}

}  // namespace

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

