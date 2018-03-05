#include <sys/resource.h>
#include <sdsl/int_vector.hpp>
#include "sdsl/bits.hpp"
#include "sdsl/util.hpp"
#include "rs_bitvector.hpp"
#include "raw_tree.hpp"
#include "bp_tree.hpp"
#include "wt_hpd.hpp"
#include "plain_tree.hpp"
#include "succinct_tree.hpp"
#include "path_query_processor.hpp"
#include "gtest/gtest.h"
#include <vector>
#include <random>
#include <cstring>
#include "tree_extraction_bitvector.hpp"
#include "bender_farach_colton.hpp"

typedef sdsl::bit_vector bit_vector;
typedef sdsl::int_vector<> int_vector;
using namespace pq_types;

typedef sdsl::int_vector<> int_vector;

namespace {

	TEST(util,RMQ) {
		int n= 1000000;
		long long *A= new long long[n];
		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0,n-1), zo(0,1);
		auto dice= std::bind(distribution,generator),
			 zod= std::bind(zo,generator);
		A[0]= dice();
		for ( auto i= 1; i < n; A[i]= A[i-1]+(zod()?1:-1), ++i ) ;
		pm_one_rmq *rmq= new pm_one_rmq(A,n);
		puts("Created the object");
		for ( auto it= 0; it < (1<<15); ++it ) {
			auto x= dice(), y= dice();
			if ( x > y ) std::swap(x,y);
			int idx= 0;
			long long mi= +(1LL<<62);
			for ( auto i= x; i <= y; ++i )
				if ( mi > A[i] )
					mi= A[i];
			ASSERT_EQ(mi,A[(*rmq)(x,y)]);
		}
		delete A;
	}
};

int main(int argc, char** argv)
{
   ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

