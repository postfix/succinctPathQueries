#include "sdsl/int_vector.hpp"
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
#include "tree_extraction.hpp"


typedef sdsl::bit_vector bit_vector;
typedef sdsl::int_vector<> int_vector;
using namespace pq_types;

typedef sdsl::int_vector<> int_vector;

namespace {

		TEST(trees,pathQueries01) {
			std::string s;
			std::cin >> s;
			// std::cout << s << std::endl;
			auto n = s.size()/2;
			std::vector<value_type> w(n);

			for ( auto l = 0; l < n; ++l ) {
				value_type ww;
				std::cin >> ww;
				w[l]= ww;
			}

			/*
			hpd_remapper *remapper = new hpd_remapper();
			auto ret = remapper->convert(s,w);
			s= std::get<0>(ret);
			w= std::get<1>(ret);
			// std::cout << s << std::endl;
			*/

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	

			path_query_processor *processor= new tree_extraction<>(s,w);
			/*
			for ( auto l= 0; l < w.size(); ++l )
				ASSERT_EQ(dynamic_cast<tree_extraction<> *>(processor)->weight_of(l),w[l]);
				*/

			/*
			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(i,(dynamic_cast<wt_hpd*>(processor))->position_in_chain(chain[i])) << "x = " << chain[i] << "\n";
			*/

			dynamic_cast<raw_tree *>(raw)->query(303,4764,w);
			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ((dynamic_cast<tree_extraction<> *>(processor))->lca(x,y),(dynamic_cast<raw_tree *>(raw))->lca(x,y));
				//auto l1= dynamic_cast<tree_extraction<> *>(processor)->lca(x,y), l2= (dynamic_cast<raw_tree *>(raw))->lca(x,y);
				//ASSERT_EQ(processor->weight_of(l1),w[l2] );
				ASSERT_EQ((dynamic_cast<raw_tree *>(raw))->query(x,y,w),processor->query(x,y));
				std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}

		TEST(trees,pathQueries02) {
			std::string s;
			std::cin >> s;
			// std::cout << s << std::endl;
			auto n = s.size()/2;
			std::vector<value_type> w(n);

			for ( auto l = 0; l < n; ++l ) {
				value_type ww;
				std::cin >> ww;
				w[l]= ww;
			}

			/*
			hpd_remapper *remapper = new hpd_remapper();
			auto ret = remapper->convert(s,w);
			s= std::get<0>(ret);
			w= std::get<1>(ret);
			// std::cout << s << std::endl;
			*/

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	

			path_query_processor *processor= new tree_extraction<>(s,w);

			/*
			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(i,(dynamic_cast<wt_hpd*>(processor))->position_in_chain(chain[i])) << "x = " << chain[i] << "\n";
			*/

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				ASSERT_EQ((dynamic_cast<tree_extraction<> *>(processor))->lca(x,y),(dynamic_cast<raw_tree *>(raw))->lca(x,y));
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}

		TEST(trees,pathQueries03) {
			std::string s;
			std::cin >> s;
			// std::cout << s << std::endl;
			auto n = s.size()/2;
			std::vector<value_type> w(n);

			for ( auto l = 0; l < n; ++l ) {
				value_type ww;
				std::cin >> ww;
				w[l]= ww;
			}

			/*
			hpd_remapper *remapper = new hpd_remapper();
			auto ret = remapper->convert(s,w);
			s= std::get<0>(ret);
			w= std::get<1>(ret);
			// std::cout << s << std::endl;
			*/

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	

			path_query_processor *processor= new tree_extraction<>(s,w);

			/*
			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(i,(dynamic_cast<wt_hpd*>(processor))->position_in_chain(chain[i])) << "x = " << chain[i] << "\n";
			*/

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				ASSERT_EQ((dynamic_cast<tree_extraction<> *>(processor))->lca(x,y),(dynamic_cast<raw_tree *>(raw))->lca(x,y));
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}

		TEST(trees,pathQueries04) {
			std::string s;
			std::cin >> s;
			// std::cout << s << std::endl;
			auto n = s.size()/2;
			std::vector<value_type> w(n);

			for ( auto l = 0; l < n; ++l ) {
				value_type ww;
				std::cin >> ww;
				w[l]= ww;
			}

			/*
			hpd_remapper *remapper = new hpd_remapper();
			auto ret = remapper->convert(s,w);
			s= std::get<0>(ret);
			w= std::get<1>(ret);
			// std::cout << s << std::endl;
			*/

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	

			path_query_processor *processor= new tree_extraction<>(s,w);

			/*
			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(i,(dynamic_cast<wt_hpd*>(processor))->position_in_chain(chain[i])) << "x = " << chain[i] << "\n";
			*/

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				ASSERT_EQ((dynamic_cast<tree_extraction<> *>(processor))->lca(x,y),(dynamic_cast<raw_tree *>(raw))->lca(x,y));
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}


};

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

