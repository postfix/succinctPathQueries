
#include <sys/resource.h>
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

typedef sdsl::bit_vector bit_vector;
typedef sdsl::int_vector<> int_vector;
using namespace pq_types;

namespace {

		TEST(trees,pathQueries01) {
			std::string s;
			std::cin>>s;
			auto n = s.size()/2;

			/*
			puts("setting up the size...");
			const rlim_t kStackSize= (1<<21)*100;
			struct rlimit r1;
			int result;
			result= getrlimit(RLIMIT_STACK,&r1);
			if ( result == 0 ) {
				if ( r1.rlim_cur < kStackSize ) {
					r1.rlim_cur= kStackSize;
					result= setrlimit(RLIMIT_STACK,&r1);
					if ( result != 0 ) {
						fprintf(stderr,"setrlimit returned result = %d\n",result);
						assert( false );
					}
				}
			}
			puts("Stack size set");
			*/

			std::vector<pq_types::value_type> w(n);

			for ( auto l = 0; l < n; std::cin >> w[l++] ) ;

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	
			succinct_tree *T   = new bp_tree(s);
			assert( T );

			hpd *H= new hpd(T);

			auto bundle= (*H)();
			auto bv= std::get<0>(bundle);
			succinct_tree *original= T,
						  *condensed= new bp_tree(&bv);
			bit_vector B= std::get<1>(bundle);
			std::vector<node_type> chain= std::get<2>(bundle);

			wt_int wt;
			
			int_vector weights= int_vector(T->size());
			for ( auto l = 0; l < T->size(); ++l ) weights[l]= w[chain[l]];
			construct_im(wt,weights);
			path_query_processor *processor= new wt_hpd(original,condensed,&wt,B);

			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(weights[i],processor->weight_of(chain[i])) << "x = " << chain[i] << "\n";

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				//std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}
		TEST(trees,pathQueries02) {
			std::string s;
			std::cin>>s;
			auto n = s.size()/2;

			/*
			puts("setting up the size...");
			const rlim_t kStackSize= (1<<21)*100;
			struct rlimit r1;
			int result;
			result= getrlimit(RLIMIT_STACK,&r1);
			if ( result == 0 ) {
				if ( r1.rlim_cur < kStackSize ) {
					r1.rlim_cur= kStackSize;
					result= setrlimit(RLIMIT_STACK,&r1);
					if ( result != 0 ) {
						fprintf(stderr,"setrlimit returned result = %d\n",result);
						assert( false );
					}
				}
			}
			puts("Stack size set");
			*/

			std::vector<pq_types::value_type> w(n);

			for ( auto l = 0; l < n; std::cin >> w[l++] ) ;

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	
			succinct_tree *T   = new bp_tree(s);
			assert( T );

			hpd *H= new hpd(T);

			auto bundle= (*H)();
			auto bv= std::get<0>(bundle);
			succinct_tree *original= T,
						  *condensed= new bp_tree(&bv);
			bit_vector B= std::get<1>(bundle);
			std::vector<node_type> chain= std::get<2>(bundle);

			wt_int wt;
			
			int_vector weights= int_vector(T->size());
			for ( auto l = 0; l < T->size(); ++l ) weights[l]= w[chain[l]];
			construct_im(wt,weights);
			path_query_processor *processor= new wt_hpd(original,condensed,&wt,B);

			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(weights[i],processor->weight_of(chain[i])) << "x = " << chain[i] << "\n";

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				//std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}
		TEST(trees,pathQueries03) {
			std::string s;
			std::cin>>s;
			auto n = s.size()/2;

			/*
			puts("setting up the size...");
			const rlim_t kStackSize= (1<<21)*100;
			struct rlimit r1;
			int result;
			result= getrlimit(RLIMIT_STACK,&r1);
			if ( result == 0 ) {
				if ( r1.rlim_cur < kStackSize ) {
					r1.rlim_cur= kStackSize;
					result= setrlimit(RLIMIT_STACK,&r1);
					if ( result != 0 ) {
						fprintf(stderr,"setrlimit returned result = %d\n",result);
						assert( false );
					}
				}
			}
			puts("Stack size set");
			*/

			std::vector<pq_types::value_type> w(n);

			for ( auto l = 0; l < n; std::cin >> w[l++] ) ;

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	
			succinct_tree *T   = new bp_tree(s);
			assert( T );

			hpd *H= new hpd(T);

			auto bundle= (*H)();
			auto bv= std::get<0>(bundle);
			succinct_tree *original= T,
						  *condensed= new bp_tree(&bv);
			bit_vector B= std::get<1>(bundle);
			std::vector<node_type> chain= std::get<2>(bundle);

			wt_int wt;
			
			int_vector weights= int_vector(T->size());
			for ( auto l = 0; l < T->size(); ++l ) weights[l]= w[chain[l]];
			construct_im(wt,weights);
			path_query_processor *processor= new wt_hpd(original,condensed,&wt,B);

			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(weights[i],processor->weight_of(chain[i])) << "x = " << chain[i] << "\n";

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				//std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}
		TEST(trees,pathQueries04) {
			std::string s;
			std::cin>>s;
			auto n = s.size()/2;

			/*
			puts("setting up the size...");
			const rlim_t kStackSize= (1<<21)*100;
			struct rlimit r1;
			int result;
			result= getrlimit(RLIMIT_STACK,&r1);
			if ( result == 0 ) {
				if ( r1.rlim_cur < kStackSize ) {
					r1.rlim_cur= kStackSize;
					result= setrlimit(RLIMIT_STACK,&r1);
					if ( result != 0 ) {
						fprintf(stderr,"setrlimit returned result = %d\n",result);
						assert( false );
					}
				}
			}
			puts("Stack size set");
			*/

			std::vector<pq_types::value_type> w(n);

			for ( auto l = 0; l < n; std::cin >> w[l++] ) ;

			std::default_random_engine generator;
			std::uniform_int_distribution<int> distribution(0,n-1);
			auto dice= std::bind(distribution,generator);

			succinct_tree *raw = new raw_tree(s);	
			succinct_tree *T   = new bp_tree(s);
			assert( T );

			hpd *H= new hpd(T);

			auto bundle= (*H)();
			auto bv= std::get<0>(bundle);
			succinct_tree *original= T,
						  *condensed= new bp_tree(&bv);
			bit_vector B= std::get<1>(bundle);
			std::vector<node_type> chain= std::get<2>(bundle);

			wt_int wt;
			
			int_vector weights= int_vector(T->size());
			for ( auto l = 0; l < T->size(); ++l ) weights[l]= w[chain[l]];
			construct_im(wt,weights);
			path_query_processor *processor= new wt_hpd(original,condensed,&wt,B);

			for ( auto i= 0; i < n; ++i )
				ASSERT_EQ(weights[i],processor->weight_of(chain[i])) << "x = " << chain[i] << "\n";

			int cs= 0;
			for ( auto it= 0; it < (1<<15); ++it ) {
				auto x= dice(), y= dice();
				//std::cout << "\n\033[1;34mquery: "<< x <<" "<<y<<":\033[0m\n";
				ASSERT_EQ(processor->query(x,y),(dynamic_cast<raw_tree *>(raw))->query(x,y,w));
				//std::cout <<"\033[1;32mPassed: "<<++cs << "\033[0;m"<<std::endl;
			}
		}


};

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

