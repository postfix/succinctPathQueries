/*
 * Input format:
 * 2 lines per test instance:
 * 	-- bp sequence of '(' and ')'
 * 	-- weights of the nodes in pre-order
 * 	The size of the tree is inferred from the bp sequence
 */
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
#include <cassert>
#include <vector>
#include <random>
#include <cstring>
#include <chrono>
#include "hpd_remapper.hpp"

typedef sdsl::bit_vector bit_vector;
typedef sdsl::int_vector<> int_vector;
typedef pq_types::wt_int wt_int;
typedef pq_types::value_type value_type;
typedef pq_types::node_type node_type;
typedef pq_types::size_type size_type;

int main( int argc, char **argv ) {
	int i,j,k,n,qr,cs = 0,oqr;
	double ax= 0,bx= 0,tmp;
	assert( argc >= 2 );
	FILE *fp = fopen(argv[1],"w");
	for ( std::string s; std::cin >> s; ++cs ) {
		auto n = s.size()/2;

		const rlim_t kStackSize= n*4*5;
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
		puts("Stack size successfully set");

		std::vector<pq_types::value_type> w(n);

		for ( auto l = 0; l < n; std::cin >> w[l++] ) ;

		hpd_remapper *remapper= new hpd_remapper();
		auto ret = remapper->convert(s,w);
		delete remapper;
		s= std::get<0>(ret);
		w= std::get<1>(ret);
		// std::cout << s << std::endl;

		puts("Remapping successful");

		std::default_random_engine generator;
		std::uniform_int_distribution<int> distribution(0,n-1);
		auto dice= std::bind(distribution,generator);

		//succinct_tree *raw = new raw_tree(s);	
		succinct_tree *T   = new bp_tree(s);
		puts("Constructed the original tree");

		hpd *H= new hpd(T);
		auto bundle= (*H)();
		delete H;
		auto bv= std::get<0>(bundle);
		succinct_tree *original= T,
					  *condensed= new bp_tree(&bv);
		bit_vector B= std::get<1>(bundle);
		std::vector<node_type> chain= std::get<2>(bundle);
		wt_int wt;
			
		int_vector weights= int_vector(T->size());
		for ( auto l = 0; l < T->size(); ++l ) weights[l]= w[chain[l]];
		construct_im(wt,weights);
		puts("Constructed the Wavelet tree");
		path_query_processor *processor= new wt_hpd(original,condensed,&wt,B);
		puts("Constructed the Processor");
		puts("Support structures built, now starting the queries");

		auto start = std::chrono::high_resolution_clock::now();
		for ( qr= (1<<9), oqr= qr; qr--; ) {
			i= dice(), j= dice();
			//printf("[%d] %d %d\n",qr,i,j);
			value_type res= processor->query(static_cast<node_type>(i),static_cast<node_type>(j));
			//printf("%d\n",(int)res);
		}
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		fprintf(fp,"%.6lf %.6lf\n",tmp=(dynamic_cast<wt_hpd*>(processor))->bits_per_node(),elapsed.count()/oqr);
		ax+= elapsed.count()/oqr, bx+= tmp;
	}
	fprintf(fp,"%.6lf %.6lf\n",bx/cs,ax/cs);
	fclose(fp);
	return 0;
}
