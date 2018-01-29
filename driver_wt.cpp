/*
 * Input format:
 * 2 lines per test instance:
 * 	-- bp sequence of '(' and ')'
 * 	-- weights of the nodes in pre-order
 * 	The size of the tree is inferred from the bp sequence
 */
#include <bits/stdc++>
#include <cstdio>
#include <chrono>
#include <vector>
#include <iostream>
#include "plain_tree.hpp"
#include "succinct_tree.hpp"
#include "wavelet_trees.hpp"
#define N (1<<22)

typedef sdsl::bit_vector bit_vector;
typedef sdsl::int_vector<value_type> int_vector;
typedef sdsl::wt_int<sdsl::rrr_vector<63>> wt_int;

sdsl::value_type w[N];

int main( int argc, char **argv ) {
	int i,j,k,n,qr,cs = 0,oqr;
	double ax = 0;
	FILE *fp = fopen(argv[1],"w");
	for ( string s; cin >> s; ) {
		succinct_tree T = bp_tree(s);
		hpd H = hpd(T);
		auto bundle = H();
		succinct_tree &original= std::get<0>(bundle),
					  condensed= bp_tree(std::get<1>(bundle));
		bit_vector B= std::get<2>(bundle);
		std::vector<node_type> chain= std::get<3>(bundle);
		wt_int wt;
		// reading the weights
		for ( auto l = 0; l < T.size(); cin >> w[l++] ) ;
		int_vector weights= int_vector<>(T.size());
		for ( auto l = 0; l < T.size(); ++l ) weights[l]= w[chain[l]];
		construct_im(wt,chain);
		path_query_processor processor = wt_hpd(original,condensed,wt,B);
		auto start = std::chrono::high_resolution_clock::now();
		for ( scanf("%d",&qr), oqr= qr; qr-- && 2 == scanf("%d %d",&i,&j); ) {
			value_type res= processor.query(static_cast<node_type>(i),static_cast<node_type>(j));
			printf("%d\n",(int)res);
		}
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		fprintf(fp,"%.6lf\n",elapsed.count()/oqr);
		ax += elapsed.count()/oqr;
	}
	fprintf(fp,"%.6lf\n",ax/cs);
	fclose(fp);
	return 0;
}
