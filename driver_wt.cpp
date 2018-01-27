#include <bits/stdc++>
#include <chrono>
#include "plain_tree.hpp"

int main( int argc, char **argv ) {
	int i,j,k,n,qr,cs = 0,oqr;
	double ax = 0;
	FILE *fp = fopen(argv[1],"w");
	for ( ;(cin >> n) && ++cs; ) {
		plain_tree T(n);
		cin >> T;
		T.normalize();
		hpd H = hpd(&T);
		std::tuple<bit_vector,bit_vector,bit_vector,int_vector<value_type>> bundle = H();
		succinct_tree original  = bp_tree(std::get<0>(bundle)),
					  condensed = bp_tree(std::get<1>(bundle));
		bit_vector<> B = std::get<2>(bundle);
		int_vector<value_type> C = std::get<3>(bundle);
		path_query_processor processor = wt_hpd(original,condensed,wt_int(C),B);
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

