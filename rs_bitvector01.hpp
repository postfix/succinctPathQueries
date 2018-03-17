#ifndef RS_BITVECTOR01_INCLUDED
#define RS_BITVECTOR01_INCLUDED

#include "pq_types.hpp"
#include <sdsl/rrr_vector.hpp>

//! Wrapper class for sd-bitvector and its rank/select supports
class rs_bitvector01 {

	sdsl::rrr_vector<> rrr_b{};
	sdsl::rrr_vector<> :: select_1_type 	rrr_sel1{&rrr_b};
	sdsl::rrr_vector<> :: select_0_type 	rrr_sel0{&rrr_b};
	sdsl::rrr_vector<> :: rank_1_type		rrr_rank1{&rrr_b};
	sdsl::rrr_vector<> :: rank_0_type		rrr_rank0{&rrr_b};
	pq_types::size_type _sz= 0;

public:

	typedef pq_types::size_type		size_type;
	typedef pq_types::value_type    value_type;

	//rs_bitvector() {};
	//rs_bitvector( const sdsl::bit_vector &bv ) : sd_b(bv), sdb_sel(&sd_b), sdb_rank(&sd_b) {}
	explicit rs_bitvector01( const sdsl::bit_vector &bv ) {
		_sz= sdsl::size_in_bytes(bv);
		rrr_b = sdsl::rrr_vector<>(bv); 
		rrr_sel1 = sdsl::rrr_vector<>::select_1_type(&rrr_b);
		rrr_rank1 = sdsl::rrr_vector<>::rank_1_type(&rrr_b);
		rrr_sel0 = sdsl::rrr_vector<>::select_0_type(&rrr_b);
		rrr_rank0 = sdsl::rrr_vector<>::rank_0_type(&rrr_b);
	}

	double size_in_bytes() const {
		return _sz+(sdsl::size_in_bytes(rrr_b)+sdsl::size_in_bytes(rrr_sel0)+sdsl::size_in_bytes(rrr_rank0)+sdsl::size_in_bytes(rrr_sel1)+sdsl::size_in_bytes(rrr_rank1));
	}

	//! return the number of 1-bits in \f$ 0 \leq j < i\f$
	/*! \param i An index i with \f0 \leq i < size()\f$
	 */
	size_type rank( size_type x, size_type i ) const {
		return i?rrr_rank1(x):rrr_rank0(x);
	}

	//! return the 0-based position of the i-th occurrence in the bit vector
	/*! \param i An index i with \f$1 \leq i \leq size()\f$
	 */
	const size_type select( size_type x, size_type i ) const {
		return i?rrr_sel1(x):rrr_sel0(x);
	}

	value_type operator[] ( size_type i ) const {
		return rrr_b[i];
	}
};

#endif
