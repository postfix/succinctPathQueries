#ifndef RS_BITVECTOR_INCLUDED
#define RS_BITVECTOR_INCLUDED

#include "pq_types.hpp"
#include <sdsl/rrr_vector.hpp>

//! Wrapper class for sd-bitvector and its rank/select supports
class rs_bitvector {

	sdsl::rrr_vector<> rrr_b{};
	sdsl::rrr_vector<> :: select_1_type 	rrr_sel{&rrr_b};
	sdsl::rrr_vector<> :: rank_1_type		rrr_rank{&rrr_b};
	pq_types::size_type _sz= 0;

public:

	typedef pq_types::size_type		size_type;
	typedef pq_types::value_type    value_type;

	//rs_bitvector() {};
	//rs_bitvector( const sdsl::bit_vector &bv ) : sd_b(bv), sdb_sel(&sd_b), sdb_rank(&sd_b) {}
	explicit rs_bitvector( const sdsl::bit_vector &bv ) {
		/*_sz+= sdsl::size_in_bytes(bv);
		std::cout << _sz << std::endl;*/
		rrr_b = sdsl::rrr_vector<>(bv); 
		rrr_sel = sdsl::rrr_vector<>::select_1_type(&rrr_b);
		rrr_rank = sdsl::rrr_vector<>::rank_1_type(&rrr_b);
	}

	double size_in_bytes() const {
		return _sz+(sdsl::size_in_bytes(rrr_b)+sdsl::size_in_bytes(rrr_sel)+sdsl::size_in_bytes(rrr_rank));
	}

	//! return the number of 1-bits in \f$ 0 \leq j < i\f$
	/*! \param i An index i with \f0 \leq i < size()\f$
	 */
	size_type rank( size_type i ) const {
		return rrr_rank(i);
	}

	//! return the 0-based position of the i-th occurrence in the bit vector
	/*! \param i An index i with \f$1 \leq i \leq size()\f$
	 */
	const size_type select( size_type i ) const {
		return rrr_sel(i);
	}

	value_type operator[] ( size_type i ) const {
		return rrr_b[i];
	}

};

#endif
