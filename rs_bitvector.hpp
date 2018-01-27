#ifndef RS_BITVECTOR_INCLUDED
#define RS_BITVECTOR_INCLUDED

#include "pq_types.h"
#include "sd_vector.hpp"

//! Wrapper class for sd-bitvector and its rank/select supports
class rs_bitvector {

	sd_vector<> sd_b;
	sd_vector<> :: select_1_type 	sdb_sel;
	sd_vector<> :: rank_1_type 		sdb_rank;

public:

	typedef pq_types::size_type		size_type;

	rs_bitvector( const bit_vector &bv ) sd_b(bv), sdb_sel(&sd_b), sdb_rank(&sd_b) {}

	//! return the number of 1-bits in \f$ 0 \leq j < i\f$
	/*! \param i An index i with \f0 \leq i < size()\f$
	 */
	size_type rank( size_type i ) const {
		return sdb_rank(i);
	}

	//! return the 0-based position of the i-th occurrence in the bit vector
	/*! \param i An index i with \f$1 \leq i \leq size()\f$
	 */
	size_type select( size_type i ) const {
		return sdb_sel(i);
	}

};

#endif
