#ifndef PQ_TYPES_INCLUDED
#define PQ_TYPES_INCLUDED

#include "sdsl/bit_vectors.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/wavelet_trees.hpp"

namespace pq_types {
	typedef sdsl::bit_vector::size_type 	    size_type;
	typedef sdsl::bit_vector::size_type 	    node_type;
	typedef sdsl::int_vector<32>::value_type	value_type;
	typedef sdsl::wt_int<sdsl::rrr_vector<63>> wt_int;
};

#endif
