#ifndef PATH_QUERY_PROCESSOR_INCLUDED
#define PATH_QUERY_PROCESSOR_INCLUDED

#include "pq_types.hpp"

class path_query_processor {
public:
	typedef pq_types::node_type		node_type;
	typedef pq_types::value_type 	value_type;

	virtual value_type weight( const node_type x ) const ;

	// median query: returns the median weight on the x-to-y path
	virtual value_type query( const node_type x, const node_type y ) const ;
};

#endif
