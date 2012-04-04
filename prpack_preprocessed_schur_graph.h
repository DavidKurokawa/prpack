#ifndef PRPACK_PREPROCESSED_SCHUR_GRAPH
#define PRPACK_PREPROCESSED_SCHUR_GRAPH
#include "prpack_graph.h"
#include "prpack_adjacency_list.h"

namespace prpack {

	class prpack_preprocessed_schur_graph : public prpack_graph {
		public:
			// instance variables
			int num_dangling_vs;
			double* ii;
			double* inv_num_outlinks;
			int* decoding;
			// constructor
			prpack_preprocessed_schur_graph(prpack_adjacency_list* al);
	};

};

#endif
