#ifndef PRPACK_PREPROCESSED_SCC_GRAPH
#define PRPACK_PREPROCESSED_SCC_GRAPH
#include "prpack_graph.h"
#include "prpack_adjacency_list.h"

namespace prpack {

	// Pre-processed graph class
	class prpack_preprocessed_scc_graph : public prpack_graph {
		public:
			// instance variables
			int num_es_inside;
			int* heads_inside;
			int* tails_inside;
			int num_es_outside;
			int* heads_outside;
			int* tails_outside;
			double* ii;
			double* inv_num_outlinks;
			int num_comps;
			int* divisions;
			int* decoding;
			// constructor
			prpack_preprocessed_scc_graph(prpack_adjacency_list* al);
	};

};

#endif
