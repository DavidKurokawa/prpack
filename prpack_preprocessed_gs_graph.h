#ifndef PRPACK_PREPROCESSED_GS_GRAPH
#define PRPACK_PREPROCESSED_GS_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_adjacency_list.h"

namespace prpack {

	// Pre-processed graph class
	class prpack_preprocessed_gs_graph : public prpack_preprocessed_graph {
		public:
			// constructor
			prpack_preprocessed_gs_graph(prpack_adjacency_list* al);
	};

};

#endif
