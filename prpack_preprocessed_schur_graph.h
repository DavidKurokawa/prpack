#ifndef PRPACK_PREPROCESSED_SCHUR_GRAPH
#define PRPACK_PREPROCESSED_SCHUR_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_base_graph.h"

namespace prpack {

	class prpack_preprocessed_schur_graph : public prpack_preprocessed_graph {
		public:
			// instance variables
			int num_dangling_vs;
			int* encoding;
			int* decoding;
			// constructor
			prpack_preprocessed_schur_graph(prpack_base_graph* bg);
	};

};

#endif
