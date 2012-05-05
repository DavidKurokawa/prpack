#ifndef PRPACK_PREPROCESSED_GS_GRAPH
#define PRPACK_PREPROCESSED_GS_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_base_graph.h"
#include "mex.h"

namespace prpack {

	// Pre-processed graph class
	class prpack_preprocessed_gs_graph : public prpack_preprocessed_graph {
		public:
			// instance variables
			int* heads;
			int* tails;
			// constructors
			prpack_preprocessed_gs_graph(prpack_base_graph* bg);
			prpack_preprocessed_gs_graph(const mxArray* a);
			// method
			mxArray* to_matlab_array() const;
	};

};

#endif
