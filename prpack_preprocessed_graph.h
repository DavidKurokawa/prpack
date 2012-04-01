#ifndef PRPACK_PREPROCESSED_GRAPH
#define PRPACK_PREPROCESSED_GRAPH
#include <list>
#include <string>
#include "prpack_graph.h"

typedef std::list<int> matrix_type;

namespace prpack {

	// Pre-processed graph class
	class prpack_preprocessed_graph : public prpack_graph {
		private:
			void convert(matrix_type* matrix, int& num_vs, int& num_es, int*& x, int*& y);
		public:
			double* ii;
			double* inv_num_outlinks;
			prpack_preprocessed_graph(const std::string& filename);
	};

};

#endif
