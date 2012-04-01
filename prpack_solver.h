#ifndef PRPACK_SOLVER
#define PRPACK_SOLVER
#include "prpack_preprocessed_graph.h"
#include "prpack_csr.h"
#include "prpack_edgelist.h"
#include "prpack_result.h"
#include <string>

namespace prpack {

	// Solver class.
	class prpack_solver {
		private:
			prpack_preprocessed_graph* g;
		public:
			prpack_solver(prpack_csr* g);
			prpack_solver(prpack_edgelist* g);
			prpack_solver(const std::string& filename);
			prpack_result* solve(double alpha, double tol);
			prpack_result* solve(double alpha, double tol, double* u, double* v);
	};

};

#endif
