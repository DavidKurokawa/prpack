#ifndef PRPACK_SOLVER
#define PRPACK_SOLVER
#include "prpack_preprocessed_graph.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include "prpack_result.h"
#include <string>

namespace prpack {

	// Solver class.
	class prpack_solver {
		private:
			// instance variables
			prpack_adjacency_list* al;
			prpack_preprocessed_gs_graph* gsg;
			prpack_preprocessed_scc_graph* sccg;
			// methods
			void initialize();
			prpack_result* solve_via_gs(double alpha, double tol, int num_vs, int num_es, double* ii, int* heads, int* tails, double* inv_num_outlinks, double* uv, int nsc, int* divisions, int* decoding);
			prpack_result* solve_via_scc_gs(double alpha, double tol, int num_vs, int num_es, double* ii, int* heads, int* tails, double* inv_num_outlinks, double* uv, int nsc, int* divisions, int* decoding);
		public:
			// constructors
			prpack_solver(prpack_csr* g);
			prpack_solver(prpack_edge_list* g);
			prpack_solver(prpack_adjacency_list* g);
			prpack_solver(const std::string& filename);
			// methods
			prpack_result* solve(double alpha, double tol);
			prpack_result* solve(double alpha, double tol, double* u, double* v);
	};

};

#endif
