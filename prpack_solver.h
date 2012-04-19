#ifndef PRPACK_SOLVER
#define PRPACK_SOLVER
#include "prpack_base_graph.h"
#include "prpack_preprocessed_gs_graph.h"
#include "prpack_preprocessed_schur_graph.h"
#include "prpack_preprocessed_scc_graph.h"
#include "prpack_result.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <string>

namespace prpack {

	// Solver class.
	class prpack_solver {
		private:
			// instance variables
			double read_time;
			prpack_base_graph* bg;
			prpack_preprocessed_gs_graph* gsg;
			prpack_preprocessed_schur_graph* sg;
			prpack_preprocessed_scc_graph* sccg;
			// methods
			void initialize();
			static prpack_result* solve_via_gs(
					double alpha,
					double tol,
					int num_vs,
					int num_es,
					int* heads,
					int* tails,
					double* ii,
					double* inv_num_outlinks,
					double* u,
					double* v);
			static prpack_result* solve_via_schur_gs(
					double alpha,
					double tol,
					int num_vs,
					int num_dangling_vs,
					int num_es,
					int* heads,
					int* tails,
					double* ii,
					double* inv_num_outlinks,
					double* uv,
					int* encoding,
					int* decoding,
					bool normalize = true);
			static prpack_result* solve_via_schur_gs_uv(
					double alpha,
					double tol,
					int num_vs,
					int num_dangling_vs,
					int num_es,
					int* heads,
					int* tails,
					double* ii,
					double* inv_num_outlinks,
					double* u,
					double* v,
					int* encoding,
					int* decoding);
			static prpack_result* solve_via_scc_gs(
					double alpha,
					double tol,
					int num_vs,
					int num_es_inside,
					int* heads_inside,
					int* tails_inside,
					int num_es_outside,
					int* heads_outside,
					int* tails_outside,
					double* ii,
					double* inv_num_outlinks,
					double* uv,
					int num_comps,
					int* divisions,
					int* encoding,
					int* decoding,
					bool normalize = true);
			static prpack_result* solve_via_scc_gs_uv(
					double alpha,
					double tol,
					int num_vs,
					int num_es_inside,
					int* heads_inside,
					int* tails_inside,
					int num_es_outside,
					int* heads_outside,
					int* tails_outside,
					double* ii,
					double* inv_num_outlinks,
					double* u,
					double* v,
					int num_comps,
					int* divisions,
					int* encoding,
					int* decoding);
			static double* permute(int length, double* a, int* coding);
			static prpack_result* combine_uv(
					int num_vs,
					double* inv_num_outlinks,
					int* encoding,
					double alpha,
					prpack_result* ret_u,
					prpack_result* ret_v);
		public:
			// constructors
			prpack_solver(prpack_csr* g);
			prpack_solver(prpack_edge_list* g);
			prpack_solver(prpack_base_graph* g);
			prpack_solver(const std::string& filename, const std::string& format);
			// methods
			prpack_result* solve(double alpha, double tol, const std::string& method);
			prpack_result* solve(double alpha, double tol, double* u, double* v, const std::string& method);
	};

};

#endif
