#ifndef PRPACK_SOLVER
#define PRPACK_SOLVER
#include "prpack_base_graph.h"
#include "prpack_csc.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include "prpack_preprocessed_ge_graph.h"
#include "prpack_preprocessed_gs_graph.h"
#include "prpack_preprocessed_scc_graph.h"
#include "prpack_preprocessed_schur_graph.h"
#include "prpack_result.h"

// TODO Make this a user configurable variable
#define PRPACK_SOLVER_MAX_ITERS 1000000

namespace prpack {

    // Solver class.
    class prpack_solver {
        private:
            // instance variables
            double read_time;
            prpack_base_graph* bg;
            prpack_preprocessed_ge_graph* geg;
            prpack_preprocessed_gs_graph* gsg;
            prpack_preprocessed_schur_graph* sg;
            prpack_preprocessed_scc_graph* sccg;
            // methods
            void initialize();
            prpack_result* solve_via_ge(
                    const double alpha,
                    const double tol,
                    const int num_vs,
                    const double* matrix,
                    const double* uv);
            prpack_result* solve_via_ge_uv(
                    const double alpha,
                    const double tol,
                    const int num_vs,
                    const double* matrix,
                    const double* d,
                    const double* u,
                    const double* v);
            static prpack_result* solve_via_gs(
                    double alpha,
                    double tol,
                    int num_vs,
                    int num_es,
                    int* heads,
                    int* tails,
                    double* vals,
                    double* ii,
                    double* d,
                    double* inv_num_outlinks,
                    double* u,
                    double* v);
            static prpack_result* solve_via_gs_err(
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
                    int num_no_in_vs,
                    int num_no_out_vs,
                    int num_es,
                    int* heads,
                    int* tails,
                    double* vals,
                    double* ii,
                    double* d,
                    double* inv_num_outlinks,
                    double* uv,
                    int* encoding,
                    int* decoding,
                    bool should_normalize = true);
            static prpack_result* solve_via_schur_gs_uv(
                    double alpha,
                    double tol,
                    int num_vs,
                    int num_no_in_vs,
                    int num_no_out_vs,
                    int num_es,
                    int* heads,
                    int* tails,
                    double* vals,
                    double* ii,
                    double* d,
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
                    double* vals_inside,
                    int num_es_outside,
                    int* heads_outside,
                    int* tails_outside,
                    double* vals_outside,
                    double* ii,
                    double* d,
                    double* inv_num_outlinks,
                    double* uv,
                    int num_comps,
                    int* divisions,
                    int* encoding,
                    int* decoding,
                    bool should_normalize = true);
            static prpack_result* solve_via_scc_gs_uv(
                    double alpha,
                    double tol,
                    int num_vs,
                    int num_es_inside,
                    int* heads_inside,
                    int* tails_inside,
                    double* vals_inside,
                    int num_es_outside,
                    int* heads_outside,
                    int* tails_outside,
                    double* vals_outside,
                    double* ii,
                    double* d,
                    double* inv_num_outlinks,
                    double* u,
                    double* v,
                    int num_comps,
                    int* divisions,
                    int* encoding,
                    int* decoding);
            static void ge(const int sz, double* A, double* b);
            static void normalize(const int length, double* x);
            static prpack_result* combine_uv(
                    int num_vs,
                    double* d,
                    double* inv_num_outlinks,
                    int* encoding,
                    double alpha,
                    prpack_result* ret_u,
                    prpack_result* ret_v);
        public:
            // constructors
            prpack_solver(prpack_csc* g);
            prpack_solver(prpack_int64_csc* g);
            prpack_solver(prpack_csr* g);
            prpack_solver(prpack_edge_list* g);
            prpack_solver(prpack_base_graph* g);
            prpack_solver(const char* filename, const char* format, bool weighted);
            // destructor
            ~prpack_solver();
            // methods
            int get_num_vs();
            prpack_result* solve(double alpha, double tol, const char* method);
            prpack_result* solve(double alpha, double tol, double* u, double* v, const char* method);
    };

};

#endif
