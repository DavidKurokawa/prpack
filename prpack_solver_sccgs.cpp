#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace prpack;
using namespace std;


/** Gauss-Seidel using strongly connected components.
 * Notes:
 *   If not weighted, then we store x[i] = "x[i]/outdegree" to 
 *   avoid additional arithmetic.  We don't do this for the weighted
 *   case because the adjustment may not be constant.
 */
prpack_result* prpack_solver::solve_via_scc_gs(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es_inside,
        const int* heads_inside,
        const int* tails_inside,
        const double* vals_inside,
        const int num_es_outside,
        const int* heads_outside,
        const int* tails_outside,
        const double* vals_outside,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* uv,
        const int num_comps,
        const int* divisions,
        const int* encoding,
        const int* decoding,
        const bool should_normalize) {
    prpack_result* ret = new prpack_result();
    const bool weighted = vals_inside != NULL;
    // initialize uv values
    const double uv_const = 1.0/num_vs;
    const int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? prpack_utils::permute(num_vs, uv, encoding) : &uv_const;
    // CHECK initialize the solution with one iteration of GS from x=0.
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        x[i] = uv[uv_exists*i]/(1 - alpha*ii[i])/((weighted) ? 1 : num_outlinks[i]);
    // create x_outside
    double* x_outside = new double[num_vs];
    // run Gauss-Seidel for (I - alpha*P)*x = uv
    ret->num_es_touched = 0;
    for (int comp_i = 0; comp_i < num_comps; ++comp_i) {
        const int start_comp = divisions[comp_i];
        const int end_comp = (comp_i + 1 != num_comps) ? divisions[comp_i + 1] : num_vs;
        const bool parallelize = end_comp - start_comp > 512;
        // initialize relevant x_outside values
        for (int i = start_comp; i < end_comp; ++i) {
            x_outside[i] = 0;
            const int start_j = tails_outside[i];
            const int end_j = (i + 1 != num_vs) ? tails_outside[i + 1] : num_es_outside;
            for (int j = start_j; j < end_j; ++j)
                x_outside[i] += x[heads_outside[j]]*((weighted) ? vals_outside[j] : 1.);
            ret->num_es_touched += end_j - start_j;
        }
        double err, c;
        do {
            int num_es_touched = 0;
            err = c = 0;
            if (parallelize) {
                // iterate through vertices
                #pragma omp parallel for firstprivate(c) reduction(+:err, num_es_touched) schedule(dynamic, 64)
                for (int i = start_comp; i < end_comp; ++i) {
                    double new_val = x_outside[i];
                    const int start_j = tails_inside[i];
                    const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
                    if (weighted) {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]]*vals_inside[j];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    } else {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]*num_outlinks[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i])/num_outlinks[i];
                    }
                    num_es_touched += end_j - start_j;
                }
            } else {
                for (int i = start_comp; i < end_comp; ++i) {
                    double new_val = x_outside[i];
                    const int start_j = tails_inside[i];
                    const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
                    if (weighted) {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]]*vals_inside[j];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    } else {
                        for (int j = start_j; j < end_j; ++j) {
                            // TODO: might want to use compensation summation for large: end_j - start_j
                            new_val += x[heads_inside[j]];
                        }
                        COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]*num_outlinks[i]), c);
                        x[i] = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i])/num_outlinks[i];
                    }
                    num_es_touched += end_j - start_j;
                }
            }
            // update iteration index
            ret->num_es_touched += num_es_touched;
        } while (err/(1 - alpha) >= tol*(end_comp - start_comp)/num_vs);
    }
    // undo num_outlinks transformation
    if (!weighted)
        for (int i = 0; i < num_vs; ++i)
            x[i] *= num_outlinks[i];
    // normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
    if (should_normalize)
        normalize(num_vs, x);
    // return results
    ret->x = prpack_utils::permute(num_vs, x, decoding);
    delete[] x;
    delete[] x_outside;
    if (uv_exists)
        delete[] uv;
    return ret;
}

prpack_result* prpack_solver::solve_via_scc_gs_uv(
        const double alpha,
        const double tol,
        const int num_vs,
        const int num_es_inside,
        const int* heads_inside,
        const int* tails_inside,
        const double* vals_inside,
        const int num_es_outside,
        const int* heads_outside,
        const int* tails_outside,
        const double* vals_outside,
        const double* ii,
        const double* d,
        const double* num_outlinks,
        const double* u,
        const double* v,
        const int num_comps,
        const int* divisions,
        const int* encoding,
        const int* decoding) {

    // TODO, check if anything is dangling!
            
    // solve uv = u
    prpack_result* ret_u = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            vals_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            vals_outside,
            ii,
            d,
            num_outlinks,
            u,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // solve uv = v
    prpack_result* ret_v = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            vals_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            vals_outside,
            ii,
            d,
            num_outlinks,
            v,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // combine u and v
    return combine_uv(num_vs, d, num_outlinks, encoding, alpha, ret_u, ret_v);
}
