#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <assert.h>
#include <vector>
using namespace prpack;
using namespace std;

/** Solve PageRank through an inner-outer iteration
 * @param alpha the value of alpha for pagerank (0 <= alpha < 1) 
 * @param tol the solution tolerance (0 < tol < 1)
 * @param num_vs the number of vertices
 * @param num_es the number of edges
 * @param heads a list of length num_es giving the list of sources
 *   for each edge (the "head" of the arrow)
 * @param tails a list of length num_vs giving the pointers into heads
 *   for the list of heads for each vertex
 * @param ii an indicator vector for self-loops,
 *    ii[i] = 0 if there are no self loops and ii[i] = 1/d[i] if there are.
 * @param degs the degree of each node (sum of weights)
 * @param vals (if not null, then a list of length num_es giving the value 
 *   associated with each edge)
 * @param valtype if 0, then assume we want to normalize by the values,
 *   if 1, then assume the values are already normalized.
 * @param u the dangling node adjustment vector
 * @param v the teleportation adjustment vector
 * @param beta the value of beta in the inner outer iteration 
 *          (0 <= beta <= alpha)
 * @param itol the inner tolerance for the inner-outer iteration 
 *          (0 < itol)
 * 
 */
prpack_result* prpack_solver::solve_via_inout(
    const double alpha,
    const double tol,
    const prpack_vtype num_vs,
    const prpack_vtype num_es,
    prpack_vtype* heads,
    prpack_vtype* tails,
    double* ii,
    double* degs,
    double* vals,
    int valtype,
    const double* u,
    const double* v,
    double beta, 
    double itol) {
        
    // verify assumptions
    assert(0 < alpha);
    assert(alpha < 1);
    assert(0 < tol);
    assert(tol < 1);
    assert(beta < alpha);
    assert(0 < beta);
    assert(0 < itol);
    
    prpack_result* ret = new prpack_result();
    
    // initialize some flags
    const bool weighted = vals != NULL;
    const prpack_vtype invtolhalf = (prpack_vtype)(1./(2.*tol)); 
        // used in checking for compensated summation
    
    // initialize u and v values
    const double u_const = 1.0/num_vs;
    const double v_const = 1.0/num_vs;
    const int u_exists = (u) ? 1 : 0;
    const int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    
    // initialize the solution vector
    double* x = new double[num_vs];
    double* y = new double[num_vs];
    double* f = new double[num_vs];
    
    for (int i = 0; i < num_vs; ++i) {
        x[i] = v[v_exists*i];
    }
    
    
    // the following analysis may not be entirely correct, but it gives
    // a general idea.
    // the idegs variable encodes out-weight information
    // if g is unweighted, then 
    //   idegs[i] = 1/deg[i] if deg[i] > 0 or -1 if deg[i] == 0
    // if g is weighted, then
    //   idegs[i] = 1/sum(weights from i) if sum(weights from i) > 0 or
    //   idegs[i] = -1 
    // if g is weighted and all sums of weights from i are less <= 1 for all i
    //   then idegs[i] <= 0 and idegs[i] = -(1-sum(weights from i))
    // thus we only ever divide by idegs if it's positive, otherwise,
    // we leave it alone.
    std::vector<double> idegs(num_vs, 0.); // allocate an array of degrees
    double dtx=0., dtx_c=0.;
    for (int i = 0; i < num_vs; ++i) {
        if (degs[i] > 0.) {
            idegs[i] = 1./degs[i];
            x[i] *= idegs[i];
        } else {
            idegs[i] = -1.; // we use this 
            COMPENSATED_SUM(dtx, x[i], dtx_c);
        }
    }

 
    ret->num_es_touched = 0;
      
    // right now, we just do the power method
    // it converges in 2*log(tol)/log(alpha) steps
    // given a starting vector of v
    double nitersd = 2*log(tol)/log(alpha);
    if (nitersd > 2.e9) { 
        ret->converged = 0;
        nitersd = 2.e9;
    } else {
        ret->converged = 1;
    }
    
    int niters = (int)nitersd;
    double delta = 2.;
    for (int iter=0; iter<niters; ++iter) {
        if (delta <= tol*(1-alpha)) {
            // we are done!
            break;
        }
        // mult y=alpha*P'*x, sumy = e'*y, dtx = d'*x
        double delta_c=0.,sumy=0.;
        double new_dtx = 0., new_dtx_c = 0.;
        delta = 0.; 
        if (weighted) {
            assert(false);
        } else { 
            #pragma omp parallel for schedule(dynamic,1000) \
                reduction(+:delta,delta_c,new_dtx,new_dtx_c,sumy)
            for (prpack_vtype i = 0; i < num_vs; ++i) {
                const prpack_vtype start_j = tails[i];
                const prpack_vtype end_j = 
                        (i + 1 != num_vs) ? tails[i + 1] : num_es;
                double new_val = 0., new_val_c = 0.;
                if (end_j - start_j >= invtolhalf) {
                    for (prpack_vtype j = start_j; j < end_j; ++j) {
                        COMPENSATED_SUM(new_val, x[heads[j]], new_val_c);
                    }
                } else {
                    for (prpack_vtype j = start_j; j < end_j; ++j) {
                        new_val += x[heads[j]];
                    }
                }
                
                y[i] = alpha*new_val + alpha*dtx*u[u_exists*i] \
                        + ((1-alpha)*v[v_exists*i]+new_val_c);
                
                sumy += y[i];
                if (idegs[i] < 0.) {
                    COMPENSATED_SUM(new_dtx, -idegs[i]*y[i], new_dtx_c);
                    { COMPENSATED_SUM(delta, fabs(y[i] - x[i]), delta_c); }
                } else {
                    // non-dangling, check for self-loops 
                    if (ii[i] > 0.) {
                        // scale by degs[i] because ii[i] is inversely scaled
                        y[i] += alpha*degs[i]*ii[i]*x[i];
                    }
                    { COMPENSATED_SUM(delta, fabs(y[i] - degs[i]*x[i]), delta_c); }
                    y[i] *= idegs[i];
                }
                
            }
            
            delta += delta_c; // add in the correction
            dtx = new_dtx + new_dtx_c;
            
            //fprintf(stderr, "iter %04i; delta = %8.1e sumy = %8.1e, dtx = %8.1e\n", 
            //   iter, delta, sumy, dtx);
        }
        ret->num_es_touched += num_es;
                
        // swap x and y;
        { double* temp = x; x = y; y = temp; }
    }
    
    // renormalize x to be a solution
    #pragma omp parallel for 
    for (prpack_vtype i = 0; i < num_vs; ++i) {
        if (degs[i] > 0.) {
            x[i] *= degs[i];
        }
    }
    
    ret->x = x;
    
    delete[] y;
    delete[] f;

    return ret;
}        
