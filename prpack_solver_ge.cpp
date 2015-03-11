#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
#include <cstdlib>
#include <cstring>
using namespace prpack;
using namespace std;


prpack_result* prpack_solver::solve_via_ge(
        const double alpha,
        const double tol,
        const int num_vs,
        const double* matrix,
        const double* uv) {
    prpack_result* ret = new prpack_result();
    // initialize uv values
    const double uv_const = 1.0/num_vs;
    const int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? uv : &uv_const;
    // create matrix A
    double* A = new double[num_vs*num_vs];
    for (int i = 0; i < num_vs*num_vs; ++i)
        A[i] = -alpha*matrix[i];
    for (int i = 0; i < num_vs*num_vs; i += num_vs + 1)
        ++A[i];
    // create vector b
    double* b = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        b[i] = uv[uv_exists*i];
    // solve and normalize
    ge(num_vs, A, b);
    normalize(num_vs, b);
    // clean up and return
    delete[] A;
    ret->num_es_touched = -1;
    ret->x = b;
    return ret;
}

prpack_result* prpack_solver::solve_via_ge_uv(
        const double alpha,
        const double tol,
        const int num_vs,
        const double* matrix,
        const double* d,
        const double* u,
        const double* v) {
    prpack_result* ret = new prpack_result();
    // initialize u and v values
    const double u_const = 1.0/num_vs;
    const double v_const = 1.0/num_vs;
    const int u_exists = (u) ? 1 : 0;
    const int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // create matrix A
    double* A = new double[num_vs*num_vs];
    for (int i = 0; i < num_vs*num_vs; ++i)
        A[i] = -alpha*matrix[i];
    for (int i = 0, inum_vs = 0; i < num_vs; ++i, inum_vs += num_vs)
        for (int j = 0; j < num_vs; ++j)
            A[inum_vs + j] -= alpha*u[u_exists*i]*d[j];
    for (int i = 0; i < num_vs*num_vs; i += num_vs + 1)
        ++A[i];
    // create vector b
    double* b = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        b[i] = (1 - alpha)*v[v_exists*i];
    // solve
    ge(num_vs, A, b);
    // clean up and return
    delete[] A;
    ret->num_es_touched = -1;
    ret->x = b;
    return ret;
}

// Run Gaussian-Elimination (note: this changes A and returns the solution in b)
void prpack_solver::ge(const int sz, double* A, double* b) {
    // put into triangular form
    for (int i = 0, isz = 0; i < sz; ++i, isz += sz)
        for (int k = 0, ksz = 0; k < i; ++k, ksz += sz)
            if (A[isz + k] != 0) {
                const double coeff = A[isz + k]/A[ksz + k];
                A[isz + k] = 0;
                for (int j = k + 1; j < sz; ++j)
                    A[isz + j] -= coeff*A[ksz + j];
                b[i] -= coeff*b[k];
            }
    // backwards substitution
    for (int i = sz - 1, isz = (sz - 1)*sz; i >= 0; --i, isz -= sz) {
        for (int j = i + 1; j < sz; ++j)
            b[i] -= A[isz + j]*b[j];
        b[i] /= A[isz + i];
    }
}
