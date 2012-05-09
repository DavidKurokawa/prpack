#include "../prpack.h"
#include "utils.h"
#include "mex.h"
#include <string>
using namespace std;
using namespace prpack;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate number of inputs and outputs
    if (nrhs != 6)
        mexErrMsgTxt("Not enough input arguments.");
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_solver_wrapper = prhs[0];
    mxArray* raw_solver = mxGetField(raw_solver_wrapper, 0, "solver");
    const mxArray* raw_alpha = prhs[1];
    const mxArray* raw_tol = prhs[2];
    const mxArray* raw_u = prhs[3];
    const mxArray* raw_v = prhs[4];
    const mxArray* raw_method = prhs[5];
    // parse variables
    prpack_solver solver(raw_solver);
    int num_vs = *(int*) mxGetData(mxGetField(mxGetField(raw_solver, 0, "bg"), 0, "num_vs")); // TODO get rid of this and just use solver.num_vs, by making it final
    double alpha = parse_alpha(raw_alpha);
    double tol = parse_tol(raw_tol);
    double* u = parse_uv(num_vs, raw_u);
    double* v = parse_uv(num_vs, raw_v);
    string method = parse_method(raw_method);
    // compute pagerank
    prpack_result* res = solver.solve(alpha, tol, u, v, method);
    solver.to_matlab_array(raw_solver);
    // return
    mxArray* ares = res->to_matlab_array();
    plhs[0] = mxGetField(ares, 0, "x");
    if (nlhs >= 2)
        plhs[1] = ares;
    delete res;
}

