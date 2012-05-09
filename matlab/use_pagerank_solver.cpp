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
	if (nlhs > 3)
		mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_prpack_solver = prhs[0];
    const mxArray* raw_alpha = prhs[1];
    const mxArray* raw_tol = prhs[2];
    const mxArray* raw_u = prhs[3];
    const mxArray* raw_v = prhs[4];
    const mxArray* raw_method = prhs[5];
    // parse variables
    prpack_solver solver(raw_prpack_solver);
    int num_vs = *(int*) mxGetData(mxGetField(mxGetField(raw_prpack_solver, 0, "bg"), 0, "num_vs")); // TODO get rid of this and just use solver.num_vs, by making it final
    double alpha = parse_alpha(raw_alpha);
    double tol = parse_tol(raw_tol);
    double* u = parse_uv(num_vs, raw_u);
    double* v = parse_uv(num_vs, raw_v);
    string method = parse_method(raw_method);
    // compute pagerank
    prpack_result* res = solver.solve(alpha, tol, u, v, method);
    // return vector
    plhs[0] = mxCreateDoubleMatrix(num_vs, 1, mxREAL);
    double* ret = mxGetPr(plhs[0]);
    for (int i = 0; i < num_vs; ++i)
        ret[i] = res->x[i];
    if (nlhs >= 2)
        plhs[1] = res->to_matlab_array();
    if (nlhs >= 3)
        plhs[2] = solver.to_matlab_array();
    delete res;
}
