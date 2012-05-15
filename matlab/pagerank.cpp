#include "../prpack.h"
#include "utils.h"
#include "mex.h"
#include <string>
using namespace std;
using namespace prpack;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate number of inputs and outputs
    if (nrhs != 8)
        mexErrMsgTxt("Not enough input arguments.");
    if (nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_num_vs = prhs[0];
    const mxArray* raw_heads = prhs[1];
    const mxArray* raw_tails = prhs[2];
    const mxArray* raw_alpha = prhs[3];
    const mxArray* raw_tol = prhs[4];
    const mxArray* raw_u = prhs[5];
    const mxArray* raw_v = prhs[6];
    const mxArray* raw_method = prhs[7];
    // parse variables
    int num_vs = parse_num_vs(raw_num_vs);
    int* heads;
    int* tails;
    int num_es = parse_heads_tails(heads, tails, raw_heads, raw_tails);
    double alpha = parse_alpha(raw_alpha);
    double tol = parse_tol(raw_tol);
    double* u = parse_uv(num_vs, raw_u);
    double* v = parse_uv(num_vs, raw_v);
    string method = parse_method(raw_method);
    // compute pagerank
    prpack_edge_list g;
    g.num_vs = num_vs;
    g.num_es = num_es;
    g.heads = heads;
    g.tails = tails;
    prpack_solver solver(&g);
    prpack_result* res = solver.solve(alpha, tol, u, v, method);
    // return
    mxArray* ares = res->to_matlab_array();
    plhs[0] = mxGetField(ares, 0, "x");
    if (nlhs >= 2)
        plhs[1] = ares;
    delete res;
}
