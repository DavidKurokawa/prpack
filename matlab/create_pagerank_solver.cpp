#include "../prpack.h"
#include "mex.h"
#include <algorithm>
#include <string>
using namespace std;
using namespace prpack;

bool is_int_scalar(const mxArray* a) {
    return (mxIsInt32(a) || mxIsInt64(a))
            && mxGetNumberOfElements(a) == 1;
}

bool is_real_scalar(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && mxGetNumberOfElements(a) == 1;
}

bool is_vector(const mxArray* a) {
    const mwSize* dims = mxGetDimensions(a);
    return mxGetNumberOfDimensions(a) == 2
            && min(dims[0], dims[1]) <= 1;
}

bool is_int_vector(const mxArray* a) {
    return (mxIsInt32(a) || mxIsInt64(a))
            && is_vector(a);
}

bool is_real_vector(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && is_vector(a);
}

bool is_string(const mxArray* a) {
    return mxIsChar(a)
            && is_vector(a);
}

void mexFunction(
		int nlhs,
		mxArray *plhs[],
		int nrhs,
		const mxArray *prhs[]) {
	// validate number of inputs and outputs
    if (nrhs != 3)
        mexErrMsgTxt("Not enough input arguments.");
	if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_num_vs = prhs[0];
    const mxArray* raw_heads = prhs[1];
    const mxArray* raw_tails = prhs[2];
    // parse num_vs
    if (!is_int_scalar(raw_num_vs))
        mexErrMsgTxt("num_vs must be an int.");
    int num_vs = *((int*) mxGetData(raw_num_vs));
    if (num_vs <= 0)
        mexErrMsgTxt("num_vs must be > 0.");
    // parse heads and tails
    if (!is_int_vector(raw_heads) || !is_int_vector(raw_tails))
        mexErrMsgTxt("heads and tails must be int vectors.");
    if (mxGetNumberOfElements(raw_heads) != mxGetNumberOfElements(raw_tails))
        mexErrMsgTxt("heads and tails must be of the same size.");
    int num_es = (int) mxGetNumberOfElements(raw_heads);
    int* heads = (int*) mxGetData(raw_heads);
    int* tails = (int*) mxGetData(raw_tails);
    // create pagerank solver
    prpack_edge_list g;
    g.num_vs = num_vs;
    g.num_es = num_es;
    g.heads = heads;
    g.tails = tails;
    prpack_solver solver(&g);
    // return the pagerank solver
    plhs[0] = solver.to_matlab_array();
}

