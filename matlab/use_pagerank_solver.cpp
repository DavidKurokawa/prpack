#include "../prpack.h"
#include "mex.h"
#include <algorithm>
#include <string>
using namespace std;
using namespace prpack;

// TODO: put these functions into prpack_utils, and have them visible for matlab

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
    // parse prpack_solver
    prpack_solver solver(raw_prpack_solver);
    int num_vs = *(int*) mxGetData(mxGetField(mxGetField(raw_prpack_solver, 0, "bg"), 0, "num_vs")); // TODO get rid of this and just use solver.num_vs, by making it final
    // parse alpha and tol
    if (!is_real_scalar(raw_alpha) || !is_real_scalar(raw_tol))
        mexErrMsgTxt("alpha and tol must be real scalars.");
    double alpha = *((double*) mxGetData(raw_alpha));
    double tol = *((double*) mxGetData(raw_tol));
    if (alpha <= 0 || 1 <= alpha)
        mexErrMsgTxt("alpha must be in (0, 1).");
    if (tol <= 0)
        mexErrMsgTxt("tol must be > 0.");
    // parse u and v
    if (!is_real_vector(raw_u) || !is_real_vector(raw_v))
        mexErrMsgTxt("u and v must be real vectors.");
    mwSize u_size = mxGetNumberOfElements(raw_u);
    mwSize v_size = mxGetNumberOfElements(raw_v);
    if ((u_size != 0 && u_size != num_vs) || (v_size != 0 && v_size != num_vs))
        mexErrMsgTxt("u and v must be the same size as the matrix, or empty.");
    double* u = (u_size == 0) ? NULL : (double*) mxGetPr(raw_u);
    double* v = (v_size == 0) ? NULL : (double*) mxGetPr(raw_v);
    // parse method
    if (!is_string(raw_method))
        mexErrMsgTxt("method must be a string");
    mwSize method_length = mxGetNumberOfElements(raw_method);
    char* s = new char[method_length + 1];
    mxGetString(raw_method, s, method_length + 1);
    string method(s);
    delete[] s;
    // compute pagerank
    prpack_result* res = solver.solve(alpha, tol, u, v, method);
    // return vector
    plhs[0] = mxCreateDoubleMatrix(num_vs, 1, mxREAL);
    double* ret = mxGetPr(plhs[0]);
    for (int i = 0; i < num_vs; ++i)
        ret[i] = res->x[i];
    plhs[1] = res->to_matlab_array();
    plhs[2] = solver.to_matlab_array();
    delete res;
    ///////////////////////////////////////////////////////////////////////
    // print all variables out
    /*
    mexPrintf("alpha = %f\n", alpha);
    mexPrintf("tol = %f\n", tol);
    mexPrintf("u =");
    for (int i = 0; i < u_size; ++i)
        mexPrintf(" %f", u[i]);
    mexPrintf("\n");
    mexPrintf("v =");
    for (int i = 0; i < v_size; ++i)
        mexPrintf(" %f", v[i]);
    mexPrintf("\n");
    mexPrintf("method = %s\n", method.c_str());
     */
    ///////////////////////////////////////////////////////////////////////
}

