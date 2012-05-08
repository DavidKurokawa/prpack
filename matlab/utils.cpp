#include <algorithm>
#include "utils.h"
using namespace std;

bool is_int_scalar(const mxArray* a) {
    return mxIsInt32(a) && mxGetNumberOfElements(a) == 1;
}

bool is_double_scalar(const mxArray* a) {
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
    return mxIsInt32(a) && is_vector(a);
}

bool is_double_vector(const mxArray* a) {
    return mxIsDouble(a)
            && !mxIsComplex(a)
            && is_vector(a);
}

bool is_string(const mxArray* a) {
    return mxIsChar(a) && is_vector(a);
}

int parse_num_vs(const mxArray* raw_num_vs) {
    if (!is_int_scalar(raw_num_vs))
        mexErrMsgTxt("num_vs must be an int.");
    int num_vs = *((int*) mxGetData(raw_num_vs));
    if (num_vs <= 0)
        mexErrMsgTxt("num_vs must be > 0.");
    return num_vs;
}

int parse_heads_tails(int*& heads, int*& tails, const mxArray* raw_heads, const mxArray* raw_tails) {
    if (!is_int_vector(raw_heads) || !is_int_vector(raw_tails))
        mexErrMsgTxt("heads and tails must be int vectors.");
    if (mxGetNumberOfElements(raw_heads) != mxGetNumberOfElements(raw_tails))
        mexErrMsgTxt("heads and tails must be of the same size.");
    heads = (int*) mxGetData(raw_heads);
    tails = (int*) mxGetData(raw_tails);
    return mxGetNumberOfElements(raw_heads);
}

double parse_alpha(const mxArray* raw_alpha) {
    if (!is_double_scalar(raw_alpha))
        mexErrMsgTxt("alpha must be a real scalar.");
    double alpha = *((double*) mxGetData(raw_alpha));
    if (alpha <= 0 || 1 <= alpha)
        mexErrMsgTxt("alpha must be in (0, 1).");
    return alpha;
}

double parse_tol(const mxArray* raw_tol) {
    if (!is_double_scalar(raw_tol))
        mexErrMsgTxt("tol must be a real scalar.");
    double tol = *((double*) mxGetData(raw_tol));
    if (tol <= 0)
        mexErrMsgTxt("tol must be > 0.");
    return tol;
}

double* parse_uv(int num_vs, const mxArray* raw_uv) {
    if (!is_double_vector(raw_uv))
        mexErrMsgTxt("u and v must be real vectors.");
    mwSize uv_size = mxGetNumberOfElements(raw_uv);
    if (uv_size != 0 && uv_size != num_vs)
        mexErrMsgTxt("u and v must be the same size as the matrix, or empty.");
    return (uv_size == 0) ? NULL : (double*) mxGetPr(raw_uv);
}

string parse_method(const mxArray* raw_method) {
    if (!is_string(raw_method))
        mexErrMsgTxt("method must be a string");
    mwSize method_length = mxGetNumberOfElements(raw_method);
    char* s = new char[method_length + 1];
    mxGetString(raw_method, s, method_length + 1);
    string method(s);
    delete[] s;
    return method;
}
