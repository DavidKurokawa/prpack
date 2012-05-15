#include "../prpack.h"
#include "utils.h"
#include "mex.h"
using namespace prpack;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate number of inputs and outputs
    if (nrhs < 1)
        mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 1)
        mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 1)
        mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_csc = prhs[0];
    // parse variables
    if (mxGetN(raw_csc) != mxGetM(raw_csc))
        mexErrMsgTxt("matrix must be square.");
    int num_vs = mxGetN(raw_csc);
    int* heads = mxGetJc(raw_csc);
    int* tails = mxGetIr(raw_csc);
    int num_es = heads[num_vs];
    // create pagerank solver
    prpack_csc g;
    g.num_vs = num_vs;
    g.num_es = num_es;
    g.heads = heads;
    g.tails = tails;
    prpack_solver* solver = new prpack_solver(&g);
    // return a pointer to the pagerank solver
    mxArray* ret = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL); // TODO handle 32/64
    *(int*) mxGetData(ret) = reinterpret_cast<int>(solver); // TODO handle 32/64
    plhs[0] = ret;
}
