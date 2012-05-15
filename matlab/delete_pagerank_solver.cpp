#include "../prpack.h"
#include "mex.h"
using namespace prpack;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    // validate number of inputs and outputs
    if (nrhs < 1)
        mexErrMsgTxt("Not enough input arguments.");
    if (nrhs > 1)
        mexErrMsgTxt("Too many input arguments.");
    if (nlhs > 0)
        mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_solver_ptr = prhs[0];
    // parse variables
    prpack_solver* solver = reinterpret_cast<prpack_solver*>(*(int*) mxGetData(raw_solver_ptr)); // TODO: handle if this is 32/64
    // delete pagerank solver
    delete solver;
}
