#include "../prpack.h"
#include "utils.h"
#include "mex.h"
using namespace std;
using namespace prpack;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
	// validate number of inputs and outputs
    if (nrhs != 3)
        mexErrMsgTxt("Not enough input arguments.");
	if (nlhs > 1)
		mexErrMsgTxt("Too many output arguments.");
    // set up raw variables
    const mxArray* raw_num_vs = prhs[0];
    const mxArray* raw_heads = prhs[1];
    const mxArray* raw_tails = prhs[2];
    // parse variables
    int num_vs = parse_num_vs(raw_num_vs);
    int* heads;
    int* tails;
    int num_es = parse_heads_tails(heads, tails, raw_heads, raw_tails);
    // create pagerank solver
    prpack_edge_list g;
    g.num_vs = num_vs;
    g.num_es = num_es;
    g.heads = heads;
    g.tails = tails;
    prpack_solver* solver = new prpack_solver(&g);
    // return the pagerank solver
    plhs[0] = solver->to_matlab_array();
}
