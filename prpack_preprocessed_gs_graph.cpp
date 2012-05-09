#include "prpack_preprocessed_gs_graph.h"
#include "prpack_utils.h"
#include <algorithm>
using namespace prpack;
using namespace std;

void prpack_preprocessed_gs_graph::initialize() {
#ifdef MATLAB_MEX_FILE
    from_matlab = false;
#endif
    heads = NULL;
    tails = NULL;
    ii = NULL;
    inv_num_outlinks = NULL;
}

prpack_preprocessed_gs_graph::prpack_preprocessed_gs_graph(prpack_base_graph* bg) {
    initialize();
	num_vs = bg->num_vs;
	num_es = bg->num_es - bg->num_self_es;
	heads = new int[num_es];
	tails = new int[num_vs];
	ii = new double[num_vs];
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		tails[tails_i] = heads_i;
		ii[tails_i] = 0;
		int start_j = bg->tails[tails_i];
		int end_j = (tails_i + 1 != num_vs) ? bg->tails[tails_i + 1]: bg->num_es;
		for (int j = start_j; j < end_j; ++j) {
			if (tails_i == bg->heads[j])
				++ii[tails_i];
			else
				heads[heads_i++] = bg->heads[j];
			++inv_num_outlinks[bg->heads[j]];
		}
	}
	for (int i = 0; i < num_vs; ++i) {
		inv_num_outlinks[i] = (inv_num_outlinks[i] == 0) ? -1 : 1/inv_num_outlinks[i];
		ii[i] *= inv_num_outlinks[i];
	}
}

#ifdef MATLAB_MEX_FILE
prpack_preprocessed_gs_graph::prpack_preprocessed_gs_graph(const mxArray* a) {
    initialize();
    from_matlab = true;
    // separate raw matlab arrays
    mxArray* raw_num_vs = mxGetField(a, 0, "num_vs");
    mxArray* raw_num_es = mxGetField(a, 0, "num_es");
    mxArray* raw_ii = mxGetField(a, 0, "ii");
    mxArray* raw_inv_num_outlinks = mxGetField(a, 0, "inv_num_outlinks");
    mxArray* raw_heads = mxGetField(a, 0, "heads");
    mxArray* raw_tails = mxGetField(a, 0, "tails");
    // initialize instance variables
    num_vs = prpack_utils::matlab_array_to_int(raw_num_vs);
    num_es = prpack_utils::matlab_array_to_int(raw_num_es);
    ii = prpack_utils::matlab_array_to_double_array(raw_ii);
    inv_num_outlinks = prpack_utils::matlab_array_to_double_array(raw_inv_num_outlinks);
    heads = prpack_utils::matlab_array_to_int_array(raw_heads);
    tails = prpack_utils::matlab_array_to_int_array(raw_tails);
}
#endif

prpack_preprocessed_gs_graph::~prpack_preprocessed_gs_graph() {
#ifdef MATLAB_MEX_FILE
    if (!from_matlab) {
        delete[] heads;
        delete[] tails;
        delete[] ii;
        delete[] inv_num_outlinks;
    }
#else
    delete[] heads;
    delete[] tails;
    delete[] ii;
    delete[] inv_num_outlinks;
#endif
}

#ifdef MATLAB_MEX_FILE
mxArray* prpack_preprocessed_gs_graph::to_matlab_array() const {
    const int num_fields = 6;
    const char* field_names[num_fields] = {"num_vs", "num_es", "ii", "inv_num_outlinks", "heads", "tails"};
    mxArray* ret = mxCreateStructMatrix(1, 1, num_fields, field_names);
    mxSetField(ret, 0, "num_vs", prpack_utils::int_to_matlab_array(num_vs));
    mxSetField(ret, 0, "num_es", prpack_utils::int_to_matlab_array(num_es));
    mxSetField(ret, 0, "ii", prpack_utils::double_array_to_matlab_array(num_vs, ii));
    mxSetField(ret, 0, "inv_num_outlinks", prpack_utils::double_array_to_matlab_array(num_vs, inv_num_outlinks));
    mxSetField(ret, 0, "heads", prpack_utils::int_array_to_matlab_array(num_es, heads));
    mxSetField(ret, 0, "tails", prpack_utils::int_array_to_matlab_array(num_vs, tails));
    return ret;
}
#endif
