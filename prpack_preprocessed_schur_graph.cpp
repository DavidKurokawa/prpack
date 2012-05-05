#include "prpack_preprocessed_schur_graph.h"
#include "prpack_utils.h"
#include <algorithm>
#include <cstring>
using namespace prpack;
using namespace std;

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(prpack_base_graph* bg) {
	// initialize instance variables
	num_vs = bg->num_vs;
	num_es = bg->num_es - bg->num_self_es;
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	for (int i = 0; i < bg->num_es; ++i)
		++inv_num_outlinks[bg->heads[i]];
	// permute no-inlink vertices to the beginning, and no-outlink vertices to the end
	encoding = new int[num_vs];
	decoding = new int[num_vs];
	num_no_in_vs = num_no_out_vs = 0;
	for (int i = 0; i < num_vs; ++i) {
		if (bg->tails[i] == ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es)) {
			decoding[encoding[i] = num_no_in_vs] = i;
			++num_no_in_vs;
		} else if (inv_num_outlinks[i] == 0) {
			decoding[encoding[i] = num_vs - 1 - num_no_out_vs] = i;
			++num_no_out_vs;
		}
	}
	// permute everything else
	for (int i = 0, p = num_no_in_vs; i < num_vs; ++i)
		if (bg->tails[i] < ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es) && inv_num_outlinks[i] > 0)
			decoding[encoding[i] = p++] = i;
	// permute inv_num_outlinks
	ii = inv_num_outlinks;
	inv_num_outlinks = new double[num_vs];
	for (int i = 0; i < num_vs; ++i)
		inv_num_outlinks[encoding[i]] = (ii[i] == 0) ? -1 : 1/ii[i];
	// convert bg to head/tail format
	tails = new int[num_vs];
	heads = new int[num_es];
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		ii[tails_i] = 0;
		tails[tails_i] = heads_i;
		int decoded = decoding[tails_i];
		int start_i = bg->tails[decoded];
		int end_i = (decoded + 1 != num_vs) ? bg->tails[decoded + 1] : bg->num_es;
		for (int i = start_i; i < end_i; ++i) {
			if (decoded == bg->heads[i])
				++ii[tails_i];
			else
				heads[heads_i++] = encoding[bg->heads[i]];
		}
		if (ii[tails_i] > 0)
			ii[tails_i] *= inv_num_outlinks[tails_i];
	}
}

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(const mxArray* a) {
	// separate raw matlab arrays
	mxArray* raw_num_vs = mxGetField(a, 0, "num_vs");
	mxArray* raw_num_es = mxGetField(a, 0, "num_es");
	mxArray* raw_ii = mxGetField(a, 0, "ii");
	mxArray* raw_inv_num_outlinks = mxGetField(a, 0, "inv_num_outlinks");
	mxArray* raw_heads = mxGetField(a, 0, "heads");
	mxArray* raw_tails = mxGetField(a, 0, "tails");
	mxArray* raw_num_no_in_vs = mxGetField(a, 0, "num_no_in_vs");
	mxArray* raw_num_no_out_vs = mxGetField(a, 0, "num_no_out_vs");
	mxArray* raw_encoding = mxGetField(a, 0, "encoding");
	mxArray* raw_decoding = mxGetField(a, 0, "decoding");
	// initialize instance variables
	num_vs = prpack_utils::matlab_array_to_int(raw_num_vs);
	num_es = prpack_utils::matlab_array_to_int(raw_num_es);
	ii = prpack_utils::matlab_array_to_double_array(raw_ii);
	inv_num_outlinks = prpack_utils::matlab_array_to_double_array(raw_inv_num_outlinks);
	heads = prpack_utils::matlab_array_to_int_array(raw_heads);
	tails = prpack_utils::matlab_array_to_int_array(raw_tails);
	num_no_in_vs = prpack_utils::matlab_array_to_int(raw_num_no_in_vs);
	num_no_out_vs = prpack_utils::matlab_array_to_int(raw_num_no_out_vs);
	encoding = prpack_utils::matlab_array_to_int_array(raw_encoding);
	decoding = prpack_utils::matlab_array_to_int_array(raw_decoding);
}

prpack_preprocessed_schur_graph::~prpack_preprocessed_schur_graph() {
	delete[] ii;
	delete[] inv_num_outlinks;
	delete[] heads;
	delete[] tails;
	delete[] encoding;
	delete[] decoding;
}

mxArray* prpack_preprocessed_schur_graph::to_matlab_array() const {
    const int num_fields = 10;
    const char* field_names[num_fields] = {"num_vs", "num_es", "ii", "inv_num_outlinks", "heads", "tails", "num_no_in_vs", "num_no_out_vs", "encoding", "decoding"};
    mxArray* ret = mxCreateStructMatrix(1, 1, num_fields, field_names);    
    mxSetField(ret, 0, "num_vs", prpack_utils::int_to_matlab_array(num_vs));
    mxSetField(ret, 0, "num_es", prpack_utils::int_to_matlab_array(num_es));
    mxSetField(ret, 0, "ii", prpack_utils::double_array_to_matlab_array(num_vs, ii));
    mxSetField(ret, 0, "inv_num_outlinks", prpack_utils::double_array_to_matlab_array(num_vs, inv_num_outlinks));
    mxSetField(ret, 0, "heads", prpack_utils::int_array_to_matlab_array(num_es, heads));
    mxSetField(ret, 0, "tails", prpack_utils::int_array_to_matlab_array(num_es, tails));
    mxSetField(ret, 0, "num_no_in_vs", prpack_utils::int_to_matlab_array(num_no_in_vs));
    mxSetField(ret, 0, "num_no_out_vs", prpack_utils::int_to_matlab_array(num_no_out_vs));
    mxSetField(ret, 0, "encoding", prpack_utils::int_array_to_matlab_array(num_vs, encoding));
    mxSetField(ret, 0, "decoding", prpack_utils::int_array_to_matlab_array(num_vs, decoding));
	return ret;
}
