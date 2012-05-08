#include "prpack_result.h"
#include "prpack_utils.h"
using namespace prpack;

#ifdef MATLAB_MEX_FILE
mxArray* prpack_result::to_matlab_array() const {    
	const int num_fields = 8;
	const char* field_names[num_fields] = {"num_vs", "num_es", "x", "read_time", "preprocess_time", "compute_time", "num_es_touched", "method"};
	mxArray* ret = mxCreateStructMatrix(1, 1, num_fields, field_names);
	mxSetField(ret, 0, "num_vs", prpack_utils::int_to_matlab_array(num_vs));
	mxSetField(ret, 0, "num_es", prpack_utils::int_to_matlab_array(num_es));
	mxSetField(ret, 0, "x", prpack_utils::double_array_to_matlab_array(num_vs, x));
	mxSetField(ret, 0, "read_time", prpack_utils::double_to_matlab_array(read_time));
	mxSetField(ret, 0, "preprocess_time", prpack_utils::double_to_matlab_array(preprocess_time));
	mxSetField(ret, 0, "compute_time", prpack_utils::double_to_matlab_array(compute_time));
	mxSetField(ret, 0, "num_es_touched", prpack_utils::ll_to_matlab_array(num_es_touched));
	mxSetField(ret, 0, "method", prpack_utils::string_to_matlab_array(method));
	return ret;
}
#endif

prpack_result::~prpack_result() {
	delete[] x;
}

