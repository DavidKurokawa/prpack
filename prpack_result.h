#ifndef PRPACK_RESULT
#define PRPACK_RESULT
#include "mex.h"
#include <string>

namespace prpack {

	// Result class.
	class prpack_result {
		public:
			// instance variables
			int num_vs;
			int num_es;
			double* x;
			double read_time;
			double preprocess_time;
			double compute_time;
			long long num_es_touched;
			std::string method;
			int converged;
			// method
			mxArray* to_matlab_array() const;
			// destructor
			~prpack_result();
	};

};

#endif
