#ifndef PRPACK_RESULT
#define PRPACK_RESULT
#include <string>

namespace prpack {

	// Result class.
	struct prpack_result {
		int num_vs;
		int num_es;
		double* x;
		double preprocess_time;
		double compute_time;
		int num_iter;
		std::string method;
	};

};

#endif
