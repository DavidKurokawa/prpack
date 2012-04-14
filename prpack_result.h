#ifndef PRPACK_RESULT
#define PRPACK_RESULT
#include <string>

namespace prpack {

	// Result class.
	class prpack_result {
		public:
			// instance variables
			int num_vs;
			int num_es;
			double* x;
			double preprocess_time;
			double compute_time;
			int num_iter;
			std::string method;
			// destructor
			~prpack_result();
	};

};

#endif
