#ifndef PRPACK_UTILS
#define PRPACK_UTILS
#include "mex.h"
#include <string>

// Computes the time taken to do X and stores it in T.
#define TIME(T, X)					\
	(T) = prpack_utils::get_time();	\
    (X);							\
    (T) = prpack_utils::get_time() - (T)

// Computes S += A using C as a carry-over.
// This is a macro over a function as it is faster this way.
#define COMPENSATED_SUM(S, A, C)						\
	double compensated_sum_y = (A) - (C);				\
    double compensated_sum_t = (S) + compensated_sum_y;	\
    (C) = compensated_sum_t - (S) - compensated_sum_y;	\
    (S) = compensated_sum_t

namespace prpack {

	class prpack_utils {
		public:
			static double get_time();
			static void validate(bool condition, const std::string& msg);
			static double* permute(int length, double* a, int* coding);
			static mxArray* int_to_matlab_array(int x);
			static int matlab_array_to_int(mxArray* a);
			static mxArray* int_array_to_matlab_array(int length, int* a);
			static int* matlab_array_to_int_array(mxArray* a);
			static mxArray* double_to_matlab_array(double x);
			static double matlab_array_to_double(mxArray* a);
			static mxArray* double_array_to_matlab_array(int length, double* a);            
			static double* matlab_array_to_double_array(mxArray* a);
			static mxArray* ll_to_matlab_array(long long x);
			static mxArray* string_to_matlab_array(const std::string& s);
			static mxArray* empty_matlab_array();
	};

};

#endif

