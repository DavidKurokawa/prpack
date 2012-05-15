#ifndef UTILS
#define UTILS
#include "../prpack.h"
#include "mex.h"
#include <string>

// validation methods
bool is_double_scalar(const mxArray* a);
bool is_vector(const mxArray* a);
bool is_double_vector(const mxArray* a);
bool is_string(const mxArray* a);

// parsing methods
prpack::prpack_solver* parse_solver(const mxArray* raw_solver_ptr);
double parse_alpha(const mxArray* raw_alpha);
double parse_tol(const mxArray* raw_tol);
double* parse_uv(int num_vs, const mxArray* raw_uv);
std::string parse_method(const mxArray* raw_method);

// matlab conversion methods
mxArray* int_to_matlab_array(int x);
mxArray* double_to_matlab_array(double x);
mxArray* double_array_to_matlab_array(int length, double* a);
mxArray* ll_to_matlab_array(long long x);
mxArray* string_to_matlab_array(const std::string& s);
mxArray* solver_to_matlab_array(prpack::prpack_solver* solver);
mxArray* result_to_matlab_array(prpack::prpack_result* res);

#endif
