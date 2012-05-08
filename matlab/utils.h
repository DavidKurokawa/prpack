#ifndef UTILS
#define UTILS
#include "mex.h"
#include <string>

// validation methods
bool is_int_scalar(const mxArray* a);
bool is_real_scalar(const mxArray* a);
bool is_vector(const mxArray* a);
bool is_int_vector(const mxArray* a);
bool is_real_vector(const mxArray* a);
bool is_string(const mxArray* a);

// parsing methods
int parse_num_vs(const mxArray* raw_num_vs);
int parse_heads_tails(int*& heads, int*& tails, const mxArray* raw_heads, const mxArray* raw_tails);
double parse_alpha(const mxArray* raw_alpha);
double parse_tol(const mxArray* raw_tol);
double* parse_uv(int num_vs, const mxArray* raw_uv);
std::string parse_method(const mxArray* raw_method);

#endif
