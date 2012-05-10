#include "prpack_utils.h"
#include <cassert>
#include <ctime>
#include <iostream>
#include <string>
using namespace prpack;
using namespace std;

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
double prpack_utils::get_time()
{
    LARGE_INTEGER t, freq;
    QueryPerformanceFrequency(&freq);
    QueryPerformanceCounter(&t);
    return (double)t.QuadPart / (double)t.freq;
}
#else
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
double prpack_utils::get_time()
{
    struct timeval t; gettimeofday(&t, 0);
    return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
}
#endif

// Fails and outputs 'msg' if 'condition' is false.
void prpack_utils::validate(bool condition, const string& msg) {
    if (!condition) {
        cerr << msg << endl;
        assert(condition);
    }
}

// Permute a vector.
double* prpack_utils::permute(int length, double* a, int* coding) {
    double* ret = new double[length];
    for (int i = 0; i < length; ++i)
        ret[coding[i]] = a[i];
    return ret;
}

#ifdef MATLAB_MEX_FILE
// Convert an int to a matlab array.
mxArray* prpack_utils::int_to_matlab_array(int x) {
    mxArray* ret = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int*) mxGetData(ret) = x;
    return ret;
}

// Convert a matlab array to an int.
int prpack_utils::matlab_array_to_int(mxArray* a) {
    return *((int*) mxGetData(a));
}

// Convert an int array to a matlab array.
mxArray* prpack_utils::int_array_to_matlab_array(int length, int* a) {
    // TODO: find a way to have this not be necessary
    mxArray* ret = mxCreateNumericMatrix(length, 1, mxINT32_CLASS, mxREAL);
    int* ret_data = (int*) mxGetData(ret);
    for (int i = 0; i < length; ++i)
        ret_data[i] = a[i];
    return ret;
}

// Convert a matlab array to an int array.
int* prpack_utils::matlab_array_to_int_array(mxArray* a) {
    return (int*) mxGetData(a);
}

// Convert a double to a matlab array.
mxArray* prpack_utils::double_to_matlab_array(double x) {
    return mxCreateDoubleScalar(x);
}

// Convert a matlab array to a double.
double prpack_utils::matlab_array_to_double(mxArray* a) {
    return *((double*) mxGetData(a));
}

// Convert a double array to a matlab array.
mxArray* prpack_utils::double_array_to_matlab_array(int length, double* a) {
    // TODO: find a way to have this not be necessary
    mxArray* ret = mxCreateDoubleMatrix(length, 1, mxREAL);
    double* ret_data = mxGetPr(ret);
    for (int i = 0; i < length; ++i)
        ret_data[i] = a[i];
    return ret;
}

// Convert a matlab array to a double array.
double* prpack_utils::matlab_array_to_double_array(mxArray* a) {
    return (double*) mxGetData(a);
}

// Convert a long long to a matlab array.
mxArray* prpack_utils::ll_to_matlab_array(long long x) {
    mxArray* ret = mxCreateNumericMatrix(1, 1, mxINT64_CLASS, mxREAL);
    *(long long*) mxGetData(ret) = x;
    return ret;
}

// Convert a string to a matlab array.
mxArray* prpack_utils::string_to_matlab_array(const string& s) {
    return mxCreateString(s.c_str());
}

// Return an empty matlab array.
mxArray* prpack_utils::empty_matlab_array() {
    return mxCreateDoubleMatrix(0, 0, mxREAL);
}
#endif
