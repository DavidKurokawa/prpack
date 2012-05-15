#include "prpack_result.h"
using namespace prpack;

prpack_result::prpack_result() {
    x = NULL;
}

prpack_result::~prpack_result() {
    delete[] x;
}

