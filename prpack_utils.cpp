#include "prpack_utils.h"
#include <sys/time.h>
#include <cassert>
#include <iostream>
#include <string>
using namespace prpack;
using namespace std;

// Returns the current time.
double prpack_utils::get_time() {
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec + t.tv_usec/double(1.0e6);
}

// Fails and outputs 'msg' if 'condition' is false.
void prpack_utils::validate(bool condition, const string& msg) {
	if (!condition) {
		cerr << msg << endl;
		assert(condition);
	}
}

