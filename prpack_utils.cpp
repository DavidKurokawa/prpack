#include "prpack_utils.h"
#include <sys/time.h>
using namespace prpack;
using namespace std;

double prpack_utils::get_time() {
	struct timeval t;
	gettimeofday(&t, 0);
	return t.tv_sec + t.tv_usec/double(1.0e6);
}

