#ifndef PRPACK_UTILS
#define PRPACK_UTILS

#define TIME(T, X) (T) = prpack_utils::get_time(); (X); (T) = prpack_utils::get_time() - (T)

namespace prpack {

	class prpack_utils {
		public:
			static double get_time();
	};

};

#endif

