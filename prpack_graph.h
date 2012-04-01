#ifndef PRPACK_GRAPH
#define PRPACK_GRAPH

namespace prpack {

	// TODO: this class should not be seeable by the users of the library.
	// Super graph class.
	class prpack_graph {
		public:
			int num_vs;
			int num_es;
			int* heads;
			int* tails;
	};

};
	
#endif
