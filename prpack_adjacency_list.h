#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <list>
#include <string>

namespace prpack {

	class prpack_adjacency_list {
		public:
			// instance variables
			int num_vs;
			int num_es;
			std::list<int>* al;
			// constructors
			prpack_adjacency_list(prpack_csr* g);
			prpack_adjacency_list(prpack_edge_list* g);
			prpack_adjacency_list(const std::string& filename);
	};

};

#endif