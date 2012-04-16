#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <string>
#include <vector>

namespace prpack {

	class prpack_adjacency_list {
		public:
			// instance variables
			int num_vs;
			int num_es;
			std::vector<int>* al;
			// constructors
			prpack_adjacency_list(prpack_csr* g);
			prpack_adjacency_list(prpack_edge_list* g);
			prpack_adjacency_list(const std::string& filename);
			prpack_adjacency_list(int num_vs_, std::vector<int>* al_);
			prpack_adjacency_list(int nverts, int nedges, 
			    std::pair<int,int>* edges);
	};

};

#endif
