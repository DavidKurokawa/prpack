#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <cstdio>
#include <string>

namespace prpack {

	class prpack_adjacency_list {
		private:
			// helper methods
			void read_smat(std::FILE* f);
			void read_edges(std::FILE* f);
			void read_ascii(std::FILE* f);
		public:
			// instance variables
			int num_vs;
			int num_es;
			int num_self_es;
			int* heads;
			int* tails;
			// constructors
			prpack_adjacency_list(prpack_csr* g);
			prpack_adjacency_list(prpack_edge_list* g);
			prpack_adjacency_list(const std::string& filename, const std::string& format);
	};

};

#endif
