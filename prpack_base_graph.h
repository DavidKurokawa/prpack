#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include "mex.h"
#include <cstdio>
#include <string>

namespace prpack {

	class prpack_base_graph {
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
			prpack_base_graph(prpack_csr* g);
			prpack_base_graph(prpack_edge_list* g);
			prpack_base_graph(const std::string& filename, const std::string& format);
			prpack_base_graph(int nverts, int nedges, 
					std::pair<int,int>* edges);
			prpack_base_graph(const mxArray* a);
			// destructor
			~prpack_base_graph();
			// method
			mxArray* to_matlab_array() const;
	};

};

#endif
