#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csc.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <cstdio>
#include <string>
#include <utility>

namespace prpack {

    class prpack_base_graph {
        private:
            // helper methods
            void initialize();
            bool read_smat(std::FILE* f);
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
            prpack_base_graph(prpack_csc* g);
            prpack_base_graph(prpack_int64_csc* g);
            prpack_base_graph(prpack_csr* g);
            prpack_base_graph(prpack_edge_list* g);
            prpack_base_graph(const std::string& filename, const std::string& format);
            prpack_base_graph(int nverts, int nedges, std::pair<int,int>* edges);
            // destructor
            ~prpack_base_graph();
    };

};

#endif
