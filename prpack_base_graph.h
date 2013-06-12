#ifndef PRPACK_ADJACENCY_LIST
#define PRPACK_ADJACENCY_LIST
#include "prpack_csc.h"
#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include <cstdio>
#include <utility>

namespace prpack {

    typedef int prpack_vtype;

    /** Convert any input graph into something prpack can process itself.
     * 
     * prpack_base_graph is an in-links based representation of a graph.
     * tails[i]:tails[i+1] is an index into the heads array for all
     * edges that end at vertex i.
     *
     * for (int i=0; i < num_vs; ++i) {
     *   int end = i==num_vs-1 ? num_es : tails[i+1]; // handle end
     *   for (int j=tails[i]; j < end; ++j) {
     *       // there is a directed edge from heads[j] -> i
     *    }
     *  }
     *
     * Thus the graph is an in-edge adjacency list. 
     *
     */
    
    class prpack_base_graph {
        private:
            // helper methods
            void initialize();
            bool read_smat(std::FILE* f, const bool weighted);
            void read_edges(std::FILE* f);
            void read_ascii(std::FILE* f);
        public:
            // instance variables
            prpack_vtype num_vs;
            prpack_vtype num_es;
            prpack_vtype num_self_es;
            prpack_vtype* heads;
            prpack_vtype* tails;
            double* vals;
            // constructors
            prpack_base_graph();    // only to support inheritance
            prpack_base_graph(const prpack_csc* g);
            prpack_base_graph(const prpack_int64_csc* g);
            prpack_base_graph(const prpack_csr* g);
            prpack_base_graph(const prpack_edge_list* g);
            prpack_base_graph(const char* filename, const char* format, const bool weighted);
            prpack_base_graph(prpack_vtype nverts, prpack_vtype nedges, 
                std::pair<prpack_vtype,prpack_vtype>* edges);
            // destructor
            ~prpack_base_graph();
            // operations
            void normalize_weights();
    };

};

#endif
