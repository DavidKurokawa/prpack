#ifndef PRPACK_IGRAPH_GRAPH
#define PRPACK_IGRAPH_GRAPH

#include "igraph_interface.h"
#include "prpack_base_graph.h"

namespace prpack {

    class prpack_igraph_graph : prpack_base_graph {

        public:
            // constructors
            prpack_igraph_graph(igraph_t* g, const igraph_vector_t* weights = 0);
    };

};

#endif
