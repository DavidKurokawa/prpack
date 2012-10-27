#include "prpack_igraph_graph.h"
#include <cstdlib>
#include <cstring>

using namespace prpack;
using namespace std;

prpack_igraph_graph::prpack_igraph_graph(igraph_t* g, const igraph_vector_t* weights) {
    const igraph_bool_t is_directed = igraph_is_directed(g);
    igraph_es_t es;
    igraph_eit_t eit;
    long int eid;
    int* p_head;
    double* p_weight;

    // Get the number of vertices and edges. For undirected graphs, we add
    // an edge in both directions.
    num_vs = igraph_vcount(g);
    num_es = igraph_ecount(g);
    if (!is_directed) {
        abort();      // TODO
        num_es *= 2;
    }

    // Allocate memory for heads and tails
    heads = new int[num_es];
    tails = new int[num_vs];
    memset(tails, 0, num_vs * sizeof(tails[0]));

    // Select all the edges and iterate over them by the target vertices
    es = igraph_ess_all(IGRAPH_EDGEORDER_TO);

    // Add the edges
    igraph_eit_create(g, es, &eit);
    p_head = heads;
    while (!IGRAPH_EIT_END(eit)) {
        eid = IGRAPH_EIT_GET(eit);
        IGRAPH_EIT_NEXT(eit);

        *p_head = IGRAPH_FROM(g, eid);
        ++p_head;
        ++tails[IGRAPH_TO(g, eid)];

        if (IGRAPH_FROM(g, eid) == IGRAPH_TO(g, eid)) {
            ++num_self_es;
        }
    }
    igraph_eit_destroy(&eit);

    // Add the weights (if any)
    if (weights != 0) {
        vals = new double[num_es];

        igraph_eit_create(g, es, &eit);
        p_weight = vals;
        while (!IGRAPH_EIT_END(eit)) {
            eid = IGRAPH_EIT_GET(eit);
            IGRAPH_EIT_NEXT(eit);

            *p_weight = VECTOR(*weights)[eid];
            ++p_weight;
        }
        igraph_eit_destroy(&eit);
    }

    // Finalize the tails vector
    for (int i = 0, sum = 0; i < num_vs; ++i) {
        const int temp = sum;
        sum += tails[i];
        tails[i] = temp;
    }
}

