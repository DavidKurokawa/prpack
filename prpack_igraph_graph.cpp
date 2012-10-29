#include "prpack_igraph_graph.h"
#include <cstdlib>
#include <cstring>

using namespace prpack;
using namespace std;

prpack_igraph_graph::prpack_igraph_graph(const igraph_t* g, const igraph_vector_t* weights,
		igraph_bool_t directed) {
    const igraph_bool_t treat_as_directed = igraph_is_directed(g) && directed;
    igraph_es_t es;
    igraph_eit_t eit;
	igraph_vector_t neis;
    long int i, j, eid, sum, temp;
    int* p_head;
    double* p_weight;

    // Get the number of vertices and edges. For undirected graphs, we add
    // an edge in both directions.
    num_vs = igraph_vcount(g);
    num_es = igraph_ecount(g);
	num_self_es = 0;
    if (!treat_as_directed) {
        num_es *= 2;
    }

    // Allocate memory for heads and tails
    p_head = heads = new int[num_es];
    tails = new int[num_vs];
    memset(tails, 0, num_vs * sizeof(tails[0]));

	if (treat_as_directed) {
		// Select all the edges and iterate over them by the target vertices
		es = igraph_ess_all(IGRAPH_EDGEORDER_TO);

		// Add the edges
		igraph_eit_create(g, es, &eit);
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
	} else {
		// Select all the edges and iterate over them by the target vertices
		igraph_vector_init(&neis, 0);

		if (weights != 0) {
			p_weight = vals = new double[num_es];
		}

		for (i = 0; i < num_vs; i++) {
			igraph_incident(g, &neis, i, IGRAPH_ALL);
			temp = igraph_vector_size(&neis);
			// TODO: should loop edges be added in both directions?
			for (j = 0; j < temp; j++) {
				*p_head = IGRAPH_OTHER(g, VECTOR(neis)[j], i);
				if (i == *p_head) {
					num_self_es++;
				}
				++p_head;
			}
			if (weights != 0) {
				for (j = 0; j < temp; j++) {
					*p_weight = VECTOR(*weights)[(long int)VECTOR(neis)[j]];
					++p_weight;
				}
			}
			tails[i] = temp;
		}

		igraph_vector_destroy(&neis);
	}

    // Finalize the tails vector
    for (i = 0, sum = 0; i < num_vs; ++i) {
        temp = sum;
        sum += tails[i];
        tails[i] = temp;
    }

	// Debug
    /*
	printf("Heads:");
	for (i = 0; i < num_es; ++i) {
		printf(" %d", heads[i]);
	}
	printf("\n");
	printf("Tails:");
	for (i = 0; i < num_vs; ++i) {
		printf(" %d", tails[i]);
	}
	printf("\n");
	if (vals) {
		printf("Vals:");
		for (i = 0; i < num_es; ++i) {
			printf(" %.4f", vals[i]);
		}
		printf("\n");
	}
	printf("===========================\n");
    */
}

