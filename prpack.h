#ifndef PRPACK
#define PRPACK

#include "prpack_csr.h"
#include "prpack_edge_list.h"
#include "prpack_base_graph.h"
#include "prpack_solver.h"
#include "prpack_result.h"

/*
HIGH LEVEL:
Want to be able to give any of:
- list of edges
- CSR for out-edges
- CSR for in-edges
which can be given in as any of (which have all supertype prpack_graph):
- prpack_uint32_csr
- prpack_uint64_csr
- prpack_uint32_edge_list
- prpack_uint64_edge_list
and return a solver (prpack_solver) which contains a prpack_preprocessed_graph (subtype of prpack_graph) that can be used for any number of:
- alpha
- tol
- u
- v
and return a result set (prpack_result) indicating:
- x
- preprocess time
- compute time
- postprocess time
- total time
- number of iterations
*/

#endif
