#ifndef PRPACK_PREPROCESSED_SCHUR_GRAPH
#define PRPACK_PREPROCESSED_SCHUR_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_base_graph.h"

namespace prpack {

    class prpack_preprocessed_schur_graph : public prpack_preprocessed_graph {
        private:
            // helper methods
            void initialize();
            void initialize_weighted(prpack_base_graph* bg);
            void initialize_unweighted(prpack_base_graph* bg);
        public:
            // instance variables
            int num_no_in_vs;
            int num_no_out_vs;
            int* heads;
            int* tails;
            double* vals;
            int* encoding;
            int* decoding;
            // constructors
            prpack_preprocessed_schur_graph(prpack_base_graph* bg);
            // destructor
            ~prpack_preprocessed_schur_graph();
    };

};

#endif
