#ifndef PRPACK_PREPROCESSED_SCHUR_GRAPH
#define PRPACK_PREPROCESSED_SCHUR_GRAPH
#include "prpack_preprocessed_graph.h"
#include "prpack_base_graph.h"
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#endif

namespace prpack {

    class prpack_preprocessed_schur_graph : public prpack_preprocessed_graph {
        private:
            // instance variable
            bool from_matlab;
            // method
            void initialize();
        public:
            // instance variables
            int num_no_in_vs;
            int num_no_out_vs;
            int* heads;
            int* tails;
            int* encoding;
            int* decoding;
            // constructors
            prpack_preprocessed_schur_graph(prpack_base_graph* bg);
#ifdef MATLAB_MEX_FILE
            prpack_preprocessed_schur_graph(const mxArray* a);
#endif
            // destructor
            ~prpack_preprocessed_schur_graph();
            // method
#ifdef MATLAB_MEX_FILE
            mxArray* to_matlab_array() const;
#endif
    };

};

#endif
