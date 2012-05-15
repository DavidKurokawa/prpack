#ifndef PRPACK_CSC
#define PRPACK_CSC

namespace prpack {

    class prpack_csc {
        public:
            int num_vs;
            int num_es;
            int* heads;
            int* tails;
    };

};

#endif
