
/** 
 * @file prpack_driver_benchmark
 * Implement a simple synthetic benchmark to evaluate scalability of 
 * the prpack codes.
 * @author David F. Gleich and Dave Kurakowa
 */

/** History
 *  :2012-04-15: Initial coding
 */
 
#include <math.h>
#include <assert.h>
#include <omp.h>
 
#include <vector>
#include <iostream>

#include "prpack_base_graph.h"
#include "prpack_solver.h"

using namespace std;

/* 
 * The following code contains a few different types of utility
 * codes I've written over the years.
 * They help with single-threaded random number generation.
 */

/* Setup using TR1 on various platforms */
#if defined(_WIN32) || defined(_WIN64)
  #pragma warning(disable:4996)
  #include <random>
  #define tr1ns std::tr1
#elif defined __GNUC__
  #define GCC_VERSION (__GNUC__ * 10000 \
                                                   + __GNUC_MINOR__ * 100 \
                                                   + __GNUC_PATCHLEVEL__)
  #if GCC_VERSION < 40600
    #include <tr1/random>
    #define tr1ns std::tr1
    #define uniform_real_distribution uniform_real
    #define uniform_int_distribution uniform_int
  #else
    #include <random>
    #define tr1ns std
  #endif
#else
  #include <random>
  #define tr1ns std  
#endif  

tr1ns::mt19937 sparfun_rand;

typedef tr1ns::mt19937 generator_t;
typedef tr1ns::uniform_real_distribution <double> distribution_t;
typedef tr1ns::variate_generator<generator_t, distribution_t> variate_t;
variate_t sparfun_rand_unif(sparfun_rand, distribution_t(0.0, 1.0));


/** Generate a uniform random number. */
double sf_rand(double min0, double max0) {
  tr1ns::uniform_real_distribution<double> dist(min0,max0);
  return dist(sparfun_rand_unif);
}


/** 
 * @param dist a cumulative probability distribution of size n
 *  where dist[n-1] = 1.0 and so dist[0] is the probability
 *  of sampling 0.
 */
size_t sf_rand_distribution(size_t n, double* dist)
{
    double rval = sf_rand(0.,1.);
    size_t i=0;
    // TODO add binary search here
    for (; i<n; i++) {
        if (dist[i]>rval) {
            break;
        }
    }
    
    return i;
}

unsigned int sf_rand_size(size_t maxval)
{
    assert(maxval > 0);
    tr1ns::uniform_int_distribution<size_t> dist(0,maxval-1);
    return dist(sparfun_rand_unif);
}

/** Compute a vector of degrees for a power-law */
template <typename VertexType>
void random_power_law_degrees(size_t n, double theta, size_t max_degree,
        VertexType* degrees)
{
    assert(theta>0);
    // first compute a distribution over degrees based on the powerlaw
    std::vector<double> dist(max_degree,0);
    double total = 0.;
    for (size_t i=0; i<max_degree; ++i) {
        dist[i] = (double)n*1./(pow(((double)i+1.),theta));
        total += dist[i];
    }
    // normalize to be a cumulative probability distribuiton
    double current_sum = 0.;
    for (size_t i=0; i<max_degree; ++i) {
        current_sum += dist[i]/total;
        dist[i] = current_sum;
    }
    
    // now sample n degrees
    for (size_t i=0; i<n; ++i) {
        degrees[i] = (VertexType)sf_rand_distribution(max_degree, &dist[0]) + 1;
    }
}



/** Generate a graph based on a really simple power-law graph model.
 *
 * TODO
 * We throw in a bit of local clustering ala Neville and collabs too...
 * http://arxiv.org/abs/1202.4805
 */
void generate_graph(int nverts, double pow, int maxdeg, 
    std::vector<std::pair<int,int> >& edges) {
    // get the power-law degree dist
    std::vector<int> degs(nverts);
    random_power_law_degrees(nverts, pow, maxdeg, &degs[0]);
    
    int nedges = 0;
    for (int i=0; i<nverts; ++i) { nedges += degs[i]; }
    
    // pick nedges with <src,dst> pairs sampled with prop proportional
    // to the degree, to do so, build a big long vector with each vertex id 
    // repeated in proportion to its degree
    std::vector<int> edgesampler(nedges);
    for (int i=0, ei=0; i<nverts; ++i) {
        for (int d=0; d<degs[i]; ++d, ++ei) {
            edgesampler[ei] = i;
        }
    }
    
    //std::vector<std::pair<int,int> > edges(nedges);
    edges.resize(nedges);
    for (int i=0; i<nedges; ++i) {
        int src = edgesampler[sf_rand_size(nedges)];
        int dst = edgesampler[sf_rand_size(nedges)];
        edges[i] = std::make_pair(src,dst);
        //cout << "picked edge (" << src << " " << dst << ")" << endl;
    }
}

/** Compute a performance benchmark on a synthetic graph. */
void benchmark() {
    int nthreads = omp_get_max_threads();
    cout << "Can use up to " << nthreads << " threads" << endl;
    
    // test sizes
    // generate a 10000 node graph
    int nverts = 200000;
    std::vector<std::pair<int,int> > edges;
    generate_graph(nverts, 1.8, 10000, edges);
    prpack::prpack_base_graph *g = new prpack::prpack_base_graph(
                                        nverts, edges.size(), &edges[0]);
    cout << "nverts = " << nverts << endl;
    cout << "nedges = " << edges.size() << endl;
	prpack::prpack_solver solver(g);
	
	{
        cout << "method = sccgs" << endl;	
        omp_set_num_threads(1);
        prpack::prpack_result* res = solver.solve(0.85, 1.e-10, NULL, NULL, "sccgs");
        cout << "  preprocess time = " << res->preprocess_time << "s" << endl;
        cout << "  1-thread compute time = " << res->compute_time << "s" << endl;
        for (int t=2; t<=nthreads; ++t) {
            omp_set_num_threads(t);
            prpack::prpack_result* res = solver.solve(0.85, 1.e-10, NULL, NULL, "sccgs");
            cout << "  " << t << "-thread compute time = " << res->compute_time << "s" 
                 << "  " << res->num_es_touched/(double)edges.size() << " eff iters" 
                 << endl;
        }
    }
	
	{
        cout << "method = gs" << endl;
        prpack::prpack_result* res = solver.solve(0.85, 1.e-10, NULL, NULL, "gs");
        cout << "  preprocess time = " << res->preprocess_time << "s" << endl;
        cout << "  " << 1 << "-thread compute time = " << res->compute_time << "s" 
                 << "  " << res->num_es_touched/(double)edges.size() << " eff iters" 
                 << endl;

    }
    
    // TODO free memory
}
