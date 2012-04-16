#include "prpack_adjacency_list.h"
#include <cassert>
#include <fstream>
#include <vector>
using namespace prpack;
using namespace std;

prpack_adjacency_list::prpack_adjacency_list(prpack_csr* g) {
	// TODO
}

prpack_adjacency_list::prpack_adjacency_list(prpack_edge_list* g) {
	// TODO
}

/** Wrap the given input after validating it. 
 * 
 * The caller is responsible for ensuring the memory hangs around throughout
 * the course of the computation.  
 */
prpack_adjacency_list::prpack_adjacency_list(int num_vs_, std::vector<int>* al_) {
    al = al_;
    num_vs = num_vs_;
    num_es = 0;
    for (int i=0; i < num_vs_; ++i) {
        num_es += (int)al[i].size();
    }
}

prpack_adjacency_list::prpack_adjacency_list(int nverts, int nedges, 
        std::pair<int,int>* edges) {
    num_vs = nverts;
    num_es = nedges;
    al = new vector<int>[num_vs];    
    for (int i = 0; i < num_es; ++i) {
	    assert(edges[i].first >= 0 && edges[i].first < num_vs);
	    assert(edges[i].second >= 0 && edges[i].second < num_vs);
	    al[edges[i].first].push_back(edges[i].second);
	}
}

prpack_adjacency_list::prpack_adjacency_list(const string& filename) {
	// TODO: handle other formats than .smat
	FILE* f = fopen(filename.c_str(), "r");
	// read in header
	float blah;
	assert(fscanf(f, "%d%f%d", &num_vs, &blah, &num_es) == 3);
	// read in all the edges
	int h, t;
	al = new vector<int>[num_vs];
	for (int i = 0; i < num_es; ++i) {
		assert(fscanf(f, "%d%d%f", &h, &t, &blah) == 3);
		al[t].push_back(h);
	}
	fclose(f);
}

