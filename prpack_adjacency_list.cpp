#include "prpack_adjacency_list.h"
#include <fstream>
using namespace prpack;
using namespace std;

prpack_adjacency_list::prpack_adjacency_list(prpack_csr* g) {
	// TODO
}

prpack_adjacency_list::prpack_adjacency_list(prpack_edge_list* g) {
	// TODO
}

prpack_adjacency_list::prpack_adjacency_list(const string& filename) {
	// TODO: handle other formats than .smat
	FILE* f = fopen(filename.c_str(), "r");
	// read in header
	float blah;
	fscanf(f, "%d%f%d", &num_vs, &blah, &num_es);
	// read in all the edges
	int h, t;
	list<int>* matrix = new list<int>[num_vs];
	for (int i = 0; i < num_es; ++i) {
		fscanf(f, "%d%d%f", &h, &t, &blah);
		matrix[t].push_back(h);
	}
}

