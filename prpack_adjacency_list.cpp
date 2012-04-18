#include "prpack_adjacency_list.h"
#include <cassert>
#include <cstring>
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
	assert(fscanf(f, "%d%f%d", &num_vs, &blah, &num_es) == 3);
	// fill in heads and tails
	num_es = 0;
	int* hs = new int[num_es];
	int* ts = new int[num_es];
	tails = new int[num_vs];
	memset(tails, 0, num_vs*sizeof(tails[0]));
	for (int i = 0; i < num_es; ++i) {
		assert(fscanf(f, "%d%d%d", &hs[i], &ts[i], &blah) == 3);
		++tails[ts[i]];
		if (hs[i] == ts[i])
			++num_self_es;
	}
	for (int i = 0, sum = 0; i < num_vs; ++i) {
		int temp = sum;
		sum += tails[i];
		tails[i] = temp;
	}
	heads = new int[num_es];
	int* osets = new int[num_vs];
	memset(osets, 0, num_vs*sizeof(osets[0]));
	for (int i = 0; i < num_nonself_es; ++i)
		heads[tails[ts[i]] + osets[ts[i]]++] = hs[i];
	// clean up
	free(hs);
	free(ts);
	free(osets);
	fclose(f);
}

