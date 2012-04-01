#include "prpack_preprocessed_graph.h"
#include <cstdio>
using namespace prpack;
using namespace std;

// Converts 'matrix' into the desired 'heads' and 'tails' representation.
void prpack_preprocessed_graph::convert(matrix_type* matrix, int& num_vs, int& num_es, int*& x, int*& y) {
	x = new int[num_vs];
	y = new int[num_es];
	for (int x_i = 0, y_i = 0; x_i < num_vs; ++x_i) {
		ii[x_i] = 0;
		x[x_i] = y_i;
		for (matrix_type::iterator it = matrix[x_i].begin(); it != matrix[x_i].end(); ++it) {
			if (x_i == *it) {
				ii[x_i] += 1;
				--num_es;
			} else {
				y[y_i++] = *it;
			}
		}
		inv_num_outlinks[x_i] = (inv_num_outlinks[x_i] == 0) ? -1 : 1/inv_num_outlinks[x_i];
	}
	// invert ii
	for (int i = 0; i < num_vs; ++i)
		if (ii[i] > 0.5)
			ii[i] *= inv_num_outlinks[i];
}

prpack_preprocessed_graph::prpack_preprocessed_graph(const string& filename) : prpack_graph() {
	// TODO: handle other formats than .smat
	FILE* f = fopen(filename.c_str(), "r");
	// read in header
	float blah;
	fscanf(f, "%d%f%d", &num_vs, &blah, &num_es);
	// set up 'inv_num_outlinks' and 'ii' now that 'num_vs' is known
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0.0);
	ii = new double[num_vs];
	// read in all the edges
	int h, t;
	matrix_type* matrix = new matrix_type[num_vs];
	for (int i = 0; i < num_es; ++i) {
		fscanf(f, "%d%d%f", &h, &t, &blah);
		matrix[t].push_back(h);
		++inv_num_outlinks[h];
	}
	// convert matrix head-tail form
	convert(matrix, num_vs, num_es, tails, heads);
}

