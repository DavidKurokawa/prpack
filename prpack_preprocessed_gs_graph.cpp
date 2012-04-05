#include "prpack_preprocessed_gs_graph.h"
#include <algorithm>
#include <list>
using namespace prpack;
using namespace std;

prpack_preprocessed_gs_graph::prpack_preprocessed_gs_graph(prpack_adjacency_list* al) {
	num_vs = al->num_vs;
	num_es = al->num_es;
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	for (int b = 0; b < num_vs; ++b)
		for (list<int>::iterator a = al->al[b].begin(); a != al->al[b].end(); ++a)
			++inv_num_outlinks[*a];
	ii = new double[num_vs];
	convert(al, tails, heads);
}

// Convert the adjacency list to heads/tails format. This method will work regardless of the inlink/outlink orientation of the adjacency list
void prpack_preprocessed_gs_graph::convert(prpack_adjacency_list* al, int*& x, int*& y) {
	x = new int[num_vs];
	y = new int[num_es];
	for (int x_i = 0, y_i = 0; x_i < num_vs; ++x_i) {
		ii[x_i] = 0;
		x[x_i] = y_i;
		for (list<int>::iterator curr = al->al[x_i].begin(); curr != al->al[x_i].end(); ++curr) {
			if (x_i == *curr) {
				ii[x_i] += 1;
				--num_es;
			} else {
				y[y_i++] = *curr;
			}
		}
		inv_num_outlinks[x_i] = (inv_num_outlinks[x_i] == 0) ? -1 : 1/inv_num_outlinks[x_i];
	}
	// invert ii
	for (int i = 0; i < num_vs; ++i)
		if (ii[i] > 0.5)
			ii[i] *= inv_num_outlinks[i];
}

