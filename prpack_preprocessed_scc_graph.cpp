#include "prpack_preprocessed_scc_graph.h"
#include <algorithm>
#include <cstring>
using namespace prpack;
using namespace std;

prpack_preprocessed_scc_graph::prpack_preprocessed_scc_graph(prpack_adjacency_list* al) {
	num_vs = al->num_vs;
	num_es = al->num_es;
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	ii = new double[num_vs];
	convert(al, tails, heads);
}

// Convert the adjacency list to heads/tails format. This method will work regardless of the inlink/outlink orientation of the adjacency list
void prpack_preprocessed_scc_graph::convert(prpack_adjacency_list* al, int*& x, int*& y) {
	x = new int[num_vs];
	y = new int[num_es];
	iterative_tarjans();
	for (int i = 0; i < num_vs; ++i)
		if (ii[i] > 0.5)
			ii[i] *= inv_num_outlinks[i];
}

void prpack_preprocessed_scc_graph::iterative_tarjans() {
	// initialize variables
	num_comps = 0;
	int mn = 0;                 // the number of vertices seen so far
	int sz = 0;					// size of st
	int vs_i = 0;
	int* scc = new int[num_vs]; // the strongly connected component this vertex is in
	int* low = new int[num_vs]; // the lowest index this vertex can reach
	int* num = new int[num_vs]; // the index of this vertex in the dfs traversal
	int* st = new int[num_vs];  // a stack for the dfs
	int* vs = x;
	memset(num, -1, num_vs*sizeof(num[0]));
	memset(scc, -1, num_vs*sizeof(scc[0]));
	int* cs1 = new int[num_vs];                                // call stack variable for dfs
	matrix::iterator* cs2 = new matrix_type::iterator[num_vs]; // call stack variable for dfs
	// run iterative Tarjan's algorithm
	for (int root = 0; root < num_vs; ++root) {
		if (num[root] != -1)
			continue;
		int csz = 1;
		cs1[0] = root;
		cs2[0] = matrix[root].begin();
		// dfs
		while (csz) {
			int p = cs1[csz - 1]; // node we're dfs-ing on
			matrix_type::iterator& it = cs2[csz - 1]; // iteration of the for loop
			if (it == matrix[p].begin()) {
				low[p] = num[p] = mn++;
				st[sz++] = p;
			} else {
				--it;
				low[p] = min(low[p], low[*it]);
				++it;
			}
			bool done = false;
			for (; it != matrix[p].end(); ++it) {
				if (scc[*it] == -1) {
					if (num[*it] == -1) {
						// dfs(*it, p);
						cs1[csz] = *it;
						cs2[csz++] = matrix[*it].begin();
						++it;
						done = true;
						break;
					}
					low[p] = min(low[p], low[*it]);
				}
			}
			if (done)
				continue;
			// if p is the first explored vertex of a scc
			if (low[p] == num[p]) {
				cs1[num_vs - 1 - num_comps] = vs_i;
				while (scc[p] != num_comps) {
					scc[st[--sz]] = num_comps;
					vs[vs_i++] = st[sz];
				}
				++num_comps;
			}
			--csz;
		}
	}
	// clean up variables
	divisions = new int[num_comps];
	int* encoding = num; // given original i, return new i
	decoding = low;      // given new i, return original i
	for (int i = 0; i < num_vs; ++i)
		encoding[decoding[i] = x[i]] = i;
	for (int x_i = 0, y_i = 0; x_i < num_vs; ++x_i) {
		ii[x_i] = 0;
		x[x_i] = y_i;
		for (int& curr : matrix[decoding[x_i]]) {
			if (x_i == encoding[curr]) {
				ii[x_i] += 1;
				--num_es;
			} else {
				y[y_i++] = encoding[curr];
			}
			++inv_num_outlinks[encoding[curr]];
		}
	}
	divisions[0] = 0;
	for (int i = 1; i < num_comps; ++i)
		divisions[i] = cs1[num_vs - 1 - i];
	for (int i = 0; i < num_vs; ++i)
		inv_num_outlinks[i] = (inv_num_outlinks[i] == 0) ? -1 : 1/inv_num_outlinks[i];
	// free memory (do not free low, as we use it for decoding)
	free(cs1);
	free(cs2);
	free(scc);
	free(num);
	free(st);
}

