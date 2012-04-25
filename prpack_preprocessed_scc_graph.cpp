#include "prpack_preprocessed_scc_graph.h"
#include <algorithm>
#include <cstdlib>
#include <cstring>
using namespace prpack;
using namespace std;

prpack_preprocessed_scc_graph::prpack_preprocessed_scc_graph(prpack_base_graph* bg) {
	// initialize instance variables
	num_vs = bg->num_vs;
	num_es = bg->num_es - bg->num_self_es;
	// initialize Tarjan's algorithm variables
	num_comps = 0;
	int mn = 0;                 // the number of vertices seen so far
	int sz = 0;					// size of st
	int decoding_i = 0;         // size of decoding currently filled in
	decoding = new int[num_vs];
	int* scc = new int[num_vs]; // the strongly connected component this vertex is in
	int* low = new int[num_vs]; // the lowest index this vertex can reach
	int* num = new int[num_vs]; // the index of this vertex in the dfs traversal
	int* st = new int[num_vs];  // a stack for the dfs
	memset(num, -1, num_vs*sizeof(num[0]));
	memset(scc, -1, num_vs*sizeof(scc[0]));
	int* cs1 = new int[num_vs]; // call stack variable for dfs
	int* cs2 = new int[num_vs]; // call stack variable for dfs
	// run iterative Tarjan's algorithm
	for (int root = 0; root < num_vs; ++root) {
		if (num[root] != -1)
			continue;
		int csz = 1;
		cs1[0] = root;
		cs2[0] = bg->tails[root];
		// dfs
		while (csz) {
			int p = cs1[csz - 1]; // node we're dfs-ing on
			int& it = cs2[csz - 1]; // iteration of the for loop
			if (it == bg->tails[p]) {
				low[p] = num[p] = mn++;
				st[sz++] = p;
			} else {
				low[p] = min(low[p], low[bg->heads[it - 1]]);
			}
			bool done = false;
			int end_it = (p + 1 != num_vs) ? bg->tails[p + 1] : bg->num_es;
			for (; it < end_it; ++it) {
				int h = bg->heads[it];
				if (scc[h] == -1) {
					if (num[h] == -1) {
						// dfs(h, p);
						cs1[csz] = h;
						cs2[csz++] = bg->tails[h];
						++it;
						done = true;
						break;
					}
					low[p] = min(low[p], low[h]);
				}
			}
			if (done)
				continue;
			// if p is the first explored vertex of a scc
			if (low[p] == num[p]) {
				cs1[num_vs - 1 - num_comps] = decoding_i;
				while (scc[p] != num_comps) {
					scc[st[--sz]] = num_comps;
					decoding[decoding_i++] = st[sz];
				}
				++num_comps;
			}
			--csz;
		}
	}
	// set up other instance variables
	divisions = new int[num_comps];
	divisions[0] = 0;
	for (int i = 1; i < num_comps; ++i)
		divisions[i] = cs1[num_vs - 1 - i];
	encoding = num;
	for (int i = 0; i < num_vs; ++i)
		encoding[decoding[i]] = i;
	// fill in inside and outside instance variables
	ii = new double[num_vs];
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	tails_inside = cs1;
	heads_inside = new int[num_es];
	tails_outside = cs2;
	heads_outside = new int[num_es];
	num_es_inside = num_es_outside = 0;
	for (int comp_i = 0; comp_i < num_comps; ++comp_i) {
		const int start_i = divisions[comp_i];
		const int end_i = (comp_i + 1 != num_comps) ? divisions[comp_i + 1] : num_vs;
		for (int i = start_i; i < end_i; ++i) {
			ii[i] = 0;
			const int decoded = decoding[i];
			const int start_j = bg->tails[decoded];
			const int end_j = (decoded + 1 != num_vs) ? bg->tails[decoded + 1] : bg->num_es;
			tails_inside[i] = num_es_inside;
			tails_outside[i] = num_es_outside;
			for (int j = start_j; j < end_j; ++j) {
				int h = encoding[bg->heads[j]];
				if (h == i) {
					++ii[i];
				} else {
					if (start_i <= h && h < end_i)
						heads_inside[num_es_inside++] = h;
					else
						heads_outside[num_es_outside++] = h;
				}
				++inv_num_outlinks[h];
			}
		}
	}
	for (int i = 0; i < num_vs; ++i) {
		inv_num_outlinks[i] = (inv_num_outlinks[i] == 0) ? -1 : 1/inv_num_outlinks[i];
		ii[i] *= inv_num_outlinks[i];
	}
	// free memory
	// do not free num <==> encoding
	// do not free cs1 <==> tails_inside
	// do not free cs2 <==> tails_outside
	free(scc);
	free(low);
	free(st);
}

