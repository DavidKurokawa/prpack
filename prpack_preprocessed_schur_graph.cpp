#include "prpack_preprocessed_schur_graph.h"
#include <algorithm>
#include <cstring>
#include <list>
using namespace prpack;
using namespace std;

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(prpack_adjacency_list* al) {
	// initialize instance variables
	num_vs = al->num_vs;
	num_es = al->num_es;
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	// permute dangling nodes to end
	int* encoding = new int[num_vs];
	memset(encoding, -1, num_vs*sizeof(encoding[0]));
	decoding = new int[num_vs];
	int seen = 0;
	for (int b = 0; b < num_vs; ++b) {
		for (list<int>::iterator a = al->al[b].begin(); a != al->al[b].end(); ++a) {
			if (encoding[*a] == -1)
				encoding[*a] = seen++;
			++inv_num_outlinks[encoding[*a]];
		}
	}
	num_dangling_vs = num_vs - seen;
	for (int i = 0; i < num_vs; ++i)
		if (encoding[i] == -1)
			encoding[i] = seen++;
	for (int i = 0; i < num_vs; ++i)
		decoding[encoding[i]] = i;
	// convert al to head/tail format
	ii = new double[num_vs];
	tails = new int[num_vs];
	heads = new int[num_es];
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		ii[tails_i] = 0;
		tails[tails_i] = heads_i;
		for (list<int>::iterator curr = al->al[decoding[tails_i]].begin(); curr != al->al[decoding[tails_i]].end(); ++curr) {
			if (decoding[tails_i] == *curr) {
				ii[tails_i] += 1;
				--num_es;
			} else {
				heads[heads_i++] = encoding[*curr];
			}
		}
		inv_num_outlinks[tails_i] = (inv_num_outlinks[tails_i] == 0) ? -1 : 1/inv_num_outlinks[tails_i];
	}
	// invert ii
	for (int i = 0; i < num_vs; ++i)
		if (ii[i] > 0.5)
			ii[i] *= inv_num_outlinks[i];
}

