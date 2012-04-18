#include "prpack_preprocessed_schur_graph.h"
#include <algorithm>
#include <cstring>
using namespace prpack;
using namespace std;

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(prpack_adjacency_list* al) {
	// initialize instance variables
	num_vs = al->num_vs;
	num_es = al->num_es - al->num_self_es;
	inv_num_outlinks = new double[num_vs];
	fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
	// permute dangling nodes to end
	encoding = new int[num_vs];
	decoding = new int[num_vs];
	memset(encoding, -1, num_vs*sizeof(encoding[0]));
	int seen = 0;
	for (int t = 0; t < num_vs; ++t) {
		int start_i = al->tails[t];
		int end_i = (t + 1 != num_vs) ? al->tails[t + 1] : al->num_es;
		for (int i = start_i; i < end_i; ++i) {
			int h = al->heads[i];
			if (encoding[h] == -1)
				decoding[encoding[h] = seen++] = h;
			++inv_num_outlinks[encoding[h]];
		}
	}
	num_dangling_vs = num_vs - seen;
	for (int i = 0; i < num_vs; ++i)
		if (encoding[i] == -1)
			decoding[encoding[i] = seen++] = i;
	// convert al to head/tail format
	ii = new double[num_vs];
	tails = new int[num_vs];
	heads = new int[num_es];
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		ii[tails_i] = 0;
		tails[tails_i] = heads_i;
		int decoded = decoding[tails_i];
		int start_i = al->tails[decoded];
		int end_i = (decoded + 1 != num_vs) ? al->tails[decoded + 1] : al->num_es;
		for (int i = start_i; i < end_i; ++i) {
			if (decoded == al->heads[i])
				++ii[tails_i];
			else
				heads[heads_i++] = encoding[al->heads[i]];
		}
		inv_num_outlinks[tails_i] = (inv_num_outlinks[tails_i] == 0) ? -1 : 1/inv_num_outlinks[tails_i];
	}
	// invert ii
	for (int i = 0; i < num_vs; ++i)
		if (ii[i] > 0.5)
			ii[i] *= inv_num_outlinks[i];
}

