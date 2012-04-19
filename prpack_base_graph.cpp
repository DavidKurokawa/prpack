#include "prpack_base_graph.h"
#include "prpack_utils.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
using namespace prpack;
using namespace std;

prpack_base_graph::prpack_base_graph(prpack_csr* g) {
	// TODO
}

prpack_base_graph::prpack_base_graph(prpack_edge_list* g) {
	// TODO
}

prpack_base_graph::prpack_base_graph(const string& filename, const string& format) {
	FILE* f = fopen(filename.c_str(), "r");
	string ext = (format == "") ? filename.substr(filename.rfind('.') + 1) : format;
	if (ext == "smat")
		read_smat(f);
	else if (ext == "edges" || ext == "eg2")
		read_edges(f);
	else if (ext == "graph-txt")
		read_ascii(f);
	else
		prpack_utils::validate(false, "Invalid graph format");
	fclose(f);
}

void prpack_base_graph::read_smat(FILE* f) {
	// read in header
	float blah;
	assert(fscanf(f, "%d %f %d", &num_vs, &blah, &num_es) == 3);
	// fill in heads and tails
	num_self_es = 0;
	int* hs = new int[num_es];
	int* ts = new int[num_es];
	tails = new int[num_vs];
	memset(tails, 0, num_vs*sizeof(tails[0]));
	for (int i = 0; i < num_es; ++i) {
		assert(fscanf(f, "%d %d %f", &hs[i], &ts[i], &blah) == 3);
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
	for (int i = 0; i < num_es; ++i)
		heads[tails[ts[i]] + osets[ts[i]]++] = hs[i];
	// clean up
	free(hs);
	free(ts);
	free(osets);
}

void prpack_base_graph::read_edges(FILE* f) {
	vector<vector<int> > al;
	int h, t;
	num_es = num_self_es = 0;
	while (fscanf(f, "%d %d", &h, &t) == 2) {
		int m = (h < t) ? t : h;
		if ((int) al.size() < m + 1)
			al.resize(m + 1);
		al[t].push_back(h);
		++num_es;
		if (h == t)
			++num_self_es;
	}
	num_vs = al.size();
	heads = new int[num_es];
	tails = new int[num_vs];
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		tails[tails_i] = heads_i;
		for (int j = 0; j < (int) al[tails_i].size(); ++j)
			heads[heads_i++] = al[tails_i][j];
	}
}

void prpack_base_graph::read_ascii(FILE* f) {
	assert(fscanf(f, "%d", &num_vs) == 1);
	while (getc(f) != '\n');
	vector<int> al[num_vs];
	num_es = num_self_es = 0;
	char s[32];
	for (int h = 0; h < num_vs; ++h) {
		bool line_ended = false;
		while (!line_ended) {
			for (int i = 0; ; ++i) {
				s[i] = getc(f);
				if ('9' < s[i] || s[i] < '0') {
					line_ended = s[i] == '\n';
					if (i != 0) {
						s[i] = '\0';
						int t = atoi(s);
						al[t].push_back(h);
						++num_es;
						if (h == t)
							++num_self_es;
					}
					break;
				}
			}
		}
	}
	heads = new int[num_es];
	tails = new int[num_vs];
	for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
		tails[tails_i] = heads_i;
		for (int j = 0; j < (int) al[tails_i].size(); ++j)
			heads[heads_i++] = al[tails_i][j];
	}
}

