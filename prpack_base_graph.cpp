#include "prpack_base_graph.h"
#include "prpack_utils.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <limits>
using namespace prpack;
using namespace std;

void prpack_base_graph::initialize() {
    heads = NULL;
    tails = NULL;
}

prpack_base_graph::prpack_base_graph(prpack_csc* g) {
    initialize();
    num_vs = g->num_vs;
    num_es = g->num_es;
    // fill in heads and tails
    num_self_es = 0;
    int* hs = g->heads;
    int* ts = g->tails;
    tails = new int[num_vs];
    memset(tails, 0, num_vs*sizeof(tails[0]));
    for (int h = 0; h < num_vs; ++h) {
        int start_ti = hs[h];
        int end_ti = (h + 1 != num_vs) ? hs[h + 1] : num_es;
        for (int ti = start_ti; ti < end_ti; ++ti) {
            int t = ts[ti];
            ++tails[t];
            if (h == t)
                ++num_self_es;
        }
    }
    for (int i = 0, sum = 0; i < num_vs; ++i) {
        int temp = sum;
        sum += tails[i];
        tails[i] = temp;
    }
    heads = new int[num_es];
    int* osets = new int[num_vs];
    memset(osets, 0, num_vs*sizeof(osets[0]));
    for (int h = 0; h < num_vs; ++h) {
        int start_ti = hs[h];
        int end_ti = (h + 1 != num_vs) ? hs[h + 1] : num_es;
        for (int ti = start_ti; ti < end_ti; ++ti) {
            int t = ts[ti];
            heads[tails[t] + osets[t]++] = h;
        }
    }
    // clean up
    delete[] osets;
}


prpack_base_graph::prpack_base_graph(prpack_int64_csc* g) {
    initialize();
    // TODO remove the assert and add better behavior
    assert(num_vs <= std::numeric_limits<int>::max());
    num_vs = (int)g->num_vs;
    num_es = (int)g->num_es;
    // fill in heads and tails
    num_self_es = 0;
    int64_t* hs = g->heads;
    int64_t* ts = g->tails;
    tails = new int[num_vs];
    memset(tails, 0, num_vs*sizeof(tails[0]));
    for (int h = 0; h < num_vs; ++h) {
        int start_ti = (int)hs[h];
        int end_ti = (h + 1 != num_vs) ? (int)hs[h + 1] : num_es;
        for (int ti = start_ti; ti < end_ti; ++ti) {
            int t = (int)ts[ti];
            ++tails[t];
            if (h == t)
                ++num_self_es;
        }
    }
    for (int i = 0, sum = 0; i < num_vs; ++i) {
        int temp = sum;
        sum += tails[i];
        tails[i] = temp;
    }
    heads = new int[num_es];
    int* osets = new int[num_vs];
    memset(osets, 0, num_vs*sizeof(osets[0]));
    for (int h = 0; h < num_vs; ++h) {
        int start_ti = (int)hs[h];
        int end_ti = (h + 1 != num_vs) ? (int)hs[h + 1] : num_es;
        for (int ti = start_ti; ti < end_ti; ++ti) {
            int t = (int)ts[ti];
            heads[tails[t] + osets[t]++] = h;
        }
    }
    // clean up
    delete[] osets;
}

prpack_base_graph::prpack_base_graph(prpack_csr* g) {
    initialize();
    assert(false);
    // TODO
}

prpack_base_graph::prpack_base_graph(prpack_edge_list* g) {
    initialize();
    num_vs = g->num_vs;
    num_es = g->num_es;
    // fill in heads and tails
    num_self_es = 0;
    int* hs = g->heads;
    int* ts = g->tails;
    tails = new int[num_vs];
    memset(tails, 0, num_vs*sizeof(tails[0]));
    for (int i = 0; i < num_es; ++i) {
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
    delete[] osets;
}

prpack_base_graph::prpack_base_graph(const string& filename, const string& format) {
    initialize();
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

prpack_base_graph::~prpack_base_graph() {
    delete[] heads;
    delete[] tails;
}

bool prpack_base_graph::read_smat(FILE* f) {
    // read in header
    int nvs2=0;
    assert(fscanf(f, "%d %d %d", &num_vs, &nvs2, &num_es) == 3);
    if (nvs2 != num_vs) {
        return false;
    }
    // fill in heads and tails
    num_self_es = 0;
    int* hs = new int[num_es];
    int* ts = new int[num_es];
    tails = new int[num_vs];
    memset(tails, 0, num_vs*sizeof(tails[0]));
    for (int i = 0; i < num_es; ++i) {
        double val;
        assert(fscanf(f, "%d %d %lf", &hs[i], &ts[i], &val) == 3);
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
    delete[] hs;
    delete[] ts;
    delete[] osets;
    return true;
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
    vector<int>* al = new vector<int>[num_vs];
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
    delete[] al;
}

prpack_base_graph::prpack_base_graph(int nverts, int nedges, 
        std::pair<int,int>* edges) {
    initialize();
    num_vs = nverts;
    num_es = nedges;

    // fill in heads and tails
    num_self_es = 0;
    int* hs = new int[num_es];
    int* ts = new int[num_es];
    tails = new int[num_vs];
    memset(tails, 0, num_vs*sizeof(tails[0]));
    for (int i = 0; i < num_es; ++i) {
        assert(edges[i].first >= 0 && edges[i].first < num_vs);
        assert(edges[i].second >= 0 && edges[i].second < num_vs);
        hs[i] = edges[i].first;
        ts[i] = edges[i].second;
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
    delete[] hs;
    delete[] ts;
    delete[] osets;
}
