#include "prpack_preprocessed_schur_graph.h"
#include "prpack_utils.h"
#include <algorithm>
#include <cstring>
using namespace prpack;
using namespace std;

void prpack_preprocessed_schur_graph::initialize() {
    ii = NULL;
    inv_num_outlinks = NULL;
    heads = NULL;
    tails = NULL;
    encoding = NULL;
    decoding = NULL;
}

prpack_preprocessed_schur_graph::prpack_preprocessed_schur_graph(prpack_base_graph* bg) {
    initialize();
    // initialize instance variables
    num_vs = bg->num_vs;
    num_es = bg->num_es - bg->num_self_es;
    inv_num_outlinks = new double[num_vs];
    fill(inv_num_outlinks, inv_num_outlinks + num_vs, 0);
    for (int i = 0; i < bg->num_es; ++i)
        ++inv_num_outlinks[bg->heads[i]];
    // permute no-inlink vertices to the beginning, and no-outlink vertices to the end
    encoding = new int[num_vs];
    decoding = new int[num_vs];
    num_no_in_vs = num_no_out_vs = 0;
    for (int i = 0; i < num_vs; ++i) {
        if (bg->tails[i] == ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es)) {
            decoding[encoding[i] = num_no_in_vs] = i;
            ++num_no_in_vs;
        } else if (inv_num_outlinks[i] == 0) {
            decoding[encoding[i] = num_vs - 1 - num_no_out_vs] = i;
            ++num_no_out_vs;
        }
    }
    // permute everything else
    for (int i = 0, p = num_no_in_vs; i < num_vs; ++i)
        if (bg->tails[i] < ((i + 1 != num_vs) ? bg->tails[i + 1] : bg->num_es) && inv_num_outlinks[i] > 0)
            decoding[encoding[i] = p++] = i;
    // permute inv_num_outlinks
    ii = inv_num_outlinks;
    inv_num_outlinks = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        inv_num_outlinks[encoding[i]] = (ii[i] == 0) ? -1 : 1/ii[i];
    // convert bg to head/tail format
    tails = new int[num_vs];
    heads = new int[num_es];
    for (int tails_i = 0, heads_i = 0; tails_i < num_vs; ++tails_i) {
        ii[tails_i] = 0;
        tails[tails_i] = heads_i;
        int decoded = decoding[tails_i];
        int start_i = bg->tails[decoded];
        int end_i = (decoded + 1 != num_vs) ? bg->tails[decoded + 1] : bg->num_es;
        for (int i = start_i; i < end_i; ++i) {
            if (decoded == bg->heads[i])
                ++ii[tails_i];
            else
                heads[heads_i++] = encoding[bg->heads[i]];
        }
        if (ii[tails_i] > 0)
            ii[tails_i] *= inv_num_outlinks[tails_i];
    }
}

prpack_preprocessed_schur_graph::~prpack_preprocessed_schur_graph() {
    delete[] ii;
    delete[] inv_num_outlinks;
    delete[] heads;
    delete[] tails;
    delete[] encoding;
    delete[] decoding;
}
