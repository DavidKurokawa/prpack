#include "prpack.h"
#include <algorithm>
#include <iostream>
#include <string>
#include <utility>
using namespace prpack;
using namespace std;

int main() {
	string data_file = "./data/web-Google_transposed.smat";
	prpack_solver solver(data_file);
	prpack_result* res = solver.solve(0.85, 1e-10);
	// output results
	cout << "---------------------------" << endl;
	cout << "graph = " << data_file << endl;
	cout << "---------------------------" << endl;
	cout << "number of vertices = " << res->num_vs << endl;
	cout << "number of edges = " << res->num_es << endl;
	cout << "number of iterations = " << res->num_iter << endl;
	cout << "---------------------------" << endl;
	pair<double, int>* ranks = new pair<double, int>[res->num_vs];
	for (int i = 0; i < res->num_vs; i++)
		ranks[i] = make_pair(res->x[i], i);
	sort(ranks, ranks + res->num_vs);
	//for (int i = res->num_vs - 1; i >= 0; i--)
	for (int i = res->num_vs - 1; i >= res->num_vs - 20; i--)
		//printf("site = %d, score = %E\n", ranks[i].second, ranks[i].first);
		cout << "site = " << ranks[i].second << ", score = " << ranks[i].first << endl;
	free(ranks);
}

