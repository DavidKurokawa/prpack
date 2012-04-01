#include "prpack.h"
#include <iostream>
#include <string>
using namespace prpack;
using namespace std;

int main() {
	string data_file = "./data/web-Google_transposed.smat";
	prpack_solver solver(data_file);
	prpack_result* res = solver.solve(0.85, 1e-10);
	cout << "number of vertices = " << res->num_vs << endl;
	cout << "number of edges = " << res->num_es << endl;
	cout << "number of iterations = " << res->num_iter << endl;
}

