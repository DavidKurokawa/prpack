#include "prpack_utils.h"
#include "prpack_solver.h"
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include <omp.h>
using namespace prpack;
using namespace std;

// in prpack_driver_benchmark.cpp
void benchmark();

// Contains all possible input parameters.
class input {
	public:
		// instance variables
		string graph;
		string format;
		double alpha;
		double tol;
		string u;
		string v;
		string method;
		string output;
		// constructor (default values can be set here)
		input() {
			graph = "";
			format = "";
			alpha = 0.85;
			tol = 1e-10;
			u = "";
			v = "";
			method = "";
			output = "";
		}
		// methods
		void parse_arg(string a, string b) {
			if (a == "-f" || a == "--format")
				format = b;
			else if (a == "-a" || a == "--alpha")
				alpha = atof(b.c_str());
			else if (a == "-t" || a == "--tol" || a == "--tolerance")
				tol = atof(b.c_str());
			else if (a == "-o" || a == "--out" || a == "--output")
				output = b;
			else if (a == "-u" || a == "--u")
				u = b;
			else if (a == "-v" || a == "--v")
				v = b;
			else if (a == "-m" || a == "--method")
				method = b;
			else
				prpack_utils::validate(false, "Error: argument '" + a + "' is not valid");
		}
};

double* read_vector(const string& filename) {
	if (filename == "")
		return NULL;
	// read into double vector
	vector<double> v;
	FILE* f = fopen(filename.c_str(), "r");
	double curr;
	while (fscanf(f, "%lf", &curr) == 1)
		v.push_back(curr);
	fclose(f);
	// convert to double array
	double* ret = new double[v.size()];
	for (int i = 0; i < (int) v.size(); ++i)
		ret[i] = v[i];
	return ret;
}

void write_vector(double *x, int n, ostream& out) {
    out.precision(16);
    out << scientific;
    for (int i=0; i<n; ++i) {
        out << x[i] << std::endl;
    }
}

int main(int argc, char** argv) {
	// parse command args
	input in;
	in.graph = string(argv[1]);
	for (int i = 2; i < argc; ++i) {
		string x(argv[i]);
		int idx = x.find("=");
		if (idx == (int) string::npos) {
			prpack_utils::validate(x.length() == 2 && x[0] == '-', "Error: argument '" + x + "' is not valid");
			prpack_utils::validate(i + 1 < argc, "Error: argument '" + x + "' does not specify value");
			in.parse_arg(x, string(argv[++i]));
		} else {
			prpack_utils::validate(x.length() > 2 && x[0] == '-' && x[1] == '-', "Error: argument '" + x + "' is not valid");
			in.parse_arg(x.substr(0, idx), x.substr(idx + 1));
		}
	}

    if (in.graph == "?") {
        benchmark();
        return 0;
    }	
	
	// solve
	prpack_solver solver(in.graph, in.format);
	double* u = read_vector(in.u);
	double* v = (in.u == in.v) ? u : read_vector(in.v);
	prpack_result* res = solver.solve(in.alpha, in.tol, u, v, in.method);
	// create output stream for text data
	ostream* out = &cout; // usually, this is cout
	if (in.output == "-") {
	    out = &cerr; 
	} 
	*out << "---------------------------" << endl;
	*out << "graph = " << in.graph << endl;
	*out << "number of vertices = " << res->num_vs << endl;
	*out << "number of edges = " << res->num_es << endl;
	*out << "---------------------------" << endl;
	*out << "method = " << res->method << endl;
	*out << "read time = " << res->read_time << "s" << endl;
	*out << "preprocess time = " << res->preprocess_time << "s" << endl;
	*out << "compute time = " << res->compute_time << "s" << endl;
	*out << "number of edges touched = " << res->num_es_touched << endl;
	*out << "converged = " << res->converged << endl;
	*out << "---------------------------" << endl;
	pair<double, int>* xval_idx = new pair<double, int>[res->num_vs];
	for (int i = 0; i < res->num_vs; ++i)
		xval_idx[i] = make_pair(res->x[i], i);
	sort(xval_idx, xval_idx + res->num_vs);
	for (int i = res->num_vs - 1; i >= 0 && i >= res->num_vs - 20; --i)
		*out << "site = " << xval_idx[i].second << ", score = " << xval_idx[i].first << endl;
	delete[] xval_idx;
	
	// write the entire vector
    if (in.output != "") {
        if (in.output == "-") {
            write_vector(res->x, res->num_vs, cout);
        } else {
            ofstream outfile(in.output.c_str());
            write_vector(res->x, res->num_vs, outfile);
        }
    }
}

