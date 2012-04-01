#include "prpack_solver.h"
#include <cmath>
using namespace prpack;
using namespace std;

prpack_solver::prpack_solver(prpack_csr* g) {
	//this->g = new prpack_preprocessed_graph(g);
}

prpack_solver::prpack_solver(prpack_edgelist* g) {
	//this->g = new prpack_preprocessed_graph(g);
}

prpack_solver::prpack_solver(const string& filename) {
	this->g = new prpack_preprocessed_graph(filename);
}

prpack_result* prpack_solver::solve(double alpha, double tol) {
	return solve(alpha, tol, NULL, NULL);
}

prpack_result* prpack_solver::solve(double alpha, double tol, double* u, double* v) {
	// set up convenience variables
	int num_vs = g->num_vs;
	int num_es = g->num_es;
	int* heads = g->heads;
	int* tails = g->tails;
	double* inv_num_outlinks = g->inv_num_outlinks;
	double* ii = g->ii;
	// set up return variable
	prpack_result* ret = new prpack_result();
	ret->num_vs = num_vs;
	ret->num_es = num_es;
	// initialize u and v values
	double u_const = 1.0/num_vs;
	double v_const = (1.0 - alpha)/num_vs;
	int u_exists = (u) ? 1 : 0;
	int v_exists = (v) ? 1 : 0;
	u = (u) ? u : &u_const;
	v = (v) ? v : &v_const;
	if (v_exists)
		for (int i = 0; i < num_vs; i++)
			v[i] *= (1.0 - alpha);
	// initialize the eigenvector (and use personalization vector)
	double* x = new double[num_vs];
	for (int i = 0; i < num_vs; ++i)
		x[i] = v[v_exists*i]/(1.0 - alpha)*inv_num_outlinks[i];
	// initialize delta
	double delta = 0;
	for (int i = 0; i < num_vs; ++i)
		if (inv_num_outlinks[i] < 0)
			delta += x[i]/inv_num_outlinks[i];
	delta *= alpha;
	// run Gauss-Seidel
	ret->num_iter = 0;
	double err, old_val, new_val, t, y, c = 0;
	do {
		// iterate through vertices
		err = 0;
		for (int i = 0; i < num_vs; i++) {
			old_val = x[i]/inv_num_outlinks[i];
			new_val = 0;
			int start_j = tails[i], end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
			for (int j = start_j; j < end_j; j++)
				// TODO: might want to use compensation summation for large: end_j - start_j
				new_val += x[heads[j]];
			new_val = alpha*new_val + v[v_exists*i];
			if (inv_num_outlinks[i] < 0) {
				delta -= alpha*old_val;
				new_val += delta*u[u_exists*i];
				new_val /= 1 - alpha*u[u_exists*i];
				delta += alpha*new_val;
			} else {
				new_val += delta*u[u_exists*i];
				new_val /= 1 - alpha*ii[i];
			}
			// use compensation summation for: err += fabs(new_val - old_val)
			y = fabs(new_val - old_val) - c;
			t = err + y;
			c = (t - err) - y;
			err = t;
			x[i] = new_val*inv_num_outlinks[i];
		}
		// update iteration index
		++ret->num_iter;
	} while (err >= tol);
	// undo inv_num_outlinks transformation
	for (int i = 0; i < num_vs; i++)
		x[i] /= inv_num_outlinks[i];
	// clean up
	if (v_exists)
		for (int i = 0; i < num_vs; i++)
			v[i] /= (1.0 - alpha);
	// return results
	ret->x = x;
	return ret;
}

