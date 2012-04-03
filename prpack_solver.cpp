#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
using namespace prpack;
using namespace std;

void prpack_solver::initialize() {
	gsg = NULL;
	sccg = NULL;
}

prpack_solver::prpack_solver(prpack_csr* g) {
	initialize();
	al = new prpack_adjacency_list(g);
}

prpack_solver::prpack_solver(prpack_edge_list* g) {
	initialize();
	al = new prpack_adjacency_list(g);
}

prpack_solver::prpack_solver(prpack_adjacency_list* g) {
	initialize();
	al = g;
}

prpack_solver::prpack_solver(const string& filename) {
	initialize();
	al = new prpack_adjacency_list(filename);
}

prpack_result* prpack_solver::solve(double alpha, double tol) {
	return solve(alpha, tol, NULL, NULL);
}

prpack_result* prpack_solver::solve(double alpha, double tol, double* u, double* v) {
	double preprocess_time = 0;
	double compute_time = 0;
	prpack_result* ret;
	if (u != v) {
		if (gsg == NULL)
			TIME(preprocess_time, gsg = new prpack_preprocessed_gs_graph(al));
		TIME(compute_time, ret = solve_via_gs(
				alpha,
				tol,
				gsg->num_vs,
				gsg->num_es,
				gsg->heads,
				gsg->tails,
				gsg->ii,
				gsg->inv_num_outlinks,
				u,
				v));
	} else {
		if (sccg == NULL)
			TIME(preprocess_time, sccg = new prpack_preprocessed_scc_graph(al));
		TIME(compute_time, ret = solve_via_scc_gs(
				alpha,
				tol,
				sccg->num_vs,
				sccg->num_es_inside,
				sccg->heads_inside,
				sccg->tails_inside,
				sccg->num_es_outside,
				sccg->heads_outside,
				sccg->tails_outside,
				sccg->ii,
				sccg->inv_num_outlinks,
				u,
				sccg->num_comps,
				sccg->divisions,
				sccg->decoding));
	}
	ret->preprocess_time = preprocess_time;
	ret->compute_time = compute_time;
	ret->num_vs = al->num_vs;
	ret->num_es = al->num_es;
	return ret;
}

// various solving methods ////////////////////////////////////////////////////////////////////////

// vanilla Gauss-Seidel
prpack_result* prpack_solver::solve_via_gs(
		double alpha,
		double tol,
		int num_vs,
		int num_es,
		int* heads,
		int* tails,
		double* ii,
		double* inv_num_outlinks,
		double* u,
		double* v) {
	prpack_result* ret = new prpack_result();
	// initialize u and v values
	double u_const = 1.0/num_vs;
	double v_const = (1.0 - alpha)/num_vs;
	int u_exists = (u) ? 1 : 0;
	int v_exists = (v) ? 1 : 0;
	u = (u) ? u : &u_const;
	v = (v) ? v : &v_const;
	if (v_exists)
		for (int i = 0; i < num_vs; ++i)
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
	double err, old_val, new_val, y, t, c = 0;
	do {
		// iterate through vertices
		err = 0;
		for (int i = 0; i < num_vs; ++i) {
			old_val = x[i]/inv_num_outlinks[i];
			new_val = 0;
			int start_j = tails[i], end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
			for (int j = start_j; j < end_j; ++j)
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
	for (int i = 0; i < num_vs; ++i)
		x[i] /= inv_num_outlinks[i];
	// clean up
	if (v_exists)
		for (int i = 0; i < num_vs; ++i)
			v[i] /= (1.0 - alpha);
	// return results
	ret->x = x;
	return ret;
}

// Gauss-Seidel using strongly connected components
prpack_result* prpack_solver::solve_via_scc_gs(
		double alpha,
		double tol,
		int num_vs,
		int num_es_inside,
		int* heads_inside,
		int* tails_inside,
		int num_es_outside,
		int* heads_outside,
		int* tails_outside,
		double* ii,
		double* inv_num_outlinks,
		double* uv,
		int num_comps,
		int* divisions,
		int* decoding) {
	prpack_result* ret = new prpack_result();
	// initialize uv values
	double uv_const = 1.0/num_vs;
	int uv_exists = (uv) ? 1 : 0;
	uv = (uv) ? uv : &uv_const;
	// initialize the eigenvector (and use personalization vector)
	double* x = new double[num_vs];
	for (int i = 0; i < num_vs; ++i)
		x[i] = uv[uv_exists*i]*inv_num_outlinks[i];
	// create x_outside
	double* x_outside = new double[num_vs];
	// run Gauss-Seidel for (I - alpha*P)*x = uv
	ret->num_iter = 0;
	for (int comp_i = 0; comp_i < num_comps; ++comp_i) {
		const int start_comp = divisions[comp_i];
		const int end_comp = (comp_i + 1 != num_comps) ? divisions[comp_i + 1] : num_vs;
		const bool parallelize = end_comp - start_comp > 512;
		// initialize relevant x_outside values
		for (int i = start_comp; i < end_comp; ++i) {
			x_outside[i] = 0;
			const int start_j = tails_outside[i];
			const int end_j = (i + 1 != num_vs) ? tails_outside[i + 1] : num_es_outside;
			for (int j = start_j; j < end_j; ++j)
				x_outside[i] += x[heads_outside[j]];
		}
		double err, c;
		do {
			if (parallelize) {
				// iterate through vertices
				#pragma omp parallel for if (parallelize) schedule(dynamic, 4)
				for (int i = start_comp; i < end_comp; ++i) {
					double new_val = x_outside[i];
					const int start_j = tails_inside[i];
					const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
					for (int j = start_j; j < end_j; ++j)
						// TODO: might want to use compensation summation for large: end_j - start_j
						new_val += x[heads_inside[j]];
					new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
					x[i] = new_val*inv_num_outlinks[i];
				}
				// compute error
				err = c = 0;
				#pragma omp parallel for if (parallelize) firstprivate(c) reduction(+:err) schedule(dynamic, 4)
				for (int i = start_comp; i < end_comp; ++i) {
					double curr = x_outside[i];
					const int start_j = tails_inside[i];
					const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
					for (int j = start_j; j < end_j; ++j)
						// TODO: might want to use compensation summation for large: end_j - start_j
						curr += x[heads_inside[j]];
					// use compensation summation for: err += fabs(uv[uv_exists*i] + alpha*curr - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i])
					double y = fabs(uv[uv_exists*i] + alpha*curr - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i]) - c;
					double t = err + y;
					c = t - err - y;
					err = t;
				}
			} else {
				// iterate through vertices
				for (int i = start_comp; i < end_comp; ++i) {
					double new_val = x_outside[i];
					const int start_j = tails_inside[i];
					const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
					for (int j = start_j; j < end_j; ++j)
						// TODO: might want to use compensation summation for large: end_j - start_j
						new_val += x[heads_inside[j]];
					new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
					x[i] = new_val*inv_num_outlinks[i];
				}
				// compute error
				err = c = 0;
				for (int i = start_comp; i < end_comp; ++i) {
					double curr = x_outside[i];
					const int start_j = tails_inside[i];
					const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
					for (int j = start_j; j < end_j; ++j)
						// TODO: might want to use compensation summation for large: end_j - start_j
						curr += x[heads_inside[j]];
					// use compensation summation for: err += fabs(uv[uv_exists*i] + alpha*curr - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i])
					double y = fabs(uv[uv_exists*i] + alpha*curr - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i]) - c;
					double t = err + y;
					c = t - err - y;
					err = t;
				}
			}
			// update iteration index
			++ret->num_iter;
		} while (err > tol*(end_comp - start_comp)/num_vs);
	}
	// undo inv_num_outlinks transformation
	for (int i = 0; i < num_vs; ++i)
		x[i] /= inv_num_outlinks[i];
	// normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
	double norm = 0;
	for (int i = 0; i < num_vs; ++i)
		norm += x[i];
	norm = 1/norm;
	for (int i = 0; i < num_vs; ++i)
		x[i] *= norm;
	// return results
	ret->x = new double[num_vs];
	for (int i = 0; i < num_vs; ++i)
		ret->x[decoding[i]] = x[i];
	free(x);
	return ret;
}

