#include "prpack_solver.h"
#include "prpack_utils.h"
#include <cmath>
#include <cstdlib>
using namespace prpack;
using namespace std;

void prpack_solver::initialize() {
    gsg = NULL;
    sg = NULL;
    sccg = NULL;
}

prpack_solver::prpack_solver(prpack_csr* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(prpack_edge_list* g) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(g));
}

prpack_solver::prpack_solver(prpack_base_graph* g) {
    initialize();
    TIME(read_time, bg = g);
}

prpack_solver::prpack_solver(const string& filename, const string& format) {
    initialize();
    TIME(read_time, bg = new prpack_base_graph(filename, format));
}

#ifdef MATLAB_MEX_FILE
prpack_solver::prpack_solver(const mxArray* a) {
    initialize();
    // separate raw matlab arrays
    mxArray* raw_read_time = mxGetField(a, 0, "read_time");
    mxArray* raw_bg = mxGetField(a, 0, "bg");
    mxArray* raw_gsg = mxGetField(a, 0, "gsg");
    mxArray* raw_sg = mxGetField(a, 0, "sg");
    mxArray* raw_sccg = mxGetField(a, 0, "sccg");
    // initialize instance variables
    read_time = prpack_utils::matlab_array_to_double(raw_read_time);
    bg = new prpack_base_graph(raw_bg);
    if (!mxIsEmpty(raw_gsg))
        gsg = new prpack_preprocessed_gs_graph(raw_gsg);
    if (!mxIsEmpty(raw_sg))
        sg = new prpack_preprocessed_schur_graph(raw_sg);
    if (!mxIsEmpty(raw_sccg))
        sccg = new prpack_preprocessed_scc_graph(raw_sccg);
}
#endif

prpack_solver::~prpack_solver() {
    delete bg;
    delete gsg;
    delete sg;
    delete sccg;
}

#ifdef MATLAB_MEX_FILE
mxArray* prpack_solver::to_matlab_array(mxArray* a) {
    if (a == NULL) {
        const int num_fields = 5;
        const char* field_names[num_fields] = {"read_time", "bg", "gsg", "sg", "sccg"};
        mxArray* ret = mxCreateStructMatrix(1, 1, num_fields, field_names);
        mxSetField(ret, 0, "read_time", prpack_utils::double_to_matlab_array(read_time));
        mxSetField(ret, 0, "bg", bg->to_matlab_array());
        mxSetField(ret, 0, "gsg", (gsg != NULL) ? gsg->to_matlab_array() : prpack_utils::empty_matlab_array());
        mxSetField(ret, 0, "sg", (sg != NULL) ? sg->to_matlab_array() : prpack_utils::empty_matlab_array());
        mxSetField(ret, 0, "sccg", (sccg != NULL) ? sccg->to_matlab_array() : prpack_utils::empty_matlab_array());
        return ret;
    } else {
        if (gsg != NULL && mxIsEmpty(mxGetField(a, 0, "gsg")))
            mxSetField(a, 0, "gsg", gsg->to_matlab_array());
        if (sg != NULL && mxIsEmpty(mxGetField(a, 0, "sg")))
            mxSetField(a, 0, "sg", sg->to_matlab_array());
        if (sccg != NULL && mxIsEmpty(mxGetField(a, 0, "sccg")))
            mxSetField(a, 0, "sccg", sccg->to_matlab_array());
    }
}
#endif

prpack_result* prpack_solver::solve(double alpha, double tol, const string& method) {
    return solve(alpha, tol, NULL, NULL, method);
}

prpack_result* prpack_solver::solve(double alpha, double tol, double* u, double* v, const string& method) {
    double preprocess_time = 0;
    double compute_time = 0;
    prpack_result* ret;
    if (method == "gs") {
        if (gsg == NULL) {
            TIME(preprocess_time, gsg = new prpack_preprocessed_gs_graph(bg));
        }
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
        ret->method = "gs";
    } else if (method == "gserr") {
        if (gsg == NULL) {
            TIME(preprocess_time, gsg = new prpack_preprocessed_gs_graph(bg));
        }
        TIME(compute_time, ret = solve_via_gs_err(
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
        ret->method = "gserr";
    } else if (method == "sgs" || (method == "" && u == v)) {
        if (sg == NULL) {
            TIME(preprocess_time, sg = new prpack_preprocessed_schur_graph(bg));
        }
        TIME(compute_time, ret = solve_via_schur_gs(
                alpha,
                tol,
                sg->num_vs,
                sg->num_no_in_vs,
                sg->num_no_out_vs,
                sg->num_es,
                sg->heads,
                sg->tails,
                sg->ii,
                sg->inv_num_outlinks,
                u,
                sg->encoding,
                sg->decoding));
        ret->method = "sgs";
    } else if (method == "sgs_uv" || (method == "" && u != v)) {
        if (sg == NULL) {
            TIME(preprocess_time, sg = new prpack_preprocessed_schur_graph(bg));
        }
        TIME(compute_time, ret = solve_via_schur_gs_uv(
                alpha,
                tol,
                sg->num_vs,
                sg->num_no_in_vs,
                sg->num_no_out_vs,
                sg->num_es,
                sg->heads,
                sg->tails,
                sg->ii,
                sg->inv_num_outlinks,
                u,
                v,
                sg->encoding,
                sg->decoding));
        ret->method = "sgs_uv";
    } else if (method == "sccgs" || (method == "" && u == v)) {
        if (sccg == NULL) {
            TIME(preprocess_time, sccg = new prpack_preprocessed_scc_graph(bg));
        }
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
                sccg->encoding,
                sccg->decoding));
        ret->method = "sccgs";
    } else {
        if (sccg == NULL) {
            TIME(preprocess_time, sccg = new prpack_preprocessed_scc_graph(bg));
        }
        TIME(compute_time, ret = solve_via_scc_gs_uv(
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
                v,
                sccg->num_comps,
                sccg->divisions,
                sccg->encoding,
                sccg->decoding));
        ret->method = "sccgs_uv";
    }
    ret->read_time = read_time;
    ret->preprocess_time = preprocess_time;
    ret->compute_time = compute_time;
    ret->num_vs = bg->num_vs;
    ret->num_es = bg->num_es;
    return ret;
}

// VARIOUS SOLVING METHODS ////////////////////////////////////////////////////////////////////////

// Vanilla Gauss-Seidel.
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
    double v_const = 1.0/num_vs;
    int u_exists = (u) ? 1 : 0;
    int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        x[i] = 0;
    // initialize delta
    double delta = 0;
    // run Gauss-Seidel
    ret->num_es_touched = 0;
    double err = 1, old_val, new_val, c = 0;
    do {
        // iterate through vertices
        for (int i = 0; i < num_vs; ++i) {
            old_val = x[i]/inv_num_outlinks[i];
            new_val = 0;
            int start_j = tails[i], end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
            for (int j = start_j; j < end_j; ++j)
                // TODO: might want to use compensation summation for large: end_j - start_j
                new_val += x[heads[j]];
            new_val = alpha*new_val + (1 - alpha)*v[v_exists*i];
            if (inv_num_outlinks[i] < 0) {
                delta -= alpha*old_val;
                new_val += delta*u[u_exists*i];
                new_val /= 1 - alpha*u[u_exists*i];
                delta += alpha*new_val;
            } else {
                new_val += delta*u[u_exists*i];
                new_val /= 1 - alpha*ii[i];
            }
            COMPENSATED_SUM(err, old_val - new_val, c);
            x[i] = new_val*inv_num_outlinks[i];
        }
        // update iteration index
        ret->num_es_touched += num_es;
    } while (err >= tol);
    // undo inv_num_outlinks transformation
    for (int i = 0; i < num_vs; ++i)
        x[i] /= inv_num_outlinks[i];
    // return results
    ret->x = x;
    return ret;
}

// Implement a gauss-seidel-like process with a strict error bound
// we return a solution with 1-norm error less than tol.
prpack_result* prpack_solver::solve_via_gs_err(
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
    double v_const = 1.0/num_vs;
    int u_exists = (u) ? 1 : 0;
    int v_exists = (v) ? 1 : 0;
    u = (u) ? u : &u_const;
    v = (v) ? v : &v_const;
    // Note to Dave, we can't rescale v because we could be running this
    // same routine from multiple threads.
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i) {
        x[i] = 0.;
    }
    // initialize delta
    double delta = 0.;
    // run Gauss-Seidel, note that we store x/deg[i] throughout this 
    // iteration.
    long long maxedges = (long long)((double)num_es*std::min(
                            log(tol)/log(alpha),
                            (double)PRPACK_SOLVER_MAX_ITERS));
    ret->num_es_touched = 0;
    double err=1., c = 0.;
    do {
        // iterate through vertices
        for (int i = 0; i < num_vs; ++i) {
            double old_val = x[i]/inv_num_outlinks[i]; // adjust back to the "true" value.
            double new_val = 0.;
            int start_j = tails[i], end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
            for (int j = start_j; j < end_j; ++j) {
                // TODO: might want to use compensation summation for large: end_j - start_j
                new_val += x[heads[j]];
            }
            new_val = alpha*new_val + alpha*ii[i]*old_val + (1.0-alpha)*v[v_exists*i];
            new_val += delta*u[u_exists*i]; // add the dangling node adjustment
            if (inv_num_outlinks[i] < 0) {
                delta += alpha*(new_val - old_val);
            } 
            // note that new_val > old_val, but the fabs is just for 
            COMPENSATED_SUM(err, -(new_val - old_val), c);
            x[i] = new_val*inv_num_outlinks[i];
        }
        // update iteration index
        ret->num_es_touched += num_es;
    } while (err >= tol && ret->num_es_touched < maxedges);
    if (err >= tol) {
        ret->converged = 0;
    } else {
        ret->converged = 1;
    }
    // undo inv_num_outlinks transformation
    for (int i = 0; i < num_vs; ++i)
        x[i] /= inv_num_outlinks[i];
    // return results
    ret->x = x;
    return ret;
}

// Gauss-Seidel using the Schur complement to separate dangling nodes.
prpack_result* prpack_solver::solve_via_schur_gs(
        double alpha,
        double tol,
        int num_vs,
        int num_no_in_vs,
        int num_no_out_vs,
        int num_es,
        int* heads,
        int* tails,
        double* ii,
        double* inv_num_outlinks,
        double* uv,
        int* encoding,
        int* decoding,
        bool normalize) {
    prpack_result* ret = new prpack_result();
    // initialize uv values
    double uv_const = 1.0/num_vs;
    int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? prpack_utils::permute(num_vs, uv, encoding) : &uv_const;
    // initialize the eigenvector (and use personalization vector)
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs - num_no_out_vs; ++i)
        x[i] = uv[uv_exists*i]*inv_num_outlinks[i]/(1 - alpha*ii[i]);
    // run Gauss-Seidel for the top left part of (I - alpha*P)*x = uv
    ret->num_es_touched = 0;
    double err, c;
    do {
        // iterate through vertices
        int num_es_touched = 0;
        err = c = 0;
        #pragma omp parallel for firstprivate(c) reduction(+:err, num_es_touched) schedule(dynamic, 64)
        for (int i = num_no_in_vs; i < num_vs - num_no_out_vs; ++i) {
            double new_val = 0;
            const int start_j = tails[i];
            const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
            for (int j = start_j; j < end_j; ++j)
                // TODO: might want to use compensation summation for large: end_j - start_j
                new_val += x[heads[j]];
            COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i]), c);
            new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
            x[i] = new_val*inv_num_outlinks[i];
            num_es_touched += end_j - start_j;
        }
        // update iteration index
        ret->num_es_touched += num_es_touched;
    } while (err/(1 - alpha) >= tol*(num_vs - num_no_out_vs)/num_vs);
    // solve for the dangling nodes
    int num_es_touched = 0;
    #pragma omp parallel for reduction(+:num_es_touched) schedule(dynamic, 64)
    for (int i = num_vs - num_no_out_vs; i < num_vs; ++i) {
        x[i] = 0;
        const int start_j = tails[i];
        const int end_j = (i + 1 != num_vs) ? tails[i + 1] : num_es;
        for (int j = start_j; j < end_j; ++j)
            x[i] += x[heads[j]];
        x[i] = (alpha*x[i] + uv[uv_exists*i])/(1 - alpha*ii[i]);
        num_es_touched += end_j - start_j;
    }
    ret->num_es_touched += num_es_touched;
    // undo inv_num_outlinks transformation
    for (int i = 0; i < num_vs - num_no_out_vs; ++i)
        x[i] /= inv_num_outlinks[i];
    // normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
    if (normalize) {
        double norm = 0;
        for (int i = 0; i < num_vs; ++i)
            norm += x[i];
        norm = 1/norm;
        for (int i = 0; i < num_vs; ++i)
            x[i] *= norm;
    }
    // return results
    ret->x = prpack_utils::permute(num_vs, x, decoding);
    delete[] x;
    if (uv_exists)
        delete[] uv;
    return ret;
}

prpack_result* prpack_solver::solve_via_schur_gs_uv(
        double alpha,
        double tol,
        int num_vs,
        int num_no_in_vs,
        int num_no_out_vs,
        int num_es,
        int* heads,
        int* tails,
        double* ii,
        double* inv_num_outlinks,
        double* u,
        double* v,
        int* encoding,
        int* decoding) {
    // solve uv = u
    prpack_result* ret_u = solve_via_schur_gs(
            alpha,
            tol,
            num_vs,
            num_no_in_vs,
            num_no_out_vs,
            num_es,
            heads,
            tails,
            ii,
            inv_num_outlinks,
            u,
            encoding,
            decoding,
            false);
    // solve uv = v
    prpack_result* ret_v = solve_via_schur_gs(
            alpha,
            tol,
            num_vs,
            num_no_in_vs,
            num_no_out_vs,
            num_es,
            heads,
            tails,
            ii,
            inv_num_outlinks,
            v,
            encoding,
            decoding,
            false);
    // combine the u and v cases
    return combine_uv(num_vs, inv_num_outlinks, encoding, alpha, ret_u, ret_v);
}

// Gauss-Seidel using strongly connected components.
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
        int* encoding,
        int* decoding,
        bool normalize) {
    prpack_result* ret = new prpack_result();
    // initialize uv values
    double uv_const = 1.0/num_vs;
    int uv_exists = (uv) ? 1 : 0;
    uv = (uv) ? prpack_utils::permute(num_vs, uv, encoding) : &uv_const;
    // initialize the eigenvector
    double* x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        x[i] = uv[uv_exists*i]*inv_num_outlinks[i]/(1 - alpha*ii[i]);
    // create x_outside
    double* x_outside = new double[num_vs];
    // run Gauss-Seidel for (I - alpha*P)*x = uv
    ret->num_es_touched = 0;
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
            ret->num_es_touched += end_j - start_j;
        }
        double err, c;
        do {
            int num_es_touched = 0;
            err = c = 0;
            if (parallelize) {
                // iterate through vertices
                #pragma omp parallel for firstprivate(c) reduction(+:err, num_es_touched) schedule(dynamic, 64)
                for (int i = start_comp; i < end_comp; ++i) {
                    double new_val = x_outside[i];
                    const int start_j = tails_inside[i];
                    const int end_j = (i + 1 != num_vs) ? tails_inside[i + 1] : num_es_inside;
                    for (int j = start_j; j < end_j; ++j)
                        // TODO: might want to use compensation summation for large: end_j - start_j
                        new_val += x[heads_inside[j]];
                    COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i]), c);
                    new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    x[i] = new_val*inv_num_outlinks[i];
                    num_es_touched += end_j - start_j;
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
                    COMPENSATED_SUM(err, fabs(uv[uv_exists*i] + alpha*new_val - (1 - alpha*ii[i])*x[i]/inv_num_outlinks[i]), c);
                    new_val = (alpha*new_val + uv[uv_exists*i])/(1 - alpha*ii[i]);
                    x[i] = new_val*inv_num_outlinks[i];
                    num_es_touched += end_j - start_j;
                }
            }
            // update iteration index
            ret->num_es_touched += num_es_touched;
        } while (err/(1 - alpha) >= tol*(end_comp - start_comp)/num_vs);
    }
    // undo inv_num_outlinks transformation
    for (int i = 0; i < num_vs; ++i)
        x[i] /= inv_num_outlinks[i];
    // normalize x to get the solution for: (I - alpha*P - alpha*u*d')*x = (1 - alpha)*v
    if (normalize) {
        double norm = 0;
        for (int i = 0; i < num_vs; ++i)
            norm += x[i];
        norm = 1/norm;
        for (int i = 0; i < num_vs; ++i)
            x[i] *= norm;
    }
    // return results
    ret->x = prpack_utils::permute(num_vs, x, decoding);
    delete[] x;
    delete[] x_outside;
    if (uv_exists)
        delete[] uv;
    return ret;
}

prpack_result* prpack_solver::solve_via_scc_gs_uv(
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
        double* u,
        double* v,
        int num_comps,
        int* divisions,
        int* encoding,
        int* decoding) {
    // solve uv = u
    prpack_result* ret_u = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            ii,
            inv_num_outlinks,
            u,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // solve uv = v
    prpack_result* ret_v = solve_via_scc_gs(
            alpha,
            tol,
            num_vs,
            num_es_inside,
            heads_inside,
            tails_inside,
            num_es_outside,
            heads_outside,
            tails_outside,
            ii,
            inv_num_outlinks,
            v,
            num_comps,
            divisions,
            encoding,
            decoding,
            false);
    // combine u and v
    return combine_uv(num_vs, inv_num_outlinks, encoding, alpha, ret_u, ret_v);
}

// VARIOUS HELPER METHODS /////////////////////////////////////////////////////////////////////////

// Combine u and v results.
prpack_result* prpack_solver::combine_uv(
        int num_vs,
        double* inv_num_outlinks,
        int* encoding,
        double alpha,
        prpack_result* ret_u,
        prpack_result* ret_v) {
    prpack_result* ret = new prpack_result();
    double delta_u = 0;
    double delta_v = 0;
    for (int i = 0; i < num_vs; ++i) {
        if (inv_num_outlinks[encoding[i]] < 0) {
            delta_u += ret_u->x[i];
            delta_v += ret_v->x[i];
        }
    }
    double s = ((1 - alpha)*alpha*delta_v)/(1 - alpha*delta_u);
    double t = 1 - alpha;
    ret->x = new double[num_vs];
    for (int i = 0; i < num_vs; ++i)
        ret->x[i] = s*ret_u->x[i] + t*ret_v->x[i];
    ret->num_es_touched = ret_u->num_es_touched + ret_v->num_es_touched;
    // clean up and return
    delete ret_u;
    delete ret_v;
    return ret;
}
