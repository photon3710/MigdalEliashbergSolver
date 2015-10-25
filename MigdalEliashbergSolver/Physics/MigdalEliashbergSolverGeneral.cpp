//
//  MigdalEliashbergSolverGeneral.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "MigdalEliashbergSolverGeneral.h"

using namespace std;

/****************************************************************************************************/

MigdalEliashbergSolverGeneral::MigdalEliashbergSolverGeneral() : report_stream(cout) {
    
    sys = new ElectronPhononSystemGeneral();
    
}

/****************************************************************************************************/

MigdalEliashbergSolverGeneral::~MigdalEliashbergSolverGeneral() {

    free_quad();
    free_ws();
    
}

/****************************************************************************************************/

ostream& MigdalEliashbergSolverGeneral::report() {
    return report_stream << "MigdalEliashbergSolverGeneral::message: ";
}

/****************************************************************************************************/

void MigdalEliashbergSolverGeneral::set_sys_params(vector<double> F_eliash_vec,
                                                   const double & omega_c_el,
                                                   const double & mu_c,
                                                   const double & ph_temp) {
    
    /* cutoff frequency for electronic integrations */
    sys->omega_c_el = omega_c_el;
    sys->mu_c = mu_c;
    this->ph_temp = ph_temp;
    sys->set_ph_spect(F_eliash_vec);
    
    report() << "Total electron-phonon coupling (alpha^2): " << sys->alpha2 << endl;
    
}

/****************************************************************************************************/

/* TODO this is temporary */
void MigdalEliashbergSolverGeneral::set_quad_params(const size_t & order_base,
                                                    const size_t & order_inf,
                                                    const double & eps_f,
                                                    const size_t & N_rep_Z_th,
                                                    const size_t & N_rep_Z_mid,
                                                    const size_t & N_rep_Z_inf,
                                                    const size_t & N_rep_Delta_th,
                                                    const size_t & N_rep_Delta_mid) {

    this->order_base = order_base;
    this->order_inf = order_inf;
    this->eps_f = eps_f;
    this->N_rep_Z_th = N_rep_Z_th;
    this->N_rep_Z_mid = N_rep_Z_mid;
    this->N_rep_Z_inf = N_rep_Z_inf;
    this->N_rep_Delta_th = N_rep_Delta_th;
    this->N_rep_Delta_mid = N_rep_Delta_mid;
    
}

void MigdalEliashbergSolverGeneral::init_el_quad() {
    
    double omega_mid, omega_th;

    double * xtab_base, * wtab_base;
    double * xtab_inf, * wtab_inf;

    vector<double> x0_vec, w0_vec;
    vector<double> x1_vec, w1_vec;
    vector<double> x2_vec, w2_vec;
    vector<double> x3_vec, w3_vec;
    vector<double> x4_vec, w4_vec;
    vector<double> x_vec_t, w_vec_t;
    size_t offset;
    
    /* free the previous quad workspace if one exists */
    if (quad_init) {
        free_quad();
        quad_init = false;
    }
    
    /* calculate the basic quadrature */
    xtab_base = new double[order_base];
    wtab_base = new double[order_base];
    kronrod_set(static_cast<int>(order_base), xtab_base, wtab_base);

    xtab_inf = new double[order_inf];
    wtab_inf = new double[order_inf];
    kronrod_set(static_cast<int>(order_inf), xtab_inf, wtab_inf);
    
    /************ quadrature for Z ***********/
    
    /* end point for the first part of integration */
    /* NOTE: omega_max_ph = 1.0 always; that's the middle term "1.0"; the last term is extra padding */
    omega_mid = sys->omega_c_el + 1.0 + (1.0/sys->N_ph);
    
    /* omega_th is the interval within which f(\omega) falls to eps_f */
    omega_th = sys->el_temp*log(1/eps_f-1);
    
    if (omega_th>omega_mid) {
        report() << "WARNING: temperature is too high and omega_c must be increased. "
            << "The result will be unreliable." << endl;
        omega_th = 0.5 * omega_mid;
    }

    /* create a composite quadrature within [0, omega_th] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_Z_th, 0.0, omega_th, x0_vec, w0_vec);
    
    /* create a composite quadrature within [omega_th, omega_mid] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_Z_mid, omega_th, omega_mid, x1_vec, w1_vec);

    /* create a semi-infinite quadrature for [omega_mid, +\infty]; since the kernel has a tail like
     $\omega'^{-2}$ decay, we use the mapping \omega' -> \omega_1 + t/(1-t), create a quadrature for
     the t integration, and map back to the original frequency variable */
    generate_composite_rule(xtab_inf, wtab_inf, order_inf, -1.0, 1.0,
                            N_rep_Z_inf, 0.0, 1.0, x_vec_t, w_vec_t);
    
    /* transform from t to \omega' */
    for (size_t j=0; j<x_vec_t.size(); j++) {
        x2_vec.push_back(omega_mid + x_vec_t[j]/(1.0-x_vec_t[j]));
        w2_vec.push_back(w_vec_t[j]/((1.0-x_vec_t[j])*(1.0-x_vec_t[j])));
    }
    
    /* put the pieces together;
     NOTE: since Gauss-Kronrod quadratures have no end-points, we don't
     worry about them; if other quadratures are used, care must be taken */

    this->N_Z = x0_vec.size() + x1_vec.size() + x2_vec.size();
    this->quad_Z_x = gsl_vector_alloc(N_Z);
    this->quad_Z_w = gsl_vector_alloc(N_Z);
    
    offset = 0;
    for (size_t j=0; j<x0_vec.size(); j++) {
        gsl_vector_set(quad_Z_x, offset + j, x0_vec[j]);
        gsl_vector_set(quad_Z_w, offset + j, w0_vec[j]);
    }
    offset += x0_vec.size();
    for (size_t j=0; j<x1_vec.size(); j++) {
        gsl_vector_set(quad_Z_x, offset + j, x1_vec[j]);
        gsl_vector_set(quad_Z_w, offset + j, w1_vec[j]);
    }
    offset += x1_vec.size();
    for (size_t j=0; j<x2_vec.size(); j++) {
        gsl_vector_set(quad_Z_x, offset + j, x2_vec[j]);
        gsl_vector_set(quad_Z_w, offset + j, w2_vec[j]);
    }
    
    /************ quadrature for \Delta ***********/
    
    /* omega_th is the interval within which f(\omega) falls to eps_f */
    omega_th = sys->el_temp*log(1/eps_f-1);
    
    if (omega_th>sys->omega_c_el) {
        report() << "WARNING: temperature is too high and omega_c must be increased. "
        << "The result will be unreliable." << endl;
        omega_th = 0.5 * sys->omega_c_el;
    }

    /* create a composite quadrature within [0, omega_th] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_Delta_th, 0.0, omega_th, x3_vec, w3_vec);
    
    /* create a composite quadrature within [omega_th, omega_c] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_Delta_mid, omega_th, sys->omega_c_el, x4_vec, w4_vec);
    
    this->N_Delta = x3_vec.size() + x4_vec.size();
    
    this->quad_Delta_x = gsl_vector_alloc(N_Delta);
    this->quad_Delta_w = gsl_vector_alloc(N_Delta);
    
    offset = 0;
    for (size_t j=0; j<x3_vec.size(); j++) {
        gsl_vector_set(quad_Delta_x, offset + j, x3_vec[j]);
        gsl_vector_set(quad_Delta_w, offset + j, w3_vec[j]);
    }
    offset += x3_vec.size();
    for (size_t j=0; j<x4_vec.size(); j++) {
        gsl_vector_set(quad_Delta_x, offset + j, x4_vec[j]);
        gsl_vector_set(quad_Delta_w, offset + j, w4_vec[j]);
    }
    
    /* set the flag */
    quad_init = true;
    
    /* free memory */
    delete[] xtab_base;
    delete[] wtab_base;
    delete[] xtab_inf;
    delete[] wtab_inf;
    
}

/****************************************************************************************************/

void MigdalEliashbergSolverGeneral::free_quad() {
    
    if (quad_init) {
        
        gsl_vector_free(quad_Z_x);
        gsl_vector_free(quad_Z_w);
        gsl_vector_free(quad_Delta_x);
        gsl_vector_free(quad_Delta_w);
        
        quad_init = false;

    }
    
}

/****************************************************************************************************/

void MigdalEliashbergSolverGeneral::alloc_ws() {

    if (ws_init) {
        report() << "Free the workspace first!" << endl;
        return;
    }
    
    solver_ws_RE_Q_Z = gsl_vector_alloc(N_Delta);
    solver_ws_IM_Q_Z = gsl_vector_alloc(N_Delta);
    
    solver_ws_kern = gsl_matrix_alloc(N_Delta, N_Delta);
    solver_ws_eig_ws = gsl_eigen_nonsymm_alloc(N_Delta);
    solver_ws_eig_vec = gsl_vector_complex_alloc(N_Delta);
    
    Delta_c = gsl_vector_complex_alloc(N_Delta);
    phi_c = gsl_vector_complex_alloc(N_Delta);
    
    ws_init = true;

}

void MigdalEliashbergSolverGeneral::free_ws() {
    
    if (ws_init) {
        
        gsl_vector_free(solver_ws_RE_Q_Z);
        gsl_vector_free(solver_ws_IM_Q_Z);
        
        gsl_matrix_free(solver_ws_kern);
        gsl_eigen_nonsymm_free(solver_ws_eig_ws);
        gsl_vector_complex_free(solver_ws_eig_vec);
        
        gsl_vector_complex_free(Delta_c);
        gsl_vector_complex_free(phi_c);

        ws_init = false;

    }
    
}

void MigdalEliashbergSolverGeneral::calc_kernel() {
    
    if (!quad_init) {
        report() << "The quadrature is not initalized yet!" << endl;
        return;
    }
    
    if (ws_init) {
        free_ws();
    }
    
    /********** memory allocation **********/
    alloc_ws();
    
    /********** evaluation of Q_Z **********/
    
    /* initialize and calculate the system workspace */
    sys->set_omega_L_vec(quad_Delta_x);
    sys->set_omega_R_vec(quad_Z_x);
    sys->update_el_dist();
    
    double c_1, c_2, omega;
    
    /* NOTE the j-loop on for \omega */
    for (size_t j=0; j<N_Delta; j++) {

        c_1 = 0;
        c_2 = 0;
        
        /* NOTE the k-loop on for \omega' integration of Z */
        for (size_t k=0; k<N_Z; k++) {
            c_1 += sys->RE_K_Z(j, k) * gsl_vector_get(quad_Z_w, k);
            c_2 += sys->IM_K_Z(j, k) * gsl_vector_get(quad_Z_w, k);
        }
        
        omega = gsl_vector_get(quad_Delta_x, j);
        
        gsl_vector_set(solver_ws_RE_Q_Z, j, c_1 / omega);
        gsl_vector_set(solver_ws_IM_Q_Z, j, c_2 / omega);
    
    }
    
    if (debug_mode && plot_debug && mpl != NULL) {
        
        mpl->mpl_xy_plot("RE_Q", helper::gsl_vector_to_STL_vector(quad_Delta_x),
                         helper::gsl_vector_to_STL_vector(solver_ws_RE_Q_Z), "\\omega", "Re[Q]");
        mpl->mpl_xy_plot("IM_Q", helper::gsl_vector_to_STL_vector(quad_Delta_x),
                         helper::gsl_vector_to_STL_vector(solver_ws_IM_Q_Z), "\\omega", "Im[Q]");
    }
    
    /********** evaluation of \Delta kernel **********/
    sys->set_omega_R_vec(quad_Delta_x);
    sys->update_el_dist();
    
    double Z1, Z2, Z_sq, K_Delta_1, K_Delta_2, kern_elem, denom, wp;
    
    for (size_t j=0; j<N_Delta; j++) {
        
        Z1   = 1.0 - gsl_vector_get(solver_ws_RE_Q_Z, j);
        Z2   = -gsl_vector_get(solver_ws_IM_Q_Z, j);
        Z_sq = Z1*Z1 + Z2*Z2;
        
        for (size_t k=0; k<N_Delta; k++) {
            
            denom     = gsl_vector_get(quad_Delta_x, k) * Z_sq;
            wp        = gsl_vector_get(quad_Delta_w, k);
            K_Delta_1 = sys->RE_K_Delta(j, k);
            K_Delta_2 = sys->IM_K_Delta(j, k);
            
            /* the \int d\omega' part */
            kern_elem = -wp * (K_Delta_1 * Z1 + K_Delta_2 * Z2) / denom
            + sys->mu_c * wp * (1-2*sys->f_vec[k]) * Z1 / denom;
            
            /* the \delta(\omega-\omega') part */
            if (j==k) kern_elem += 1.0;
            
            gsl_matrix_set(solver_ws_kern, j, k, kern_elem);
            
        }
    }
    
}

/****************************************************************************************************/

gsl_complex MigdalEliashbergSolverGeneral::eval_min_eig(const double & el_temp, const double & ph_temp) {

    /* set the temperatures */
    sys->set_el_temp(el_temp);
    sys->set_ph_temp(ph_temp);
    sys->update_ph_dist();
    sys->update_interp_coeffs();
    
    /* init the quadrature */
    init_el_quad();
    
    /* calculate the kernel matrix elements */
    calc_kernel();
    
    /* find the eigenvalues of the kernel and locate its smallest eigenvalue */
    gsl_eigen_nonsymm(solver_ws_kern, solver_ws_eig_vec, solver_ws_eig_ws);
    
    gsl_complex min_eig = gsl_vector_complex_get(solver_ws_eig_vec, 0);
    
    for (size_t j=1; j<N_Delta; j++) {
        gsl_complex c_eig = gsl_vector_complex_get(solver_ws_eig_vec, j);
        if (GSL_REAL(c_eig) < GSL_REAL(min_eig)) {
            min_eig = c_eig;
        }
    }

    return min_eig;

}

/****************************************************************************************************/

double gsl_wrapper_general_get_min_eig(double T, void * p) {
    
    gsl_MigdalEliashbergSolverGeneral_func_params * par =
    (gsl_MigdalEliashbergSolverGeneral_func_params*) p;
    
    gsl_complex eig;
    if (par->ph_temp>=0) {
        eig = par->ptr_ME_solver->eval_min_eig(T, par->ph_temp);
    } else {
        eig = par->ptr_ME_solver->eval_min_eig(T, T);
    }
    par->IM_eig = GSL_IMAG(eig);
    
    return GSL_REAL(eig);

}

int MigdalEliashbergSolverGeneral::find_Tc(const double & Tc_lower, const double & Tc_upper,
                                             const long & max_iter, const double & epsabs) {
    
    int status, exit_status = EXIT_SUCCESS;
    int iter = 0;
    
    gsl_MigdalEliashbergSolverGeneral_func_params p;
    gsl_function gsl_func_get_min_eig;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    p.ptr_ME_solver = this;
    p.ph_temp = ph_temp;
    gsl_func_get_min_eig.params = &p;
    gsl_func_get_min_eig.function = &gsl_wrapper_general_get_min_eig;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &gsl_func_get_min_eig, Tc_lower, Tc_upper);
    
    ios::fmtflags old_settings = report_stream.flags();
    streamsize old_precision = report_stream.precision();
    
    report_stream.width(12);
    report_stream.precision(12);
    report_stream.setf(ios::fixed, ios::floatfield);

    clock_t t0 = clock();

    report() << "starting root finding..." << endl;
    
    do {
        
        iter++;
        status = gsl_root_fsolver_iterate(s);
        
        double r    = gsl_root_fsolver_root(s);
        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);
        
        status = gsl_root_test_interval(x_lo, x_hi, 0.0, epsabs);
        
        report() << "iter: " << iter << ", Tc_lower: " << x_lo << ", Tc_upper: " << x_hi
        << ", Tc: " << r << endl;
        
    }
    
    while (status == GSL_CONTINUE && iter < max_iter);
    
    if (status == GSL_SUCCESS) {
        report() << "Root found!" << endl;
        if (abs(p.IM_eig)<epsabs) {
            report() << "Acceptable solution: the imaginary part of the eigenvalue is " << p.IM_eig << endl;
            exit_status = EXIT_SUCCESS;
        } else {
            report() << "Unacceptable solution: the imaginary part of the eigenvalue is " << p.IM_eig << endl;
            exit_status = EXIT_FAILURE;
        }
        Tc = gsl_root_fsolver_root(s);
    } else {
        report() << "Root could not be found!" << endl;
        Tc = 0;
        exit_status = EXIT_FAILURE;
    }
    
    report_stream.flags(old_settings);
    report_stream.precision(old_precision);
    
    /* calculate the critical gap function */
    if (status == GSL_SUCCESS) {
        
        report() << "Calculating the gap function..." << endl;
        sys->set_el_temp(Tc);
        sys->set_ph_temp(Tc);
        sys->update_ph_dist();
        sys->update_interp_coeffs();
        init_el_quad();
        calc_kernel();

        /* workspace for eigenvector calculation */
        gsl_eigen_nonsymmv_workspace * solver_ws_eig_v_ws = gsl_eigen_nonsymmv_alloc(N_Delta);
        gsl_matrix_complex * solver_ws_evec = gsl_matrix_complex_alloc(N_Delta, N_Delta);
        
        gsl_eigen_nonsymmv(solver_ws_kern, solver_ws_eig_vec, solver_ws_evec, solver_ws_eig_v_ws);
        
        size_t min_eig_idx = 0;
        gsl_complex min_eig = gsl_vector_complex_get(solver_ws_eig_vec, min_eig_idx);
        
        for (size_t j=1; j<N_Delta; j++) {
            gsl_complex c_eig = gsl_vector_complex_get(solver_ws_eig_vec, j);
            if (GSL_REAL(c_eig) < GSL_REAL(min_eig)) {
                min_eig_idx = j;
                min_eig = c_eig;
            }
        }
        
        /* store the eigenvalue */
        for (size_t j=0; j<N_Delta; j++) {
            gsl_complex c_Delta = gsl_matrix_complex_get(solver_ws_evec, j, min_eig_idx);
            gsl_complex c_Z;
            GSL_SET_COMPLEX(&c_Z, 1.0 - gsl_vector_get(solver_ws_RE_Q_Z, j),
                            -gsl_vector_get(solver_ws_IM_Q_Z, j));
            gsl_vector_complex_set(Delta_c, j, c_Delta);
            gsl_vector_complex_set(phi_c, j, gsl_complex_mul(c_Delta, c_Z));
        }
        
        /* free the memory */
        gsl_eigen_nonsymmv_free(solver_ws_eig_v_ws);
        gsl_matrix_complex_free(solver_ws_evec);
        
    } else {
        
        gsl_vector_complex_set_zero(Delta_c);
        gsl_vector_complex_set_zero(phi_c);
        
    }
    
    clock_t t1 = clock();
    
    report() << "Finding Tc took " << (t1-t0)/(double)CLOCKS_PER_SEC << " seconds." << endl;
    
    /* free up the memory */
    gsl_root_fsolver_free(s);
    
    solver_status = exit_status;
    
    return status;
    
}

/****************************************************************************************************/

void MigdalEliashbergSolverGeneral::save_data() {
    
    /* make a random ID */
    size_t id_len = 12;
    srand(time(NULL));
    
    auto randchar = []() -> char {
        const char charset[] = "0123456789abcdefghijklmnopqrstuvwxyz";
        const size_t max_index = sizeof(charset) - 1;
        return charset[rand() % max_index];
    };
    
    string id(id_len, 0);
    std::generate_n(id.begin(), id_len, randchar);
    
    // string c_output_path = output_path + "/" + id;
    string c_output_path = output_path;
    
    if (!boost::filesystem::create_directories(boost::filesystem::path(c_output_path))) {
        report() << "can not make the output directory!" << endl;
        return;
    }
    
    report() << "saving the results..." << endl;
    
    ofstream f;
    f.precision(12);
    
    /* write the params */
    f.open(c_output_path + "/params.csv", std::ofstream::out);
    
    f << "mu_c, " << sys->mu_c << endl;
    f << "omega_c, " << sys->omega_c_el << endl;
    f << "N_ph, " << sys->N_ph << endl;
    f << "alpha2, " << sys->alpha2 << endl;
    
    f.close();
    
    if (solver_status == EXIT_SUCCESS) {
        
        f.open(c_output_path + "/omega.csv", std::ofstream::out);
        helper::write_gsl_vector(f, quad_Delta_x);
        f.close();
        
        f.open(c_output_path + "/gap.csv", std::ofstream::out);
        helper::write_gsl_vector_complex(f, Delta_c);
        f.close();
        
        f.open(c_output_path + "/Omega_ph.csv", std::ofstream::out);
        helper::write_STL_vector(f, sys->Omega_ph_vec);
        f.close();

        f.open(c_output_path + "/eliash.csv", std::ofstream::out);
        helper::write_STL_vector(f, sys->F_eliash_vec);
        f.close();

        f.open(c_output_path + "/sol.csv", std::ofstream::out);
        f << "Tc, " << Tc << endl;
        
    }
    
}





/****************************************************************************************************/
/******************************************* tests **************************************************/
/****************************************************************************************************/

int MigdalEliashbergSolverGeneral::test_quad() {
    
    double tol = 1e-4;
    int status = EXIT_SUCCESS;
    
    size_t prec = report_stream.precision();
    report_stream.precision(12);
    
    report() << "The quadrature for Z has " << quad_Z_x->size << " points" << endl;

    if (debug_mode && plot_debug) {
        vector<double> idx_vec;
        for (size_t j=0; j<quad_Z_x->size; j++) {
            idx_vec.push_back(static_cast<double>(j));
        }
        if (mpl != NULL) {
            mpl->mpl_xy_plot("quad_Z", idx_vec, helper::gsl_vector_to_STL_vector(quad_Z_x),
                             "index", "quad_Z grid points");
        }
    }

    double sum, req;
    
    /* calculate the integral of f(\omega') = 1/(1+x^2) */
    sum = 0;
    req = 0.5*M_PI;
    for (size_t j=0; j<quad_Z_x->size; j++) {
        double x = gsl_vector_get(quad_Z_x, j);
        double w = gsl_vector_get(quad_Z_w, j);
        sum += w/(1+x*x);
    }
    
    report() << "calculated: " << sum << ", required: " << req << ", error: "
    << abs(sum-req) << endl;
    
    if (abs(sum-req)>tol) { status = EXIT_FAILURE; }
    
    /* calculate the integral of f(\omega') = 1/(1+x^2) */
    sum = 0;
    req = 2*M_PI/(3*sqrt(3));
    for (size_t j=0; j<quad_Z_x->size; j++) {
        double x = gsl_vector_get(quad_Z_x, j);
        double w = gsl_vector_get(quad_Z_w, j);
        sum += w*x/(1+x*x*x);
    }
    
    report() << "calculated: " << sum << ", required: " << req << ", error: "
    << abs(sum-req) << endl;
    
    if (abs(sum-req)>tol) { status = EXIT_FAILURE; }
    
    report() << "The quadrature for Delta has " << quad_Delta_x->size << " points" << endl;
    
    if (debug_mode && plot_debug) {
        vector<double> idx_vec;
        for (size_t j=0; j<quad_Delta_x->size; j++) {
            idx_vec.push_back(static_cast<double>(j));
        }
        if (mpl != NULL) {
            mpl->mpl_xy_plot("quad_Delta", idx_vec, helper::gsl_vector_to_STL_vector(quad_Delta_x),
                             "index", "quad_Delta grid points");
        }
    }
    
    /* calculate the integral of f(\omega') = 1 */
    sum = 0;
    req = sys->omega_c_el;
    for (size_t j=0; j<quad_Delta_x->size; j++) {
        double x = gsl_vector_get(quad_Delta_x, j);
        double w = gsl_vector_get(quad_Delta_w, j);
        sum += w;
    }
    
    report() << "calculated: " << sum << ", required: " << req << ", error: "
    << abs(sum-req) << endl;
    
    if (abs(sum-req)>tol) { status = EXIT_FAILURE; }
    
    if (status==EXIT_SUCCESS) { report() << "Quadratures pass the test!" << endl; }
    
    report_stream.precision(prec);
    
    return status;
}






