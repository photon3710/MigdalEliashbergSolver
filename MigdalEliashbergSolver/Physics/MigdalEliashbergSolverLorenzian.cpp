//
//  MigdalEliashbergSolverLorenzian.cpp
//  MigdalEliashbergSolverLorenzian
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "MigdalEliashbergSolverLorenzian.h"

using namespace std;

/****************************************************************************************************/

double gsl_wrapper_get_min_eig(double alpha2, void * p) {
    
    gsl_ME_func_params * par = (gsl_ME_func_params*) p;
    
    gsl_complex eig = par->ptr_ME_solver->get_kern_min_eig(alpha2);
    par->IM_eig = GSL_IMAG(eig);
    
    return GSL_REAL(eig);
}

/****************************************************************************************************/

MigdalEliashbergSolverLorenzian::MigdalEliashbergSolverLorenzian() : report_stream(cout) {

    sys = new ElectronPhononSystemLorenzian();
    
}


MigdalEliashbergSolverLorenzian::~MigdalEliashbergSolverLorenzian() {

    /* free up the memory */
    free_memory();
    
    delete sys;
    
}

ostream& MigdalEliashbergSolverLorenzian::report() {
    return report_stream << "MigdalEliashbergSolverLorenzian::message: ";
}

void MigdalEliashbergSolverLorenzian::init_el_ph_system(const double & el_temp, const double & n0,
                                               const double & Omega_0, const double & Omega_c,
                                               const double & gamma, const double & mu_c) {
    
    /* instantiate the electron phonon system */
    sys->set_system_params(el_temp, n0, Omega_0, Omega_c, gamma, mu_c);
    
    /* set the member variables */
    this->omega_c = 10*Omega_0 + Omega_c;

}

/****************************************************************************************************/

double MigdalEliashbergSolverLorenzian::RE_K_Z_truncated_Lorenzian(const double & omega,
                                                          const double & omegap) {
    
    double n = sys->ph_dist(0);
    double f = sys->el_dist(omegap);
    double fb = 1.0 - f;
    double Omega_0 = sys->Omega_0;
    
    return (n+fb)*(sys->RE_truncated_Lorenzian(omega - omegap - Omega_0)
                   + sys->RE_truncated_Lorenzian(omega + omegap + Omega_0)) +
        (n+f)*(sys->RE_truncated_Lorenzian(omega - omegap + Omega_0)
               + sys->RE_truncated_Lorenzian(omega + omegap - Omega_0));

}


double MigdalEliashbergSolverLorenzian::IM_K_Z_truncated_Lorenzian(const double & omega,
                                                          const double & omegap) {
    
    double n = sys->ph_dist(0);
    double f = sys->el_dist(omegap);
    double fb = 1.0 - f;
    double Omega_0 = sys->Omega_0;
    
    return (n+fb)*(sys->IM_truncated_Lorenzian(omega - omegap - Omega_0)
                   + sys->IM_truncated_Lorenzian(omega + omegap + Omega_0)) +
    (n+f)*(sys->IM_truncated_Lorenzian(omega - omegap + Omega_0)
           + sys->IM_truncated_Lorenzian(omega + omegap - Omega_0));
    
}


double MigdalEliashbergSolverLorenzian::RE_K_Delta_truncated_Lorenzian(const double & omega,
                                                              const double & omegap) {
    
    double n = sys->ph_dist(0);
    double f = sys->el_dist(omegap);
    double fb = 1.0 - f;
    double Omega_0 = sys->Omega_0;
    
    return (n+fb)*(sys->RE_truncated_Lorenzian(omega + omegap + Omega_0)
                   - sys->RE_truncated_Lorenzian(omega - omegap - Omega_0)) +
    (n+f)*(sys->RE_truncated_Lorenzian(omega + omegap - Omega_0)
           - sys->RE_truncated_Lorenzian(omega - omegap + Omega_0));
    
}


double MigdalEliashbergSolverLorenzian::IM_K_Delta_truncated_Lorenzian(const double & omega,
                                                              const double & omegap) {
    
    double n = sys->ph_dist(0);
    double f = sys->el_dist(omegap);
    double fb = 1.0 - f;
    double Omega_0 = sys->Omega_0;
    
    return (n+fb)*(sys->IM_truncated_Lorenzian(omega + omegap + Omega_0)
                   - sys->IM_truncated_Lorenzian(omega - omegap - Omega_0)) +
    (n+f)*(sys->IM_truncated_Lorenzian(omega + omegap - Omega_0)
           - sys->IM_truncated_Lorenzian(omega - omegap + Omega_0));
    
}

/******** memory management *************************************************************************/

void MigdalEliashbergSolverLorenzian::free_memory() {
    
    if (quad_initialized) {
        
        gsl_vector_free(quad_Z_x);
        gsl_vector_free(quad_Z_w);
        gsl_vector_free(quad_Delta_x);
        gsl_vector_free(quad_Delta_w);

    }
    
    quad_initialized = false;
    
    free_ws_normal();
    
}


int MigdalEliashbergSolverLorenzian::alloc_ws_normal() {
    
    if (!quad_initialized) {
        report() << "The quadrature is not initalized yet!" << endl;
        return EXIT_FAILURE;
    }
    
    if (ws_normal_alloc) {
        report() << "The workspace is already alloced!" << endl;
        return EXIT_FAILURE;
    }
    
    __solver_ws_RE_Q_Z = gsl_vector_alloc(N_Delta);
    __solver_ws_IM_Q_Z = gsl_vector_alloc(N_Delta);
    
    __solver_ws_RE_K_Delta = gsl_matrix_alloc(N_Delta, N_Delta);
    __solver_ws_IM_K_Delta = gsl_matrix_alloc(N_Delta, N_Delta);
    
    __solver_ws_el_fact = gsl_vector_alloc(N_Delta);
    __solver_ws_kern = gsl_matrix_alloc(N_Delta, N_Delta);
    __solver_ws_eig_ws = gsl_eigen_nonsymm_alloc(N_Delta);
    __solver_ws_eig_vec = gsl_vector_complex_alloc(N_Delta);
    
    Delta_crit = gsl_vector_complex_alloc(N_Delta);
    phi_crit = gsl_vector_complex_alloc(N_Delta);
    
    ws_normal_alloc = true;
    
    return EXIT_SUCCESS;
    
}

/****************************************************************************************************/


int MigdalEliashbergSolverLorenzian::free_ws_normal() {
    
    if (ws_normal_alloc) {
        
        gsl_vector_free(__solver_ws_RE_Q_Z);
        gsl_vector_free(__solver_ws_IM_Q_Z);
        
        gsl_matrix_free(__solver_ws_RE_K_Delta);
        gsl_matrix_free(__solver_ws_IM_K_Delta);
        
        gsl_vector_free(__solver_ws_el_fact);
        
        gsl_matrix_free(__solver_ws_kern);
        gsl_eigen_nonsymm_free(__solver_ws_eig_ws);
        gsl_vector_complex_free(__solver_ws_eig_vec);
        
        gsl_vector_complex_free(Delta_crit);
        gsl_vector_complex_free(phi_crit);

        ws_normal_alloc = false;
        ws_normal_calc = false;
        solver_status = EXIT_FAILURE;

    }
    
    return EXIT_SUCCESS;
    
}


int MigdalEliashbergSolverLorenzian::init_quad_normal(const size_t & order_base, const size_t & order_inf,
                                             const size_t & N_rep_inf, const double & eps_f,
                                             const double & eta, const size_t & N_th,
                                             const size_t & N_rep_Z_def,
                                             const size_t & N_rep_Delta_def) {
    
    /* create a quadrature for $\omega'$ integration in $[0, \infty]$
     for calculating $Q(\Omega), and in $[0, \omega_c]$ for \Delta */
    
    double quad_dx, omega_1, req_dx;
    int N_rep_0;
    double * xtab_base, * wtab_base;
    double * xtab_inf, * wtab_inf;
    vector<double> x0_vec, w0_vec;
    vector<double> x1_vec, w1_vec;
    vector<double> x2_vec, w2_vec;
    vector<double> x_vec_t, w_vec_t;
    
    /* calculate the base quadrature */
    xtab_base = new double[order_base];
    wtab_base = new double[order_base];
    kronrod_set(static_cast<int>(order_base), xtab_base, wtab_base);
    
    /************ quadrature for Z ***********/

    /* end point for the first part of integration */
    omega_1 = omega_c + sys->Omega_max + 2*sys->gamma;

    if (N_rep_Z_def > 0) {
        N_rep_0 = static_cast<int>(N_rep_Z_def);
    }
    else {
        /* find the largest spacing between the absicca of the quadrature */
        quad_dx = 0;
        for (size_t j=0; j<order_base-1; j++) {
            quad_dx = max(quad_dx, xtab_base[j+1]-xtab_base[j]);
        }
        
        quad_dx = 2*quad_dx*omega_1; /* 2 comes from the range of quadrature */
        
        /* require a resolution of eta\gamma and N_th points in the effective support of f(\omega') */
        req_dx = min(eta*sys->gamma, sys->el_temp*log(1/eps_f-1)/N_th);
        N_rep_0 = static_cast<int>(ceil(max(1.0,quad_dx/req_dx)));
    }

    
    /* create a GK quadrature for [0, omega_1] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_0, 0, omega_1, x0_vec, w0_vec);
    
    /* create a semi-infinite quadrature for [omega_1, +\infty]; since the kernel has a tail like
       $\omega'^{-2}$ decay, we use the mapping \omega' -> \omega_1 + t/(1-t), create a GK quadrature for the t integration, and finally map back to the original frequency variable */
    xtab_inf = new double[order_inf];
    wtab_inf = new double[order_inf];
    kronrod_set(static_cast<int>(order_inf), xtab_inf, wtab_inf);
    generate_composite_rule(xtab_inf, wtab_inf, order_inf, -1.0, 1.0,
                            N_rep_inf, 0.0, 1.0, x_vec_t, w_vec_t);
    
    /* transform from t to \omega' */
    for (size_t j=0; j<x_vec_t.size(); j++) {
        x1_vec.push_back(omega_1 + x_vec_t[j]/(1.0-x_vec_t[j]));
        w1_vec.push_back(w_vec_t[j]/((1.0-x_vec_t[j])*(1.0-x_vec_t[j])));
    }
    
    /* put the pieces together; NOTE: since GK quadratures have no end-points, we don't worry about them; if other quadratures are used, care must be taken */
    this->N_Z = x0_vec.size() + x1_vec.size();
    this->quad_Z_x = gsl_vector_alloc(N_Z);
    this->quad_Z_w = gsl_vector_alloc(N_Z);
    for (size_t j=0; j<x0_vec.size(); j++) {
        gsl_vector_set(quad_Z_x, j, x0_vec[j]);
        gsl_vector_set(quad_Z_w, j, w0_vec[j]);
    }
    for (size_t j=0; j<x1_vec.size(); j++) {
        gsl_vector_set(quad_Z_x, x0_vec.size()+j, x1_vec[j]);
        gsl_vector_set(quad_Z_w, x0_vec.size()+j, w1_vec[j]);
    }
    
    /************ quadrature for \Delta ***********/
    
    if (N_rep_Delta_def > 0) {
        N_rep_0 = static_cast<int>(N_rep_Delta_def);
    }
    else {
        /* find the largest spacing between the absicca of the quadrature */
        quad_dx = 0;
        for (size_t j=0; j<order_base-1; j++) {
            quad_dx = max(quad_dx, xtab_base[j+1]-xtab_base[j]);
        }
        quad_dx = 2*quad_dx*omega_c; /* 2 comes from the range of quadrature */
        
        /* require a resolution of 0.1\gamma and 20 points in the effective support of f(\omega') */
        req_dx = min(eta*sys->gamma, sys->el_temp*log(1/eps_f-1)/N_th);
        N_rep_0 = static_cast<int>(ceil(max(1.0,quad_dx/req_dx)));
    }
    
    /* create a GK quadrature for [0, omega_c] */
    generate_composite_rule(xtab_base, wtab_base, order_base, -1.0, 1.0,
                            N_rep_0, 0, omega_c, x2_vec, w2_vec);
    
    this->N_Delta = x2_vec.size();
    this->quad_Delta_x = gsl_vector_alloc(N_Delta);
    this->quad_Delta_w = gsl_vector_alloc(N_Delta);
    for (size_t j=0; j<x2_vec.size(); j++) {
        gsl_vector_set(quad_Delta_x, j, x2_vec[j]);
        gsl_vector_set(quad_Delta_w, j, w2_vec[j]);
    }

    quad_initialized = true;
    
    /* free memory */
    delete[] xtab_base;
    delete[] wtab_base;
    delete[] xtab_inf;
    delete[] wtab_inf;
    
    return EXIT_SUCCESS;
    
}


int MigdalEliashbergSolverLorenzian::calc_ws_normal() {
    
    if (!quad_initialized) {
        report() << "The quadrature is not initalized yet!" << endl;
        return EXIT_FAILURE;
    }
    
    if (!ws_normal_alloc) {
        report() << "The workspace is not alloced yet!" << endl;
        return EXIT_FAILURE;
    }
    
    double c_1, c_2, omega, omegap, wp;
    clock_t t0, t1;
    
    /****** calculate Q_Z(\omega) *******/

    if (debug_mode) {
        t0 = clock();
    }
    
    /* NOTE the j-loop on for \omega */
    for (size_t j=0; j<N_Delta; j++) {
        c_1 = 0; c_2 = 0;
        omega = gsl_vector_get(quad_Delta_x, j);
        /* NOTE the k-loop on for \omega' integration of Z */
        for (size_t k=0; k<N_Z; k++) {
            omegap = gsl_vector_get(quad_Z_x, k);
            wp     = gsl_vector_get(quad_Z_w, k);
            c_1 += wp * RE_K_Z_truncated_Lorenzian(omega, omegap);
            c_2 += wp * IM_K_Z_truncated_Lorenzian(omega, omegap);
        }
        gsl_vector_set(__solver_ws_RE_Q_Z, j, c_1/omega);
        gsl_vector_set(__solver_ws_IM_Q_Z, j, c_2/omega);
    }
    
    if (debug_mode) {
        
        t1 = clock();
        
        report() << "calculation of Q_Z took " << (t1-t0)/(double)CLOCKS_PER_SEC
        << " seconds." << endl;
        
        //        mpl->mpl_xy_plot("RE_Q", MPLPyPlotter::gsl_vector_to_STL_vector(quad_Delta_x),
        //                         MPLPyPlotter::gsl_vector_to_STL_vector(__solver_ws_RE_Q_Z), "\\omega", "Re[Q]");
        //        mpl->mpl_xy_plot("IM_Q", MPLPyPlotter::gsl_vector_to_STL_vector(quad_Delta_x),
        //                         MPLPyPlotter::gsl_vector_to_STL_vector(__solver_ws_IM_Q_Z), "\\omega", "Im[Q]");
        //        mpl->show();
        
    }
    
    /****** calculate K_\Delta(\omega, \omega') *******/

    if (debug_mode) {
        t0 = clock();
    }
    
    for (size_t j=0; j<N_Delta; j++) {
        omega = gsl_vector_get(quad_Delta_x, j);
        for (size_t k=0; k<N_Delta; k++) {
            omegap = gsl_vector_get(quad_Delta_x, k);
            gsl_matrix_set(__solver_ws_RE_K_Delta, j, k,
                           RE_K_Delta_truncated_Lorenzian(omega, omegap));
            gsl_matrix_set(__solver_ws_IM_K_Delta, j, k,
                           IM_K_Delta_truncated_Lorenzian(omega, omegap));
        }
    }
    
    if (debug_mode) {
        
        t1 = clock();
        
        report() << "calculation of K_Delta took " << (t1-t0)/(double)CLOCKS_PER_SEC
        << " seconds." << endl;
        
    }

    /****** calculate 1-2*f(\omega') *******/

    if (debug_mode) {
        t0 = clock();
    }

    
    for (size_t j=0; j<N_Delta; j++) {
        omegap = gsl_vector_get(quad_Delta_x, j);
        gsl_vector_set(__solver_ws_el_fact, j, 1-2*sys->el_dist(omegap));
    }
    
    if (debug_mode) {
        
        t1 = clock();
        
        report() << "calculation of 1-2*f took " << (t1-t0)/(double)CLOCKS_PER_SEC
        << " seconds." << endl;
        
    }
    
    ws_normal_calc = true;

    return EXIT_SUCCESS;
    
}


void MigdalEliashbergSolverLorenzian::__solver_eval_kernel(const double & alpha2) {
    
    double Z1, Z2, Z_sq, kern_elem, denom, wp;
    
    for (size_t j=0; j<N_Delta; j++) {
        
        Z1   = 1.0 - alpha2*gsl_vector_get(__solver_ws_RE_Q_Z, j);
        Z2   = -alpha2*gsl_vector_get(__solver_ws_IM_Q_Z, j);
        Z_sq = Z1*Z1 + Z2*Z2;
        
        for (size_t k=0; k<N_Delta; k++) {
            
            denom = gsl_vector_get(quad_Delta_x, k) * Z_sq;
            wp    = gsl_vector_get(quad_Delta_w, k);
                                   
            /* the \int d\omega' part */
            kern_elem = -alpha2 * wp * (gsl_matrix_get(__solver_ws_RE_K_Delta, j ,k) * Z1
                                   + gsl_matrix_get(__solver_ws_IM_K_Delta, j, k) * Z2) / denom
            + sys->mu_c * wp * gsl_vector_get(__solver_ws_el_fact, k) * Z1 / denom;
            
            /* the \delta(\omega-\omega') part */
            if (j==k) kern_elem += 1.0;
            
            gsl_matrix_set(__solver_ws_kern, j, k, kern_elem);
            
        }
    }
    
}


gsl_complex MigdalEliashbergSolverLorenzian::get_kern_min_eig(const double & alpha2) {
    
    /* evaluate the kernel */
    __solver_eval_kernel(alpha2);
    
    /* find the eigenvalues of the kernel and locate its smallest eigenvalue */
    gsl_eigen_nonsymm(__solver_ws_kern, __solver_ws_eig_vec, __solver_ws_eig_ws);
    
    gsl_complex min_eig = gsl_vector_complex_get(__solver_ws_eig_vec, 0);

    for (size_t j=1; j<N_Delta; j++) {
        gsl_complex c_eig = gsl_vector_complex_get(__solver_ws_eig_vec, j);
        if (GSL_REAL(c_eig) < GSL_REAL(min_eig)) {
            min_eig = c_eig;
        }
    }
    
    return min_eig;
    
}


int MigdalEliashbergSolverLorenzian::find_alpha2_critical(const double & alpha2_lower,
                                                 const double & alpha2_upper,
                                                 const long & max_iter,
                                                 const double & epsabs) {
    
    if (!quad_initialized) {
        report() << "The quadrature is not initalized yet!" << endl;
        return EXIT_FAILURE;
    }
    
    if (!ws_normal_calc) {
        report() << "The workspace is not calculated yet!" << endl;
        return EXIT_FAILURE;
    }
    
    int status, exit_status = EXIT_SUCCESS;
    int iter = 0;
    
    gsl_ME_func_params p;
    gsl_function gsl_func_get_min_eig;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    
    p.ptr_ME_solver = this;
    gsl_func_get_min_eig.params = &p;
    gsl_func_get_min_eig.function = &gsl_wrapper_get_min_eig;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc(T);
    gsl_root_fsolver_set(s, &gsl_func_get_min_eig, alpha2_lower, alpha2_upper);
    
    ios::fmtflags old_settings = report_stream.flags();
    streamsize old_precision = report_stream.precision();
    
    report_stream.width(12);
    report_stream.precision(12);
    report_stream.setf(ios::fixed, ios::floatfield);
    
    report() << "starting root finding..." << endl;
    
    do {
        
        iter++;
        status = gsl_root_fsolver_iterate(s);
        
        double r    = gsl_root_fsolver_root(s);
        double x_lo = gsl_root_fsolver_x_lower(s);
        double x_hi = gsl_root_fsolver_x_upper(s);
        
        status = gsl_root_test_interval(x_lo, x_hi, 0.0, epsabs);
        
        report() << "iter: " << iter << ", alpha2_lower: " << x_lo << ", alpha2_upper: " << x_hi
                 << ", alpha2_crit: " << r << endl;
        
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
        alpha2_crit = gsl_root_fsolver_root(s);
    } else {
        report() << "Root could not be found!" << endl;
        alpha2_crit = 0;
        exit_status = EXIT_FAILURE;
    }
    
    report_stream.flags(old_settings);
    report_stream.precision(old_precision);
    
    /* calculate the critical gap function */
    if (status == GSL_SUCCESS) {
        
        report() << "Calculating the gap function..." << endl;
        __solver_eval_kernel(alpha2_crit);
        
        /* workspace for eigenvector calculation */
        gsl_eigen_nonsymmv_workspace * __solver_ws_eig_v_ws = gsl_eigen_nonsymmv_alloc(N_Delta);
        gsl_matrix_complex * __solver_ws_evec = gsl_matrix_complex_alloc(N_Delta, N_Delta);

        gsl_eigen_nonsymmv(__solver_ws_kern, __solver_ws_eig_vec, __solver_ws_evec, __solver_ws_eig_v_ws);
        
        size_t min_eig_idx = 0;
        gsl_complex min_eig = gsl_vector_complex_get(__solver_ws_eig_vec, min_eig_idx);
        
        for (size_t j=1; j<N_Delta; j++) {
            gsl_complex c_eig = gsl_vector_complex_get(__solver_ws_eig_vec, j);
            if (GSL_REAL(c_eig) < GSL_REAL(min_eig)) {
                min_eig_idx = j;
                min_eig = c_eig;
            }
        }
        
        /* store the eigenvalue */
        for (size_t j=0; j<N_Delta; j++) {
            gsl_complex c_Delta = gsl_matrix_complex_get(__solver_ws_evec, j, min_eig_idx);
            gsl_complex c_Z;
            GSL_SET_COMPLEX(&c_Z, 1.0 - alpha2_crit*gsl_vector_get(__solver_ws_RE_Q_Z, j),
                            -alpha2_crit*gsl_vector_get(__solver_ws_IM_Q_Z, j));
            gsl_vector_complex_set(Delta_crit, j, c_Delta);
            gsl_vector_complex_set(phi_crit, j, gsl_complex_mul(c_Delta, c_Z));
        }
        
        /* free the memory */
        gsl_eigen_nonsymmv_free(__solver_ws_eig_v_ws);
        gsl_matrix_complex_free(__solver_ws_evec);
        
    } else {
        
        gsl_vector_complex_set_zero(Delta_crit);
        gsl_vector_complex_set_zero(phi_crit);
    
    }
    
    /* free up the memory */
    gsl_root_fsolver_free(s);
    
    solver_status = exit_status;
    
    return status;
    
}

/****************************************************************************************************/

void MigdalEliashbergSolverLorenzian::save_data() {
    
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
    
    string c_output_path = output_path + "/" + id;

    if (!boost::filesystem::create_directories(boost::filesystem::path(c_output_path))) {
        report() << "can not make the output directory!" << endl;
        return;
    }
    
    report() << "saving the results..." << endl;
    
    ofstream f;
    f.precision(12);
    
    /* write the params */
    f.open(c_output_path + "/params.csv", std::ofstream::out);

    f << "T, " << sys->el_temp << endl;
    f << "n0, " << sys->n0 << endl;
    f << "Omega_0, " << sys->Omega_0 << endl;
    f << "Omega_c, " << sys->Omega_c << endl;
    f << "gamma, " << sys->gamma << endl;
    f << "mu_c, " << sys->mu_c << endl;
    f << "omega_c, " << omega_c << endl;
    
    f.close();
    
    if (solver_status == EXIT_SUCCESS) {
    
        f.open(c_output_path + "/omega.csv", std::ofstream::out);
        helper::write_gsl_vector(f, quad_Delta_x);
        f.close();
        
        f.open(c_output_path + "/gap.csv", std::ofstream::out);
        helper::write_gsl_vector_complex(f, Delta_crit);
        f.close();
        
        f.open(c_output_path + "/sol.csv", std::ofstream::out);
        f << "alpha2_crit, " << alpha2_crit << endl;
        
    }
    
}


/****************************************************************************************************/
/******************************************* tests **************************************************/
/****************************************************************************************************/

int MigdalEliashbergSolverLorenzian::test_quad() {

    double tol = 1e-4;
    int status = EXIT_SUCCESS;
    
    if (!quad_initialized) {
        report() << "The quadrature is not initialized!" << endl;
        status = EXIT_FAILURE;
        return status;
    }
    
    size_t prec = report_stream.precision();
    report_stream.precision(12);
    
    report() << "The quadrature for Z has " << quad_Z_x->size << " points" << endl;
    
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
    
    /* calculate the integral of f(\omega') = 1 */
    sum = 0;
    req = omega_c;
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

