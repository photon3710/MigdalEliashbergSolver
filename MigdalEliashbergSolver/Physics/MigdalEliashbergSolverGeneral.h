//
//  MigdalEliashbergSolverGeneral.h
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef __MigdalEliashbergSolver__MigdalEliashbergSolverGeneral__
#define __MigdalEliashbergSolver__MigdalEliashbergSolverGeneral__

#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include <functional>
#include <algorithm>

#include <boost/filesystem.hpp>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_linalg.h>

#include "Quadrature.h"
#include "ElectronPhononSystemGeneral.h"
#include "MPLPyPlotter.h"

#include <Accelerate/Accelerate.h>

struct gsl_MigdalEliashbergSolverGeneral_func_params;

class MigdalEliashbergSolverGeneral {
//private:
public:
    
    /* flags */
    bool quad_init = false;
    bool ws_init = false;
    
    /* quadrature */
    gsl_vector * quad_Z_x;
    gsl_vector * quad_Z_w;
    gsl_vector * quad_Delta_x;
    gsl_vector * quad_Delta_w;
    size_t N_Z, N_Delta;

    /* workspace */
    gsl_vector * solver_ws_RE_Q_Z;
    gsl_vector * solver_ws_IM_Q_Z;
    
    gsl_matrix * solver_ws_kern;

    gsl_eigen_nonsymm_workspace * solver_ws_eig_ws;
    gsl_vector_complex * solver_ws_eig_vec;
    
    /* subroutines */
    /* create a quadrature for $\omega'$ integration in $[0, \infty]$ for calculating $Q(\Omega)$, and in
     $[0, \omega_c]$ for calculting $\Delta(\omega)$ */
    void init_el_quad();
    void free_quad();
    
    void alloc_ws();
    void free_ws();
    void calc_kernel();

public:
    
    bool debug_mode = true;
    bool plot_debug = false;
    MPLPyPlotter * mpl = NULL;
    std::string output_path;
    std::ostream & report_stream;
    
    /* the system */
    ElectronPhononSystemGeneral * sys;
    double ph_temp;
    
    /* quadrature params */
    size_t order_base;
    size_t order_inf;
    double eps_f;
    
    size_t N_rep_Z_th;
    size_t N_rep_Z_mid;
    size_t N_rep_Z_inf;
    size_t N_rep_Delta_th;
    size_t N_rep_Delta_mid;

    /* solution containers */
    bool solver_status = EXIT_FAILURE;
    double Tc;
    gsl_vector_complex * Delta_c; /* the found gap function at criticality */
    gsl_vector_complex * phi_c; /* phi = Z*\Delta */
    
    /* methods */
    void debug_on() { debug_mode = true; }
    void set_mpl_plotter(MPLPyPlotter * mpl) { this->mpl = mpl; };
    void set_output_path(const std::string & path) { this->output_path = path; }
    void save_data();
    
    std::ostream& report();
    
    /* system setup */
    void set_sys_params(std::vector<double> F_eliash_vec, const double & omega_c_el,
                        const double & mu_c, const double & ph_temp);

    void set_quad_params(const size_t & order_base = 15,
                         const size_t & order_inf = 15,
                         const double & eps_f = 1e-2,
                         const size_t & N_rep_Z_th = 20,
                         const size_t & N_rep_Z_mid = 20,
                         const size_t & N_rep_Z_inf = 20,
                         const size_t & N_rep_Delta_th = 20,
                         const size_t & N_rep_Delta_mid = 20);
    
    gsl_complex eval_min_eig(const double & el_temp, const double & ph_temp);
    
    int find_Tc(const double & Tc_lower, const double & Tc_upper,
                const long & max_iter, const double & epsabs);
    
    MigdalEliashbergSolverGeneral();
    ~MigdalEliashbergSolverGeneral();
    
    /************* test functions ***************/
    int test_quad();
    
};

struct gsl_MigdalEliashbergSolverGeneral_func_params {
    MigdalEliashbergSolverGeneral * ptr_ME_solver;
    double IM_eig, ph_temp;
};


#endif /* defined(__MigdalEliashbergSolver__MigdalEliashbergSolverGeneral__) */
