//
//  MigdalEliashbergSolverLorenzian.h
//  MigdalEliashbergSolverLorenzian
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

// TODO: setup an integration grid for Z in the normal phase
// TODO: calculate Z in the normal phase
//

#ifndef __MigdalEliashbergSolverLorenzian__MigdalEliashbergSolverLorenzian__
#define __MigdalEliashbergSolverLorenzian__MigdalEliashbergSolverLorenzian__

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

#include "HelperModule.h"
#include "Quadrature.h"
#include "ElectronPhononSystem.h"
#include "ElectronPhononSystemLorenzian.h"
#include "MPLPyPlotter.h"

#include <Accelerate/Accelerate.h>

struct gsl_ME_func_params;

class MigdalEliashbergSolverLorenzian {
private:
    
    bool quad_initialized = false;
    bool ws_normal_calc = false;
    bool ws_normal_alloc = false;
    bool solver_status = EXIT_FAILURE;
    
    gsl_vector * __solver_ws_RE_Q_Z;
    gsl_vector * __solver_ws_IM_Q_Z;
    
    gsl_matrix * __solver_ws_RE_K_Delta;
    gsl_matrix * __solver_ws_IM_K_Delta;
    gsl_vector * __solver_ws_el_fact; /* 1 - 2*f(\omega') */
    
    gsl_matrix * __solver_ws_kern;
    gsl_eigen_nonsymm_workspace * __solver_ws_eig_ws;
    gsl_vector_complex * __solver_ws_eig_vec;

    void __solver_eval_kernel(const double & alpha2);
    
public:
    
    /* workspace */

    bool debug_mode = false;
    MPLPyPlotter * mpl;
    std::string output_path;
    std::ostream & report_stream;
    
    ElectronPhononSystemLorenzian * sys;
    
    double omega_c; /* $\omega_c$: cutoff frequency for gap integrals */

    /* quadrature */
    gsl_vector * quad_Z_x;
    gsl_vector * quad_Z_w;
    gsl_vector * quad_Delta_x;
    gsl_vector * quad_Delta_w;
    size_t N_Z, N_Delta;
    
    /* solution containers */
    double alpha2_crit; /* storage for the found critical alpha^2 */
    gsl_vector_complex * Delta_crit; /* the found gap function at criticality */
    gsl_vector_complex * phi_crit; /* phi = Z*\Delta */
    
    /* methods */
    void debug_on() { debug_mode = true; }
    void set_mpl_plotter(MPLPyPlotter * mpl) { this->mpl = mpl; };
    void set_output_path(const std::string & path) { this->output_path = path; }
    void save_data();
    
    std::ostream& report();
    
    void init_el_ph_system(const double & el_temp, const double & n0,
                           const double & Omega_0, const double & Omega_c,
                           const double & gamma, const double & mu_c);
    
    /* ME kernal functions for a truncated Lorenzian */
    /* Re[K_Z(\omega,\omega') */
    double RE_K_Z_truncated_Lorenzian(const double & omega, const double & omegap);
    double IM_K_Z_truncated_Lorenzian(const double & omega, const double & omegap);
    double RE_K_Delta_truncated_Lorenzian(const double & omega, const double & omegap);
    double IM_K_Delta_truncated_Lorenzian(const double & omega, const double & omegap);
    
    /* routines for calculating T_c */
    int init_quad_normal(const size_t & order_base, const size_t & order_inf,
                         const size_t & N_rep_inf, const double & eps_f,
                         const double & eta, const size_t & N_th,
                         const size_t & N_rep_Z_def = 0,
                         const size_t & N_rep_Delta_def = 0);
    
    int alloc_ws_normal();
    int calc_ws_normal();
    int free_ws_normal();
    
    gsl_complex get_kern_min_eig(const double & alpha2);
    double get_kern_det(const double & alpha2);
    int find_alpha2_critical(const double & alpha2_lower,
                             const double & alpha2_upper,
                             const long & max_iter = 100,
                             const double & epsabs = 1e-4);
    
    /* constructor and destructors */
    void free_memory();
    
    MigdalEliashbergSolverLorenzian();
    ~MigdalEliashbergSolverLorenzian();
    
    /************* test functions ***************/
    int test_quad();
    
};

struct gsl_ME_func_params {
    MigdalEliashbergSolverLorenzian * ptr_ME_solver;
    double IM_eig;
};

void write_gsl_vector(std::ofstream & of, gsl_vector * vec);
void write_gsl_vector_complex(std::ofstream & of, gsl_vector_complex * vec);


#endif /* defined(__MigdalEliashbergSolverLorenzian__MigdalEliashbergSolverLorenzian__) */
