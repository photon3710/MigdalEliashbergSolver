//
//  ElectronPhononSystemGeneral.h
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//
//
//  NOTE:
//
//  The Eliashberg function $\alpha^2 F(\Omega)$ samples F_eliash_vec[0 ... N]
//  must have its two first elements and the last element set to 0. We simply
//  ignore these elements in the implementation.
//
//  This restriction is due to the pathological behavior of the
//  low energy part of the phonon spectrum, leading to $\log(z)$ terms
//  that require a careful treatment of numerical integrations. We would like
//  to avoid this complication. Out convention of having the first two samples
//  set to zero effectively implies an IR cut-off on the Eliashberg function.
//  For large enough $N$, this restriction is immaterial.
//

#ifndef __MigdalEliashbergSolver__ElectronPhononSystemGeneral__
#define __MigdalEliashbergSolver__ElectronPhononSystemGeneral__

#include <cmath>
#include <vector>
#include <stdio.h>

#include "gsl/gsl_vector.h"

#include "ElectronPhononSystem.h"
#include "HelperModule.h"

class ElectronPhononSystemGeneral : public ElectronPhononSystem {
public:

    double el_temp = 0; /* temperature of electrons */
    double ph_temp = 0; /* temperature of phohons */
    double mu_c = 0; /* Coulomb repulsion */
    double omega_c_el = 10.0; /* default: 10 x max phonon freq */
    double alpha2;
    
    size_t N_ph = 0; /* number of phonon spectrum subdivisions */
    
    std::vector<double> omega_L_vec; /* omega_L = omega */
    std::vector<double> omega_R_vec; /* omega_R = omega' */
    std::vector<double> Omega_ph_vec; /* phonon freq */
    
    std::vector<double> f_vec; /* electron dist */
    std::vector<double> n_vec; /* phonon dist */
    
    std::vector<double> F_eliash_vec; /* Eiashberg samples */

    std::vector<double> F_a_vec; /* slopes of the linear interpolant of $F(\Omega)$ */
    std::vector<double> F_b_vec; /* intercepts of the linear interpolant of $F(\Omega)$ */
    std::vector<double> F_alpha_vec; /* coefficients of the log terms */

    std::vector<double> nF_a_vec; /* slopes of the linear interpolant of $n(\Omega) F(\Omega)$ */
    std::vector<double> nF_b_vec; /* intercepts of the linear interpolant of $n(\Omega) F(\Omega)$*/
    std::vector<double> nF_alpha_vec; /* coefficients of the log terms */
    
    /* the virtual functions */
    double eliash_func(const double &); /* Eliashberg function: F(\Omega) */
    double el_dist(const double &); /* Electron distribution function: f(\omega) */
    double ph_dist(const double &); /* Phonon distribution function: n(\Omega) */
    
    /* set physical quantities */
    void set_el_temp(const double & el_temp) { this->el_temp = el_temp; }
    void set_ph_temp(const double & ph_temp) { this->ph_temp = ph_temp; }
    
    void set_ph_spect(const std::vector<double> & F_eliash_vec);
    void oversample_ph_spect(const size_t & n_os);
    
    void set_omega_L_vec(gsl_vector * omega_L_vec);
    void set_omega_R_vec(gsl_vector * omega_R_vec);
    
    /* update the distributions */
    void update_el_dist();
    void update_ph_dist();
    void update_interp_coeffs();
    void update_F_interp_coeffs();
    void update_nF_interp_coeffs();
    
    /* Cauchy transform of a trapezoidal piece of Eliashberg function */
    static double xlog(const double & x);
    double RE_I(const double & a, const double & b,
                const double & Omega_0, const double & Omega_1, const double & z);
    double IM_I(const double & a, const double & b,
                const double & Omega_0, const double & Omega_1, const double & z);
    double RE_J(const double & a, const double & b,
                const double & Omega_0, const double & Omega_1, const double & z);
    double IM_J(const double & a, const double & b,
                const double & Omega_0, const double & Omega_1, const double & z);
    
    /* The RE and IM parts of the Cauchy trasform of $\alpha^2 F(\Omega)$ */
    double RE_I_F_full(const double & z);
    double IM_I_F_full(const double & z);
    double RE_J_F_full(const double & z);
    double IM_J_F_full(const double & z);

    /* The RE and IM parts of the Cauchy trasform of $n(\Omega) \alpha^2 F(\Omega)$ */
    double RE_I_nF_full(const double & z);
    double IM_I_nF_full(const double & z);
    double RE_J_nF_full(const double & z);
    double IM_J_nF_full(const double & z);

    /* the Migdal-Eliashberg kernels */
    double RE_K_Z(const size_t & L_idx, const size_t & R_idx);
    double IM_K_Z(const size_t & L_idx, const size_t & R_idx);
    double RE_K_Delta(const size_t & L_idx, const size_t & R_idx);
    double IM_K_Delta(const size_t & L_idx, const size_t & R_idx);
    
    /* constructor and deconstructor */
    ElectronPhononSystemGeneral();
    ~ElectronPhononSystemGeneral();
    
    /************* test functions ***************/
    void calculate_total_eliash_weight();
    
};

#endif /* defined(__MigdalEliashbergSolver__ElectronPhononSystemGeneral__) */
