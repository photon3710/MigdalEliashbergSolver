//
//  ElectronPhononSystemLorenzian.h
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef __MigdalEliashbergSolver__ElectronPhononSystemLorenzian__
#define __MigdalEliashbergSolver__ElectronPhononSystemLorenzian__

#include <cmath>
#include <vector>
#include <stdio.h>
#include "ElectronPhononSystem.h"
#include "Quadrature.h"

class ElectronPhononSystemLorenzian : public ElectronPhononSystem {
public:
    
    double el_temp; /* temperature of electrons */
    double n0; /* constant phonon distribution function */
    
    /* Eliashberg function parameters */
    double A_norm, Omega_0, gamma, Omega_c, gamma2, Omega_c2, c0;
    
    double Omega_max; /* maximum phonon frequency */
    double mu_c; /* Coulomb repulsion */
    
    /* the virtual functions */
    double eliash_func(const double &); /* Eliashberg function: F(\Omega) */
    double el_dist(const double &); /* Electron distribution function: f(\omega) */
    double ph_dist(const double &); /* Phonon distribution function: n(\Omega) */
    
    /* Cauchy transform of the truncated Lorenzian */
    double RE_truncated_Lorenzian(const double &);
    double IM_truncated_Lorenzian(const double &);
        
    /* constructor and deconstructor */
    ElectronPhononSystemLorenzian();
    ~ElectronPhononSystemLorenzian();
    
    void set_system_params(const double & el_temp, const double & n0,
                           const double & Omega_0, const double & Omega_c,
                           const double & gamma, const double & mu_c);
    
    /************* test functions ***************/
    void test_calculate_total_eliash_weight();
    
};

#endif /* defined(__MigdalEliashbergSolver__ElectronPhononSystemLorenzian__) */
