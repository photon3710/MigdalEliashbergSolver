//
//  ElectronPhononSystem.h
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef __MigdalEliashbergSolver__ElectronPhononSystem__
#define __MigdalEliashbergSolver__ElectronPhononSystem__

#include <stdio.h>

/* the electron phonon class definition */
class ElectronPhononSystem {
public:
    
    /* there are pure abstract functions */
    virtual double eliash_func(const double &) = 0; /* Eliashberg function: F(\Omega) */
    virtual double el_dist(const double &) = 0; /* Electron distribution function: f(\omega) */
    virtual double ph_dist(const double &) = 0; /* Phonon distribution function: n(\Omega) */
    
    ElectronPhononSystem() {}
    
    ~ElectronPhononSystem() {}

};

#endif /* defined(__MigdalEliashbergSolver__ElectronPhononSystem__) */
