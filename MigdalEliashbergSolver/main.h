//
//  main.h
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef MigdalEliashbergSolver_main_h
#define MigdalEliashbergSolver_main_h

#include <iostream>
#include <iomanip>
#include <cmath>
#include <ctime>
#include <cstdio>
#include <vector>
#include <algorithm>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_linalg.h>

#include <libconfig.h++>

#include "HelperModule.h"
#include "MPLPyPlotter.h"

#include "ElectronPhononSystem.h"

#include "MigdalEliashbergSolverLorenzian.h"
#include "ElectronPhononSystemLorenzian.h"

#include "MigdalEliashbergSolverGeneral.h"
#include "ElectronPhononSystemGeneral.h"


#include <Accelerate/Accelerate.h>

#endif
