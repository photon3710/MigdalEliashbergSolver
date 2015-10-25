//
//  MPLPyPlotter.h
//  A c++ interface to MPLPyPlotter python code
//
//  Created by Mehrtash Babadi on 9/26/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef MPLPyPlotter_h
#define MPLPyPlotter_h

#include <iostream>
#include <vector>
#include <cstdio>

#include <gsl/gsl_vector.h>

#include "HelperModule.h"

class MPLPyPlotter {
public:
    
    MPLPyPlotter(const std::string & python_code, const std::string & python_exec) :
    python_code(python_code), python_exec(python_exec), process_running(false) {}
    
    ~MPLPyPlotter();
    
    int start();
    int stop();
    int push(const std::string & cmd);
    int flush();
    int show();
    
    int mpl_xy_plot(const std::string & plot_name, const std::vector<double> & x,
                    const std::vector<double> & y, const std::string & x_label = "x",
                    const std::string & y_label = "y");
    
    static std::vector<double> gsl_vector_to_STL_vector(gsl_vector * gsl_vec);

private:
    
    std::string python_code;
    std::string python_exec;
    bool process_running;
    FILE * proc;

};

#endif
