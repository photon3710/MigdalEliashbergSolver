//
//  MPLPyPlotter.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/26/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "MPLPyPlotter.h"

using namespace std;

int MPLPyPlotter::start() {
    
    std::cout << "starting MPLPyPlotter..." << std::endl;

    std::string cmd = python_exec + " -i " + python_code;
    proc = popen(cmd.c_str(), "w");
    process_running = true;
    
    return EXIT_SUCCESS;
}

int MPLPyPlotter::stop() {
    
    if (process_running==true) {
        std::cout << "stopping MPLPyPlotter..." << std::endl;
        fprintf(proc, "exit\n");
        pclose(proc);
        process_running = false;
    }
    
    return EXIT_SUCCESS;
}

MPLPyPlotter::~MPLPyPlotter() {
    stop();
}

int MPLPyPlotter::push(const std::string & cmd) {
    fprintf(proc, (cmd + "\n").c_str());
    return EXIT_SUCCESS;
}

int MPLPyPlotter::flush() {
    fflush(proc);
    return EXIT_SUCCESS;
}

int MPLPyPlotter::show() {
    push("show");
    flush();
    return EXIT_SUCCESS;
}

// make an xy plot
int MPLPyPlotter::mpl_xy_plot(const std::string & plot_name, const std::vector<double> & x,
                              const std::vector<double> & y, const std::string & x_label,
                              const std::string & y_label) {
    
    if (process_running == false) {
        return EXIT_FAILURE;
    }
    
    if (x.size() < 1 || x.size() != y.size() || plot_name.length() < 1) {
        return EXIT_FAILURE;
    }
    
    // make a new figure
    std::string cmd = "new_figure " + plot_name;
    push(cmd);
    
    // add in the xy data
    for (size_t idx = 0; idx < x.size(); idx++) {
        double c_x = x[idx];
        double c_y = y[idx];
        push("add " + std::to_string(c_x) + " " + std::to_string(c_y));
    }
    
    // add in labels
    push("set_x_label " + x_label);
    push("set_y_label " + y_label);

    // make the figure
    push("make_figure");
    flush();
    
    return EXIT_SUCCESS;
    
}


