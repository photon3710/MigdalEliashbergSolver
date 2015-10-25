//
//  main.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/26/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//
//
//  TODO:
//
//  [DONE] Implemeny oversampling: for any given test Eliashberg function, oversampling will increase the
//   accuracy of phonon spectral summations. This is due to the divergence of $n(\Omega)$ diverges
//   for small $\Omega$ and the necessity of using a fine grid for phonon spectral summations.
//
//   For an unbiased comparison between different test Eliashberg functions with different number of
//   points, we must oversample the coarser test functions to match the finer ones.
//
//

#include "main.h"

using namespace std;
using namespace libconfig;

int job_runner(const string & config_file);

void help() {
    
    cout << "Real-axis Migdal-Eliashberg solver tool v0   (c) 2015 Mehrtash Babadi" << endl
         << "---------------------------------------------------------------------" << endl << endl
         << "options: " << endl << endl
         << " -h or --help:   show this message" << endl
         << " <config_file>:  provide a job configuration file" << endl;
    
}

int main(int argc, const char * argv[]) {

//#warning remove this before the final build
//    /* <DEBUG CODE> */
//    argc = 2;
//    argv[1] = "/Users/mehrtashbabadi/Codes/cpp/MigdalEliashbergSolver/Output/GA alpha2=0.5/"
//    "MESolverBase.cfg";
//    /* </DEBUG CODE> */
    
    string config_file;
    
    if (argc > 1) {
        
        string arg1 = argv[1];
        
        if(arg1.compare("--help") == 0 || arg1.compare("-h") == 0) {
            help();
            return EXIT_SUCCESS;
        }
        
        config_file = arg1;
        
    } else {
        
        help();
        return EXIT_FAILURE;
    
    }
    
    return job_runner(config_file);
    
}

int job_runner(const string & config_file) {
    
    Config cfg;
    cfg.setAutoConvert(true);

    try
    {
        cfg.readFile(config_file.c_str());
    }
    catch(const FileIOException &fioex)
    {
        std::cerr << "I/O error while reading file." << std::endl;
        return EXIT_FAILURE;
    }
    catch(const ParseException &pex)
    {
        std::cerr << "Parse error at " << pex.getFile() << ":" << pex.getLine()
        << " - " << pex.getError() << std::endl;
        return EXIT_FAILURE;
    }

    const Setting & root = cfg.getRoot();
    
    string out_path_base, out_path_ext;
    double mu_c, omega_c, ph_temp;
    unsigned order_base, order_inf, N_rep_Z_th, N_rep_Z_mid, N_rep_Z_inf, N_rep_Delta_th, N_rep_Delta_mid;
    double eps_f;
    double Tc_lo, Tc_hi, tol_Tc_rel;
    unsigned max_iter;
    vector<double> F_eliash_vec;
    unsigned n_os;

    try
    {
        /* output */
        out_path_base = (const char*)(root["output"].lookup("out_path_base"));
        out_path_ext = (const char*)(root["output"].lookup("out_path_ext"));
        
        /* system */
        mu_c = (double)(root["system"].lookup("mu_c"));
        omega_c = (double)(root["system"].lookup("omega_c"));
        ph_temp = (double)(root["system"].lookup("ph_temp"));
        n_os = (unsigned)(root["system"].lookup("n_os"));

        const Setting & eliash_array = root["system"]["F_eliash_vec"];
        for (libconfig::Setting::const_iterator it = eliash_array.begin(); it != eliash_array.end(); ++it) {
            F_eliash_vec.push_back((double)(*it));
        }
        
        /* quadrature */
        order_base = (unsigned)(root["quadrature"].lookup("order_base"));
        order_inf = (unsigned)(root["quadrature"].lookup("order_inf"));
        N_rep_Z_th = (unsigned)(root["quadrature"].lookup("N_rep_Z_th"));
        N_rep_Z_mid = (unsigned)(root["quadrature"].lookup("N_rep_Z_mid"));
        N_rep_Z_inf = (unsigned)(root["quadrature"].lookup("N_rep_Z_inf"));
        N_rep_Delta_th = (unsigned)(root["quadrature"].lookup("N_rep_Delta_th"));
        N_rep_Delta_mid = (unsigned)(root["quadrature"].lookup("N_rep_Delta_mid"));
        eps_f = (double)(root["quadrature"].lookup("eps_f"));
        
        /* solver */
        Tc_lo = (double)(root["solver"].lookup("Tc_lo"));
        Tc_hi = (double)(root["solver"].lookup("Tc_hi"));
        tol_Tc_rel = (double)(root["solver"].lookup("tol_Tc_rel"));
        max_iter = (unsigned)(root["solver"].lookup("max_iter"));
        
    }
    catch (const SettingNotFoundException & nfex)
    {
        cout << "Value not provided for " << nfex.getPath() << endl;
        return EXIT_FAILURE;
    }
    
    MigdalEliashbergSolverGeneral solver;
    
    /* set the system parameters */
    solver.set_output_path(out_path_base + "/" + out_path_ext);
    solver.set_sys_params(F_eliash_vec, omega_c, mu_c, ph_temp);
    solver.sys->update_F_interp_coeffs();
    solver.sys->oversample_ph_spect(n_os);
    solver.report() << "oversampled to " << solver.sys->N_ph << " spectral points." << endl;
    
    solver.set_quad_params(order_base, order_inf, eps_f, N_rep_Z_th, N_rep_Z_mid, N_rep_Z_inf, N_rep_Delta_th, N_rep_Delta_mid);
    
    gsl_set_error_handler_off();
    
    solver.find_Tc(Tc_lo, Tc_hi, max_iter, tol_Tc_rel);

    solver.save_data();
    
    return EXIT_SUCCESS;

}


