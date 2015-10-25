//
//  ElectronPhononSystemLorenzian.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 9/27/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "ElectronPhononSystemLorenzian.h"

using namespace std;

/****************************************************************************************************/

ElectronPhononSystemLorenzian::ElectronPhononSystemLorenzian() : ElectronPhononSystem() { }


ElectronPhononSystemLorenzian::~ElectronPhononSystemLorenzian() { }


void ElectronPhononSystemLorenzian::set_system_params(const double & el_temp, const double & n0,
                                              const double & Omega_0, const double & Omega_c,
                                              const double & gamma, const double & mu_c) {

    this->el_temp = el_temp;
    this->n0 = n0;
    this->Omega_0 = Omega_0;
    this->Omega_c = Omega_c;
    this->gamma = gamma;
    this->mu_c = mu_c;
    
    /* aux */
    this->gamma2 = gamma*gamma;
    this->Omega_c2 = Omega_c*Omega_c;
    this->c0 = (gamma2+Omega_c2)*std::atan(Omega_c/gamma);
    
    /* calculate the truncated Lorenzian normalization constant */
    A_norm = (0.5*gamma)/(std::atan(Omega_c/gamma) - gamma*Omega_c/(gamma2 + Omega_c2));
    
    /* set the maximum phonon spectrum frequecny */
    Omega_max = Omega_0 + Omega_c;

}

double ElectronPhononSystemLorenzian::el_dist(const double & omega) {

    return 1.0 / (exp(omega/el_temp) + 1.0);
    
}

double ElectronPhononSystemLorenzian::ph_dist(const double & Omega) {
    
    return n0;
}

double ElectronPhononSystemLorenzian::eliash_func(const double & Omega) {
    
    return IM_truncated_Lorenzian(Omega - Omega_0)/M_PI;
    
}

double ElectronPhononSystemLorenzian::RE_truncated_Lorenzian(const double & x) {
    
    double num = 2*c0*x + gamma*(x*x-Omega_c2)*std::log(std::abs((x-Omega_c)/(x+Omega_c)));
    double den = 2*(x*x+gamma2)*(c0 - gamma*Omega_c);
    
    return num/den;
    
}

double ElectronPhononSystemLorenzian::IM_truncated_Lorenzian(const double & x) {
    
    return (std::abs(x) < Omega_c) ?
    M_PI*A_norm*(1.0/(x*x+gamma2) - 1.0/(Omega_c2+gamma2)) : 0;
    
}

/****************************************************************************************************/
/******************************************* tests **************************************************/
/****************************************************************************************************/


void ElectronPhononSystemLorenzian::test_calculate_total_eliash_weight() {
    
    /* generate GK15 weights */
    int order = 15;
    double xtab[order], weight[order];
    kronrod_set(order, xtab, weight);

    /* generate a composite rule */
    vector<double> xtab_vec, weight_vec;
    generate_composite_rule(xtab, weight, order, -1.0, 1.0, 20, 0, this->Omega_max, xtab_vec,weight_vec);
    
    double sum = 0;
    for (size_t i=0; i<xtab_vec.size(); i++) {
        sum += weight_vec[i]*eliash_func(xtab_vec[i]);
    }
    
    cout << "Total weight of $F(\\Omega)$ is: " << sum << endl;
    
}

