//
//  ElectronPhononSystemGeneral.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "ElectronPhononSystemGeneral.h"

using namespace std;

/****************************************************************************************************/

ElectronPhononSystemGeneral::ElectronPhononSystemGeneral() {

}

ElectronPhononSystemGeneral::~ElectronPhononSystemGeneral() {
    
    /* the STL vectors are automatically destructed, so nothing to do here */

}

/****************************************************************************************************/

void ElectronPhononSystemGeneral::set_omega_L_vec(gsl_vector * omega_L_vec) {
    
    this->omega_L_vec.clear();
    
    for (size_t j=0; j<omega_L_vec->size; j++) {
        this->omega_L_vec.push_back(gsl_vector_get(omega_L_vec, j));
    }
    
}

void ElectronPhononSystemGeneral::set_omega_R_vec(gsl_vector * omega_R_vec) {
    
    this->omega_R_vec.clear();
    
    for (size_t j=0; j<omega_R_vec->size; j++) {
        this->omega_R_vec.push_back(gsl_vector_get(omega_R_vec, j));
    }
    
}

/****************************************************************************************************/

void ElectronPhononSystemGeneral::set_ph_spect(const vector<double> & F_eliash_vec) {
    
    /* N_ph is the number of Omega points */
    this->N_ph = F_eliash_vec.size() + 3;

    /* set Omega_ph_vec */
    this->Omega_ph_vec.clear();
    for (size_t j=0; j<N_ph; j++) {
        this->Omega_ph_vec.push_back(static_cast<double>(j)/(N_ph-1));
    }

    /* pad with 2 leading zeroes */
    this->F_eliash_vec.clear();
    this->F_eliash_vec.push_back(0);
    this->F_eliash_vec.push_back(0);
    for (vector<double>::const_iterator j = F_eliash_vec.begin(); j != F_eliash_vec.end(); ++j) {
        this->F_eliash_vec.push_back(*j);
    }
    /* pad with 1 trailing zero */
    this->F_eliash_vec.push_back(0);
    
    /* calculate the total weight */
    alpha2 = 0;
    for (size_t j=0; j<(N_ph-1); j++) {
        alpha2 += 0.5 * (Omega_ph_vec[j+1] - Omega_ph_vec[j])
        * (this->F_eliash_vec[j] + this->F_eliash_vec[j+1]);
    }
}

void ElectronPhononSystemGeneral::oversample_ph_spect(const size_t & n_os) {
    
    size_t new_N_ph = n_os*(N_ph-1)+1;
    vector<double> new_Omega_ph_vec;
    vector<double> new_F_eliash_vec;
    
    new_Omega_ph_vec.resize(new_N_ph);
    new_F_eliash_vec.resize(new_N_ph);
    size_t j=0;
    generate_n(new_Omega_ph_vec.begin(), new_N_ph,
               [&]() -> double {return static_cast<double>(j++)/(new_N_ph-1);});
    j=0;
    generate_n(new_F_eliash_vec.begin(), new_N_ph,
               [&]() -> double {return eliash_func(new_Omega_ph_vec[j++]);});
    
    //cout << F_eliash_vec.back() << endl;
    //cout << eliash_func(0.9999) << endl;
    //helper::cout_vector(new_Omega_ph_vec);
    //helper::cout_vector(new_F_eliash_vec);
    
    new_F_eliash_vec.erase(new_F_eliash_vec.begin());
    new_F_eliash_vec.erase(new_F_eliash_vec.begin());
    new_F_eliash_vec.pop_back();
    set_ph_spect(new_F_eliash_vec);
    
    //cout << alpha2 << endl;
    
}

/* NOTE: calcualate the electron distribution function on omega_R_vec */
void ElectronPhononSystemGeneral::update_el_dist() {
    
    this->f_vec.clear();
    for (size_t j=0; j<omega_R_vec.size(); j++) {
        f_vec.push_back(el_dist(omega_R_vec[j]));
    }
    
}

void ElectronPhononSystemGeneral::update_ph_dist() {
    
    this->n_vec.clear();
    
    /* the first two elements are irrelevant; see the notes in the header file */
    n_vec.push_back(0.0);
    for (size_t j=1; j<N_ph; j++) {
        n_vec.push_back(ph_dist(Omega_ph_vec[j]));
    }
    
}

void ElectronPhononSystemGeneral::update_interp_coeffs() {

    update_F_interp_coeffs();
    update_nF_interp_coeffs();
    
}

void ElectronPhononSystemGeneral::update_F_interp_coeffs() {
    
    /* the coefficients for F(\Omega) */
    this->F_a_vec.clear();
    this->F_b_vec.clear();
    this->F_alpha_vec.clear();
    
    for (size_t j=0; j<(N_ph-1); j++) {
        double dOmega = Omega_ph_vec[j+1] - Omega_ph_vec[j];
        double a = (F_eliash_vec[j+1] - F_eliash_vec[j]) / dOmega;
        double b = (F_eliash_vec[j] * Omega_ph_vec[j+1] - F_eliash_vec[j+1] * Omega_ph_vec[j]) / dOmega;
        F_a_vec.push_back(a);
        F_b_vec.push_back(b);
    }
    
    this->F_alpha_vec.push_back(F_a_vec[0]);
    for (size_t j=1; j<(N_ph-1); j++) {
        this->F_alpha_vec.push_back(F_a_vec[j]-F_a_vec[j-1]);
    }
    this->F_alpha_vec.push_back(-F_a_vec[N_ph-2]);
    
}

void ElectronPhononSystemGeneral::update_nF_interp_coeffs() {

    /* the coefficients for n(\Omega) F(\Omega) */
    this->nF_a_vec.clear();
    this->nF_b_vec.clear();
    this->nF_alpha_vec.clear();
    
    for (size_t j=0; j<(N_ph-1); j++) {
        double dOmega = Omega_ph_vec[j+1] - Omega_ph_vec[j];
        double a = (n_vec[j+1] * F_eliash_vec[j+1] - n_vec[j] * F_eliash_vec[j]) / dOmega;
        double b = (n_vec[j] * F_eliash_vec[j] * Omega_ph_vec[j+1]
                    - n_vec[j+1] * F_eliash_vec[j+1] * Omega_ph_vec[j]) / dOmega;
        nF_a_vec.push_back(a);
        nF_b_vec.push_back(b);
    }
    
    this->nF_alpha_vec.push_back(nF_a_vec[0]);
    for (size_t j=1; j<(N_ph-1); j++) {
        this->nF_alpha_vec.push_back(nF_a_vec[j]-nF_a_vec[j-1]);
    }
    this->nF_alpha_vec.push_back(-nF_a_vec[N_ph-2]);

}

/****************************************************************************************************/

double ElectronPhononSystemGeneral::el_dist(const double & omega) {
    
    return (el_temp==0) ? 0.0 : 1.0 / (exp(omega/el_temp) + 1.0);
    
}

double ElectronPhononSystemGeneral::ph_dist(const double & Omega) {
    
    return (ph_temp==0) ? 0.0 : 1.0 / (exp(Omega/ph_temp) - 1.0);
    
}

double ElectronPhononSystemGeneral::eliash_func(const double & Omega) {
    
    if (Omega<0 || Omega>1) {
        return 0.0;
    } else {
        size_t bin = min(static_cast<size_t>(floor((N_ph-1)*Omega)),N_ph-2);
        return F_a_vec[bin]*Omega + F_b_vec[bin];
    }
    
}

/****************************************************************************************************/

double ElectronPhononSystemGeneral::xlog(const double & x) {

    return (abs(x)>1e-14) ? x*log(abs(x)) : 0.0;
    
}

/****************************************************************************************************/

double ElectronPhononSystemGeneral::RE_I_F_full(const double & z) {

    double out = 0.0;
    for (size_t j=0; j<N_ph; j++) {
        out += F_alpha_vec[j] * xlog(z-Omega_ph_vec[j]);
    }
    return out;
    
}

double ElectronPhononSystemGeneral::IM_I_F_full(const double & z) {

    if (z<0 || z>1) {
        return 0.0;
    } else {
        size_t bin = min(static_cast<size_t>(floor((N_ph-1)*z)),N_ph-2);
        return -M_PI*(F_a_vec[bin]*z + F_b_vec[bin]);
    }

}

double ElectronPhononSystemGeneral::RE_J_F_full(const double & z) {
    return -RE_I_F_full(-z);
}

double ElectronPhononSystemGeneral::IM_J_F_full(const double & z) {
    return +IM_I_F_full(-z);
}

/****************************************************************************************************/

double ElectronPhononSystemGeneral::RE_I_nF_full(const double & z) {
    
    double out = 0.0;
    for (size_t j=0; j<N_ph; j++) {
        out += nF_alpha_vec[j]*xlog(abs(z-Omega_ph_vec[j]));
    }
    return out;
    
}

double ElectronPhononSystemGeneral::IM_I_nF_full(const double & z) {
    
    if (z<0 || z>1) {
        return 0.0;
    } else {
        size_t bin = min(static_cast<size_t>(floor((N_ph-1)*z)), N_ph-2);
        return -M_PI*(nF_a_vec[bin]*z + nF_b_vec[bin]);
    }
    
}

double ElectronPhononSystemGeneral::RE_J_nF_full(const double & z) {
    return -RE_I_nF_full(-z);
}

double ElectronPhononSystemGeneral::IM_J_nF_full(const double & z) {
    return +IM_I_nF_full(-z);
}

/****************************************************************************************************/
#warning check these!
double ElectronPhononSystemGeneral::RE_K_Z(const size_t & L_idx, const size_t & R_idx) {

    double omega = omega_L_vec[L_idx];
    double omegap = omega_R_vec[R_idx];
    double f = f_vec[R_idx];
    double fb = 1 - f;

    return
    + fb * RE_I_F_full(omega - omegap) + RE_I_nF_full(omega - omegap)
    + fb * RE_J_F_full(omega + omegap) + RE_J_nF_full(omega + omegap)
    + f * RE_J_F_full(omega - omegap) + RE_J_nF_full(omega - omegap)
    + f * RE_I_F_full(omega + omegap) + RE_I_nF_full(omega + omegap);

}

double ElectronPhononSystemGeneral::IM_K_Z(const size_t & L_idx, const size_t & R_idx) {

    double omega = omega_L_vec[L_idx];
    double omegap = omega_R_vec[R_idx];
    double f = f_vec[R_idx];
    double fb = 1 - f;
    
    return
    + fb * IM_I_F_full(omega - omegap) + IM_I_nF_full(omega - omegap)
    + fb * IM_J_F_full(omega + omegap) + IM_J_nF_full(omega + omegap)
    + f * IM_J_F_full(omega - omegap) + IM_J_nF_full(omega - omegap)
    + f * IM_I_F_full(omega + omegap) + IM_I_nF_full(omega + omegap);

}

double ElectronPhononSystemGeneral::RE_K_Delta(const size_t & L_idx, const size_t & R_idx) {

    double omega = omega_L_vec[L_idx];
    double omegap = omega_R_vec[R_idx];
    double f = f_vec[R_idx];
    double fb = 1 - f;
    
    return
    + fb * RE_J_F_full(omega + omegap) + RE_J_nF_full(omega + omegap)
    - fb * RE_I_F_full(omega - omegap) - RE_I_nF_full(omega - omegap)
    + f * RE_I_F_full(omega + omegap) + RE_I_nF_full(omega + omegap)
    - f * RE_J_F_full(omega - omegap) - RE_J_nF_full(omega - omegap);

}

double ElectronPhononSystemGeneral::IM_K_Delta(const size_t & L_idx, const size_t & R_idx) {

    double omega = omega_L_vec[L_idx];
    double omegap = omega_R_vec[R_idx];
    double f = f_vec[R_idx];
    double fb = 1 - f;
    
    return
    + fb * IM_J_F_full(omega + omegap) + IM_J_nF_full(omega + omegap)
    - fb * IM_I_F_full(omega - omegap) - IM_I_nF_full(omega - omegap)
    + f * IM_I_F_full(omega + omegap) + IM_I_nF_full(omega + omegap)
    - f * IM_J_F_full(omega - omegap) - IM_J_nF_full(omega - omegap);

}


/****************************************************************************************************/
/****** deprecated **********************************************************************************/
/****************************************************************************************************/

//
//double ElectronPhononSystemGeneral::RE_I(const double & a, const double & b,
//                                         const double & Omega_0, const double & Omega_1,
//                                         const double & z) {
//    double val = (a*z+b)*log(abs((z-Omega_0)/(z-Omega_1))) + a*(Omega_0-Omega_1);
//    return (isnan(val) || abs(val)>100) ? 0.0 : val;
//
//}
//
//
//double ElectronPhononSystemGeneral::IM_I(const double & a, const double & b,
//                                         const double & Omega_0, const double & Omega_1,
//                                         const double & z) {
//
//    return (z>Omega_0 && z<Omega_1) ? (-M_PI*(a*z+b)) : 0.0;
//
//}
//
//double ElectronPhononSystemGeneral::RE_J(const double & a, const double & b,
//                                         const double & Omega_0, const double & Omega_1,
//                                         const double & z) {
//
//    return -RE_I(a, b, Omega_0, Omega_1, -z);
//
//}
//
//double ElectronPhononSystemGeneral::IM_J(const double & a, const double & b,
//                                         const double & Omega_0, const double & Omega_1,
//                                         const double & z) {
//
//    return IM_I(a, b, Omega_0, Omega_1, -z);
//
//}
//
//double ElectronPhononSystemGeneral::RE_K_Z(const size_t & L_idx, const size_t & R_idx) {
//
//    double omega = omega_L_vec[L_idx];
//    double omegap = omega_R_vec[R_idx];
//    double f = f_vec[R_idx];
//    double fb = 1 - f;
//    
//    double val = 0; /* output */
//    
//    /* NOTE: we assume n_vec[0] = 0.0; we treat the first interval properly afterwards */
//    for (size_t j=0; j<N_ph; j++) {
//        val += (n_vec[j] + fb)
//        * (+ RE_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap)
//           + RE_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap))
//        + (n_vec[j] + f)
//        * (+ RE_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap)
//           + RE_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap));
//    }
//    
//    /* contribution from the n-terms in the first interval [0, \Omega_1] */
//    val += (F_eliash_vec[1]*ph_temp/Omega_ph_vec[1])
//    * (+ RE_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + RE_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       + RE_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + RE_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap));
//    
//    return val;
//}
//
//double ElectronPhononSystemGeneral::IM_K_Z(const size_t & L_idx, const size_t & R_idx) {
//    
//    double omega = omega_L_vec[L_idx];
//    double omegap = omega_R_vec[R_idx];
//    double f = f_vec[R_idx];
//    double fb = 1 - f;
//    
//    double val = 0; /* output */
//    
//    /* NOTE: we assume n_vec[0] = 0.0; we treat the first interval properly afterwards */
//    for (size_t j=0; j<N_ph; j++) {
//        val += (n_vec[j] + fb)
//        * (+ IM_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap)
//           + IM_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap))
//        + (n_vec[j] + f)
//        * (+ IM_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap)
//           + IM_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap));
//    }
//    
//    /* contribution from the n-terms in the first interval [0, \Omega_1] */
//    val += (F_eliash_vec[1]*ph_temp/Omega_ph_vec[1])
//    * (+ IM_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + IM_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       + IM_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + IM_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap));
//    
//    return val;
//}
//
//double ElectronPhononSystemGeneral::RE_K_Delta(const size_t & L_idx, const size_t & R_idx) {
//    
//    double omega = omega_L_vec[L_idx];
//    double omegap = omega_R_vec[R_idx];
//    double f = f_vec[R_idx];
//    double fb = 1 - f;
//    
//    double val = 0; /* output */
//    
//    /* NOTE: we assume n_vec[0] = 0.0; we treat the first interval properly afterwards */
//    for (size_t j=0; j<N_ph; j++) {
//        val += (n_vec[j] + fb)
//        * (+ RE_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap)
//           - RE_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap))
//        + (n_vec[j] + f)
//        * (+ RE_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap)
//           - RE_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap));
//    }
//    
//    /* contribution from the n-terms in the first interval [0, \Omega_1] */
//    val += (F_eliash_vec[1]*ph_temp/Omega_ph_vec[1])
//    * (+ RE_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       - RE_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + RE_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       - RE_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap));
//    
//    return val;
//}
//
//double ElectronPhononSystemGeneral::IM_K_Delta(const size_t & L_idx, const size_t & R_idx) {
//    
//    double omega = omega_L_vec[L_idx];
//    double omegap = omega_R_vec[R_idx];
//    double f = f_vec[R_idx];
//    double fb = 1 - f;
//    
//    double val = 0; /* output */
//    
//    /* NOTE: we assume n_vec[0] = 0.0; we treat the first interval properly afterwards */
//    for (size_t j=0; j<N_ph; j++) {
//        val += (n_vec[j] + fb)
//        * (+ IM_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap)
//           - IM_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap))
//        + (n_vec[j] + f)
//        * (+ IM_I(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega + omegap)
//           - IM_J(a_vec[j], b_vec[j], Omega_ph_vec[j], Omega_ph_vec[j+1], omega - omegap));
//    }
//    
//    /* contribution from the n-terms in the first interval [0, \Omega_1] */
//    val += (F_eliash_vec[1]*ph_temp/Omega_ph_vec[1])
//    * (+ IM_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       - IM_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap)
//       + IM_I(0.0, 1.0, 0.0, Omega_ph_vec[1], omega + omegap)
//       - IM_J(0.0, 1.0, 0.0, Omega_ph_vec[1], omega - omegap));
//    
//    return val;
//}
//


