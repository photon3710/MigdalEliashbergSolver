//
//  HelperModule.cpp
//  MigdalEliashbergSolver
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#include "HelperModule.h"

using namespace std;

void helper::write_gsl_vector(std::ofstream & of, gsl_vector * vec) {
    
    for (size_t j=0; j<vec->size; j++) {
        of << gsl_vector_get(vec, j) << endl;
    }
    
}

void helper::write_gsl_vector_complex(std::ofstream & of, gsl_vector_complex * vec) {
    
    for (size_t j=0; j<vec->size; j++) {
        of << GSL_VECTOR_REAL(vec, j) << ", " << GSL_VECTOR_IMAG(vec, j) << endl;
    }
    
}

void helper::write_STL_vector(std::ofstream & of, const vector<double> & vec) {
    
    for (size_t j=0; j<vec.size(); j++) {
        of << vec[j] << endl;
    }
    
}

vector<double> helper::gsl_vector_to_STL_vector(gsl_vector * gsl_vec) {
    vector<double> vec;
    vec.reserve(gsl_vec->size);
    for (size_t j=0; j<gsl_vec->size; j++) {
        vec.push_back(gsl_vector_get(gsl_vec, j));
    }
    return vec;
}

void helper::cout_vector(std::vector<double> vec) {
    for (vector<double>::iterator it = vec.begin(); it != vec.end(); it++) {
        cout << *it;
        if (it!=(--vec.end())) {
            cout << ", ";
        }
    }
    cout << endl;
}

vector<double> helper::linspace(const double & a, const double & b, const size_t & N) {
    vector<double> vec;
    for (size_t j=0; j<N; j++) {
        vec.push_back(a + (b-a)*j/(N-1));
    }
    return vec;
}

