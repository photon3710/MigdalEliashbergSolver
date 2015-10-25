//
//  HelperModule.h
//
//  Created by Mehrtash Babadi on 10/6/15.
//  Copyright (c) 2015 mbabadi. All rights reserved.
//

#ifndef __HelperModule__
#define __HelperModule__

#include <stdio.h>
#include <iostream>
#include <vector>
#include <fstream>

#include <gsl/gsl_vector.h>

namespace helper {
    
    void write_gsl_vector(std::ofstream & of, gsl_vector * vec);
    
    void write_gsl_vector_complex(std::ofstream & of, gsl_vector_complex * vec);

    void write_STL_vector(std::ofstream & of, const std::vector<double> & vec);

    std::vector<double> gsl_vector_to_STL_vector(gsl_vector * gsl_vec);
    
    std::vector<double> linspace(const double & a, const double & b, const size_t & N);

    void cout_vector(std::vector<double> vec);
}

#endif /* defined(__HelperModule__) */
