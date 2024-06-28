#ifndef BIOCONSTANTS_H
#define BIOCONSTANTS_H

#include <limits>
#include <cmath>
#include <iomanip>
#include "bioTypes.h"

namespace constants {

  // Function to calculate and cache the sqrt_max_float at runtime
  inline bioReal get_sqrt_max_float() {
    static const bioReal sqrt_max_float = std::sqrt(std::numeric_limits<bioReal>::max());
    return sqrt_max_float;
  }

  inline bioReal get_sqrt_machine_epsilon() {
    static const bioReal sqrt_machine_epsilon = std::sqrt(std::numeric_limits<bioReal>::epsilon());
    return sqrt_machine_epsilon;
  }

  inline bioReal log_sqrt_max_float() {
    return std::log(get_sqrt_max_float());
  }

}

inline void validate(bioReal alpha) {
bioReal sqrt_max_float = constants::get_sqrt_max_float();
if (alpha > sqrt_max_float || alpha < -sqrt_max_float) {
    std::stringstream str ;
    str << "The input value must be between  ["
    << std::setprecision(3) << std::scientific << -sqrt_max_float
    << " and "
    << std::setprecision(3) << std::scientific << sqrt_max_float
    << "]. Invalid value: "
    << alpha << std::endl ;
      throw bioExceptions(__FILE__,__LINE__,str.str()) ;
  }
}

inline bioReal project(bioReal alpha) {
    bioReal sqrt_max_float = constants::get_sqrt_max_float();
    if (alpha >= sqrt_max_float) {
        return sqrt_max_float;
    } else if (alpha <= -sqrt_max_float) {
        return -sqrt_max_float;
    } else {
        return alpha;
    }
}
#endif // BIOCONSTANTS_H