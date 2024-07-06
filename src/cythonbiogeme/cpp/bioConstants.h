#ifndef BIOCONSTANTS_H
#define BIOCONSTANTS_H

#include <limits>
#include <cmath>
#include <iomanip>
#include "bioTypes.h"
#include "bioExceptions.h"

namespace constants {

  // Function to calculate and cache the upper bound on floating values
  inline bioReal get_upper_bound() {
    static const bioReal sqrt_max_float = std::sqrt(std::numeric_limits<bioReal>::max());
    return sqrt_max_float;
  }

  //inline bioReal get_machine_epsilon() {
  //  return std::sqrt(std::numeric_limits<bioReal>::epsilon());
  //}

  inline bioReal get_almost_zero() {
    static const bioReal machine_epsilon = std::numeric_limits<bioReal>::epsilon();
    return machine_epsilon;
  }

}

inline void validate(bioReal alpha) {
bioReal upper_bound = constants::get_upper_bound();
if (alpha > upper_bound || alpha < -upper_bound) {
    std::stringstream str ;
    str << "The input value must be between  ["
    << std::setprecision(3) << std::scientific << -upper_bound
    << " and "
    << std::setprecision(3) << std::scientific << upper_bound
    << "]. Invalid value: "
    << alpha << std::endl ;
      throw bioExceptions(__FILE__,__LINE__,str.str()) ;
  }
}

inline bioReal project(bioReal alpha) {
    bioReal upper_bound = constants::get_upper_bound();
    if (alpha >= upper_bound) {
        return upper_bound;
    } else if (alpha <= -upper_bound) {
        return -upper_bound;
    } else {
        return alpha;
    }
}
#endif // BIOCONSTANTS_H