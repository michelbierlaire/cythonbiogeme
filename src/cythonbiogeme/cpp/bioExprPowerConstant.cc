//-*-c++-*------------------------------------------------------------
//
// File name : bioExprPowerConstant.cc
// Date:   Wed Jul 3 10:20:56 2024
// Author: Michel Bierlaire
//
//--------------------------------------------------------------------

#include "bioExprPowerConstant.h"
#include "bioExceptions.h"
#include <cmath>
#include <sstream>
#include "bioDebug.h"
#include "bioConstants.h"

bioExprPowerConstant::bioExprPowerConstant(bioExpression* l, bioReal an_exponent) :
  child(l), exponent(an_exponent) {
  listOfChildren.push_back(l) ;
}

bioExprPowerConstant::~bioExprPowerConstant() {

}

const bioDerivatives* bioExprPowerConstant::getValueAndDerivatives(std::vector<bioUInt> literalIds,
						     bioBoolean gradient,
						     bioBoolean hessian) {

  static const bioReal upper_bound = constants::get_upper_bound();
  static const bioReal almost_zero = constants::get_almost_zero();

  theDerivatives.with_g = gradient ;
  theDerivatives.with_h = hessian ;

  bioUInt n = literalIds.size() ;
  theDerivatives.resize(n) ;

  const bioDerivatives* childResult = child->getValueAndDerivatives(literalIds,gradient,hessian) ;

  bool is_exponent_integer = (std::trunc(exponent) == exponent);
  // Check for domain errors that might be raised by std::pow
  if (childResult->f < 0 && !is_exponent_integer) {
    throw bioExceptions(__FILE__,__LINE__,"Negative base with a non-integer exponent is invalid in the real number domain.");
  }

  if (exponent == 0) {
    theDerivatives.setDerivativesToZero() ;
    theDerivatives.f = 1.0 ;
    return &theDerivatives ;
  }

  if (childResult->f >= almost_zero || exponent > 2 || is_exponent_integer) {
    theDerivatives.f = std::pow(childResult->f, exponent) ;
    if (gradient) {
      bioReal to_k_minus_1 = exponent * std::pow(childResult->f, exponent-1.0);
	  for (bioUInt i = 0 ; i < n ; ++i) {
	    if (childResult->g[i] != 0.0) {
	      theDerivatives.g[i] = to_k_minus_1 * childResult->g[i] ;
	    }
	    else {
	      theDerivatives.g[i] = 0.0 ;
	    }
	    if (hessian) {
          for (bioUInt j = i ; j < n ; ++j) {
            theDerivatives.h[i][j] = 0.0 ;
            if (childResult->g[i] != 0.0 && childResult->g[i] != 0.0) {
              theDerivatives.h[i][j] +=
                exponent *
                (exponent - 1.0) *
                std::pow(childResult->f, exponent-2.0) *
                childResult->g[i] *
                childResult->g[j] ;
            }
            if (childResult->h[i][j] != 0.0) {
              theDerivatives.h[i][j] += to_k_minus_1 * childResult->h[i][j] ;
            }
            if (i != j) {
              theDerivatives.h[j][i] = theDerivatives.h[i][j] ;
            }
          }
        }
      }
	}
    theDerivatives.dealWithNumericalIssues() ;
    return &theDerivatives ;
  }
  if (childResult->f >= 0 && exponent == 2.0) {
    theDerivatives.f = childResult->f * childResult->f ;
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
	    if (childResult->g[i] != 0.0) {
	      theDerivatives.g[i] = 2.0 * childResult->f  * childResult->g[i] ;
	    }
	    else {
	      theDerivatives.g[i] = 0.0 ;
	    }
	    if (hessian) {
          for (bioUInt j = i ; j < n ; ++j) {
            theDerivatives.h[i][j] = 0.0 ;
            if (childResult->g[i] != 0.0 && childResult->g[i] != 0.0) {
              theDerivatives.h[i][j] += 2.0 *
                childResult->g[i] *
                childResult->g[j] ;
            }
            if (childResult->h[i][j] != 0.0) {
              theDerivatives.h[i][j] += 2.0 * childResult->h[i][j] ;
            }
            if (i != j) {
              theDerivatives.h[j][i] = theDerivatives.h[i][j] ;
            }
          }
        }
      }
	}
    theDerivatives.dealWithNumericalIssues() ;
    return &theDerivatives ;
  }
  if (childResult->f >= 0 && exponent > 0 && exponent <= 2.0) {
    bioReal the_power = std::pow(almost_zero, exponent-1) ;
    theDerivatives.f = the_power * childResult->f ;
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
	    if (childResult->g[i] != 0.0) {
	      theDerivatives.g[i] = the_power * childResult->g[i] ;
	    }
	    else {
	      theDerivatives.g[i] = 0.0 ;
	    }
	    if (hessian) {
          for (bioUInt j = i ; j < n ; ++j) {
            if (childResult->h[i][j] == 0) {
              theDerivatives.h[i][j] = 0.0 ;
            }
            else {
               theDerivatives.h[i][j] = the_power * childResult->h[i][j] ;
            }
            if (i != j) {
              theDerivatives.h[j][i] = theDerivatives.h[i][j] ;
            }
          }
        }
      }
	}
    theDerivatives.dealWithNumericalIssues() ;
    return &theDerivatives ;
  }
  if (childResult->f >= 0 && exponent < 0) {
    theDerivatives.f =
      std::pow(almost_zero, exponent-1) * childResult->f +
      upper_bound * (1.0 - childResult->f / almost_zero);
    if (gradient) {
      bioReal the_power = std::pow(almost_zero, exponent) ;
      for (bioUInt i = 0 ; i < n ; ++i) {
	    if (childResult->g[i] != 0.0) {
	      theDerivatives.g[i] = (the_power - upper_bound) * childResult->g[i] / almost_zero;
	    }
	    else {
	      theDerivatives.g[i] = 0.0 ;
	    }
	    if (hessian) {
          for (bioUInt j = i ; j < n ; ++j) {
            if (childResult->h[i][j] == 0) {
              theDerivatives.h[i][j] = 0.0 ;
            }
            else {
              theDerivatives.h[i][j] = (the_power - upper_bound) * childResult->h[i][j] / almost_zero;
            }
            if (i != j) {
              theDerivatives.h[j][i] = theDerivatives.h[i][j] ;
            }
          }
        }
      }
	}
    theDerivatives.dealWithNumericalIssues() ;
    return &theDerivatives ;
  }
  std::stringstream str ;
  str << "This should not be reached, as all conditions have been enumerated. exponent = " << exponent << " base = " << childResult->f ;
  throw bioExceptions(__FILE__,__LINE__,str.str());
}

bioString bioExprPowerConstant::print(bioBoolean hp) const {
  std::stringstream str ;
  if (hp) {
    str << "^(" << child->print(hp) << "*" << exponent << ")" ;

  }
  else {
    str << "(" << child->print(hp) << "^" << exponent << ")" ;
  } 
  return str.str() ;
}
