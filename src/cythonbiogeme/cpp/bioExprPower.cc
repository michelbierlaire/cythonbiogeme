//-*-c++-*------------------------------------------------------------
//
// File name : bioExprPower.cc
// @date   Fri Apr 13 12:20:46 2018
// @author Michel Bierlaire
// @version Revision 1.0
//
//--------------------------------------------------------------------

#include "bioExprPower.h"
#include <cmath>
#include <sstream>
#include "bioDebug.h"
#include "bioConstants.h"
#include "bioExceptions.h"

bioExprPower::bioExprPower(bioExpression* l, bioExpression* r) :
  left(l), right(r) {
  listOfChildren.push_back(l) ;
  listOfChildren.push_back(r) ;

}

bioExprPower::~bioExprPower() {

}

const bioDerivatives* bioExprPower::getValueAndDerivatives(std::vector<bioUInt> literalIds,
						     bioBoolean gradient,
						     bioBoolean hessian) {

  static const bioReal upper_bound = constants::get_upper_bound();
  static const bioReal almost_zero = constants::get_almost_zero();
  static const bioReal log_almost_zero = std::log(almost_zero) ;

  // If the exponent is a constant, we call the corresponding code.
  if (!right->containsLiterals(literalIds)) {
    return getValueAndDerivativesConstant(literalIds, gradient, hessian) ;
  }


  theDerivatives.with_g = gradient ;
  theDerivatives.with_h = hessian ;

  bioUInt n = literalIds.size() ;
  theDerivatives.resize(n) ;

  const bioDerivatives* leftResult = left->getValueAndDerivatives(literalIds,gradient,hessian) ;

  const bioDerivatives* rightResult = right->getValueAndDerivatives(literalIds,gradient,hessian) ;


  // As the exponent is not constant, the base must be non negative
  if (leftResult->f < 0) {
    std::stringstream str ;
    str << "Current values of the literals: " << std::endl ;
    std::map<bioString,bioReal> m = getAllLiteralValues() ;
    for (std::map<bioString,bioReal>::iterator i = m.begin() ;
	  i != m.end() ;
	  ++i) {
	  str << i->first << " = " << i->second << std::endl ;
    }
    if (rowIndex != NULL) {
	  str << "row number: " << *rowIndex << ", ";
    }
    str << "Cannot raise a non positive number to a non integer [" << leftResult->f << "]" << std::endl ;
      throw bioExceptions(__FILE__,__LINE__,str.str()) ;
  }

  if (leftResult->f >= almost_zero) {
    theDerivatives.f = pow(leftResult->f,rightResult->f) ;
    if (gradient) {
      std::vector<bioReal> G(n,0.0) ;
      for (bioUInt i = 0 ; i < n ; ++i) {
        theDerivatives.g[i] = 0.0 ;
        if (theDerivatives.f != 0.0) {
	      if (leftResult->g[i] != 0.0 && rightResult->f != 0.0) {
	        bioReal term = leftResult->g[i] * rightResult->f / leftResult->f ;
	        G[i] += term ;
	      }
	      if (rightResult->g[i] != 0.0) {
	        G[i] += rightResult->g[i] * log(leftResult->f) ;
	      }
	      theDerivatives.g[i] = theDerivatives.f * G[i] ;
        }
      }
      if (hessian) {
        for (bioUInt i = 0 ; i < n ; ++i) {
	      for (bioUInt j = i ; j < n ; ++j) {
	        bioReal v = G[i] * theDerivatives.g[j] ;
	        if (theDerivatives.f != 0 && rightResult != NULL) {
	          bioReal term(0.0) ;
	          bioReal hright = rightResult->h[i][j] ;
	          if (hright != 0.0) {
	            term += hright * log(leftResult->f) ;
	          }
	          if (leftResult->g[j] != 0.0 && rightResult->g[i] != 0.0) {
	            term += leftResult->g[j] * rightResult->g[i] / leftResult->f ;
	          }
	          if (leftResult->g[i] != 0.0 && rightResult->g[j] != 0.0) {
	            term += leftResult->g[i] * rightResult->g[j] / leftResult->f ;
	          }
	          if (leftResult->g[i] != 0.0 && leftResult->g[j] != 0.0) {
	            bioReal asquare = leftResult->f * leftResult->f ;
	            term -= leftResult->g[i] * leftResult->g[j] * rightResult->f / asquare ;
	          }
	          bioReal hleft = leftResult->h[i][j] ;
	          if (hleft != 0.0) {
	            term += hleft * rightResult->f / leftResult->f ;
	          }
	          if (term != 0.0) {
	            v += term * theDerivatives.f ;
	          }
	        }
	        if (i == j) {
	        theDerivatives.h[i][i] = v ;
	        }
	        else {
	          theDerivatives.h[i][j] = theDerivatives.h[j][i] = v ;
	        }
	      }
        }
      }
    }
    theDerivatives.dealWithNumericalIssues() ;
    return &theDerivatives ;
  }


  // Now we treat the case where the base is close to zero
  // We start py precomputing various quantities.

  bioReal the_power = std::pow(almost_zero, rightResult->f - 1.0) ;
  bioReal f_plus = leftResult->f * the_power ;
  std::vector<bioReal> gradient_f_plus(n) ;
  std::vector<bioReal> gradient_right_times_log_xi(n) ;
  std::vector<std::vector<bioReal> > hessian_f_plus(n,std::vector<bioReal>(n, 0.0)) ;
  if (gradient) {
    for (bioUInt i = 0 ; i < n ; ++i) {
      gradient_right_times_log_xi[i] = log_almost_zero * rightResult->g[i] ;
      gradient_f_plus[i] =
        f_plus * gradient_right_times_log_xi[i] +
        the_power * leftResult->g[i] ;
    }
    if (hessian) {
      for (bioUInt i = 0 ; i < n ; ++i) {
        for (bioUInt j = i ; j < n ; ++j) {
          hessian_f_plus[i][j] =
            gradient_f_plus[j] * leftResult->f * gradient_right_times_log_xi[i]
            + the_power * leftResult->g[j] * gradient_right_times_log_xi[i]
            + f_plus * log_almost_zero * rightResult->h[i][j]
            + gradient_f_plus[j] * leftResult->g[i]
            + the_power * leftResult->h[i][j] ;
          if (i != j) {
            hessian_f_plus[j][i] = hessian_f_plus[i][j] ;
          }
        }
      }
    }
  }
  // Case where the exponent is non negative.

  if (rightResult->f >= 0.0) {
    theDerivatives.f = f_plus ;
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
        theDerivatives.g[i] = gradient_f_plus[i] ;
        if (hessian) {
          for (bioUInt j = 0 ; j < n ; ++j) {
            theDerivatives.h[i][j] = hessian_f_plus[i][j] ;
          }
        }
      }
    }
    return &theDerivatives ;
  }

  // Case where the exponent is negative.

  if (rightResult->f < 0.0) {
    theDerivatives.f = f_plus + upper_bound * (1.0 - leftResult->f / almost_zero);
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
        theDerivatives.g[i] = gradient_f_plus[i] - leftResult->g[i] * upper_bound / almost_zero;
        if (hessian) {
          for (bioUInt j = 0 ; j < n ; ++j) {
            theDerivatives.h[i][j] = hessian_f_plus[i][j] - leftResult->h[i][j] * upper_bound / almost_zero ;
          }
        }
      }
    }
    return &theDerivatives ;
  }
  std::stringstream str ;
  str << "This should not be reached, as all conditions have been enumerated. exponent = " << leftResult->f << " base = " << rightResult->f ;
  throw bioExceptions(__FILE__,__LINE__,str.str());
}

const bioDerivatives* bioExprPower::getValueAndDerivativesConstant(std::vector<bioUInt> literalIds,
						     bioBoolean gradient,
						     bioBoolean hessian) {

  // This code is copied from the bioExprPowerConstant.cc file.
  static const bioReal upper_bound = constants::get_upper_bound();
  static const bioReal almost_zero = constants::get_almost_zero();

  theDerivatives.with_g = gradient ;
  theDerivatives.with_h = hessian ;

  bioUInt n = literalIds.size() ;
  theDerivatives.resize(n) ;

  const bioDerivatives* childResult = left->getValueAndDerivatives(literalIds,gradient,hessian) ;

  bioReal exponent = right->getValue() ;
  bool is_exponent_integer = (std::trunc(exponent) == exponent);
  // Check for domain errors that might be raised by std::pow
  if (childResult->f < 0 && !is_exponent_integer) {
    std::stringstream str ;
    str << "Current values of the literals: " << std::endl ;
    std::map<bioString,bioReal> m = getAllLiteralValues() ;
    for (std::map<bioString,bioReal>::iterator i = m.begin() ;
	  i != m.end() ;
	  ++i) {
	  str << i->first << " = " << i->second << std::endl ;
    }
    if (rowIndex != NULL) {
	  str << "row number: " << *rowIndex << ", ";
    }
    str << "Negative base with a non-integer exponent is invalid in the real number domain: " << childResult->f << " ** " << exponent << std::endl ;
    throw bioExceptions(__FILE__,__LINE__,str.str()) ;
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




bioString bioExprPower::print(bioBoolean hp) const {
  std::stringstream str ;
  if (hp) {
    str << "^(" << left->print(hp) << "*" << right->print(hp) << ")" ;

  }
  else {
    str << "(" << left->print(hp) << "^" << right->print(hp) << ")" ;
  } 
  return str.str() ;
}
