//-*-c++-*------------------------------------------------------------
//
// File name : bioExprLogzero.cc
// @date   Mon Oct 24 09:48:09 2022
// @author Michel Bierlaire
// @version Revision 1.0
//
//--------------------------------------------------------------------

#include "bioExprLogzero.h"
#include "bioDebug.h"
#include "bioExceptions.h"
#include "bioConstants.h"

#include <cmath>
#include <sstream>

bioExprLogzero::bioExprLogzero(bioExpression* c) :
  child(c) {
  listOfChildren.push_back(c) ;
}

bioExprLogzero::~bioExprLogzero() {
  
}

const bioDerivatives* bioExprLogzero::getValueAndDerivatives(std::vector<bioUInt> literalIds,
							     bioBoolean gradient,
							     bioBoolean hessian) {
  
theDerivatives.with_g = gradient ;
  theDerivatives.with_h = hessian ;

  bioUInt n = literalIds.size() ;
  theDerivatives.resize(n) ;

  static const bioReal upper_bound = constants::get_upper_bound();
  static const bioReal almost_zero = constants::get_almost_zero();
  static const bioReal slope = (log(almost_zero) + upper_bound) / almost_zero;


  const bioDerivatives* childResult = child->getValueAndDerivatives(literalIds,gradient,hessian) ;
  bioReal cf = childResult->f ;
  validate(cf) ;
  if (cf == 0.0) {
    theDerivatives.f = 0.0 ;
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
        theDerivatives.g[i] = 0.0 ;
        if (hessian) {
          theDerivatives.h[i][i] = 0 ;
	      for (bioUInt j = 0 ; j < n ; ++j) {
	        theDerivatives.h[i][j] = theDerivatives.h[j][i] = 0.0 ;
	      }
	    }
	  }
	}
	return &theDerivatives ;
  }
  if (cf < 0) {
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
    str << "Cannot take the log of a non positive number [" << childResult->f << "]" << std::endl ;
      throw bioExceptions(__FILE__,__LINE__,str.str()) ;
  }


  if (cf < almost_zero) {
    theDerivatives.f = log(almost_zero) * cf / almost_zero
      - (1 - cf/almost_zero) * upper_bound ;
    if (gradient) {
      for (bioUInt i = 0 ; i < n ; ++i) {
        theDerivatives.g[i] = slope * childResult->g[i] ;
        if (hessian) {
          theDerivatives.h[i][i] = 0 ;
	      for (bioUInt j = i ; j < n ; ++j) {
	        theDerivatives.h[i][j] = slope * childResult->h[i][j];
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

  theDerivatives.f = log(cf) ;
  if (gradient) {
    for (bioUInt i = 0 ; i < n ; ++i) {
      theDerivatives.g[i] = childResult->g[i] / cf ;
      if (hessian) {
	    for (bioUInt j = i ; j < n ; ++j) {
	      bioReal fsquare = cf * cf ;
	      theDerivatives.h[i][j] = childResult->h[i][j] / cf - childResult->g[i] *  childResult->g[j] / fsquare ;
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

bioString bioExprLogzero::print(bioBoolean hp) const {
  std::stringstream str ; 
  str << "logzero(" << child->print(hp) << ")";
  return str.str() ;
}

