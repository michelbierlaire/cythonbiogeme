//-*-c++-*------------------------------------------------------------
//
// File name : bioMemoryManagement.h
// @date   Sat Sep 26 12:19:14 2020
// @author Michel Bierlaire
//
//--------------------------------------------------------------------

#ifndef bioMemoryManagement_h
#define bioMemoryManagement_h

#include <vector>
#include <set>
#include <unordered_map>
#include "bioTypes.h"
#include "bioString.h"
#include "bioExprLinearUtility.h"

class bioExpression ;
class bioExprFreeParameter ;
class bioExprFixedParameter ;
class bioExprVariable ;
class bioExprDraws ;
class bioExprRandomVariable ;
class bioExprNumeric ;
class bioExprPlus ;
class bioExprMinus ;
class bioExprTimes ;
class bioExprDivide ;
class bioExprPower ;
class bioExprPowerConstant ;
class bioExprAnd ;
class bioExprOr ;
class bioExprEqual ;
class bioExprNotEqual ;
class bioExprLess ;
class bioExprLessOrEqual ;
class bioExprGreater ;
class bioExprGreaterOrEqual ;
class bioExprMin ;
class bioExprMax ;
class bioExprUnaryMinus ;
class bioExprMontecarlo ;
class bioExprNormalCdf ;
class bioExprPanelTrajectory ;
class bioExprExp ;
class bioExprSin ;
class bioExprCos ;
class bioExprLog ;
class bioExprLogzero ;
class bioExprBelongsTo ;
class bioExprDerive ;
class bioExprIntegrate ;
class bioExprLogLogit ;
class bioExprLogLogitFullChoiceSet ;
class bioExprMultSum ;
class bioExprConditionalSum ;
class bioExprElem ;

class bioSeveralExpressions ;

class bioMemoryManagement {

public:
  static bioMemoryManagement* the() ;
  ~bioMemoryManagement() ;
  void releaseAllMemory() ;
  bioExprFreeParameter* get_bioExprFreeParameter(bioUInt literalId,
						 bioUInt parameterId,
						 bioString name) ;
  bioExprFixedParameter* get_bioExprFixedParameter(bioUInt literalId,
						   bioUInt parameterId,
						   bioString name) ; 
  bioExprVariable* get_bioExprVariable(bioUInt literalId,
				       bioUInt variableId,
				       bioString name) ; 
  bioExprDraws* get_bioExprDraws(bioUInt uniqueId,
				 bioUInt drawId,
				 bioString name) ;
  bioExprRandomVariable* get_bioExprRandomVariable(bioUInt literalId,
						   bioUInt id,
						   bioString name) ;
  bioExprNumeric* get_bioExprNumeric(bioReal v) ;
  bioExprPlus* get_bioExprPlus(bioExpression* ell, bioExpression* r) ;
  bioExprMinus* get_bioExprMinus(bioExpression* ell, bioExpression* r) ;
  bioExprTimes* get_bioExprTimes(bioExpression* ell, bioExpression* r) ;
  bioExprDivide* get_bioExprDivide(bioExpression* ell, bioExpression* r) ;
  bioExprPower* get_bioExprPower(bioExpression* ell, bioExpression* r) ;
  bioExprPowerConstant* get_bioExprPowerConstant(bioExpression* ell, bioReal exponent) ;
  bioExprAnd* get_bioExprAnd(bioExpression* ell, bioExpression* r) ;
  bioExprOr* get_bioExprOr(bioExpression* ell, bioExpression* r) ;
  bioExprEqual* get_bioExprEqual(bioExpression* ell, bioExpression* r) ;
  bioExprNotEqual* get_bioExprNotEqual(bioExpression* ell, bioExpression* r) ;
  bioExprLess* get_bioExprLess(bioExpression* ell, bioExpression* r) ;
  bioExprLessOrEqual* get_bioExprLessOrEqual(bioExpression* ell, bioExpression* r) ;
  bioExprGreater* get_bioExprGreater(bioExpression* ell, bioExpression* r) ;
  bioExprGreaterOrEqual* get_bioExprGreaterOrEqual(bioExpression* ell, bioExpression* r) ;
  bioExprMin* get_bioExprMin(bioExpression* ell, bioExpression* r) ;
  bioExprMax* get_bioExprMax(bioExpression* ell, bioExpression* r) ;
  bioExprUnaryMinus* get_bioExprUnaryMinus(bioExpression* ell) ;
  bioExprMontecarlo* get_bioExprMontecarlo(bioExpression* ell) ;
  bioExprNormalCdf* get_bioExprNormalCdf(bioExpression* ell) ;
  bioExprPanelTrajectory* get_bioExprPanelTrajectory(bioExpression* ell) ;
  bioExprExp* get_bioExprExp(bioExpression* ell) ;
  bioExprSin* get_bioExprSin(bioExpression* ell) ;
  bioExprCos* get_bioExprCos(bioExpression* ell) ;
  bioExprLog* get_bioExprLog(bioExpression* ell) ;
  bioExprLogzero* get_bioExprLogzero(bioExpression* ell) ;
  bioExprDerive* get_bioExprDerive(bioExpression* c, bioUInt lid) ;
  bioExprBelongsTo* get_bioExprBelongsTo(bioExpression* c, const std::set<bioReal>& the_set) ;
  bioExprIntegrate* get_bioExprIntegrate(bioExpression* c, bioUInt lid) ;
  bioExprLinearUtility* get_bioExprLinearUtility(std::vector<bioLinearTerm> t) ;
  bioExprLogLogit* get_bioExprLogLogit(bioExpression* c,
				       std::map<bioUInt,bioExpression*> u,
				       std::map<bioUInt,bioExpression*> a) ;
  bioExprLogLogitFullChoiceSet* get_bioExprLogLogitFullChoiceSet(bioExpression* c,
								 std::map<bioUInt,bioExpression*> u) ;
  bioExprConditionalSum* get_bioExprConditionalSum(std::unordered_map<bioExpression*, bioExpression*> d) ;
  bioExprMultSum* get_bioExprMultSum(std::vector<bioExpression*> e) ;
  bioExprElem* get_bioExprElem(bioExpression* k, std::map<bioUInt,bioExpression*> d) ;
  bioSeveralExpressions* get_bioSeveralExpressions(std::vector<bioExpression*> exprs) ;
private:
  bioMemoryManagement() ;
  std::vector<std::unique_ptr<bioExprFreeParameter> > a_bioExprFreeParameter ;
  std::vector<std::unique_ptr<bioExprFixedParameter> > a_bioExprFixedParameter ;
  std::vector<std::unique_ptr<bioExprVariable> > a_bioExprVariable ;
  std::vector<std::unique_ptr<bioExprDraws> > a_bioExprDraws ;
  std::vector<std::unique_ptr<bioExprRandomVariable> > a_bioExprRandomVariable ;
  std::vector<std::unique_ptr<bioExprNumeric> > a_bioExprNumeric ;
  std::vector<std::unique_ptr<bioExprPlus> > a_bioExprPlus ;
  std::vector<std::unique_ptr<bioExprMinus> > a_bioExprMinus ;
  std::vector<std::unique_ptr<bioExprTimes> > a_bioExprTimes ;
  std::vector<std::unique_ptr<bioExprDivide> > a_bioExprDivide ;
  std::vector<std::unique_ptr<bioExprPower> > a_bioExprPower ;
  std::vector<std::unique_ptr<bioExprPowerConstant> > a_bioExprPowerConstant ;
  std::vector<std::unique_ptr<bioExprAnd> > a_bioExprAnd ;
  std::vector<std::unique_ptr<bioExprOr> > a_bioExprOr ;
  std::vector<std::unique_ptr<bioExprEqual> > a_bioExprEqual ;
  std::vector<std::unique_ptr<bioExprNotEqual> > a_bioExprNotEqual ;
  std::vector<std::unique_ptr<bioExprLess> > a_bioExprLess ;
  std::vector<std::unique_ptr<bioExprLessOrEqual> > a_bioExprLessOrEqual ;
  std::vector<std::unique_ptr<bioExprGreater> > a_bioExprGreater ;
  std::vector<std::unique_ptr<bioExprGreaterOrEqual> > a_bioExprGreaterOrEqual ;
  std::vector<std::unique_ptr<bioExprMin> > a_bioExprMin ;
  std::vector<std::unique_ptr<bioExprMax> > a_bioExprMax ;
  std::vector<std::unique_ptr<bioExprUnaryMinus> > a_bioExprUnaryMinus ;
  std::vector<std::unique_ptr<bioExprMontecarlo> > a_bioExprMontecarlo ;
  std::vector<std::unique_ptr<bioExprNormalCdf> > a_bioExprNormalCdf ;
  std::vector<std::unique_ptr<bioExprPanelTrajectory> > a_bioExprPanelTrajectory ;
  std::vector<std::unique_ptr<bioExprExp> > a_bioExprExp ;
  std::vector<std::unique_ptr<bioExprSin> > a_bioExprSin ;
  std::vector<std::unique_ptr<bioExprCos> > a_bioExprCos ;
  std::vector<std::unique_ptr<bioExprLog> > a_bioExprLog ;
  std::vector<std::unique_ptr<bioExprLogzero> > a_bioExprLogzero ;
  std::vector<std::unique_ptr<bioExprDerive> > a_bioExprDerive ;
  std::vector<std::unique_ptr<bioExprBelongsTo> > a_bioExprBelongsTo ;
  std::vector<std::unique_ptr<bioExprIntegrate> > a_bioExprIntegrate ;
  std::vector<std::unique_ptr<bioExprLinearUtility> > a_bioExprLinearUtility ;
  std::vector<std::unique_ptr<bioExprLogLogit> > a_bioExprLogLogit ;
  std::vector<std::unique_ptr<bioExprLogLogitFullChoiceSet> > a_bioExprLogLogitFullChoiceSet ;
  std::vector<std::unique_ptr<bioExprMultSum> > a_bioExprMultSum ;
  std::vector<std::unique_ptr<bioExprConditionalSum> > a_bioExprConditionalSum ;
  std::vector<std::unique_ptr<bioExprElem> > a_bioExprElem ;
  std::vector<std::unique_ptr<bioSeveralExpressions> > a_bioSeveralExpressions ;
};
#endif
