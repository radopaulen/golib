// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FPRGUROBI_HPP
#define MC__FPRGUROBI_HPP

#include <map>
#include "gurobi_c++.h"

#include "fprelax.hpp"

extern "C"{
  #include <fenv.h>
  int fedisableexcept( int );
}

#define MC__FPRGUROBI_RELAX_AUXILIARY_BOUNDS 0.

namespace mc
{

//! @brief C++ template class for the definition and linear relaxation of factorable programs, and the solution of these relaxations using Gurobi
////////////////////////////////////////////////////////////////////////
//! mc::FPRGurobi<T> is a C++ class derived from FPRelax<T> that allows
//! definition of factorable programs, generation of linear relaxations
//! for such programs, and solution of these relaxations using Gurobi
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPRGurobi: public FPRelax<T>
////////////////////////////////////////////////////////////////////////
{
private:

  GRBEnv* _env;
  GRBModel* _model;

  typedef std::map<long,GRBVar> t_GRBVar;
  typedef typename t_GRBVar::iterator it_GRBVar;
  typedef typename t_GRBVar::const_iterator cit_GRBVar;
  t_GRBVar _var;
  t_GRBVar _aux;

  typedef std::multimap< const FPOp<T>*, GRBConstr, lt_FPOp<T> > t_GRBCut;
  typedef typename t_GRBCut::iterator it_GRBCut;
  typedef typename t_GRBCut::const_iterator cit_GRBCut;
  t_GRBCut _cut;

public:

  /** @ingroup FP
   *  @{
   */
  //! @brief Default Constructor
  FPRGurobi()
    {
      _env = new GRBEnv();
      _model = new GRBModel(*_env);
    }

  //! @brief Default Constructor
  ~FPRGurobi()
    {
      delete _model;
      delete _env;
    }

  //! @brief Reset the model (all variables and operations)
  void reset()
    {
      FPRelax<T>::reset();
      delete _model; _model = new GRBModel( *_env );
      _var.clear(); _aux.clear(); _cut.clear();
    }

  //! @brief Update bounds for a given variable
  void update_bounds
    ( FPVar<T>&pVar, const T&Bnd )
    { _update_var( &pVar, Bnd ); }

  //! @brief Set objective function in factorable program
  void set_objective
    ( const typename FPRGurobi<T>::OBJTYPE type, const FPVar<T>&objVar );

  //! @brief Remove constraint <a>constr</a> from factorable program
  bool remove_constraint
    ( FPOp<T>* constr )
    { std::pair<it_GRBCut,it_GRBCut> pitcut = _cut.equal_range( constr );
      _model->update();
      for( it_GRBCut itcutp=pitcut.first; itcutp!=pitcut.second; ){
        it_GRBCut itcut = itcutp;
        itcutp++;
        _model->remove( (*itcut).second );
        _cut.erase( itcut );
      }
      //return( pitcut.first != pitcut.second );
      return FPRelax<T>::remove_constraint( constr );
    }

  // //! @brief Solve (linear) relaxation of factorable program -- return value is status
  //void generate_cuts
  //  ( bool reset=false )
  //  {     
  //    FPRelax<T>::generate_cuts( reset );
  //    _model->update();
  //  }

  //! @brief Solve (linear) relaxation of factorable program -- return value is status
  int solve
    ( const int FPMethod=-1 )
    {
      _model->getEnv().set( GRB_IntParam_Method,            FPMethod );
      _model->getEnv().set( GRB_IntParam_OutputFlag,        FPRelax<T>::options.SOLVER_DISPLAY );
      _model->getEnv().set( GRB_DoubleParam_FeasibilityTol, FPRelax<T>::options.LP_FEASTOL );
      _model->getEnv().set( GRB_DoubleParam_OptimalityTol,  FPRelax<T>::options.LP_OPTIMTOL );
      _model->getEnv().set( GRB_DoubleParam_MIPGap,         FPRelax<T>::options.MILP_RELGAP );
      _model->getEnv().set( GRB_DoubleParam_MIPGapAbs,      FPRelax<T>::options.MILP_ABSGAP );
      _model->getEnv().set( GRB_IntParam_Presolve,          FPRelax<T>::options.SOLVER_PRESOLVE?-1:0  );
      if( FPRelax<T>::options.SOLVER_OUTPUT_FILE != "" ){
        _model->update();
        _model->write( FPRelax<T>::options.SOLVER_OUTPUT_FILE );
      }
      fedisableexcept(FE_ALL_EXCEPT);
      _model->optimize();
      return _model->get(GRB_IntAttr_Status);
    }

  //! @brief Get dual variable in relaxation of constraint <a>constr</a>
  double get_dual
    ( const GRBConstr& constr ) const
    {
      return constr.get(GRB_DoubleAttr_Pi);
    }

  //! @brief Get slack variable in relaxation of constraint <a>constr</a>
  double get_slack
    ( const GRBConstr& constr ) const
    {
      return constr.get(GRB_DoubleAttr_Slack);
    }

  //! @brief Get value of variable w/ index <a>ix</a> in relaxation
  double get_variable
    ( const long ix ) const
    {
      cit_GRBVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return (*it).second.get(GRB_DoubleAttr_X);
    }

  //! @brief Get reduced cost of variable w/ index <a>ix</a> in relaxation
  double get_reduced_cost
    ( const long ix ) const
    {
      cit_GRBVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return (*it).second.get(GRB_DoubleAttr_RC);
    }

  //! @brief Get basis status of variable w/ index <a>ix</a> in relaxation
  int get_basis_status
    ( const long ix ) const
    {
      cit_GRBVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return (*it).second.get(GRB_IntAttr_VBasis);
    }

  //! @brief Get optimal objective value in relaxation
  double get_objective() const
    {
      return _model->get( GRB_DoubleAttr_ObjVal );
    }

  //! @brief Get runtime for relaxation solution
  double get_runtime() const
    {
      return _model->get( GRB_DoubleAttr_Runtime );
    }

  //! @brief Retreive pointer to Gurobi model
  const GRBModel* get_model() const
    { return _model; }
  /** @} */

private:

  //! @brief Appends new auxiliary variable
  virtual void _append_aux
    ( FPVar<T>*pAux );

  //! @brief Appends new original variable
  virtual void _append_var
    ( FPVar<T>*pVar );

  //! @brief Update the bounds of the variable <a>pVar</a> as <a>Bnd</a> 
  virtual void _update_var
    ( FPVar<T>*pVar, const T&Bnd );

  //! @brief Appends new relaxation cut w/ 1 variable
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1 )
    {
      FPCut<T>*pCut = FPRelax<T>::_append_cut( op, type, b, id1, a1 );
      _append_cut( pCut );
      return pCut;
    }
  //! @brief Appends new relaxation cut w/ 2 variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1,
      const typename FPVar<T>::pt_idVar id2, const double a2 )
    {
      FPCut<T>*pCut = FPRelax<T>::_append_cut( op, type, b, id1, a1, id2, a2 );
      _append_cut( pCut );
      return pCut;
    }
  //! @brief Appends new relaxation cut w/ 3 variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1,
      const typename FPVar<T>::pt_idVar id2, const double a2,
      const typename FPVar<T>::pt_idVar id3, const double a3 )
    {
      FPCut<T>*pCut = FPRelax<T>::_append_cut( op, type, b, id1, a1, id2, a2,
        id3, a3 );
      _append_cut( pCut );
      return pCut;
    }
  //! @brief Appends new relaxation cut w/ <a>n</a> variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type, const double b,
      const unsigned int n, const typename FPVar<T>::pt_idVar* id,
      const double*a )
    {
      FPCut<T>*pCut = FPRelax<T>::_append_cut( op, type, b, n, id, a );
      _append_cut( pCut );
      return pCut;
    }
  //! @brief Appends new relaxation cut <a>pCut</a>
  void _append_cut
    ( const FPCut<T>* pCut );
};

template <typename T> inline void
FPRGurobi<T>::_append_cut
( const FPCut<T>* pCut )
{
  GRBVar* VarSOS = 0;
  double* WeiSOS = 0;
  int TypSOS = 0;
  if( pCut->type() == FPCut<T>::SOS1 || pCut->type() == FPCut<T>::SOS2 ){
    VarSOS = new GRBVar[pCut->nvar()];
    WeiSOS = new double[pCut->nvar()];
    TypSOS = ( pCut->type() == FPCut<T>::SOS1? GRB_SOS_TYPE1: GRB_SOS_TYPE2 );
  }

  GRBLinExpr lhs;
  for( unsigned int k=0; k<pCut->nvar(); k++ ){
    switch( pCut->idvar()[k].first ){
      case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
      { it_GRBVar ivar = _var.find( pCut->idvar()[k].second );
        if( ivar==_var.end() )
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
        GRBVar&Var = (*ivar).second;
        lhs += GRBLinExpr( Var, pCut->coef()[k] );
        if( pCut->type() == FPCut<T>::SOS1 || pCut->type() == FPCut<T>::SOS2 ){
          VarSOS[k] = Var; WeiSOS[k] = (double)k/(double)pCut->nvar();
        }
        break; }
      default:
      { it_GRBVar iaux = _aux.find( pCut->idvar()[k].second );
        if( iaux==_aux.end() )
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
        GRBVar&Aux = (*iaux).second;
        lhs += GRBLinExpr( Aux, pCut->coef()[k] );
        if( pCut->type() == FPCut<T>::SOS1 || pCut->type() == FPCut<T>::SOS2 ){
          VarSOS[k] = Aux; WeiSOS[k] = (double)k/(double)pCut->nvar();
        }
        break; }
    }
  }

  switch( pCut->type() ){
    case FPCut<T>::SOS1:
    case FPCut<T>::SOS2:
      _model->addSOS( VarSOS, WeiSOS, pCut->nvar(), TypSOS );
      delete [] VarSOS;
      delete [] WeiSOS;
    case FPCut<T>::EQ:
      _cut.insert( std::make_pair( pCut->op(), _model->addConstr( lhs,
        GRB_EQUAL, pCut->rhs() ) ) );
      break;
    case FPCut<T>::LE:
      _cut.insert( std::make_pair( pCut->op(), _model->addConstr( lhs,
        GRB_LESS_EQUAL, pCut->rhs() ) ) );
      break;
    case FPCut<T>::GE:
      _cut.insert( std::make_pair( pCut->op(), _model->addConstr( lhs,
        GRB_GREATER_EQUAL, pCut->rhs() ) ) );
      break;
  }
}

template <typename T> inline void
FPRGurobi<T>::_append_aux
( FPVar<T>*pAux )
{
  FPRelax<T>::_append_aux( pAux );

  switch( pAux->id().first ){
    case FPVar<T>::AUXCONT:
    { it_GRBVar itv = _aux.find( pAux->id().second );
      if( itv != _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      _aux.insert( std::make_pair( pAux->id().second, _model->addVar(
        Op<T>::l(pAux->num().I)-MC__FPRGUROBI_RELAX_AUXILIARY_BOUNDS,
        Op<T>::u(pAux->num().I)+MC__FPRGUROBI_RELAX_AUXILIARY_BOUNDS,
        0.0, GRB_CONTINUOUS, pAux->name() ) ) );
      _model->update();
      break;
    }
    
    case FPVar<T>::AUXBIN:
    { it_GRBVar itv = _aux.find( pAux->id().second );
      if( itv != _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      _aux.insert( std::make_pair( pAux->id().second, _model->addVar(
        Op<T>::l(pAux->num().I), Op<T>::u(pAux->num().I), 0.0,
        GRB_BINARY, pAux->name() ) ) );
      _model->update();
      break;
    }
    
    case FPVar<T>::AUXINT: case FPVar<T>::AUXREAL:
      break;

    default:
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }
}

template <typename T> inline void
FPRGurobi<T>::_append_var
( FPVar<T>*pVar )
{
  FPRelax<T>::_append_var( pVar );
  switch( pVar->id().first ){
    case FPVar<T>::VARCONT:
    { it_GRBVar itv = _var.find( pVar->id().second );
      if( itv != _var.end() ){
        _model->remove( itv->second );
        _var.erase( itv );
      }
      GRBVar newVar = _model->addVar( Op<T>::l(pVar->num().I),
        Op<T>::u(pVar->num().I), 0.0, GRB_CONTINUOUS, pVar->name() );
      _var.insert( std::make_pair( pVar->id().second, newVar ) );
      _model->update();
      break;
    }

    case FPVar<T>::VARBIN:
    { it_GRBVar itv = _var.find( pVar->id().second );
      if( itv != _var.end() ){
        _model->remove( itv->second );
        _var.erase( itv );
      }
      _var.insert( std::make_pair( pVar->id().second, _model->addVar(
        Op<T>::l(pVar->num().I), Op<T>::u(pVar->num().I), 0.0,
        GRB_BINARY, pVar->name() ) ) );
      _model->update();
      break;
    }

    default:
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }
}

template <typename T> inline void
FPRGurobi<T>::_update_var
( FPVar<T>*pVar, const T&Bnd )
{
  switch( pVar->id().first ){
    case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
    { it_GRBVar itv = _var.find( pVar->id().second );
      if( itv == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      (*itv).second.set( GRB_DoubleAttr_LB, Op<T>::l(Bnd) ); 
      (*itv).second.set( GRB_DoubleAttr_UB, Op<T>::u(Bnd) );
      break;
    }
    case FPVar<T>::AUXCONT: case FPVar<T>::AUXBIN:
    { it_GRBVar itv = _aux.find( pVar->id().second );
      if( itv == _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      (*itv).second.set( GRB_DoubleAttr_LB, Op<T>::l(Bnd) ); 
      (*itv).second.set( GRB_DoubleAttr_UB, Op<T>::u(Bnd) );
      break;
    }
    default:;
  }
  return FPRelax<T>::_update_var( pVar, Bnd );
}

template <typename T> inline void
FPRGurobi<T>::set_objective
( const typename FPRGurobi<T>::OBJTYPE type, const FPVar<T>&objVar )
{
  // Declare objective in LP model
  FPRelax<T>::set_objective( type, objVar );

  double vobj = 0.;
  it_GRBVar iobjVar;
  switch( objVar.id().first ){
    case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
      iobjVar = _var.find( objVar.id().second );
      if( iobjVar==_var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      break;
    default:
      iobjVar = _aux.find( objVar.id().second );
      if( iobjVar==_aux.end() )
        vobj = ( objVar.num().t==FPVarNum<T>::INT? objVar.num().n: objVar.num().x );
      break;
  }
  
  switch( type ){
    case FPRGurobi<T>::MIN:
      _model->set( GRB_IntAttr_ModelSense, 1 ); break;
    case FPRGurobi<T>::MAX:
      _model->set( GRB_IntAttr_ModelSense, -1 ); break;
  }
  if( iobjVar == _aux.end() )
    _model->setObjective( GRBLinExpr( vobj ) );
  else
    //(*iobjVar).second.set( GRB_DoubleAttr_Obj,  1. );
    _model->setObjective( GRBLinExpr( (*iobjVar).second,  1. ) );      
}

} // namespace mc

#endif
