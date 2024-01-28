// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FPRCPLEX_HPP
#define MC__FPRCPLEX_HPP

#include <map>
#include "ilcplex/ilocplex.h"

#include "fprelax.hpp"

#define MC__FPRCPLEX_RELAX_AUXILIARY_BOUNDS 0.

namespace mc
{

//! @brief C++ template class for the definition and linear relaxation of factorable programs, and the solution of these relaxations using Cplex
////////////////////////////////////////////////////////////////////////
//! mc::FPRCplex<T> is a C++ class derived from FPRelax<T> that allows
//! definition of factorable programs, generation of linear relaxations
//! for such programs, and solution of these relaxations using Cplex
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPRCplex: public FPRelax<T>
////////////////////////////////////////////////////////////////////////
{
private:

  IloEnv* _env;
  IloModel* _model;
  IloCplex* _cplex;
  IloNum _time;
  IloObjective _curObj;
  bool _setObj;
  
  typedef std::map<long,IloNumVar> t_IloNumVar;
  typedef typename t_IloNumVar::iterator it_IloNumVar;
  typedef typename t_IloNumVar::const_iterator cit_IloNumVar;
  t_IloNumVar _var;
  t_IloNumVar _aux;

  typedef std::multimap< const FPOp<T>*, IloRange, lt_FPOp<T> > t_IloRange;
  typedef typename t_IloRange::iterator it_IloRange;
  typedef typename t_IloRange::const_iterator cit_IloRange;
  t_IloRange _cut;

public:

  /** @ingroup FP
   *  @{
   */
  //! @brief Default Constructor
  FPRCplex()
    {
      _env = new IloEnv;
      _model = new IloModel(*_env);
      _cplex = new IloCplex(*_env);
      _setObj = false;
    }
 
  //! @brief Default Constructor
  ~FPRCplex()
    {
      delete _model;
      delete _cplex;
      _env->end();
      delete _env;
    }

  //! @brief Reset the model (all variables and operations)
  void reset()
    {
      FPRelax<T>::reset();
      delete _model; delete _cplex;
      _env->end();
      delete _env;
      _env = new IloEnv;
      _model = new IloModel( *_env );
      _cplex = new IloCplex( *_env );
      _setObj = false;
      _var.clear(); _aux.clear(); _cut.clear();
    }

  //! @brief Update bounds for a given variable
  void update_bounds
    ( FPVar<T>&pVar, const T&Bnd )
    { _update_var( &pVar, Bnd ); }

  //! @brief Set objective function in factorable program
  void set_objective
    ( const typename FPRCplex<T>::OBJTYPE type, const FPVar<T>&objVar );

  //! @brief Remove constraint <a>constr</a> from factorable program
  bool remove_constraint
    ( FPOp<T>* constr )
    { std::pair<it_IloRange,it_IloRange> pitcut = _cut.equal_range( constr );
      for( it_IloRange itcutp=pitcut.first; itcutp!=pitcut.second; ){
        it_IloRange itcut = itcutp;
        itcutp++;
        _model->remove( (*itcut).second );
        _cut.erase( itcut );
      }
      //return( pitcut.first != pitcut.second );
      return FPRelax<T>::remove_constraint( constr );
    }

  //! @brief Solve (linear) relaxation of factorable program -- return value is status
  typename IloAlgorithm::Status solve
    (  const int FPMethod=0, const bool outputFlag=false,
       const std::string modelFileName="" )
    {
      _cplex->extract(*_model);
      _cplex->setWarning( outputFlag? std::cout: _env->getNullStream() );
      _cplex->setOut( outputFlag? std::cout: _env->getNullStream() );
      _cplex->setParam( IloCplex::RootAlg, FPMethod );
      _cplex->setParam( IloCplex::EpOpt, 1e-9 );
      _cplex->setParam( IloCplex::EpRHS, 1e-9 );
      if( modelFileName != "" ) _cplex->exportModel( modelFileName.c_str() );
      _time = -_cplex->getCplexTime();
      _cplex->solve();
      _time += _cplex->getCplexTime();
      return _cplex->getStatus();
    }

  //! @brief Get dual variable in relaxation of constraint <a>constr</a>
  IloNum get_dual
    ( const IloRange& constr ) const
    {
      return _cplex->getDual( constr );
    }

  //! @brief Get slack variable in relaxation of constraint <a>constr</a>
  IloNum get_slack
    ( const IloRange& constr ) const
    {
      return _cplex->getSlack( constr );
    }

  //! @brief Get value of variable w/ index <a>ix</a> in relaxation
  IloNum get_variable
    ( const long ix ) const
    {
      cit_IloNumVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return _cplex->getValue( (*it).second );
    }

  //! @brief Get reduced cost of variable w/ index <a>ix</a> in relaxation
  IloNum get_reduced_cost
    ( const long ix ) const
    {
      cit_IloNumVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return _cplex->getReducedCost( (*it).second );
    }

  //! @brief Get basis status of variable w/ index <a>ix</a> in relaxation
  typename IloCplex::BasisStatus get_basis_status
    ( const long ix ) const
    {
      cit_IloNumVar it = _var.find( ix );
      if( it == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
      return _cplex->getBasisStatus( (*it).second );
    }

  //! @brief Get LP objective value
  IloNum get_objective() const
    {
      return _cplex->getObjValue();
    }

  //! @brief Get LP runtime
  double get_runtime() const
    {
      return _time;
    }

  //! @brief Retreive pointer to Cplex model
  const IloModel* get_model() const
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
FPRCplex<T>::_append_cut
( const FPCut<T>* pCut )
{
  bool isSOS = ( pCut->type() == FPCut<T>::SOS1
              || pCut->type() == FPCut<T>::SOS2 );
  IloNumVarArray VarSOS( *_env, (isSOS? pCut->nvar(): 0) );
  IloNumArray WeiSOS( *_env, (isSOS? pCut->nvar(): 0) );

  IloExpr lhs( *_env );
  for( unsigned int k=0; k<pCut->nvar(); k++ ){
    switch( pCut->idvar()[k].first ){
      case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
      { it_IloNumVar ivar = _var.find( pCut->idvar()[k].second );
        if( ivar==_var.end() )
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
        lhs += pCut->coef()[k] * (*ivar).second;
        if( isSOS ){ VarSOS[k] = (*ivar).second;
                     WeiSOS[k] = (double)k/(double)pCut->nvar(); };
        break; }
      default:
      { it_IloNumVar iaux = _aux.find( pCut->idvar()[k].second );
        if( iaux==_aux.end() )
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
        lhs += pCut->coef()[k] * (*iaux).second;
        if( isSOS ){ VarSOS[k] = (*iaux).second;
                     WeiSOS[k] = (double)k/(double)pCut->nvar(); };
        break; }
    }
  }
  IloRange constr;
  switch( pCut->type() ){
    case FPCut<T>::SOS1:
      _model->add( IloSOS1( *_env, VarSOS, WeiSOS ) );
      constr = ( lhs == pCut->rhs() ); break;
    case FPCut<T>::SOS2:
      _model->add( IloSOS2( *_env, VarSOS, WeiSOS ) );
      constr = ( lhs == pCut->rhs() ); break;
    case FPCut<T>::EQ:
      constr = ( lhs == pCut->rhs() ); break;
    case FPCut<T>::LE:
      constr = ( lhs <= pCut->rhs() ); break;
    case FPCut<T>::GE:
      constr = ( lhs >= pCut->rhs() ); break;
  }
  _model->add( constr );
  _cut.insert( std::make_pair( pCut->op(), constr ) );
}

template <typename T> inline void
FPRCplex<T>::_append_aux
( FPVar<T>*pAux )
{
  FPRelax<T>::_append_aux( pAux );

  switch( pAux->id().first ){
    case FPVar<T>::AUXCONT:
    { it_IloNumVar itv = _aux.find( pAux->id().second );
      if( itv != _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      IloNumVar newAux( *_env,
        Op<T>::l(pAux->num().I)-MC__FPRCPLEX_RELAX_AUXILIARY_BOUNDS,
        Op<T>::u(pAux->num().I)+MC__FPRCPLEX_RELAX_AUXILIARY_BOUNDS,
        ILOFLOAT, pAux->name().c_str() );
      _model->add( newAux );
      _aux.insert( std::make_pair( pAux->id().second, newAux ) );
      break;
    }

    case FPVar<T>::AUXBIN:
    { it_IloNumVar itv = _aux.find( pAux->id().second );
      if( itv != _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      IloNumVar newAux( *_env, Op<T>::l(pAux->num().I),
        Op<T>::u(pAux->num().I), ILOBOOL, pAux->name().c_str() );
      _model->add( newAux );
      _aux.insert( std::make_pair( pAux->id().second, newAux ) );
      break;
    }

    case FPVar<T>::AUXINT: case FPVar<T>::AUXREAL:
      break;

    default:
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }
}

template <typename T> inline void
FPRCplex<T>::_append_var
( FPVar<T>*pVar )
{
  FPRelax<T>::_append_var( pVar );

  switch( pVar->id().first ){
    case FPVar<T>::VARCONT:
    { it_IloNumVar itv = _var.find( pVar->id().second );
      if( itv != _var.end() ){
        _model->remove( itv->second );
        _var.erase( itv );
      }
      //char vname[pVar->name().size()+2];
      //strcpy( vname, pVar->name().c_str() );
      IloNumVar newVar( *_env, Op<T>::l(pVar->num().I),
        Op<T>::u(pVar->num().I), ILOFLOAT, pVar->name().c_str() );
      _model->add( newVar );
      _var.insert( std::make_pair( pVar->id().second, newVar ) );
      break;
    }

    case FPVar<T>::VARBIN:
    { it_IloNumVar itv = _var.find( pVar->id().second );
      if( itv != _var.end() ){
        _model->remove( itv->second );
        _var.erase( itv );
      }
      IloNumVar newVar( *_env, Op<T>::l(pVar->num().I),
        Op<T>::u(pVar->num().I), ILOBOOL, pVar->name().c_str() );
      _model->add( newVar );
      break;
    }

    default:
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }
}

template <typename T> inline void
FPRCplex<T>::_update_var
( FPVar<T>*pVar, const T&Bnd )
{
  switch( pVar->id().first ){
    case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
    { it_IloNumVar itv = _var.find( pVar->id().second );
      if( itv == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      (*itv).second.setBounds( Op<T>::l(Bnd), Op<T>::u(Bnd) );
      break;
    }
    case FPVar<T>::AUXCONT: case FPVar<T>::AUXBIN:
    { it_IloNumVar itv = _aux.find( pVar->id().second );
      if( itv == _aux.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      (*itv).second.setBounds( Op<T>::l(Bnd), Op<T>::u(Bnd) );
      break;
    }
    default:;
  }
  return FPRelax<T>::_update_var( pVar, Bnd );
}

template <typename T> inline void
FPRCplex<T>::set_objective
( const typename FPRCplex<T>::OBJTYPE type, const FPVar<T>&objVar )
{
  // Declare objective in LP model
  FPRelax<T>::set_objective( type, objVar );

  IloNum vobj = 0;
  it_IloNumVar iobjVar;
  switch( objVar.id().first ){
    case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
      iobjVar = _var.find( objVar.id().second );
      if( iobjVar == _var.end() )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      break;
    default:
      iobjVar = _aux.find( objVar.id().second );
      if( iobjVar == _aux.end() )
        vobj = ( objVar.num().t==FPVarNum<T>::INT? objVar.num().n: objVar.num().x );
      break;
  }
  
  if( _setObj ) _model->remove( _curObj );
  switch( type ){
    case FPRCplex<T>::MIN:
      if( iobjVar == _aux.end() )
        _curObj = IloMinimize( *_env, vobj );
      else
        _curObj = IloMinimize( *_env, (*iobjVar).second );
      break;
    case FPRCplex<T>::MAX:
      if( iobjVar == _aux.end() )
        _curObj = IloMaximize( *_env, vobj );
      else
        _curObj = IloMaximize( *_env, (*iobjVar).second );
      break;
  }
  _model->add( _curObj );
  _setObj = true;
  
}

} // namespace mc

#endif
