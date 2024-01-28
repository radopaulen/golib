// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_NLPSBB Deterministic Global Optimization of Nonlinear Programs using Spatial Branch-and-Bound
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

Consider a nonlinear program (NLP) of the form
\f{align*}

  \min_{{\bf x}\in X}\ & f({\bf x})\\
  \text{s.t.} \ & g_j({\bf x})\{\leq,=,\geq\} 0,\  j=1,\ldots,n_g\\
  & x_i^\ell\leq x_i\leq x_i^u,\ i=1,\ldots,n_x
\f}
where \f$X\subset\mathbb{R}^{n_x}\f$, and \f$f, g_1, \ldots,g_{n_g}\f$ are potentially nonconvex functions that are factorable and twice continuously differentiable on \f$X\f$. The class mc::NLPSBB in MC++ uses branch-and-bound search to solve such NLPs to global optimality. An LP relaxation of the NLP is constructed using the class mc::FPRelax and solved using either <A href="http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/">CPLEX</A> or <A href="http://www.gurobi.com/">GUROBI</A>, while local solution are obtained based on the class mc::IPOPT which calls the software package <A href="https://projects.coin-or.org/Ipopt">IPOPT</A>.

mc::NLPSBB is templated in two variable types T and NLP:
- T is the class that performs the underlying interval arithmetics. Both verified and non-verified types are supported. For verified computations, the libraries <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> and <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> are supported; the non-verified interval type mc::Interval of MC++ can also be used. More generally, any interval type can be used provided that the templated structure mc::Op is instantiated accordingly -- see the header files <tt>mcprofil.h</tt> and <tt>mcfilib.h</tt> for examples.
- NLP is the class that defines the optimization problem to be solved, which has to be a class derived from mc::NLPSTRUCT. For example, suppose we want to solve the following NLP:
\f{align*}
  \max_{\bf x}\ & x_1+x_2 \\
  \text{s.t.} \ & x_1\,x_2 \leq 4 \\
  & 3 \leq x_1 \leq 6\\
  & 0 \leq x_2 \leq 4.
\f}
The following class is defined:
\code
  #include "nlpsbb.h"
  
  const unsigned int NX = 2;
  const unsigned int NG = 1;

  class NLP : public virtual mc::NLPSTRUCT
  {
  public:
    NLP(): mc::NLPSTRUCT( NX, NG )
      {}

    template <typename U>
    std::pair<U,t_OBJ> OBJ
      ( const U*x )
      {
        return std::make_pair( x[0]+x[1], MAX );
      }

    template <typename U>
    std::pair<U,t_CTR> CTR
      ( const unsigned int ig, const U*x )
      {
        assert( ig < ng() );
        switch( ig ){
          case 0: return std::make_pair( x[0]*x[1]-4, LE );
          default: throw std::runtime_error("invalid size");
        }
      }
  };
\endcode
where NX and NG stand for the number of parameters and constraints, respectively; and the templated functions OBJ and CTR define the objective function and left-hand side of the constraints, respectively.

\section sec_NLPBB_solve How to Solve an NLP Model using mc::NLPSBB?

Start by defining the parameter set \f$X\f$ and initial guess \f$x_0\f$:

\code
  std::pair<double,double> X[NX];
  X[0] = std::make_pair( 0., 6. );
  X[1] = std::make_pair( 0., 4. );
  double x0[NX];
  x0[0] = 1.;
  x0[1] = 2.;
\endcode

Then, an mc::NLPSBB object is defined---the default interval type mc::Interval is used here:

\code
  #include "interval.h"
  typedef mc::Interval I;
  mc::NLPSBB<I,NLP> pSBB = new mc::NLPSBB<I,NLP>;
\endcode

The NLP model is optimized as follows:

\code
  std::cout << pSBB;
  std::pair<double,const double*> optimum = pSBB->solve( X, x0 );
\endcode

where the arguments of the NLPSBB::solve method are the variable bounds and the initial guess used by the local optimizer, respectively. The return value is of a pair, whose elements are the best solution found (incumbent) and the point at which it is attained. If the problem is infeasible, the first argument turns out to be equal to mc::SBB::Options::INF (for a minimize problem and -mc::SBB::Options::INF for a maximize problem) and the second one is a 0 pointer.

The following result is obtained for this problem:

\verbatim
______________________________________________________________________

           GLOBAL NONLINEAR (NLP) OPTIMIZER IN MC++
______________________________________________________________________

SPATIAL BRANCH-AND-BOUND OPTIONS:

  ABSOLUTE CONVERGENCE TOLERANCE                            1.0e-03
  RELATIVE CONVERGENCE TOLERANCE                            1.0e-03
  MAXIMUM NUMBER OF ITERATIONS (PARTITIONS)                 NO LIMIT
  MAXIMUM CPU TIME (SEC)                                    1.0e+06
  DISPLAY LEVEL                                             2
  INTERNAL VALUE FOR INFINITY                               1.0e+20
  UPDATE ENTIRE TREE AFTER AN INCUMBENT IMPROVEMENT?        Y
  BRANCHING STRATEGY FOR NODE PARTITIONING                  OMEGA
  BRANCHING TOLERANCE FOR VARIABLE AT BOUND                 1.00e-03
  BRANCHING STRATEGY FOR VARIABLE SELECTION                 RGREL
  INITIALIZATION THRESHOLD IN RELIABILITY BRANCHING         1
  SCORE PARAMETER IN RELIABILITY BRANCHING                  0.17
  USE CONSTRAINT PROPAGATION?                               N
  USE REDUCTION CONSTRAINTS?                                N
  USE OPTIMIZATION-BASED DOMAIN REDUCTION?                  Y
  OPTIMIZATION-BASED REDUCTION MAX LOOPS                    4
  LOOP THRESHOLD FOR OPTIMIZATION-BASED REDUCTION           20%
  CUTS FROM TAYLOR MODELS                                   NONE
  ORDER OF TAYLOR MODEL                                     4
 _______________________________________________________________________

 INDEX STACK    CUMUL TIME        RELAX         INC   PARENT        UBD           LBD         ACTION  
     1     1  3.200200e-02  1.000000e+20 -1.000000e+20     0  6.666667e+00  5.000000e+00        FATHOM

#  NORMAL TERMINATION: 3.304207 CPU SEC  (LBD:64.8%  UBD:33.8%)
#  INCUMBENT VALUE:  6.666667e+00
#  INCUMBENT POINT:  6.000000e+00  6.666667e-01
#  INCUMBENT FOUND AT NODE: 1
#  TOTAL NUMBER OF NODES:   1
#  MAXIMUM NODES IN STACK:  0
\endverbatim

The spatial B&B method solves this problem in just 1 node. It turns out that the lower-bounding problem (local NLP solver) is a local solution with objective value equal to 5. However, the solution point of the upper-bounding problem (LP relaxation) turns out to be feasible with objective value equal to 6.667. In other word, we are lucky here that the solution of the relaxation at the root node is the global optimum. It follows that the root node is fathomed and the B&B algorithm terminates.

*/

#ifndef MC__NLPSBB_HPP
#define MC__NLPSBB_HPP

#undef USE_CPLEX

#include <stdexcept>
#include <cassert>
#include "sbb.hpp"
#include "mcipopt.hpp"
#ifdef USE_CPLEX
  #include "fprcplex.hpp"
#else
  #include "fprgurobi.hpp"
#endif

#undef MC__NLPSBB_DEBUG
#undef MC__NLPSBB_SHOW_FAILURE
#undef MC__NLPSBB_SHOW_REDUCTION
#undef MC__NLPSBB_SHOW_FPREL
#undef MC__NLPSBB_TRACE

/* TO DO:
- Documentation: NLPSBB options
- Extend to encompass MINLP models
- Implement Taylor model cuts
- Implement edge-concave envelopes for multilinear terms (via FPRelax)

BUGS:
*/

namespace mc
{
//! @brief C++ class for global solution of nonconvex NLPs using spatial branch-and-bound
////////////////////////////////////////////////////////////////////////
//! mc::NLPSBB is a C++ class for solving nonconvex NLPs to global
//! optimality using spatial branch-and-bound
////////////////////////////////////////////////////////////////////////
template <typename T, typename NLP>
class NLPSBB: protected SBB<T>, public virtual NLPSTRUCT
{
  // Overloading stdout operator
  template <typename U, typename DU> friend std::ostream& operator<<
    ( std::ostream&os, const NLPSBB<U,DU>& );

public:
  typedef McCormick<T> MCT;
  typedef TModel<T>    TMT;
  typedef TModel<MCT>  TMMCT;
  typedef TVar<T>      TVT;
  typedef TVar<MCT>    TVMCT;
#ifdef USE_CPLEX
  typedef FPRCplex<T>  FPRT;
#else
  typedef FPRGurobi<T> FPRT;
#endif
  typedef FPVar<T>     FPVT;
  typedef FPOp<T>      FPCT;
  typedef TModel<FPVT> TMFPVT;

  /** @defgroup NLPSBB Nonconvex NLP using Spatial Branch-and-Bound
   *  @{
   */
  //! @brief Constructor
  NLPSBB();

  //! @brief Destructor
  virtual ~NLPSBB();

  //! @brief NLPSBB options
  struct Options
  {
    //! @brief Constructor
    Options():
      USE_CONSTRAINT_PROPAGATION(false), USE_REDUCTION_CONSTRAINTS(false),
      USE_DOMAIN_REDUCTION(true), DOMAIN_REDUCTION_MAXLOOPS(4),
      DOMAIN_REDUCTION_THRESHOLD(0.2), TAYLOR_MODEL_CUTS(NONE),
      TAYLOR_MODEL_ORDER(4)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching strategy
    enum TMCUTS{
      NONE=0,	//!< No cuts
      TM,	//!< Taylor model cuts
      TMMC	//!< McCormick-Taylor model cuts
    };
    //! @brief Whether to pre-process nodes using constraint propagation
    bool USE_CONSTRAINT_PROPAGATION;
    //! @brief Whether to pre-process nodes using reduction constraints
    bool USE_REDUCTION_CONSTRAINTS;
    //! @brief Whether to post-process nodes using optimality-based domain reduction
    bool USE_DOMAIN_REDUCTION;
    //! @brief Maximum number of domain reduction loops
    unsigned int DOMAIN_REDUCTION_MAXLOOPS;
    //! @brief Threshold for repeating domain reduction (minimum domain reduction ratio)
    double DOMAIN_REDUCTION_THRESHOLD;
    //! @brief Bounding strategy
    TMCUTS TAYLOR_MODEL_CUTS;
    //! @brief Order of Taylor/McCormick-Taylor model for IVP bounding
    unsigned int TAYLOR_MODEL_ORDER;
  } options;

  //! @brief NLPSBB exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLPSBB exception handling
    enum TYPE{
      REDUC=0,	//!< Error during variable domain reduction
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    std::string what(){
      switch( _ierr ){
      case REDUC:
        return "Error during variable domain reduction";
      case UNDEF: default:
        return "Error due to calling a feature not yet implemented in NLPSBB";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    double cumul_NLPLOC;
    double cumul_FPREL;
    double cumul_NLPREL;
    double cumul_REDUC;
    double cumul_SBB;
  } stats;

  //! @brief Reference to SBB options
  typename SBB<T>::Options& options_SBB()
    { return SBB<T>::options; }
  //! @brief Reference to IPOPT options
  typename IPOPT<NLP>::Options& options_IPOPT()
    { return _pIPOPT->options; }
  //! @brief Reference to FPRelax options
  typename FPRelax<T>::Options& options_FPRelax()
    { return _pFPR->options; }
  //! @brief Reference to TModel options
  typename TModel<T>::Options& options_TModel()
    { return *_pTMOpt; }

  //! @brief Solve DO model -- return value is status
  std::pair<double,const double*> solve
    ( const T*P, const double*p0, std::ostream&os=std::cout );
  /** @} */

private:
  //! @brief Local NLP solver
  Ipopt::SmartPtr< IPOPT<NLP> > _pIPOPT;
  //! @brief FPRelax environment for polyhedral relaxation
#ifdef USE_CPLEX
  FPRCplex<T>* _pFPR;
#else
  FPRGurobi<T>* _pFPR;
#endif
  //! @brief FPRelax variable array for polyhedral relaxation
  FPVar<T>* _pFPV;
  //! @brief TM environment for NLP bounding
  TModel<T>* _pTM;
  //! @brief TMMC environment for NLP bounding
  TModel<MCT>* _pTMMC;
  //! @brief TMFP environment for NLP polyhedral relaxation
  TModel<FPVT>* _pTMFP;
  //! @brief Array holding variable bounds
  std::pair<double,double>* _pP;
  //! @brief TM options
  typename TModel<T>::Options* _pTMOpt;

  //! @brief Variable holding NLP problem structure
  Structure _NLPstruct;

  //! @brief Structure holding current optimization data
  struct Data{
    const double*p0;
    const T*P;
  } _data;

  //! @brief User-function to subproblems
  typename SBB<T>::STATUS subproblems
    ( const typename SBB<T>::TASK task, const unsigned int np, T*P,
      double*p, double&f, const double INC );
  //! @brief Subproblem for local optimization
  typename SBB<T>::STATUS _subproblem_local
    ( const unsigned int np, const T*P, double*p, double&f );
  //! @brief Subproblem for feasibility test
  typename SBB<T>::STATUS _subproblem_feasibility
    ( const unsigned int np, const double*p, double&f );
  //! @brief Subproblem for relaxation
  typename SBB<T>::STATUS _subproblem_relaxed
    ( const unsigned int np, T*P, double*p, double&f, const double INC );

  //! @brief Set NLP problem structure
  typename SBB<T>::TYPE _set_NLP_structure();
  //! @brief Initialize NLP relaxation
  void _init_NLP_relax();

  //! @brief Create NLP polyhedral relaxations via McCormick-Taylor models
  void _bound_NLP_TMMC
    ( const double*p, FPVT*FPVfg );
  //! @brief Create NLP polyhedral relaxations via Taylor models
  void _bound_NLP_TM
    ( const double*p, FPVT*FPVfg );

  //! @brief Private methods to block default compiler methods
  NLPSBB(const NLPSBB&);
  NLPSBB& operator=(const NLPSBB&);
};

template <typename T, typename NLP>
inline
NLPSBB<T,NLP>::NLPSBB()
: NLPSTRUCT(NLP())
{
  _pIPOPT  = new IPOPT<NLP>;
  _pFPR    = new FPRT;
  _pFPV    = new FPVT[_np];
  _pTM = 0;
  _pTMMC = 0;
  _pTMFP = 0;
  _pP = new std::pair<double,double>[_np];
  _pTMOpt = new typename TModel<T>::Options;
}

template <typename T, typename NLP>
inline
NLPSBB<T,NLP>::~NLPSBB()
{
  // do NOT erase pIPOPT (SmartPtr used)
  delete _pFPR;
  delete[] _pFPV;
  delete _pTM;
  delete _pTMMC;
  delete _pTMFP;
  delete[] _pP;
  delete _pTMOpt;
}

template <typename T, typename NLP>
inline typename SBB<T>::TYPE
NLPSBB<T,NLP>::_set_NLP_structure
()
{
  // Structure of objective function and constraints (parameter dependence)
  Structure *pS = new Structure[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    pS[ip].indep(ip);   

  std::pair<Structure,NLPSTRUCT::t_OBJ> pobj = NLP().OBJ( pS );
  typename SBB<T>::TYPE pb;
  switch( pobj.second ){
    case NLPSTRUCT::MIN: pb = SBB<T>::MINIMIZE; break;
    default:             pb = SBB<T>::MAXIMIZE; break;
  }
  _NLPstruct = pobj.first;
  for( unsigned int ic=0; ic<_nc; ic++ )
    _NLPstruct += NLP().CTR( ic, pS ).first;
#ifdef MC__NLPSBB_DEBUG
  std::cout << "NLP <- " << _NLPstruct << std::endl;
  int dum; std::cin >> dum;
#endif

  delete[] pS;

  // Variables to exclude from branching
  SBB<T>::_exclude_vars.clear();
  std::map<int,bool>& NLPdep = _NLPstruct.dep();
  for( unsigned int ip=0; ip<_np; ip++ ){
    std::map<int,bool>::iterator it = NLPdep.find(ip);
    if( it==NLPdep.end() || it->second ) SBB<T>::_exclude_vars.insert(ip);
  }

  return pb;
}

template <typename T, typename NLP>
inline void
NLPSBB<T,NLP>::_init_NLP_relax
()
{
  delete _pTM;
  delete _pTMMC;
  delete _pTMFP;

  switch( options.TAYLOR_MODEL_CUTS ){
  // Set McCormick-Taylor model environments
  case Options::TMMC:
    _pTMMC = new TMMCT( _np, options.TAYLOR_MODEL_ORDER );
    _pTMFP = new TMFPVT( _np, options.TAYLOR_MODEL_ORDER );
    _pTMFP->options = _pTMMC->options = *_pTMOpt;
    break;
  // Set Taylor model environments
  case Options::TM:
    _pTM = new TMT( _np, options.TAYLOR_MODEL_ORDER );
    _pTMFP = new TMFPVT( _np, options.TAYLOR_MODEL_ORDER );
    _pTMFP->options = _pTM->options = *_pTMOpt;
    break;
  default:
    break;
  }
}

template <typename T, typename NLP>
inline std::pair<double,const double*>
NLPSBB<T,NLP>::solve
( const T*P, const double*p0, std::ostream&os )
{
  // Keep track of bounds and initial guess
  _data.p0 = p0;
  _data.P  = P;

  // Keep track of execution times
  stats.cumul_FPREL = stats.cumul_REDUC = stats.cumul_NLPREL
  = stats.cumul_NLPLOC = 0.;

  // Determine NLP problem structure
  typename SBB<T>::TYPE pb = _set_NLP_structure();
  _init_NLP_relax();

  // Run SBB solver
  std::set<unsigned int>&exclude = SBB<T>::_exclude_vars;
  SBB<T>::variables( _np, P, &exclude );
  stats.cumul_SBB = -time();
  const std::pair<double,const double*>& opt = SBB<T>::solve( pb, p0, os );
  stats.cumul_SBB += time();
  return opt;
}

template <typename T, typename NLP>
inline typename SBB<T>::STATUS
NLPSBB<T,NLP>::_subproblem_local
( const unsigned int np, const T*P, double*p, double&f )
{
  // Set current variable range
  for( unsigned int ip=0; ip<_np; ip++ )
    _pP[ip] = std::make_pair(Op<T>::l(P[ip]),Op<T>::u(P[ip]));

  // Compute local optimum in current variable range
  Ipopt::ApplicationReturnStatus status;
  stats.cumul_NLPLOC -= time();
  status = _pIPOPT->solve( _pP, p );
  stats.cumul_NLPLOC += time();
  if( status == Ipopt::Solve_Succeeded
   || status == Ipopt::Solved_To_Acceptable_Level
   || status == Ipopt::Feasible_Point_Found ){
    // Copy solution value and solution point
    f = _pIPOPT->solution().f;
#ifdef MC__NLPSBB_DEBUG
    std::cout << "  floc = " << f << std::endl;
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = _pIPOPT->solution().p[ip];
#ifdef MC__NLPSBB_DEBUG
      std::cout << std::scientific << std::setprecision(4)
                << "  ploc(" << ip << ") = " << p[ip] << std::endl;
#endif
    }
    return SBB<T>::NORMAL;
  }
  else 
    return SBB<T>::FAILURE;
}

template <typename T, typename NLP>
inline typename SBB<T>::STATUS
NLPSBB<T,NLP>::_subproblem_feasibility
( const unsigned int np, const double*p, double&f )
{
  // Check current point feasibility
  double maxinfeas = 0.;
  std::pair<double,NLPSTRUCT::t_CTR> g;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    g = NLP().CTR( ic, p );
    switch( g.second ){
    case NLPSTRUCT::EQ:
      maxinfeas = std::max(maxinfeas,std::fabs(g.first)); break;
    case NLPSTRUCT::LE:
      maxinfeas = std::max(maxinfeas,g.first); break;
    case NLPSTRUCT::GE:
      maxinfeas = std::max(maxinfeas,-g.first); break;
    }
  }

  // Compute objective function value
  f = NLP().OBJ( p ).first;

  return( maxinfeas < _pIPOPT->options.PRIMALTOL?
          SBB<T>::NORMAL: SBB<T>::INFEASIBLE );
}

// template <typename T, typename NLP>
// inline void
// NLPSBB<T,NLP>::_bound_IVP_TM
// ( const double*p )
// {
//   // Set Taylor model for parameters
//   stats.cumul_ODEBND -= time();
//   TVT TVp[_np];
//   for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
//     if( _IVPdepend.find(ip) != _IVPdepend.end() ){
//       TVp[ip] = TVT( _pTM, isub, _pFPV[ip].I() );
//       isub++;
//     }
//     else
//       TVp[ip] = _pFPV[ip].I();
//   }
//     
//   // Set Taylor model for states at stage times
//   TVT* TVxk[_data.ns+1];
//   for( unsigned int is=0; is<=_data.ns; is++ )
//     TVxk[is] = new TVT[_nx];
//   
//   // Compute state bounds at stage times
//   if( _pODEBND->bounds( _data.ns, _data.tk, TVp, TVxk )
//       != ODEBND<T,NLP>::NORMAL ){
//     for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];
//     stats.cumul_ODEBND += time();
//     return SBB<T>::FAILURE;
//   }
// #ifdef MC__NLPSBB_DEBUG
//   std::cout << *_pODEBND;
// #endif
//   stats.cumul_ODEBND += time();
//   
//   // Define an FPVar for the state variables participating in FPRelax environment
//   stats.cumul_FPREL -= time();
//   _pFPR->TModelFP( _pTMFP, _pTM, _IVPparam );
//   std::map<int,bool>& DOdep = _DOstruct[1].dep();
//   bool failed = false;
//   for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
//     for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
//       if( Op<T>::diam(TVxk[is][ix].B()) >= SBB<T>::options.INF ){
//         failed = true; break;
//       }
//       std::map<int,bool>::iterator it = DOdep.find(ixs);
//       if( it!=DOdep.end() )
//         FPVxk[is][ix] = FPVT( _pFPR, TVxk[is][ix], _pTMFP );
//     }
//   stats.cumul_FPREL += time();
//   
//   // Clean-up
//   stats.cumul_ODEBND -= time();
//   for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];
//   stats.cumul_ODEBND += time();
//   return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
// }
// 
// template <typename T, typename NLP>
// inline typename SBB<T>::STATUS
// NLPSBB<T,NLP>::_bound_IVP_TMMC
// ( const double*p, FPVT**FPVxk )
// {
//   // Set McCormick-Taylor model for parameters
//   stats.cumul_ODEBND -= time();
//   TVMCT TVMCp[_np];
//   const unsigned int nsub = _IVPdepend.size();
//   double MCpref[nsub];
//   for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
//     if( _IVPdepend.find(ip) != _IVPdepend.end() ){
//       TVMCp[ip] = TVMCT( _pTMMC, isub, MCT(_pFPV[ip].I(),p[ip]).sub(nsub,isub) );
//       MCpref[isub++] = p[ip];
//     }
//     else
//       TVMCp[ip] = MCT(_pFPV[ip].I(),p[ip]);
//   }
//     
//   // Set McCormick-Taylor model for states at stage times
//   TVMCT* TVMCxk[_data.ns+1];
//   for( unsigned int is=0; is<=_data.ns; is++ )
//     TVMCxk[is] = new TVMCT[_nx];
//   
//   // Compute state bounds at stage times
//   if( _pODEBND->bounds( _data.ns, _data.tk, TVMCp, TVMCxk )
//       != ODEBND<T,NLP>::NORMAL ){
//     for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVMCxk[is];
//     stats.cumul_ODEBND += time();
//     return SBB<T>::FAILURE;
//   }
// #ifdef MC__NLPSBB_DEBUG
//   std::cout << *_pODEBND;
// #endif
//   stats.cumul_ODEBND += time();
// 
//   // Define an FPVar for the state variables participating in FPRelax environment
//   stats.cumul_FPREL -= time();
//   _pFPR->TModelFP( _pTMFP, _pTMMC, _IVPparam );
//   std::map<int,bool>& DOdep = _DOstruct[1].dep();
//   bool failed = false;
//   for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
//     for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
//       if( Op<T>::diam(TVMCxk[is][ix].B().I()) >= SBB<T>::options.INF ){
//         failed = true; break;
//       }
//       std::map<int,bool>::iterator it = DOdep.find(ixs);
//       if( it!=DOdep.end() )
//         FPVxk[is][ix] = FPVT( _pFPR, TVMCxk[is][ix], _pTMFP, MCpref, _IVPparam );
//     }
//   stats.cumul_FPREL += time();
//   
//   // Clean-up
//   stats.cumul_ODEBND -= time();
//   for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVMCxk[is];
//   stats.cumul_ODEBND += time();
//   return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
// }

template <typename T, typename NLP>
inline typename SBB<T>::STATUS
NLPSBB<T,NLP>::_subproblem_relaxed
( const unsigned int np, T*P, double*p, double&f, const double INC )
{
  // Main loop for relaxation and domain reduction
  for( unsigned int ired=0; !ired || ired < options.DOMAIN_REDUCTION_MAXLOOPS; ired++ ){

    // Set current variables
    stats.cumul_FPREL -= time();
    _pFPR->reset();
    for( unsigned int ip=0; ip<_np; ip++ )
      _pFPV[ip] = FPVT( _pFPR, ip, P[ip] );
#ifdef MC__NLPSBB_DEBUG
    std::cout << "FPR:\n" << *_pFPR;
#endif

    // Objective and constraint definition
    std::pair<FPVT,NLPSTRUCT::t_OBJ> pobj = NLP().OBJ( _pFPV );
    for( unsigned int ic=0; ic<_nc; ic++ ){
      std::pair<FPVT,NLPSTRUCT::t_CTR> pctr = NLP().CTR( ic, _pFPV );
      switch( pctr.second ){
        case NLPSTRUCT::EQ: _pFPR->add_constraint( pctr.first, FPRT::EQ, 0. ); break;
        case NLPSTRUCT::LE: _pFPR->add_constraint( pctr.first, FPRT::LE, 0. ); break;
        default:            _pFPR->add_constraint( pctr.first, FPRT::GE, 0. ); break;
      }
    }

//     // Additional cuts based on polyhedral relaxation of Taylor models
//     switch( options.TAYLOR_MODEL_CUTS ){
//     // State polyhedral relaxations via Taylor models
//     case Options::TM:
//       _bound_NLP_TM( p ); break;
//     // State polyhedral relaxations via McCormick-Taylor models
//     case Options::TMMC:
//       _bound_NLP_TMMC( p ); break;
//     default:
//       break;
//     }

    // Incumbent cuts -- ALWAYS since this also defines the objective function!
    //if( options.USE_CONSTRAINT_PROPAGATION || options.USE_DOMAIN_REDUCTION ){
      FPCT* pCut_inc = 0;
      switch( pobj.second ){
        case NLPSTRUCT::MIN:
          pCut_inc = _pFPR->add_constraint( pobj.first, FPRT::LE, INC ); break;
        default:
          pCut_inc = _pFPR->add_constraint( pobj.first, FPRT::GE, INC ); break;
      }
    //}
    stats.cumul_FPREL += time();

    // Node processing: Reduction constraints
    stats.cumul_REDUC -= time();
    if ( options.USE_REDUCTION_CONSTRAINTS )
      _pFPR->generate_reduction_constraints();

    // Node processing: Constraint propagation
    if ( options.USE_CONSTRAINT_PROPAGATION ){
      if ( !_pFPR->propagate_constraints() ){
        stats.cumul_REDUC += time();
        return SBB<T>::INFEASIBLE;
      }
      for( unsigned int ip=0; ip<_np; ip++ )
        if( !Op<T>::inter(  P[ip], P[ip],_pFPV[ip].I() ) )
          throw Exceptions( Exceptions::REDUC );
    }
    stats.cumul_REDUC += time();

    // Polyhedral relaxations
    stats.cumul_FPREL -= time();
    _pFPR->generate_cuts();
    stats.cumul_FPREL += time();
#ifdef MC__NLPSBB_SHOW_FPREL
    std::cout << "FPR:\n" << *_pFPR;
    { int dum; std::cin >> dum; }
#endif

    // Node processing: Domain reduction
    if( options.USE_DOMAIN_REDUCTION && ired < options.DOMAIN_REDUCTION_MAXLOOPS ){
      stats.cumul_REDUC -= time();
      double maxred = 0.;
      for( unsigned int ip=0; ip<_np; ip++ ){
        T P0 = P[ip];
        for( unsigned int ipb=0; ipb<2; ipb++ ){
          _pFPR->set_objective( (ipb? FPRT::MIN: FPRT::MAX), _pFPV[ip] );
#ifdef MC__NLPSBB_DEBUG
          switch( _pFPR->solve( 1, true ) ){
#else
          switch( _pFPR->solve() ){
#endif
#ifdef USE_CPLEX
          case IloAlgorithm::Optimal:
#else
          case GRB_OPTIMAL:
#endif
            // Update variable bounds
            switch( ipb ){
            case 0: // range upper bound
              if( !Op<T>::inter(  P[ip], P[ip],
                T(-SBB<T>::options.INF,_pFPR->get_objective()) ) ){
                P[ip] = Op<T>::l( P[ip] );
                //throw Exceptions( Exceptions::REDUC );
                //return SBB<T>::INFEASIBLE;
              }
              break;
            case 1: // range lower bound
              if( !Op<T>::inter(  P[ip], P[ip],
                  T(_pFPR->get_objective(),SBB<T>::options.INF) ) ){
                P[ip] = Op<T>::u( P[ip] );
                //throw Exceptions( Exceptions::REDUC );
                //return SBB<T>::INFEASIBLE;
              }
              break;
            }
            _pFPR->update_bounds( _pFPV[ip], P[ip] );
#ifdef MC__NLPSBB_DEBUG
            std::cout << "  P(" << ip << ") = " << P[ip] << std::endl;
#endif
            break;
#ifdef USE_CPLEX
          case IloAlgorithm::Infeasible:
          case IloAlgorithm::InfeasibleOrUnbounded:
#else
          case GRB_INFEASIBLE:
          case GRB_INF_OR_UNBD:
#endif
            stats.cumul_REDUC += time();
            return SBB<T>::INFEASIBLE;
          default:
            stats.cumul_REDUC += time();
#ifdef MC__NLPSBB_SHOW_FAILURE
            std::cout << "FAILURE IN DOMAIN REDUCTION LP\n"; 
#endif
            return SBB<T>::FAILURE;
          }
        }
        maxred = std::max( 1.-Op<T>::diam(P[ip])/Op<T>::diam(P0), maxred );
      }
      stats.cumul_REDUC += time();

      // Reconstruct lower bounding problem for tightened bounds
#ifdef MC__NLPSBB_SHOW_REDUCTION
      std::cerr << "max. domain reduction: " << maxred << std::endl;
#endif
      if( maxred > options.DOMAIN_REDUCTION_THRESHOLD
       && ired+1 < options.DOMAIN_REDUCTION_MAXLOOPS ) continue;
    }

    // Compute lower bound
    stats.cumul_NLPREL -= time();
    if( pCut_inc ) _pFPR->remove_constraint( pCut_inc );
    switch( pobj.second ){
      case NLPSTRUCT::MIN: _pFPR->set_objective( FPRT::MIN, pobj.first ); break;
      default:		   _pFPR->set_objective( FPRT::MAX, pobj.first ); break;
    }
#ifdef MC__NLPSBB_DEBUG
    switch( _pFPR->solve( 1, true ) ){
#else
    switch( _pFPR->solve() ){
#endif
#ifdef USE_CPLEX
    case IloAlgorithm::Optimal:
#else
    case GRB_OPTIMAL:
#endif
      // Copy solution value and solution point
      f = _pFPR->get_objective();
#ifdef MC__NLPSBB_DEBUG
      std::cout << "  frel = " << f << std::endl;
#endif
      for( unsigned int ip=0; ip<_np; ip++ ){
        p[ip] = _pFPR->get_variable(ip);
#ifdef MC__NLPSBB_DEBUG
        std::cout << std::scientific << std::setprecision(4)
                  << "  prel(" << ip << ") = " << p[ip] << std::endl;
#endif
      }
      stats.cumul_NLPREL += time();
      return SBB<T>::NORMAL;
#ifdef USE_CPLEX
    case IloAlgorithm::Infeasible:
    case IloAlgorithm::InfeasibleOrUnbounded:
#else
    case GRB_INFEASIBLE:
    case GRB_INF_OR_UNBD:
#endif
      stats.cumul_NLPREL += time();
      return SBB<T>::INFEASIBLE;
    default:
#ifdef MC__NLPSBB_SHOW_FAILURE
      std::cout << "FAILURE IN LOWER BOUNDING LP\n"; 
#endif
      stats.cumul_NLPREL += time();
      return SBB<T>::FAILURE;
    }
  }
#ifdef MC__NLPSBB_SHOW_FAILURE
  std::cout << "FAILURE AT EXIT\n"; 
#endif
  stats.cumul_REDUC += time();
  return SBB<T>::FAILURE;
}

template <typename T, typename NLP>
inline typename SBB<T>::STATUS
NLPSBB<T,NLP>::subproblems
( const typename SBB<T>::TASK task, const unsigned int np, T*P,
  double*p, double&f, const double INC )
{
  switch( task ){

  // UPPER BOUND
  case SBB<T>::UPPERBD:
    switch( SBB<T>::_pb ){
    case SBB<T>::MINIMIZE:
      return _subproblem_local( np, P, p, f );
    default:
      return _subproblem_relaxed( np, P, p, f, INC );
    }

  // LOWER BOUND
  case SBB<T>::LOWERBD:
    switch( SBB<T>::_pb ){
    case SBB<T>::MAXIMIZE:
      return _subproblem_local( np, P, p, f );
    default:
      return _subproblem_relaxed( np, P, p, f, INC );
    }

  // FEASIBILITY TEST
  case SBB<T>::FEASTEST:
    return _subproblem_feasibility( np, p, f );

  // PRE/POST-PROCESSING
  case SBB<T>::PREPROC:
  case SBB<T>::POSTPROC:
    return SBB<T>::NORMAL;

  default:
    return SBB<T>::FAILURE;
  }
}

template <typename T, typename NLP>
inline void
NLPSBB<T,NLP>::Options::display
( std::ostream&out ) const
{
  // Display NLPSBB Options
  out << std::setw(60) << "  USE CONSTRAINT PROPAGATION?"
      << (USE_CONSTRAINT_PROPAGATION?"Y\n":"N\n");
  out << std::setw(60) << "  USE REDUCTION CONSTRAINTS?"
      << (USE_REDUCTION_CONSTRAINTS?"Y\n":"N\n");
  out << std::setw(60) << "  USE OPTIMIZATION-BASED DOMAIN REDUCTION?"
      << (USE_DOMAIN_REDUCTION?"Y\n":"N\n");
  out << std::setw(60) << "  OPTIMIZATION-BASED REDUCTION MAX LOOPS"
      << DOMAIN_REDUCTION_MAXLOOPS << std::endl;
  out << std::setw(60) << "  LOOP THRESHOLD FOR OPTIMIZATION-BASED REDUCTION"
      << std::fixed << std::setprecision(0)
      << DOMAIN_REDUCTION_THRESHOLD*1e2 << "%\n";
  out << std::setw(60) << "  CUTS FROM TAYLOR MODELS";
  switch( TAYLOR_MODEL_CUTS ){
  case NONE: out << "NONE\n"; break;
  case TM:   out << "TM\n";   break;
  case TMMC: out << "TMMC\n"; break;
  }
  out << std::setw(60) << "  ORDER OF TAYLOR MODEL"
      << TAYLOR_MODEL_ORDER << std::endl;
}

template <typename T, typename NLP>
inline std::ostream&
operator <<
( std::ostream&out, const NLPSBB<T,NLP>&pDO )
{
  out << std::endl
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ')
      << std::setw(52) << "GLOBAL NONLINEAR (NLP) OPTIMIZER IN MC++\n"
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ');

  // Display NLPSBB Options
  out << std::left << "SPATIAL BRANCH-AND-BOUND OPTIONS:\n\n";
  pDO.SBB<T>::options.display( out );
  pDO.NLPSBB<T,NLP>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::setfill(' ')
      << std::endl << std::endl;
  return out;
}

} // end namescape mc

#endif
