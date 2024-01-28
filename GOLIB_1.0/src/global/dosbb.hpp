// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_DOSBB Deterministic Global Optimization of Parametric IVPs in ODEs using Spatial Branch-and-Bound
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

Consider a dynamic optimization (DO) problem of the form
\f{align*}
  \min_{{\bf p}\in P}\ & \phi_0({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_0({\bf p}, {\bf x}(t))\, dt\\
  \text{s.t.} \ & \phi_k({\bf p}, {\bf x}(t_0),\ldots,{\bf x}(t_N)) + \int_{t_0}^{t_N} \psi_k({\bf p}, {\bf x}(t))\, dt \{\geq,=,\leq\} 0, \quad k=1,\ldots,n_c\\
  & \dot{{\bf x}}(t)={\bf f}({\bf p},{\bf x}(t)),\ t\in(t_0,t_{N}]; \quad {\bf x}(t_0)={\bf h}({\bf p})
\f}
where \f${\bf p}\in P\subset\mathbb{R}^{n_p}\f$, \f${\bf x}\in\mathbb{R}^{n_x}\f$, and the real-valued functions \f${\bf f}\f$, \f${\bf h}\f$, \f$\phi_k\f$ and \f$\psi_k\f$ are twice continuously differentiable in all their arguments and factorable. The class mc::DOSBB in MC++ uses branch-and-bound search to solve such DO problems to global optimality. A relaxation of the DO problem is constructed using the class mc::ODEBND that computes bounds on the solutions of parametric IVPs, and local solution are computed based on the class mc::DOSEQ which implements a sequential approach using numerical integration and local NLP search.

mc::DOSBB is templated in two variable types T and DO:
- T is the class that performs the underlying interval arithmetics. Both verified and non-verified types are supported. For verified computations, the libraries <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> and <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> are supported; the non-verified interval type mc::Interval of MC++ can also be used. More generally, any interval type can be used provided that the templated structure mc::Op is instantiated accordingly -- see the header files <tt>mcprofil.h</tt> and <tt>mcfilib.h</tt> for examples.
- DO is the class that defines the DO problem to be solved, which has to be a class derived from mc::DOSTRUCT. For example, suppose we want to solve the following dynamic optimization problem:
\f{align*}
  \max_{T,x_1^0,u_1,\ldots,u_N}\ & T \\
  \text{s.t.} \ & x_1(1) = 1\\
  & x_2(1) = 0\\
  & \left.\begin{array}{l}
    \dot{x_1}(t) = x_2(t)\\
    \dot{x_2}(t) = (u_k-x_1(t)-2*x_2(t))\theta(t)\\
    \dot{\theta}(t) = T
    \end{array}\right\},\ t\in({\textstyle\frac{k-1}{N},\frac{k}{N}}],\ k=1,\ldots,N\\
  & x_1(0) = x_1^0,\ x_2(0) = 0,\ \theta(0) = 0
\f}
with \f$N=5\f$ stages. The following class is defined:
\code
      #include "dosbb.h"

      // Number of time stages
      const unsigned int NS = 5;
      // Number of decision variables
      const unsigned int NP = 2+NS;
      // Number of state variables
      const unsigned int NX = 3;
      // Number of constraints
      const unsigned int NC = 2;

      class DO : public virtual mc::DOSTRUCT
      {
      public:
        DO()
          : mc::DOSTRUCT( NP, NX, NC )
          {}

        // ODEs right-hand side
	template <typename TX, typename TP>
        TX RHS
          ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
          {
            assert( ix < nx() );
            switch( ix ){
              case 0: return p[0]*x[1];
              case 1: return p[0]*(p[2+is]-x[0]-2*x[1])*x[2];
              case 2: return p[0];
              default: throw std::runtime_error("invalid size");
            }
          }
      
        // Initial conditions
        template <typename T>
        T IC
          ( const unsigned int ix, const T*p )
          {
            assert( ix < nx() );
            switch( ix ){
              case 0: return p[1];
              case 1: return 0.;
              case 2: return 0.;
              default: throw std::runtime_error("invalid size");
            }
          }

        // Objective functional
        template <typename T>
        std::pair<T,t_OBJ> OBJ
          ( const T*p, T* const*x, const unsigned int ns )
          {
            return std::make_pair( p[0], MIN );
          }

        // Constraint functionals
        template <typename T>
        std::pair<T,t_CTR> CTR
          ( const unsigned int ic, const T*p, T* const*x, const unsigned int ns )
          {
            switch( ic ){
              case 0: return std::make_pair( x[ns][0]-1., EQ );
              case 1: return std::make_pair( x[ns][1]-0., EQ );
              default: throw std::runtime_error("invalid size");
            }
          }
      };
\endcode

\section sec_DOSBB_solve How to Solve a DO Model using mc::DOSBB?

Start by defining the parameter set \f$P\f$ and initial guess \f$p_0\f$:

\code
  double p0[NP], tk[NS+1];
  std::pair<double,double> P[NP];
  p0[0] = 1.;   P[0] = std::make_pair( 0.1, 2. );
  p0[1] = 0.5;  P[1] = std::make_pair( -1., 1. );
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[2+is]  = std::make_pair( -0.5, 1.5 );
    p0[2+is] = 0.;
    tk[1+is] = tk[is]+1./(double)NS;
  }
\endcode

Then, define an mc::DOSBB object:

\code
  mc::DOSBB<T,DO> pDO = new mc::DOSBB<T,DO>;
\endcode

The DO model is optimized as follows:

\code
  std::pair<double,const double*> optimum = pDO->solve( NS, tk, P, p0 );
\endcode

where the first and second arguments are the number and value of the time stages, the third one is the variable bounds, and the last one is the initial guess used by the optimizer (IPOPT). The return value is of a pair, whereby the first element is the best solution found (incumbent) and the point at which it is attained. If the problem is infeasible, the first argument turns out to be equal to mc::SBB::Options::INF and the second one is a 0 pointer.

The following result is obtained for this problem:

\verbatim
 INDEX STACK   CUMUL TIME      RELAX           INC     PARENT      LBD           UBD       ACTION  
     1     1  0.000000e+00 -1.000000e+20  1.000000e+20     0  1.000000e-01  1.563176e+00   BRANCH0
     2     2  8.480530e-01  1.000000e-01  1.563176e+00     1  1.266758e+00  1.563176e+00   BRANCH2
     3     3  5.204326e+00  1.000000e-01  1.563176e+00     1  1.563176e+00       SKIPPED    FATHOM
     4     2  7.656479e+00  1.266758e+00  1.563176e+00     2    INFEASIBLE       SKIPPED    FATHOM
     5     1  9.016564e+00  1.266758e+00  1.563176e+00     2  1.563175e+00       SKIPPED    FATHOM

#  NORMAL TERMINATION: 11.512720 CPU SEC  (LBD:99.1%  UBD:0.9%)
#  INCUMBENT VALUE:  1.563176e+00
#  INCUMBENT POINT:  1.563176e+00  1.000000e+00  1.500000e+00  1.500000e+00
#                    1.500000e+00  2.808999e-02 -5.000000e-01
#  INCUMBENT FOUND AT NODE: 1
#  TOTAL NUMBER OF NODES:   5
#  MAXIMUM NODES IN STACK:  3
\endverbatim

Finally, computational statistics for the SBB algorithm can be obtained from the <a>stat</a> public member of mc::DOSBB, a structure with the following fields:
- <a>stats.cumul_SBB</a>    total time taken by the spatial B&B algorithm
- <a>stats.cumul_DOLOC</a>  cumulated time spent for solving the local DO subproblems
- <a>stats.cumul_ODEBND</a> cumulated time spent for bounding the ODE solutions
- <a>stats.cumul_FPREL</a>  cumulated time spent for generating the LP relaxations of the DO subproblems
- <a>stats.cumul_REDUC</a>  cumulated time spent for tightening variable bounds (domain reduction)
- <a>stats.cumul_DOLOC</a>  cumulated time spent for solving the relaxed DO subproblems
\code
    std::cout << std::fixed << std::setprecision(1)
      << "\n  Total Time: " << GDO.stats.cumul_SBB
      << "\n  DOLOC: "  << GDO.stats.cumul_DOLOC/GDO.stats.cumul_SBB*1e2  << "%"
      << "    ODEBND: " << GDO.stats.cumul_ODEBND/GDO.stats.cumul_SBB*1e2 << "%"
      << "    LPREL: "  << GDO.stats.cumul_FPREL/GDO.stats.cumul_SBB*1e2  << "%"
      << "    REDUC: "  << GDO.stats.cumul_REDUC/GDO.stats.cumul_SBB*1e2  << "%"
      << "    DOREL: "  << GDO.stats.cumul_DOREL/GDO.stats.cumul_SBB*1e2  << "%"
      << std::endl;
\endcode

The following output is obtained for the illustrative example, indicating that about 85% of the computational time is used for tightening the varaible bounds and for bounding the parametric ODE solutions:

\verbatim
  DOLOC: 0.9%    ODEBND: 31.5%    FPREL: 13.0%    REDUC: 53.3%    DOREL: 1.3%
\endverbatim

*/

#ifndef MC__DOSBB_HPP
#define MC__DOSBB_HPP

#define USE_CPLEX

#include <stdexcept>
#include <cassert>
#include "sbb.hpp"
#include "doseq.hpp"
//#include "odebnd.h"
#include "odebnd_val.hpp"
#include "odebnd_gsl.hpp"
#ifdef USE_CPLEX
  #include "fprcplex.hpp"
#else
  #include "fprgurobi.hpp"
#endif

#undef MC__DOSBB_DEBUG
#undef MC__DOSBB_SHOW_FAILURE
#undef MC__DOSBB_SHOW_REDUCTION
#undef MC__DOSBB_TRACE

/* TO DO:
- Describe options in documentation

BUGS:
- [OK] some issues with constraint propagation that becomes infeasible... this is caused by rounding when using a nonvalidated interval type such as mc::Interval.
- [FIXED] HYDRID option as range bounder for Taylor model gives erroneous results...
- [FIXED] repeated solution with CPLEX becomes slower and slower... sth not properly reinitialized!
*/

namespace mc
{
//! @brief C++ class for global solution of dynamic optimization using spatial branch-and-bound
////////////////////////////////////////////////////////////////////////
//! mc::DOSBB is a C++ class for solving dynamic optimization problems
//! to global optimality using spatial branch-and-bound
////////////////////////////////////////////////////////////////////////
template <typename T, typename DO>
class DOSBB: protected SBB<T>, public DOSTRUCT
{
  // Overloading stdout operator
  template <typename U, typename DU> friend std::ostream& operator<<
    ( std::ostream&os, const DOSBB<U,DU>& );

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
  typedef TModel<FPVT> TMFPVT;

  /** @defgroup DOSBB Global Dynamic Optimization using Spatial Branch-and-Bound
   *  @{
   */
  //! @brief Constructor
  DOSBB();

  //! @brief Destructor
  virtual ~DOSBB();

  //! @brief DOSBB options
  struct Options
  {
    //! @brief Constructor
    Options():
      IVP_BOUNDING_TYPE(TM), IVP_BOUNDING_PROPAGATION(VERIF),
      TAYLOR_MODEL_ORDER(4), USE_CONSTRAINT_PROPAGATION(false),
      USE_DOMAIN_REDUCTION(true), USE_RLT_TMODEL(true),
      DOMAIN_REDUCTION_MAXLOOPS(4), DOMAIN_REDUCTION_THRESHOLD(0.2)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief Branching strategy
    enum IVPBOUND{
      IA=0,	//!< Interval analysis
      MC,	//!< McCormick relaxation
      TM,	//!< Taylor model
      TMMC	//!< McCormick-Taylor model
    };
    //! @brief Branching strategy
    enum IVPMETH{
      DINEQ=0,	//!< Differential methods not accouting for truncation errors
      VERIF	//!< Verified integration accounting for truncation errors
    };
    //! @brief IVP bounder type
    IVPBOUND IVP_BOUNDING_TYPE;
    //! @brief IVP bounding method
    IVPMETH IVP_BOUNDING_PROPAGATION;
    //! @brief Order of Taylor/McCormick-Taylor model for IVP bounding
    unsigned int TAYLOR_MODEL_ORDER;
    //! @brief Whether to pre-process nodes using constraint propagation
    bool USE_CONSTRAINT_PROPAGATION;
    //! @brief Whether to post-process nodes using optimality-based domain reduction
    bool USE_DOMAIN_REDUCTION;
     //! @brief Whether to tighten Taylor model relaxations using RLT constraints
    bool USE_RLT_TMODEL;
   //! @brief Maximum number of domain reduction loops
    unsigned int DOMAIN_REDUCTION_MAXLOOPS;
    //! @brief Threshold for repeating domain reduction (minimum domain reduction ratio)
    double DOMAIN_REDUCTION_THRESHOLD;
  } options;

  //! @brief DOSBB exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for DOSBB exception handling
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
        return "Error due to calling a feature not yet implemented in DOSBB";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    double cumul_DOLOC;
    double cumul_ODEBND;
    double cumul_FPREL;
    double cumul_DOREL;
    double cumul_REDUC;
    double cumul_SBB;
  } stats;

  //! @brief Reference to SBB options
  typename SBB<T>::Options& options_SBB()
    { return SBB<T>::options; }
  //! @brief Reference to DOSEQ options
  typename DOSEQ<DO>::Options& options_DOSEQ()
    { return _pDOSEQ->options; }
  //! @brief Reference to CVODES options
  typename CVODES<DO>::Options& options_CVODES()
    { return _pDOSEQ->CVODES<DO>::options; }
  //! @brief Reference to ODEBND options
  typename ODEBND<T,DO>::Options& options_ODEBND()
    { return _pODEBND->options; }
  //! @brief Reference to ODEBND_GSL options
  typename ODEBND_GSL<T,DO>::Options& options_ODEBND_GSL()
    { return _pODEDINEQ->options; }
  //! @brief Reference to FPRelax options
  typename FPRelax<T>::Options& options_FPRelax()
    { return _pFPR->options; }
  //! @brief Reference to Taylor model options
  typename TModel<T>::Options& options_TModel()
    { return *_pTMOpt; }

  //! @brief Solve DO model -- return value is status
  std::pair<double,const double*> solve
    ( const unsigned int ns, const double*tk, const T*P, const double*p0,
      std::ostream&os=std::cout );
  /** @} */

private:
  //! @brief ODE Bounder
  ODEBND<T,DO>* _pODEBND;
  //! @brief ODE Bounder
  ODEBND_GSL<T,DO>* _pODEDINEQ;
  //! @brief Local DO solver
  Ipopt::SmartPtr< DOSEQ<DO> > _pDOSEQ;
  //! @brief FPRelax environment for polyhedral relaxation
#ifdef USE_CPLEX
  FPRCplex<T>* _pFPR;
#else
  FPRGurobi<T>* _pFPR;
#endif
  //! @brief FPRelax variable array for polyhedral relaxation
  FPVar<T>* _pFPV;
  //! @brief TM environment for IVP bounding
  TModel<T>* _pTM;
  //! @brief TMMC environment for IVP bounding
  TModel<MCT>* _pTMMC;
  //! @brief TMFP environment for IVP polyhedral relaxation
  TModel<FPVT>* _pTMFP;
  //! @brief Array holding variable bounds
  std::pair<double,double>* _pP;
  //! @brief TM options
  typename TModel<T>::Options* _pTMOpt;

  //! @brief Variable holding DO problem structure
  Structure _DOstruct[2];
  //! @brief Set of variables participating in IVP
  std::set<unsigned int> _IVPdepend;
  //! @brief Array holding IVP parameter indices
  unsigned int *_IVPparam;

  //! @brief Structure holding current optimization data
  struct Data{
    const double*p0;
    const T*P;
    unsigned int ns;
    const double *tk;  
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

  //! @brief Set DO problem structure
  typename SBB<T>::TYPE _set_DO_structure();
  //! @brief Initialize IVP relaxation
  void _init_IVP_relax();

  //! @brief Create state polyhedral relaxations via McCormick relaxations
  typename SBB<T>::STATUS _bound_IVP_MC
    ( const double*p, FPVT**FPVxk );
  //! @brief Create state polyhedral relaxations via Intervals
  typename SBB<T>::STATUS _bound_IVP_IA
    ( const double*p, FPVT**FPVxk, typename Options::IVPMETH meth );
  //! @brief Create state polyhedral relaxations via McCormick-Taylor models
  typename SBB<T>::STATUS _bound_IVP_TMMC
    ( const double*p, FPVT**FPVxk );
  //! @brief Create state polyhedral relaxations via Taylor models
  typename SBB<T>::STATUS _bound_IVP_TM
    ( const double*p, FPVT**FPVxk, typename Options::IVPMETH meth );
//   //! @brief Create state polyhedral relaxations via Taylor models using GSL
//   typename SBB<T>::STATUS _bound_IVP_TM_dineq
//     ( const double*p, FPVT**FPVxk );
  //! @brief Create cost and constraint polyhedral relaxations
  typename SBB<T>::STATUS _relax_cost_constraints
    ( FPVT**FPVxk, const double INC, std::pair<FPVT,DOSTRUCT::t_OBJ>&DOobj );
  //! @brief Perform optimization-based domain reduction
  typename SBB<T>::STATUS _domain_reduction
    ( FPVT**FPVxk, T*P, double &maxred );
  //! @brief Compute relaxation bound
  typename SBB<T>::STATUS _relaxation_solution
    ( const std::pair<FPVT,DOSTRUCT::t_OBJ>&DOobj, double*p, double&f );

  //! @brief Private methods to block default compiler methods
  DOSBB(const DOSBB&);
  DOSBB& operator=(const DOSBB&);
};

template <typename T, typename DO>
inline
DOSBB<T,DO>::DOSBB()
: DOSTRUCT(DO())
{
  _pDOSEQ  = new DOSEQ<DO>;
  _pODEBND = new ODEBND<T,DO>;
  _pODEDINEQ  = new ODEBND_GSL<T,DO>;
  _pFPR    = new FPRT;
  _pFPV    = new FPVT[_np];
  _pTM = 0;
  _pTMMC = 0;
  _pTMFP = 0;
  _IVPparam = 0;
  _pP = new std::pair<double,double>[_np];
  _pTMOpt = new typename TModel<T>::Options;
}

template <typename T, typename DO>
inline
DOSBB<T,DO>::~DOSBB()
{
  delete _pODEBND;
  delete _pODEDINEQ;
  //delete _pDOSEQ;
  delete _pFPR;
  delete[] _pFPV;
  delete _pTM;
  delete _pTMMC;
  delete _pTMFP;
  delete[] _IVPparam;
  delete[] _pP;
  delete _pTMOpt;
}

template <typename T, typename DO>
inline typename SBB<T>::TYPE
DOSBB<T,DO>::_set_DO_structure
()
{
  Structure *pS = new Structure[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    pS[ip].indep(ip);   
  Structure **xkS = new Structure*[_data.ns+1];

  // Structure of initial conditions in IVP
  xkS[0] = new Structure[_nx];
  for( unsigned int ix=0; ix<_nx; ix++ ){
    xkS[0][ix] = DO().IC( ix, pS );
#ifdef MC__DOSBB_DEBUG
    std::cout << "xkS[0][" << ix << "] <- " << xkS[0][ix]
              << std::endl;
#endif
  }

  // Structure of right-hand side in each stage in IVP
  for( unsigned int is=0; is<_data.ns; is++ ){
    xkS[is+1] = new Structure[_nx];
    for( unsigned int ix=0; ix<_nx; ix++ )
      xkS[is+1][ix] = DO().RHS( ix, pS, xkS[is], is ) + xkS[is][ix];
    bool iterate = true;
    while( iterate ){
      iterate = false;
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Structure xiSnew = DO().RHS( ix, pS, xkS[is+1], is ) + xkS[is+1][ix];
        if( xkS[is+1][ix] != xiSnew ){
	  xkS[is+1][ix] = xiSnew; iterate = true;
	}
      }
    }
#ifdef MC__DOSBB_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xkS[" << is+1 << "][" << ix << "] <- " << xkS[is+1][ix]
                << std::endl;
#endif
  }

  // Variable dependence in IVP
  _IVPdepend.clear();
  for( unsigned int ip=0; ip<_np; ip++ ){
    bool keepgoing = true;
    for( unsigned int ix=0; keepgoing && ix<_nx; ix++ ){
      std::map<int,bool>& IVPdep = xkS[_data.ns][ix].dep();
      std::map<int,bool>::iterator it = IVPdep.find(ip);
      if( it != IVPdep.end() ){
        _IVPdepend.insert(ip);
        keepgoing = false;
      }
    }
  }

  // Structure of objective function and constraints (parameter dependence)
  std::pair<Structure,DOSTRUCT::t_OBJ> pobj = DO().OBJ( pS, xkS, _data.ns );
  typename SBB<T>::TYPE pb;
  switch( pobj.second ){
    case DOSTRUCT::MIN: pb = SBB<T>::MINIMIZE; break;
    default:            pb = SBB<T>::MAXIMIZE; break;
  }
  _DOstruct[0] = pobj.first;
  for( unsigned int ic=0; ic<_nc; ic++ )
    _DOstruct[0] += DO().CTR( ic, pS, xkS, _data.ns ).first;
#ifdef MC__DOSBB_DEBUG
  std::cout << "DO <- " << _DOstruct[0] << std::endl;
  int dum; std::cin >> dum;
#endif

  // Variables to exclude from branching
  SBB<T>::_exclude_vars.clear();
  std::map<int,bool>& DOdep = _DOstruct[0].dep();
  for( unsigned int ip=0; ip<_np; ip++ ){
    std::map<int,bool>::iterator it = DOdep.find(ip);
    if( it==DOdep.end() || it->second ) SBB<T>::_exclude_vars.insert(ip);
  }

  // Structure of objective function and constraints (parameter & state dependence)
  for( unsigned int is=0, ixs=_np; is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++)
      xkS[is][ix].indep(ixs);
  _DOstruct[1] = DO().OBJ( pS, xkS, _data.ns ).first;
  for( unsigned int ic=0; ic<_nc; ic++ )
    _DOstruct[1] += DO().CTR( ic, pS, xkS, _data.ns ).first;
#ifdef MC__DOSBB_DEBUG
  std::cout << "DO <- " << _DOstruct[1] << std::endl;
  std::cin >> dum;
#endif

  // clean up
  delete[] pS;
  for( unsigned int is=0; is<=_data.ns; is++ )
    delete[] xkS[is];
  delete[] xkS;

  return pb;
}

template <typename T, typename DO>
inline void
DOSBB<T,DO>::_init_IVP_relax
()
{
  delete _pTM;   _pTM   = 0;
  delete _pTMMC; _pTMMC = 0;
  delete _pTMFP; _pTMFP = 0;

  // Match IVP parameters to DO problem
  delete[] _IVPparam;
  _IVPparam = new unsigned int[_IVPdepend.size()];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ )
    if( _IVPdepend.find(ip) != _IVPdepend.end() ) _IVPparam[isub++] = ip;

  switch( options.IVP_BOUNDING_TYPE ){
  // Set McCormick-Taylor model environments
  case Options::TMMC:
    _pTMMC = new TMMCT( _IVPdepend.size(), options.TAYLOR_MODEL_ORDER );
    _pTMFP = new TMFPVT( _IVPdepend.size(), options.TAYLOR_MODEL_ORDER );
    _pTMFP->options = _pTMMC->options = *_pTMOpt;
    break;
  // Set Taylor model environments
  case Options::TM:
    _pTM = new TMT( _IVPdepend.size(), options.TAYLOR_MODEL_ORDER );
    _pTMFP = new TMFPVT( _IVPdepend.size(), options.TAYLOR_MODEL_ORDER );
    _pTMFP->options = _pTM->options = *_pTMOpt;
    break;
  default:
    break;
  }
}

template <typename T, typename DO>
inline std::pair<double,const double*>
DOSBB<T,DO>::solve
( const unsigned int ns, const double*tk, const T*P, const double*p0,
  std::ostream&os )
{
  // Keep track of bounds, initial guess and time stages
  _data.p0 = p0;
  _data.P  = P;
  _data.tk = tk;
  _data.ns = ns;

  // Keep track of execution times
  stats.cumul_ODEBND = stats.cumul_FPREL = stats.cumul_REDUC = stats.cumul_DOREL
    = stats.cumul_DOLOC = 0.;

  // Determine DO problem structure
  typename SBB<T>::TYPE pb = _set_DO_structure();
  _init_IVP_relax();

  // Run SBB solver
  std::set<unsigned int>&exclude = SBB<T>::_exclude_vars;
  SBB<T>::variables( _np, P, &exclude );
  stats.cumul_SBB = -time();
  const std::pair<double,const double*>& opt = SBB<T>::solve( pb, p0, os );
  stats.cumul_SBB += time();
  return opt;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_subproblem_local
( const unsigned int np, const T*P, double*p, double&f )
{
  // Set current variable range
  for( unsigned int ip=0; ip<_np; ip++ )
    _pP[ip] = std::make_pair(Op<T>::l(P[ip]),Op<T>::u(P[ip]));

  // Compute local optimum in current variable range
  Ipopt::ApplicationReturnStatus status;
  stats.cumul_DOLOC -= time();
  status = _pDOSEQ->solve( _data.ns, _data.tk, _pP, p );
  stats.cumul_DOLOC += time();
  if( status == Ipopt::Solve_Succeeded
   || status == Ipopt::Solved_To_Acceptable_Level
   || status == Ipopt::Feasible_Point_Found ){
    // Copy solution value and solution point
    f = _pDOSEQ->solution().f;
#ifdef MC__DOSBB_DEBUG
    std::cout << std::scientific << std::setprecision(5);
    std::cout << "  floc = " << f << std::endl;
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = _pDOSEQ->solution().p[ip];
#ifdef MC__DOSBB_DEBUG
      std::cout << "  ploc(" << ip << ") = " << p[ip] << std::endl;
#endif
    }
    return SBB<T>::NORMAL;
  }
  else 
    return SBB<T>::FAILURE;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_subproblem_feasibility
( const unsigned int np, const double*p, double&f )
{
  // Check current point feasibility
  double maxinfeas;
  switch( _pDOSEQ->feasibility_test( _data.ns, _data.tk, p, f, maxinfeas ) ){
  case CVODES<DO>::NORMAL:
    return( maxinfeas < _pDOSEQ->options.PRIMALTOL?
            SBB<T>::NORMAL: SBB<T>::INFEASIBLE );
  default:
    return SBB<T>::FAILURE;
  }
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_bound_IVP_IA
( const double*p, FPVT**FPVxk, typename Options::IVPMETH meth )
{
  // Set interval bounds for parameters
  stats.cumul_ODEBND -= time();
  T Ip[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    Ip[ip] = _pFPV[ip].I();
    
  // Set interval bounds for states at stage times
  T* Ixk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    Ixk[is] = new T[_nx];
  
  // Compute state bounds at stage times
  bool failed = true;
  switch( meth ){
  case Options::VERIF:
    if( _pODEBND->bounds( _data.ns, _data.tk, Ip, Ixk )
        == ODEBND<T,DO>::NORMAL ) failed = false; break;
  case Options::DINEQ: default:
    if( _pODEDINEQ->bounds( _data.ns, _data.tk, Ip, Ixk )
        == ODEBND_GSL<T,DO>::NORMAL ) failed = false; break;
  }
  if( failed ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] Ixk[is];
    stats.cumul_ODEBND += time();
    return SBB<T>::FAILURE;
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << *_pODEBND;
#endif
  stats.cumul_ODEBND += time();
  
  // Define an FPVar for the state variables participating in FPRelax environment
  stats.cumul_FPREL -= time();
  std::map<int,bool>& DOdep = _DOstruct[1].dep();
  for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
      if( Op<T>::diam(Ixk[is][ix]) >= SBB<T>::options.INF ){
        failed = true; break;
      }
      std::map<int,bool>::iterator it = DOdep.find(ixs);
      if( it!=DOdep.end() )
        FPVxk[is][ix] = Ixk[is][ix];
    }
  stats.cumul_FPREL += time();
  
  // Clean-up
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] Ixk[is];

  return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_bound_IVP_MC
( const double*p, FPVT**FPVxk )
{
  // Set McCormick relaxations for parameters
  stats.cumul_ODEBND -= time();
  MCT MCp[_np];
  const unsigned int nsub = _IVPdepend.size();
  double MCpref[nsub];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
    if( _IVPdepend.find(ip) != _IVPdepend.end() ){
      MCp[ip] = MCT(_pFPV[ip].I(),p[ip]).sub(nsub,isub);
      MCpref[isub++] = p[ip];
    }
    else
      MCp[ip] = MCT(_pFPV[ip].I(),p[ip]);
  }
    
  // Set McCormick relaxations for states at stage times
  MCT* MCxk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    MCxk[is] = new MCT[_nx];

  // Compute state bounds at stage times
  if( _pODEBND->bounds( _data.ns, _data.tk, MCp, MCxk )
      != ODEBND<T,DO>::NORMAL ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] MCxk[is];
    stats.cumul_ODEBND += time();
    return SBB<T>::FAILURE;
  }
#ifdef MC__DOSBB_DEBUG
  std::cout << *_pODEBND;
#endif
  stats.cumul_ODEBND += time();

  // Define an FPVar for the state variables participating in FPRelax environment
  stats.cumul_FPREL -= time();
  std::map<int,bool>& DOdep = _DOstruct[1].dep();
  bool failed = false;
  for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
      if( Op<T>::diam(MCxk[is][ix].I()) >= SBB<T>::options.INF ){
        failed = true; break;
      }
      std::map<int,bool>::iterator it = DOdep.find(ixs);
      if( it!=DOdep.end() )
        FPVxk[is][ix] = FPVT( _pFPR, MCxk[is][ix], MCpref, _IVPparam );
    }
  stats.cumul_FPREL += time();
  
  // Clean-up
  stats.cumul_ODEBND -= time();
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] MCxk[is];
  stats.cumul_ODEBND += time();
  return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_bound_IVP_TM
( const double*p, FPVT**FPVxk, typename Options::IVPMETH meth )
{
  // Set Taylor model for parameters
  stats.cumul_ODEBND -= time();
  TVT TVp[_np];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
    if( _IVPdepend.find(ip) != _IVPdepend.end() ){
      TVp[ip] = TVT( _pTM, isub, _pFPV[ip].I() );
      isub++;
    }
    else
      TVp[ip] = _pFPV[ip].I();
  }
    
  // Set Taylor model for states at stage times
  TVT* TVxk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    TVxk[is] = new TVT[_nx];
  
  // Compute state bounds at stage times
  bool failed = true;
  switch( meth ){
  case Options::VERIF:
    if( _pODEBND->bounds( _data.ns, _data.tk, TVp, TVxk )
        == ODEBND<T,DO>::NORMAL ) failed = false; break;
  case Options::DINEQ: default:
    if( _pODEDINEQ->bounds( _data.ns, _data.tk, TVp, TVxk )
        == ODEBND_GSL<T,DO>::NORMAL ) failed = false; break;
  }
  if( failed ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];
    stats.cumul_ODEBND += time();
    return SBB<T>::FAILURE;
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << *_pODEBND;
#endif
  stats.cumul_ODEBND += time();
  
  // Define an FPVar for the state variables participating in FPRelax environment
  stats.cumul_FPREL -= time();
  _pFPR->update_TModel( _pTMFP, _pTM, _IVPparam );
  std::map<int,bool>& DOdep = _DOstruct[1].dep();
  for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
      if( Op<T>::diam(TVxk[is][ix].B()) >= SBB<T>::options.INF ){
        failed = true; break;
      }
      std::map<int,bool>::iterator it = DOdep.find(ixs);
      if( it!=DOdep.end() )
        FPVxk[is][ix] = FPVT( _pFPR, TVxk[is][ix], _pTMFP );
    }
  stats.cumul_FPREL += time();
  
  // Clean-up
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];

  return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
}
// 
// template <typename T, typename DO>
// inline typename SBB<T>::STATUS
// DOSBB<T,DO>::_bound_IVP_TM_dineq
// ( const double*p, FPVT**FPVxk )
// {
//   // Set Taylor model for parameters
//   stats.cumul_ODEBND -= time();
//   T Ip[_np];
//   TVT TVp[_np];
//   for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
//     if( _IVPdepend.find(ip) != _IVPdepend.end() ){
//       Ip[ip]  = _pFPV[ip].I();
//       TVp[ip] = TVT( _pTM, isub, Ip[ip] );
//       isub++;
//     }
//     else{
//       Ip[ip]  = _pFPV[ip].I();
//       TVp[ip] = _pFPV[ip].I();
//     }
//   }
//     
//   // Set Taylor model for states at stage times
//   T* Ixk[_data.ns+1];
//   TVT* TVxk[_data.ns+1];
//   for( unsigned int is=0; is<=_data.ns; is++ ){
//     Ixk[is]  = new T[_nx];
//     TVxk[is] = new TVT[_nx];
//   }
// 
// #ifdef MC__DOSBB_DEBUG
//   _pODEDINEQ->options.DISPLAY = 2;
// #endif
//   
//   // Compute state bounds at stage times
//   typename ODEBND_GSL<T,DO>::STATUS flag_TM; //flag_IA, 
//   //flag_IA = _pODEDINEQ->bounds( _data.ns, _data.tk, Ip,  Ixk );
//   flag_TM = _pODEDINEQ->bounds( _data.ns, _data.tk, TVp, TVxk );
//   if( //flag_IA != ODEBND_GSL<T,DO>::NORMAL &&
//       flag_TM != ODEBND_GSL<T,DO>::NORMAL ){
// #ifdef MC__DOSBB_SHOW_FAILURE
//     std::cout << "FAILURE IN IVP BOUNDING: tf=" << _pODEDINEQ->final_time()
//               << std::endl; 
// #endif
//     for( unsigned int is=0; is<=_data.ns; is++ ){
//       delete[] Ixk[is]; delete[] TVxk[is];
//     }
//     stats.cumul_ODEBND += time();
//     return SBB<T>::FAILURE;
//   }
//   stats.cumul_ODEBND += time();
//   
//   // Define an FPVar for the state variables participating in FPRelax environment
//   stats.cumul_FPREL -= time();
//   _pFPR->update_TModel( _pTMFP, _pTM, _IVPparam );
//   std::map<int,bool>& DOdep = _DOstruct[1].dep();
//   bool failed = false;
//   for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ ){
//     for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
//       //if( Op<T>::diam(TVxk[is][ix].B()) >= SBB<T>::options.INF ){
//       //  failed = true; break;
//       //}
//       std::map<int,bool>::iterator it = DOdep.find( ixs );
//       if( it == DOdep.end() ) continue;
//       //if( flag_IA == ODEBND_GSL<T,DO>::NORMAL )
//       //  FPVxk[is][ix] = FPVT( _pFPR, Ixk[is][ix], _pTMFP );
//       //if( flag_TM == ODEBND_GSL<T,DO>::NORMAL )
//         FPVxk[is][ix] = FPVT( _pFPR, TVxk[is][ix], _pTMFP );
//     }
//   }
//   stats.cumul_FPREL += time();
//   
//   // Clean-up
//   stats.cumul_ODEBND -= time();
//   for( unsigned int is=0; is<=_data.ns; is++ ){
//     delete[] Ixk[is]; delete[] TVxk[is];
//   }
//   stats.cumul_ODEBND += time();
// 
//   return( failed? SBB<T>::FAILURE: SBB<T>::NORMAL );
// }

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_bound_IVP_TMMC
( const double*p, FPVT**FPVxk )
{
  // Set McCormick-Taylor model for parameters
  stats.cumul_ODEBND -= time();
  TVMCT TVMCp[_np];
  const unsigned int nsub = _IVPdepend.size();
  double MCpref[nsub];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
    if( _IVPdepend.find(ip) != _IVPdepend.end() ){
      TVMCp[ip] = TVMCT( _pTMMC, isub, MCT(_pFPV[ip].I(),p[ip]).sub(nsub,isub) );
      MCpref[isub++] = p[ip];
    }
    else
      TVMCp[ip] = MCT(_pFPV[ip].I(),p[ip]);
  }
    
  // Set McCormick-Taylor model for states at stage times
  TVMCT* TVMCxk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    TVMCxk[is] = new TVMCT[_nx];
  
  // Compute state bounds at stage times
  if( _pODEBND->bounds( _data.ns, _data.tk, TVMCp, TVMCxk )
      != ODEBND<T,DO>::NORMAL ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVMCxk[is];
    stats.cumul_ODEBND += time();
    return SBB<T>::FAILURE;
  }
#ifdef MC__DOSBB_DEBUG
  std::cout << *_pODEBND;
#endif
  stats.cumul_ODEBND += time();

  // Define an FPVar for the state variables participating in FPRelax environment
  stats.cumul_FPREL -= time();
  _pFPR->update_TModel( _pTMFP, _pTMMC, _IVPparam );
  std::map<int,bool>& DOdep = _DOstruct[1].dep();
  bool failed = false;
  for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
      if( Op<T>::diam(TVMCxk[is][ix].B().I()) >= SBB<T>::options.INF ){
        failed = true; break;
      }
      std::map<int,bool>::iterator it = DOdep.find(ixs);
      if( it!=DOdep.end() )
        FPVxk[is][ix] = FPVT( _pFPR, TVMCxk[is][ix], _pTMFP, MCpref, _IVPparam );
    }
  stats.cumul_FPREL += time();
  
  // Clean-up
  stats.cumul_ODEBND -= time();
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVMCxk[is];
  stats.cumul_ODEBND += time();
  return failed? SBB<T>::FAILURE: SBB<T>::NORMAL;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_relax_cost_constraints
( FPVT**FPVxk, const double INC, std::pair<FPVT,DOSTRUCT::t_OBJ>&DOobj )
{
  // Constraint definition
  stats.cumul_FPREL -= time();
  for( unsigned int ic=0; ic<_nc; ic++ ){
    std::pair<FPVT,DOSTRUCT::t_CTR> pctr = DO().CTR( ic, _pFPV, FPVxk, _data.ns );
    switch( pctr.second ){
      case DOSTRUCT::EQ: _pFPR->add_constraint( pctr.first, FPRT::EQ, 0. ); break;
      case DOSTRUCT::LE: _pFPR->add_constraint( pctr.first, FPRT::LE, 0. ); break;
      default:           _pFPR->add_constraint( pctr.first, FPRT::GE, 0. ); break;
    }
  }

  // Incumbent bound -- ALWAYS since this also defines the objective function!
   DOobj = DO().OBJ( _pFPV, FPVxk, _data.ns );
   switch( DOobj.second ){
      case DOSTRUCT::MIN: _pFPR->add_constraint( DOobj.first, FPRT::LE, INC ); break;
      default:            _pFPR->add_constraint( DOobj.first, FPRT::GE, INC ); break;
    }
  stats.cumul_FPREL += time();

  // Node processing: Constraint propagation
  stats.cumul_REDUC -= time();
  if ( options.USE_CONSTRAINT_PROPAGATION && !_pFPR->propagate_constraints() ){
    stats.cumul_REDUC += time();
    return SBB<T>::INFEASIBLE;
  }
  stats.cumul_REDUC += time();

  // Polyhedral relaxations
  stats.cumul_FPREL -= time();
  _pFPR->generate_cuts();
  if( options.USE_RLT_TMODEL ) _pFPR->generate_RLT_TModel( _pTMFP );
  stats.cumul_FPREL += time();

  return SBB<T>::NORMAL;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_domain_reduction
( FPVT**FPVxk, T*P, double&maxred )
{
  maxred = 0.;
  for( unsigned int ip=0; ip<_np; ip++ ){
    T P0 = P[ip];
    for( unsigned int ipb=0; ipb<2; ipb++ ){
      _pFPR->set_objective( (ipb? FPRT::MIN: FPRT::MAX), _pFPV[ip] );
#ifdef MC__DOSBB_DEBUG
      _pFPR->options.SOLVER_DISPLAY = true;
#endif
      switch( _pFPR->solve() ){
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
            //for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
            //throw Exceptions( Exceptions::REDUC );
            //return SBB<T>::INFEASIBLE;
          }
          break;
        case 1: // range lower bound
          if( !Op<T>::inter(  P[ip], P[ip],
              T(_pFPR->get_objective(),SBB<T>::options.INF) ) ){
            P[ip] = Op<T>::u( P[ip] );
            //for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
            //throw Exceptions( Exceptions::REDUC );
            //return SBB<T>::INFEASIBLE;
          }
          break;
        }
        _pFPR->update_bounds( _pFPV[ip], P[ip] );
#ifdef MC__DOSBB_DEBUG
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
        return SBB<T>::INFEASIBLE;
      default:
#ifdef MC__DOSBB_SHOW_FAILURE
        std::cout << "FAILURE IN DOMAIN REDUCTION LP\n"; 
#endif
        return SBB<T>::FAILURE;
      }
    }
    maxred = std::max( 1.-Op<T>::diam(P[ip])/Op<T>::diam(P0), maxred );
  }

  return SBB<T>::NORMAL;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_relaxation_solution
( const std::pair<FPVT,DOSTRUCT::t_OBJ>&DOobj, double*p, double&f )
{
  // Compute lower bound
  switch( DOobj.second ){
    case DOSTRUCT::MIN: _pFPR->set_objective( FPRT::MIN, DOobj.first ); break;
    default:            _pFPR->set_objective( FPRT::MAX, DOobj.first ); break;
  }
#ifdef MC__DOSBB_DEBUG
  _pFPR->options.SOLVER_DISPLAY = true;
#endif
  switch( _pFPR->solve() ){
#ifdef USE_CPLEX
  case IloAlgorithm::Optimal:
#else
  case GRB_OPTIMAL:
#endif
    // Copy solution value and solution point
    f = _pFPR->get_objective();
#ifdef MC__DOSBB_DEBUG
    std::cout << std::scientific << std::setprecision(5);
    std::cout << "  frel = " << f << std::endl;
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = _pFPR->get_variable(ip);
#ifdef MC__DOSBB_DEBUG
      std::cout << "  prel(" << ip << ") = " << p[ip] << std::endl;
#endif
    }
    return SBB<T>::NORMAL;
#ifdef USE_CPLEX
  case IloAlgorithm::Infeasible:
  case IloAlgorithm::InfeasibleOrUnbounded:
#else
  case GRB_INFEASIBLE:
  case GRB_INF_OR_UNBD:
#endif
    return SBB<T>::INFEASIBLE;
  default:
#ifdef MC__DOSBB_SHOW_FAILURE
    std::cout << "FAILURE IN LOWER BOUNDING LP\n"; 
#endif
      return SBB<T>::FAILURE;
  }
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::_subproblem_relaxed
( const unsigned int np, T*P, double*p, double&f, const double INC )
{
  stats.cumul_FPREL -= time();
  FPVT* FPVxk[_data.ns];
  for( unsigned int is=0; is<=_data.ns; is++ )
    FPVxk[is] = new FPVT[_nx];
  stats.cumul_FPREL += time();

  // Main loop for relaxation and domain reduction
  for( unsigned int ired=0; !ired || ired < options.DOMAIN_REDUCTION_MAXLOOPS; ired++ ){

    // Set current variables
    stats.cumul_FPREL -= time();
    _pFPR->reset();
    for( unsigned int ip=0; ip<_np; ip++ )
      _pFPV[ip] = FPVT( _pFPR, ip, P[ip] );
#ifdef MC__DOSBB_DEBUG
    std::cout << "FPR:\n" << *_pFPR;
#endif
    stats.cumul_FPREL += time();

    // State polyhedral relaxation
    typename SBB<T>::STATUS status;
    try{
      switch( options.IVP_BOUNDING_TYPE ){
      // State polyhedral relaxations via Intervals
      case Options::IA:
        status = _bound_IVP_IA( p, FPVxk, options.IVP_BOUNDING_PROPAGATION ); break;
      // State polyhedral relaxations via McCormick relaxations
      case Options::MC:
        status = _bound_IVP_MC( p, FPVxk ); break;
      // State polyhedral relaxations via Taylor models
      case Options::TM:
        status = _bound_IVP_TM( p, FPVxk, options.IVP_BOUNDING_PROPAGATION ); break;
      // State polyhedral relaxations via McCormick-Taylor models
      case Options::TMMC:
        status = _bound_IVP_TMMC( p, FPVxk ); break;
      default:
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        throw Exceptions( Exceptions::UNDEF );
      }

     // State bounding/relaxation unsuccessful
     if( status != SBB<T>::NORMAL ){
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
#ifdef MC__DOSBB_SHOW_FAILURE
        std::cout << "FAILURE IN ODE BOUNDING: t=" << _pODEBND->final_time()
                  << std::endl;
#endif
        return SBB<T>::FAILURE;
      }
    }
    catch(...){
#ifdef MC__DOSBB_SHOW_FAILURE
      std::cout << "FAILURE IN ODE BOUNDING: t=" << _pODEBND->final_time()
                << std::endl;
#endif
      return SBB<T>::FAILURE;
    }

    // Cost and constraint polyhedral relaxation
    std::pair<FPVT,DOSTRUCT::t_OBJ> DOobj;
    status = _relax_cost_constraints( FPVxk, INC, DOobj );
    if( status != SBB<T>::NORMAL ){
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        return status;
    }

    // Node processing: Domain reduction
    if( options.USE_DOMAIN_REDUCTION && ired < options.DOMAIN_REDUCTION_MAXLOOPS ){
      double maxred;
      stats.cumul_REDUC -= time();
      status = _domain_reduction( FPVxk, P, maxred );
      stats.cumul_REDUC += time();
      if( status != SBB<T>::NORMAL ){
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        return status;
      }

      // Reconstruct lower bounding problem when sufficient tightening was achieved
#ifdef MC__DOSBB_SHOW_REDUCTION
      std::cerr << "max. domain reduction: " << maxred << std::endl;
#endif
      if( maxred > options.DOMAIN_REDUCTION_THRESHOLD
       && ired+1 < options.DOMAIN_REDUCTION_MAXLOOPS ) continue;
    }

    // Solve relaxed problem
    stats.cumul_DOREL -= time();
    status = _relaxation_solution( DOobj, p, f );
    stats.cumul_DOREL += time();
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
    return status;
  }

#ifdef MC__DOSBB_SHOW_FAILURE
  std::cout << "FAILURE AT EXIT\n"; 
#endif
  return SBB<T>::FAILURE;
}

template <typename T, typename DO>
inline typename SBB<T>::STATUS
DOSBB<T,DO>::subproblems
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

template <typename T, typename DO>
inline void
DOSBB<T,DO>::Options::display
( std::ostream&out ) const
{
  // Display DOSBB Options
  out << std::setw(60) << "  BOUNDING STRATEGY FOR PARAMETRIC IVP";
  switch( IVP_BOUNDING_TYPE ){
  case IA:   out << "IA\n";   break;
  case MC:   out << "MC\n";   break;
  case TM:   out << "TM\n";   break;
  case TMMC: out << "TMMC\n"; break;
  }
  out << std::setw(60) << "  ORDER OF TAYLOR MODEL"
      << TAYLOR_MODEL_ORDER << std::endl;
  out << std::setw(60) << "  USE CONSTRAINT PROPAGATION?"
      << (USE_CONSTRAINT_PROPAGATION?"Y\n":"N\n");
  out << std::setw(60) << "  USE OPTIMIZATION-BASED DOMAIN REDUCTION?"
      << (USE_DOMAIN_REDUCTION?"Y\n":"N\n");
  out << std::setw(60) << "  USE RLT TO TIGHTEN TAYLOR MODEL RELAXATIONS?"
      << (USE_RLT_TMODEL?"Y\n":"N\n");
  out << std::setw(60) << "  OPTIMIZATION-BASED REDUCTION MAX LOOPS"
      << DOMAIN_REDUCTION_MAXLOOPS << std::endl;
  out << std::setw(60) << "  LOOP THRESHOLD FOR OPTIMIZATION-BASED REDUCTION"
      << std::fixed << std::setprecision(0)
      << DOMAIN_REDUCTION_THRESHOLD*1e2 << "%\n";
}

template <typename T, typename DO>
inline std::ostream&
operator <<
( std::ostream&out, const DOSBB<T,DO>&pDO )
{
  out << std::endl
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ')
      << std::setw(52) << "GLOBAL DYNAMIC OPTIMIZER IN MC++\n"
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ');

  // Display DOSBB Options
  out << std::left << "SPATIAL BRANCH-AND-BOUND OPTIONS:\n\n";
  pDO.SBB<T>::options.display( out );
  pDO.DOSBB<T,DO>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::setfill(' ')
      << std::endl << std::endl;
  return out;
}

} // end namescape mc

#endif
