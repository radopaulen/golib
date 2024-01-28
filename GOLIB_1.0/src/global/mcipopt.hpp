// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_MCIPOPT Local (Continuous) Optimization using IPOPT and FADBAD++
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

Consider the mathematical program (NLP) \f$\mathcal{P}:\min\{f({\bf x}):g_j({\bf x})\{\leq,=,\geq\} 0,  j=1,\ldots,n_g; x_i^\ell\leq x_i\leq x_i^u, i=1,\ldots,n_x\}\f$, where \f$f, g_1, \ldots,g_{n_g}\f$ are factorable, potentially nonlinear, real-valued functions. The class mc::IPOPT in MC++ solves such NLP problems using the software package <A href="https://projects.coin-or.org/Ipopt">IPOPT</A>, which implements a local solution method (interior point). IPOPT requires the first and second derivatives as well as the sparsity pattern of the objective and constraint functions in the NLP model. This information is generated automatically in mc::IPOPT using <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for the differentiations and mc::Structure in MC++ for the sparsity pattern.

mc::IPOPT is templated in the NLP problem to be solved, which has to be a class derived from mc::NLPSTRUCT. For example, suppose we want to solve the following NLP:
\f{align*}
  \max_{\bf p}\ & p_1+p_2 \\
  \text{s.t.} \ & p_1\,p_2 \leq 4 \\
  & 3 \leq p_1 \leq 6\\
  & 0 \leq p_2 \leq 4.
\f}
The following class is defined:
\code
  #include "mcipopt.h"
  
  const unsigned int NP = 2;
  const unsigned int NC = 1;

  class NLP : public virtual mc::NLPSTRUCT
  {
  public:
    NLP(): mc::NLPSTRUCT( NP, NC )
      {}

    template <typename U>
    std::pair<U,t_OBJ> OBJ
      ( const U*p )
      {
        return std::make_pair( p[0]+p[1], MAX );
      }

    template <typename U>
    std::pair<U,t_CTR> CTR
      ( const unsigned int ic, const U*p )
      {
        switch( ic ){
          case 0: return std::make_pair( p[0]*p[1]-4, LE );
          default: throw std::runtime_error("invalid size");
        }
      }
  };
\endcode
where NP and NC stand for the number of parameters and constraints, respectively; and the templated functions OBJ and CTR define the objective function and left-hand side of the constraints, respectively.

\section sec_MCIPOPT_solve How to Solve an NLP Model using mc::IPOPT?

We start by defining the parameter set \f$P\f$ and initial guess \f$p_0\f$:

\code
  std::pair<double,double> P[NP];
  P[0] = std::make_pair( 3., 6. );
  P[1] = std::make_pair( 0., 4. );
  double p0[NP];
  p0[0] = 1.;
  p0[1] = 2.;
\endcode

An mc::IPOPT object is then defined:

\code
  Ipopt::SmartPtr< mc::IPOPT<NLP> > pNLP = new mc::IPOPT<NLP>;
\endcode

Note that the class <A href="http://www.coin-or.org/Doxygen/Ipopt/class_ipopt_1_1_smart_ptr.html">Ipopt::SmartPtr</A> is used here to comply with the IPOPT guidelines. 

The NLP model is optimized as follows:

\code
  Ipopt::ApplicationReturnStatus status = pNLP->solve( P, p0 );
\endcode

where the first argument is the variable bounds and the second one, the initial guess used by the optimizer. The return value is of the enumeration type <A href="http://www.coin-or.org/Doxygen/Ipopt/namespace_ipopt.html#efa0497854479cde8b0994cdf132c982">Ipopt::ApplicationReturnStatus</A>. Moreover, the output level, maximum number of iterations, tolerance, and maximum CPU time can all be modified through the subsequent arguments of mc::IPOPT::solve. 

Finally, the optimization results can be retrieved as follows:

\code
  #include <iostream>
  #include <iomanip>

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << pNLP->get_objective() << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pNLP->get_variable(ip) << std::endl;
  }
\endcode

The following result is displayed:

\verbatim
NLP (LOCAL) SOLUTION: 
  f* = 6.66667
  p*(0) = 6
  p*(1) = 0.666667
\endverbatim
*/

#ifndef MC__IPOPT_HPP
#define MC__IPOPT_HPP

#include <stdexcept>
#include <assert.h>
#include "fadiff.h"
#include "badiff.h"
#include "coin/IpTNLP.hpp"
#include "coin/IpIpoptApplication.hpp"

#include "nlpstruct.hpp"
#include "structure.hpp"

#undef MC__IPOPT_DEBUG
#undef MC__IPOPT_TRACE

/* TO DO:
- [DONE] documentation
- [DONE] set options
- [DONE] retreive solution information
- [DONE] implement sparse Hessian
   [DONE] . make Sparse compatible with FADBAD++
   [DONE] . extend Sparse to determine linearity
   [DONE] . rename Sparse -> Structure
- [DONE] revise implementation to pass NLP model as template argument
- overload << operator for NLP solution
*/

namespace mc
{

//! @brief C++ derived class for NLP solution using IPOPT and FADBAD++
////////////////////////////////////////////////////////////////////////
//! mc::IPOPT is a C++ derived class for solving NLP problems using
//! IPOPT and FADBAD++
////////////////////////////////////////////////////////////////////////
template <typename NLP>
class IPOPT: public Ipopt::TNLP, public virtual NLPSTRUCT
{
  // Overloading stdout operator
  template <typename T> friend std::ostream& operator<<
    ( std::ostream&os, const IPOPT<T>& );

private:
  // AD types in FADBAD++
  typedef fadbad::F< double > t_F;
  typedef fadbad::B< double > t_B;
  typedef fadbad::B< fadbad::F< double > > t_BF;

  //! @brief Pointer to constraint structure
  std::map<int,bool>* _gstruct;
  //! @brief Pointer to Lagrangian structure
  std::map<int,bool>* _Lstruct;
  //! @brief Internal scaling factor for the objective function (+1: minimize; -1: maximize)
  double _scaling;

  //! @brief Structure holding solution information
  struct DATA{
    const double*p0;
    const std::pair<double,double>*P;
  } _data;

  //! @brief Structure holding function evaluation results
  struct EVALUATION{
    std::map<int,bool>* gstruct;
    std::map<int,bool>* Lstruct;
    double scaling;
  } _evaluation;

  //! @brief Structure holding solution information
  struct SOLUTION{
    Ipopt::SolverReturn status;
    Ipopt::Index n;
    Ipopt::Number*p;
    Ipopt::Number*upL;
    Ipopt::Number*upU;
    Ipopt::Index m;
    Ipopt::Number*g;
    Ipopt::Number*ug;
    Ipopt::Number f;
  } _solution;

public:
  /** @defgroup MCIPOPT Local (Continuous) Optimization using IPOPT and FADBAD++
   *  @{
   */
  //! @brief Constructor
  IPOPT()
    : NLPSTRUCT(NLP())
    {
      _evaluation.gstruct = _evaluation.Lstruct = 0;
      _solution.n = _solution.m = 0;
    }

  //! @brief Destructor
  virtual ~IPOPT()
    { delete[] _evaluation.gstruct; delete _evaluation.Lstruct;
      if( _solution.n ){
        delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
      }
      if( _solution.m ){
        delete[] _solution.g; delete[] _solution.ug;
      }
    }

  //! @brief Structure for setting up storing the solver options
  struct Options
  {
    //! @brief Constructor
    Options():
      CVTOL(1e-8), PRIMALTOL(1e-8), DUALTOL(1e-4), COMPLTOL(1e-7),
      MAXITER(100), MAXCPU(1e6), SPARSE(true), GRADIENT(FORWARD),
      HESSIAN(EXACT), DISPLAY(0)
      {} 
    //! @brief Enumeration type for Hessian strategy
    enum HESSIAN_STRATEGY{
      EXACT=0, 	//!< Use exact second derivatives
      LBFGS,	//!< Perform a limited-memory quasi-Newton approximation
    };
    //! @brief Enumeration type for gradient strategy
    enum GRADIENT_STRATEGY{
      FORWARD=0,	//!< Forward AD
      BACKWARD		//!< Backward AD
    };
    //! @brief Convergence tolerance
    double CVTOL;
    //! @brief Tolerance on primal feasibility
    double PRIMALTOL;
    //! @brief Tolerance on dual feasibility
    double DUALTOL;
     //! @brief Tolerance on complementarity conditions
    double COMPLTOL;
   //! @brief Maximum number of iterations
    int MAXITER;
    //! @brief Maximum run time (seconds)
    double MAXCPU;
    //! @brief Whether or not to account for the sparsity of the NLP model?
    bool SPARSE;
    //! @brief Strategy for gradient computation
    GRADIENT_STRATEGY GRADIENT;
    //! @brief Strategy for Hessian computation during linesearch
    HESSIAN_STRATEGY HESSIAN;
    //! @brief Display level in IPOPT
    int DISPLAY;
  } options;

  //! @brief Solve NLP model -- return value is status
  Ipopt::ApplicationReturnStatus solve
    ( const std::pair<double,double>*P, const double*p0 )
    {
      Ipopt::SmartPtr<Ipopt::IpoptApplication> IpoptApp
        = new Ipopt::IpoptApplication();

      // Set (a few) IPOPT options
      set_options( IpoptApp );

      // Keep track of bounds and initial guess
      _data.p0 = p0;
      _data.P  = P;

      // Run NLP solver
      Ipopt::ApplicationReturnStatus status = IpoptApp->Initialize();
      if( status == Ipopt::Solve_Succeeded ){
        status = IpoptApp->OptimizeTNLP( this );
      }

      return status;
    }

  //! @brief Get IPOPT internal scaling value
  double get_scaling()
    {
      Structure *pS = new Structure[_np];
      std::pair<Structure,t_OBJ> f = NLP().OBJ( pS );
      // set scaling factor. 1: minimize; -1: maximize
      switch( f.second ){
        case MIN: _evaluation.scaling =  1.; break;
        case MAX: _evaluation.scaling = -1.; break;
      }

      delete[] pS;
      return _evaluation.scaling;
    }

  //! @brief Get IPOPT solution info
  const SOLUTION& solution() const
    {
      return _solution;
    }
// 
//   //! @brief Get IPOPT objective value
//   double get_objective() const
//     {
//       return _solution.f;
//     }
//   //! @brief Get IPOPT ip variable value
//   double get_variable
//     ( const unsigned int ix ) const
//     {
//       assert( ix < _solution.n );
//       return _solution.p[ix];
//     }
//   //! @brief Get IPOPT ig constraint value
//   double get_constraint
//     ( const unsigned int ig ) const
//     {
//       assert( ig < _solution.m );
//       return _solution.g[ig];
//     }
  /** @} */

protected:
  //! @brief Method to return some info about the NLP
  virtual bool get_nlp_info
    ( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
      Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style );

  //! @brief Method to return the variable and constraint bounds
  virtual bool get_bounds_info
    ( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
      Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u );

  //! @brief Method to return the initial point for the NLP solver
  virtual bool get_starting_point
    ( Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
      Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m,
      bool init_lambda, Ipopt::Number* lambda );

  //! @brief Method to return the objective function value
  virtual bool eval_f
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number& obj_value );

  //! @brief Method to return the objective function gradient
  virtual bool eval_grad_f
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number* grad_f );

  //! @brief Method to return the constraint residuals
  virtual bool eval_g
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Index m, Ipopt::Number* g );

  //! @brief Method to return the structure of the jacobian (if "values" is NULL) and the values of the jacobian (otherwise)
  virtual bool eval_jac_g
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
      Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
      Ipopt::Number* values );

  //! @brief Method to return the structure of the hessian of the lagrangian (if "values" is NULL) and the values of the hessian of the lagrangian (otherwise)
  virtual bool eval_h
    ( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
      Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
      bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
      Ipopt::Index* jCol, Ipopt::Number* values );

  //! @brief Method called when the algorithm is complete
  virtual void finalize_solution
    ( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x,
      const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index m,
      const Ipopt::Number* g, const Ipopt::Number* lambda,
      Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq );

  //! @brief Method called at the end of each iteration
  virtual bool intermediate_callback
    ( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
      Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
      Ipopt::Number d_norm, Ipopt::Number regularization_size,
      Ipopt::Number alpha_du, Ipopt::Number alpha_pr,
      Ipopt::Index ls_trials, const Ipopt::IpoptData* ip_data,
      Ipopt::IpoptCalculatedQuantities* ip_cq );

private:
  //! @brief Set IPOPT options
  void set_options
    ( Ipopt::SmartPtr<Ipopt::IpoptApplication>&IpoptApp )
    {
      IpoptApp->Options()->SetNumericValue( "tol", options.CVTOL<0.?
        1e-12: options.CVTOL );
      IpoptApp->Options()->SetNumericValue( "dual_inf_tol", options.DUALTOL<=0.?
        1e-12: options.DUALTOL );
      IpoptApp->Options()->SetNumericValue( "constr_viol_tol", options.PRIMALTOL<=0.?
        1e-12: options.PRIMALTOL );
      IpoptApp->Options()->SetNumericValue( "compl_inf_tol", options.COMPLTOL<=0.?
        1e-12: options.COMPLTOL );
      IpoptApp->Options()->SetIntegerValue( "max_iter", options.MAXITER );
      IpoptApp->Options()->SetNumericValue( "max_cpu_time", options.MAXCPU);
      IpoptApp->Options()->SetNumericValue( "obj_scaling_factor", get_scaling() );
      switch( options.HESSIAN ){
      case Options::EXACT:
        IpoptApp->Options()->SetStringValue( "hessian_approximation", "exact" );
	break;
      case Options::LBFGS:
        IpoptApp->Options()->SetStringValue( "hessian_approximation", "limited-memory" );
	break;
      }
      IpoptApp->Options()->SetIntegerValue( "print_level",
        (options.DISPLAY<0? 0: (options.DISPLAY>12? 12: options.DISPLAY )) );
    }

  //! @brief Private methods to block default compiler methods
  IPOPT(const IPOPT&);
  IPOPT& operator=(const IPOPT&);
};


template <typename NLP> inline
bool
IPOPT<NLP>::get_nlp_info
( Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
  Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::get_nlp_info\n";
#endif
  // set size
  n = _np; m = _nc;

  // determine constraint sparsity
  delete[] _evaluation.gstruct;
  _evaluation.gstruct = new std::map<int,bool>[_nc];
  Structure *pS = new Structure[_np];
  for( unsigned int ip=0; ip<_np; ip++ ) pS[ip].indep(ip);   
  nnz_jac_g = 0;
  for( unsigned int ic=0; ic<_nc; ic++ ){
    std::pair<Structure,t_CTR> g = NLP().CTR( ic, pS );
    _evaluation.gstruct[ic] = g.first.dep();
    nnz_jac_g += _evaluation.gstruct[ic].size();
  }

  // determine Lagrangian sparsity
  delete _evaluation.Lstruct;
  std::pair<Structure,t_OBJ> f = NLP().OBJ( pS );
  Structure L = f.first;
#ifdef MC__IPOPT_DEBUG
  std::cout << "h_lag:" << L << std::endl;
#endif
  for( unsigned int ic=0; ic<_nc; ic++ ){
    std::pair<Structure,t_CTR> g = NLP().CTR( ic, pS );
    L += g.first;
#ifdef MC__IPOPT_DEBUG
    std::cout << "h_lag:" << L << std::endl;
#endif
  }
  _evaluation.Lstruct = new std::map<int,bool>( L.dep() ); 
  if( !options.SPARSE )
    nnz_h_lag = (n*(n+1))/2;
  else{
    int nnz_jac_lag = 0;
    std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
    for( ; cit != _evaluation.Lstruct->end(); ++cit )
      if( !cit->second ) nnz_jac_lag++;
    nnz_h_lag = (nnz_jac_lag*(nnz_jac_lag+1))/2;
  }

  // use the C style indexing (0-based)
  index_style = Ipopt::TNLP::C_STYLE;

#ifdef MC__IPOPT_DEBUG
  std::cout << "n:" << n << std::endl;
  std::cout << "m:" << m << std::endl;
  std::cout << "nnz_jac_g:" << nnz_jac_g << std::endl;
  std::cout << "nnz_h_lag:" << nnz_h_lag << std::endl;
  std::cout << "h_lag:" << L << std::endl;
#endif

  // clean up
  delete[] pS;

  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::get_bounds_info
( Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
  Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::get_bounds_info\n";
#endif
  // set variable bounds
  for( unsigned int ip=0; ip<_np; ip++ ){   
    x_l[ip] = ( _data.P? _data.P[ip].first: -INF );
    x_u[ip] = ( _data.P? _data.P[ip].second: INF );
#ifdef MC__IPOPT_DEBUG
    std::cout << "  x_l[" << ip << "] = " << x_l[ip]
              << "  x_u[" << ip << "] = " << x_u[ip] << std::endl;
#endif
  }

  // set constraint bounds
  Structure *pS = new Structure[_np];
  for( unsigned int ic=0; ic<_nc; ic++ ){
    switch( NLP().CTR( ic, pS ).second ){
      case EQ: g_l[ic] = g_u[ic] = 0.; break;
      case LE: g_l[ic] = -INF; g_u[ic] = 0.; break;
      case GE: g_l[ic] = 0.; g_u[ic] = INF; break;
    }
#ifdef MC__IPOPT_DEBUG
    std::cout << "  g_l[" << ic << "] = " << g_l[ic]
              << "  g_u[" << ic << "] = " << g_u[ic] << std::endl;
#endif
  }
  delete[] pS;

  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::get_starting_point
( Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z,
  Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m,
  bool init_lambda, Ipopt::Number* lambda )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::get_starting_point  "
              << init_x << init_z << init_lambda << std::endl;
#endif
  // Here, we assume we only have starting values for x, if you code
  // your own NLP, you can provide starting values for the dual variables
  // if you wish
  if( !init_x || init_z || init_lambda ) return false;

  // initialize to the given starting point
  if( !_data.p0 ) return false;
  for( unsigned int ip=0; ip<_np; ip++ ){
    x[ip] = _data.p0[ip];
#ifdef MC__IPOPT_DEBUG
    std::cout << "  x_0[" << ip << "] = " << x[ip] << std::endl;
#endif
  }
  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::eval_f
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number& f )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::eval_f  " << new_x << std::endl;
#endif
  // evaluate objective function
  f = NLP().OBJ( x ).first;
#ifdef MC__IPOPT_DEBUG
  std::cout << "  f = " << f << std::endl;
#endif

  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::eval_grad_f
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number* df )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::eval_grad_f  " << new_x << std::endl;
#endif
  switch( options.GRADIENT ){
   // evaluate objective function gradient (Forward AD)
   case Options::FORWARD:{
    t_F*dx = new t_F[_np];
    for( unsigned int ip=0; ip<_np; ip++ ){
      dx[ip] = x[ip];
      dx[ip].diff(ip,_np);
    }
    t_F dfobj = NLP().OBJ( dx ).first;
    for( unsigned int ip=0; ip<_np; ip++ ){
      df[ip] = dfobj.d(ip);  
#ifdef MC__IPOPT_DEBUG
      std::cout << "  df[" << ip << "] = " << df[ip] << std::endl;
#endif
    }
    delete[] dx; break;
   }

   // evaluate objective function gradient (Backward AD)
   case Options::BACKWARD:{
    t_B*dx = new t_B[_np];
    for( unsigned int ip=0; ip<_np; ip++ )
      dx[ip] = x[ip];
    t_B dfobj = NLP().OBJ( dx ).first;
    dfobj.diff(0,1);
    for( unsigned int ip=0; ip<_np; ip++ ){
      df[ip] = dx[ip].d(0);
#ifdef MC__IPOPT_DEBUG
      std::cout << "  df[" << ip << "] = " << df[ip] << std::endl;
#endif
    }
    delete[] dx; break;
   }
  }
  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::eval_g
( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
  Ipopt::Number* g )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::eval_g  " << new_x << std::endl;
#endif
  // evaluate constraint functions
  for( unsigned int ic=0; ic<_nc; ic++ ){
    g[ic] = NLP().CTR( ic, x ).first;
#ifdef MC__IPOPT_DEBUG
    std::cout << "  g[" << ic << "] = " << g[ic] << std::endl;
#endif
  }
  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::eval_jac_g
( Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m,
  Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
  Ipopt::Number* dg )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::eval_jac_g  " << new_x << std::endl;
#endif
  // return the constraint Jacobian structure
  if( !dg ){
    unsigned int ie = 0;
    for( unsigned int ic=0; ic<_nc; ic++ ){
      std::map<int,bool>::const_iterator it = _evaluation.gstruct[ic].begin();
      for( ; it != _evaluation.gstruct[ic].end(); ++it, ++ie ){
        iRow[ie] = ic;
        jCol[ie] = it->first;
#ifdef MC__IPOPT_DEBUG
        std::cout << "  dg[" << ic << ", " << it->first << "]" << std::endl;
#endif
      }
    }
    assert( ie == (unsigned int)nele_jac );
    return true;
  }
  
  switch( options.GRADIENT ){
   // evaluate constraint Jacobian (Forward AD)
   case Options::FORWARD:{
    t_F*dx = new t_F[_np];
    for( unsigned int ip=0; ip<_np; ip++ ){
      dx[ip] = x[ip];
      dx[ip].diff(ip,_np);
    }
    unsigned int ie = 0;
    for( unsigned int ic=0; ic<_nc; ic++ ){
      t_F dgctr = NLP().CTR( ic, dx ).first;
      std::map<int,bool>::const_iterator it = _evaluation.gstruct[ic].begin();
      for( ; it != _evaluation.gstruct[ic].end(); ++it, ie++ ){
        dg[ie] = dgctr.d(it->first);  
#ifdef MC__IPOPT_DEBUG
        std::cout << "  dg[" << ic << ", " << it->first << "] = " << dg[ie]
                  << std::endl;
#endif
      }
    }
    delete[] dx; break;
   }

   // evaluate constraint Jacobian (Forward AD)
   case Options::BACKWARD:{
    t_B*dx = new t_B[_np];
    unsigned int ie = 0;
    for( unsigned int ic=0; ic<_nc; ic++ ){
      for( unsigned int ip=0; ip<_np; ip++ )
        dx[ip] = x[ip];
      t_B dgctr = NLP().CTR( ic, dx ).first;
      dgctr.diff(0,1);
      std::map<int,bool>::const_iterator it = _evaluation.gstruct[ic].begin();
      for( ; it != _evaluation.gstruct[ic].end(); ++it, ie++ ){
        dg[ie] = dx[it->first].d(0);
#ifdef MC__IPOPT_DEBUG
        std::cout << "  dg[" << ic << ", " << it->first << "] = " << dg[ie]
                  << std::endl;
#endif
      }
    }
    delete[] dx; break;
   }
  }

#ifdef MC__IPOPT_DEBUG
  int dum; std::cin >> dum;
#endif

  return true;
}

template <typename NLP> inline
bool
IPOPT<NLP>::eval_h
( Ipopt::Index n, const Ipopt::Number* x, bool new_x,
  Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
  bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
  Ipopt::Index* jCol, Ipopt::Number* d2L )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::eval_h  " << new_x  << new_lambda << std::endl;
#endif
  // return the Lagangian Hessian structure (dense)
  if( !d2L ){
    unsigned int ie = 0;

    // Dense Hessian structure
    if( !options.SPARSE ){
      for( unsigned int ip=0; ip<_np; ip++ ){
        for( unsigned int jp=ip; jp<_np; jp++, ie++ ){
          iRow[ie] = ip;
          jCol[ie] = jp;
#ifdef MC__IPOPT_DEBUG
          std::cout << "  d2L[" << ip << ", " << jp << "]" << std::endl;
#endif
        }
      }
    }

    // Sparse Hessian structure
    else{
      std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
      for( ; cit != _evaluation.Lstruct->end(); ++cit ){
        if( cit->second ) continue;
        std::map<int,bool>::const_iterator cjt = cit;
        for( ; cjt != _evaluation.Lstruct->end(); ++cjt ){
          if( cjt->second ) continue;
          iRow[ie] = cit->first;
          jCol[ie] = cjt->first;
#ifdef MC__IPOPT_DEBUG
          std::cout << "  d2L[" << cit->first << ", " << cjt->first << "]" << std::endl;
#endif
          ie++;
        }
      }
    }

    assert( ie == (unsigned int)nele_hess );
    return true;
  }
  
  // return the Lagrangian Hessian values
  t_F*dx = new t_F[_np];
  for( unsigned int ip=0; ip<_np; ip++ ){
    dx[ip] = x[ip];
    dx[ip].diff(ip,_np);
  }
  t_BF*d2x = new t_BF[_np];
  for( unsigned int jp=0; jp<_np; jp++ ){
    d2x[jp] = dx[jp];
  }
  t_BF BFL = NLP().OBJ( d2x ).first;
  for( unsigned int ic=0; ic<_nc; ic++ )
    BFL += lambda[ic] * NLP().CTR( ic, d2x ).first;
  BFL.diff(0,1);
  unsigned int ie = 0;

  // Dense Hessian structure
  if( !options.SPARSE ){
    for( unsigned int ip=0; ip<_np; ip++ ){
      for( unsigned int jp=ip; jp<_np; jp++, ie++ ){
        d2L[ie] = (d2x[ip].d(0)).d(jp);
#ifdef MC__IPOPT_DEBUG
        std::cout << "  d2L[" << ip << ", " << jp << "] = " << d2L[ie]
                  << std::endl;
#endif
      }
    }
  }

  // Sparse Hessian structure
  else{
    std::map<int,bool>::const_iterator cit = _evaluation.Lstruct->begin();
    for( ; cit != _evaluation.Lstruct->end(); ++cit ){
      if( cit->second ) continue;
      std::map<int,bool>::const_iterator cjt = cit;
      for( ; cjt != _evaluation.Lstruct->end(); ++cjt ){
        if( cjt->second ) continue;
        d2L[ie] = (d2x[cit->first].d(0)).d(cjt->first);
#ifdef MC__IPOPT_DEBUG
        std::cout << "  d2L[" << cit->first << ", " << cjt->first << "] = " << d2L[ie]
                  << std::endl;
#endif
        ie++;
      }
    }
  }

  delete[] dx;
  delete[] d2x;

  return true;
}

template <typename NLP> inline
void
IPOPT<NLP>::finalize_solution
( Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* p,
  const Ipopt::Number* upL, const Ipopt::Number* upU, Ipopt::Index m,
  const Ipopt::Number* g, const Ipopt::Number* ug, Ipopt::Number f,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::finalize_solution\n";
#endif
  _solution.status = status;

  // Successful (or near-successful) completion
  //if( status == Ipopt::SUCCESS || status == Ipopt::STOP_AT_ACCEPTABLE_POINT
  // || status == Ipopt::STOP_AT_TINY_STEP ){
    // resize solution arrays
    if( _solution.n != n ){
      if( _solution.n ){
        delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
      }
      _solution.n = n;
      _solution.p = new Ipopt::Number[n];
      _solution.upL = new Ipopt::Number[n];
      _solution.upU = new Ipopt::Number[n];
    }
    if( _solution.m != m ){
      if( _solution.m ){
        delete[] _solution.g; delete[] _solution.ug;
      }
      _solution.m = m;
      _solution.g = new Ipopt::Number[m];
      _solution.ug = new Ipopt::Number[m];
    }
    // copy solution values into _solution
    _solution.f = f;
    for( Ipopt::Index i=0; i<n; i++ ){
      _solution.p[i] = p[i]; _solution.upL[i] = upL[i]; _solution.upU[i] = upU[i];
    }
    for( Ipopt::Index j=0; j<m; j++ ){
      _solution.g[j] = g[j]; _solution.ug[j] = ug[j];
    }
  //}

  // Failure
  //else{
  //  if( _solution.n ){
  //    delete[] _solution.p; delete[] _solution.upL; delete[] _solution.upU;
  //    _solution.n = 0;
  //  }
  //  if( _solution.m ){
  //    delete[] _solution.g; delete[] _solution.ug;
  //    _solution.m = 0;
  //  }
  //}
  return;
}

template <typename NLP> inline
bool
IPOPT<NLP>::intermediate_callback
( Ipopt::AlgorithmMode mode, Ipopt::Index iter, Ipopt::Number obj_value,
  Ipopt::Number inf_pr, Ipopt::Number inf_du, Ipopt::Number mu,
  Ipopt::Number d_norm, Ipopt::Number regularization_size,
  Ipopt::Number alpha_du, Ipopt::Number alpha_pr, Ipopt::Index ls_trials,
  const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq )
{
#ifdef MC__IPOPT_TRACE
    std::cout << "  IPOPT<NLP>::intermediate_callback\n";
#endif
  return true;
}

} // end namescape mc

#endif
