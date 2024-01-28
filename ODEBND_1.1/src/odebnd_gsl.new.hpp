// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_GSL_HPP
#define MC__ODEBND_GSL_HPP

#undef  MC__ODEBND_GSL_TM_DEBUG
#undef  MC__ODEBND_GSL_SAMPLE_DEBUG
#undef  MC__ODEBND_GSL_LINTR_DEBUG
#define MC__ODEBND_GSL_LINTR_TMJAC
#undef  MC__ODEBND_GSL_ELLTR_SCALED
#undef  MC__ODEBND_GSL_FI_USE

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "odestruct.hpp"
#include "ellipsoid_cpplapack.hpp"
#include "odeslv_gsl.hpp"

#include "mcop.hpp"
#include "tmodel.hpp"
#ifdef MC__ODEBND_GSL_FI_USE
  #include "finterval.hpp"
#endif

#include "mcfadbad.hpp"
//#include "fadiff.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using GSL and MC++.
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_GSL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using GSL and MC++. It implements the methods of differential 
//! inequalities and ellipsoidal bounding, as well as their combination
//! with Taylor models to enable high-order convergence.
////////////////////////////////////////////////////////////////////////
template <typename T, typename IVP>
class ODEBND_GSL: public BASE_GSL, public ODESTRUCT
{
private:

  typedef fadbad::F< double > F;
  typedef fadbad::F< T > FT;
  typedef Ellipsoid E;
  typedef TModel<T> TMT;
  typedef TVar<T> TVT;
  typedef fadbad::F< TVT > FTVT;
#ifdef MC__ODEBND_GSL_FI_USE
  typedef FI< T > FIT;
#endif

  //! @brief GSL data type for ODE system (auxiliary ODEs: differential inequalities)
  gsl_odeiv2_system _sys_dineqI;

  //! @brief GSL driver for ODE system (auxiliary ODEs: differential inequalities)
  gsl_odeiv2_driver *_driver_dineqI;

  //! @brief GSL data type for ODE system (auixiliary ODEs: differential inequalities w/ Taylor models)
  gsl_odeiv2_system _sys_dineqTM;

  //! @brief GSL driver for ODE system (auixiliary ODEs: differential inequalities w/ Taylor models)
  gsl_odeiv2_driver *_driver_dineqTM;

  //! @brief ODE solver based on GSL
  ODESLV_GSL<T,IVP>* _ODESLV_GSL;

  //! @brief state interval bounds
  T *_Ix;

  //! @brief parameter interval bounds -- do *NOT* free!
  const T *_Ip;

#ifdef MC__ODEBND_GSL_FI_USE
  //! @brief state bounds in forward AD
  FIT *_FIx;
#endif

  //! @brief mc::TModel environment
  TMT *_TMenv;

  //! @brief state Taylor model
  TVT *_TMx;

  //! @brief state derivative Taylor model
  TVT *_TMxdot;

  //! @brief state TM remainder lower bound time derivatives
  double *_RxLdot;

  //! @brief state TM remainder upper bound time derivatives
  double *_RxUdot;

  //! @brief parameter Taylor model -- do *NOT* free!
  const TVT *_TMp;

  //! @brief parameter reference
  double *_pref;

  //! @brief state reference
  double *_xref;

  //! @brief linear transformation A matrix
  double *_A;

  //! @brief linear transformation Z=A^{-1} matrix inverse
  double *_Z;

  //! @brief linear transformation B matrix
  double *_B;

  //! @brief linear transformed state interval bounds
  T *_Id;

  //! @brief linear transformed state ellipsoidal bounds
  E _Ed;

  //! @brief shape matrix (lower triangular) in ellipsoidal bounds
  double *_Q;

  //! @brief state reference time derivatives
  double *_xrefdot;

  //! @brief linear transformation A matrix time derivatives
  double *_Adot;

  //! @brief linear transformation Z=A^{-1} matrix inverse time derivatives
  double *_Zdot;

  //! @brief linear transformation B matrix time derivatives
  double *_Bdot;

  //! @brief rotated state interval bounds time derivatives
  T *_Iddot;

  //! @brief rotated state lower bound time derivatives
  double *_dLdot;

  //! @brief rotated state upper bound time derivatives
  double *_dUdot;

  //! @brief shape matrix time directives in ellipsoidal bounds
  double *_Qdot;

  //! @brief mc::TModel environment for mean-value theorem in (X,P)
  TMT *_MVXPenv;

  //! @brief rotated state Taylor model for mean-value theorem
  TVT *_MVXPd;

  //! @brief state Taylor model for mean-value theorem
  TVT *_MVXPx;

  //! @brief parameter Taylor model for mean-value theorem
  TVT *_MVXPp;

  //! @brief mc::TModel environment for linear transformation of intial value
  TMT *_MVPenv;

  //! @brief parameter Taylor model for linear transformation of initial value
  TVT *_MVPp;

public:
  /** @ingroup ODEBND
   *  @{
   */
  //! @brief Default constructor
  ODEBND_GSL
    ()
    : BASE_GSL(), ODESTRUCT(IVP())
    {
      // Initalize state/parameter
      _Ix = 0;
      _Ip = 0;
#ifdef MC__ODEBND_GSL_FI_USE
      _FIx = 0;
#endif
      _TMp = _TMx = _TMxdot = 0;
      _RxLdot = _RxUdot = 0;
      _TMenv = 0;

      // Initialize linear transformation
      _pref = _xref = _xrefdot = 0;
      _A = _Adot = _Z = _Zdot = _B = _Bdot = _Q = _Qdot = 0;
      _Id = _Iddot = 0;
      _dLdot = _dUdot = 0;

      // Initialize Taylor model
      _MVXPenv = new TMT( _nx+_np, options.ORDMIT );
      _MVXPenv->options.BOUNDER_TYPE = TMT::Options::HYBRID;
      _MVXPd = _MVXPx = _MVXPp = 0;
      _MVPenv = new TMT( _np, 1 );
      _MVPenv->options.BOUNDER_TYPE = TMT::Options::HYBRID;
      _MVPp = 0;

      // Initialize ellipsoidal calculus
      E::options.PSDCHK = false;
      
      // Initalize GSL
      _driver_dineqI = 0;
      _driver_dineqTM = 0;
      _ODESLV_GSL = new ODESLV_GSL<T,IVP>();
    }

  //! @brief Default destructor
  virtual ~ODEBND_GSL()
    {
      // Free state arrays -- Do *NOT* delete _p _Ip _TMp _TMenv
      delete[] _Ix;
#ifdef MC__ODEBND_GSL_FI_USE
      delete[] _FIx;
#endif
      delete[] _TMx;
      delete[] _TMxdot;
      delete[] _RxLdot;
      delete[] _RxUdot;
      
      // Free linear transformation arrays
      delete[] _pref;
      delete[] _xref;
      delete[] _A;
      delete[] _Z;
      delete[] _B;
      delete[] _Id;
      delete[] _Q;
      delete[] _xrefdot;
      delete[] _Adot;
      delete[] _Zdot;
      delete[] _Bdot;
      delete[] _Iddot;
      delete[] _dLdot;
      delete[] _dUdot;
      delete[] _Qdot;
      delete[] _MVXPd;
      delete[] _MVXPx;
      delete[] _MVXPp;
      delete _MVXPenv;
      delete[] _MVPp;
      delete _MVPenv;
      
      // Free GSL arrays
      if( _driver_dineqI )  gsl_odeiv2_driver_free( _driver_dineqI );
      if( _driver_dineqTM ) gsl_odeiv2_driver_free( _driver_dineqTM );
      delete _ODESLV_GSL;
    }

  //! @brief Integration results at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( const double tk, const unsigned int nxk, const T*Ixk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = Ixk[ix]; }
    Results
      ( const double tk, const unsigned int nxk, const TVT*TMxk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = TMxk[ix].B(); }
    Results
      ( const Results&res ):
      t( res.t ), nx( res.nx )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = res.X[ix]; }
    //! @brief Destructor
    ~Results()
      { delete[] X; }
    //! @brief Time instant
    double t;
    //! @brief Solution dimension
    unsigned int nx;
    //! @brief Solution bounds
    T* X;
  };

  //! @brief Integrator options
  struct Options: public BASE_GSL::Options
  {
    //! @brief Constructor
    Options():
      BASE_GSL::Options(), WRAPMIT(ELLIPS), ORDMIT(2), QTOL(machprec()), DISPLAY(0),
      RESRECORD(false)
      {}
    //! @brief Enumeration of wrapping mitigation strategies
    enum WRAPPING_STRATEGY{
      NONE=0,		//!< No wrapping mitigation
      DINEQ,		//!< Differential inequality contractor
      LINTRANS, 	//!< Linear preconditioning
      DITRANS,		//!< Differential inequality with linear preconditioning
      ELLIPS		//!< Ellipsoidal contractor with linear preconditioning [Default]
    };
    //! @brief Wrapping mitigation strategy
    WRAPPING_STRATEGY WRAPMIT;
    //! @brief Order of wrapping mitigation strategy (Default: 2)
    unsigned int ORDMIT;
    //! @brief Tolerance when dividing by trace of shape matrix in ellipsoidal bounds (Default: machprec())
    double QTOL;
    //! @brief Display level
    int DISPLAY;
    //! @brief Whether or not to record results (default: false)
    bool RESRECORD;
  } options;

  //! @brief Structure for setting up storing the solver exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for NLPSBB exception handling
    enum TYPE{
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    //! @brief Inline function returning the error description
    std::string what(){
      switch( _ierr ){
      case UNDEF: default:
        return "ODEBND_GSL::Exceptions  Error due to calling a not yet implemented feature";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Computes interval enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
      std::ostream&os=std::cout );

  //! @brief Computes Taylor model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned int ns, const double*tk, const TVT*TMp, TVT**TMxk,
      std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  STATUS hausdorff
    ( const unsigned int ns, const double*tk, const T*Ip, double**Hxk,
      const unsigned int nsamp, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distances between Taylor model enclosure and actual reachable set and between Taylor model remainder bound and actual remainder function range, using parameter sampling
  STATUS hausdorff
    ( const unsigned int ns, const double*tk, const TVT*TMp, double**HBxk,
      double**HRxk, const unsigned int nsamp, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned int iprec=5 ) const;

  //! @brief Return value of final time reached
  double final_time() const
    { return _t; };
  /** @} */

  //! @brief static pointer to class
  static ODEBND_GSL<T,IVP> *pODEBND_GSL;

private:

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > _results_IA;

  //! @brief Function to finalize statistics for GSL
  void _GSL_term();

  //! @brief Function to initialize GSL driver
  void _GSL_init
    ( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver );

  //! @brief Function to bound the remainder function relative to a Taylor model at sampling points -- return value is status
  STATUS _remainders
    ( const unsigned int ns, const double*tk, const T*Ip, const TVT*const*TMxk,
      T**Rxk, const unsigned int nsamp, std::ostream&os=std::cout );

  //! @brief Recrusive function computing bounds on errors between solutions of IVP in ODEs and truncated Taylor model using sampling
  STATUS _remainders
    ( const unsigned int ns, const double*tk, const T*Ip, const TVT*const*TMxk,
      T**Rxk, const unsigned int nsamp, unsigned int* vsamp,
      const unsigned int ip, double*p, double**xk, std::ostream&os );

  //! @brief Function to initialize GSL for state interval bounds
  void _GSL_init
    ( const T*Ip );

  //! @brief Function to initialize state interval bounds
  void _init
    ( T*Ix0 );

  //! @brief Function converting interval bounds to GSL array
  void _IA2vec
    ( const double*dL, const double*dU, const double*xref, const double*A,
      const double*Z, const double*B, double*vec );

  //! @brief Function converting interval bounds to GSL array
  void _IA2vec
    ( const T*Id, const double*xref, const double*A, const double*Z,
      const double*B, double*vec );

  //! @brief Function converting GSL array to interval bounds
  bool _vec2var
    ( const double*vec, T*Ix );

  //! @brief Function converting rotated variables into original coordinates
  template<typename U> void _dp2x
    ( const U*d, const U*p, U*x ) const;

  //! @brief Function converting rotated variables into original coordinates
  template<typename U> void _ep2x
    ( const double*scaling, const U*d, const U*p, U*x ) const;

  //! @brief Function converting ellipsoids into interval bounds
  bool _ep2x
    ( const double*Q, E&Ed, T*Id, const T*Ip, T*Ix ) const;

  //! @brief Static wrapper to function to calculate the DINEQs RHS values
  static int MC_GSLRHSI__
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Function to calculate the DINEQs RHS values
  int _ODEBND_GSL_RHSI
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Static wrapper to function to calculate the DINEQs RHS derivatives
  static int MC_GSLJACI__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to calculate the DINEQs RHS derivatives
  int _ODEBND_GSL_JACI
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to initialize GSL for state Taylor models
  void _GSL_init
    ( const TVT*TMp );

  //! @brief Function to initialize state Taylor model
  void _init
    ( TVT*TMx0 );

  //! @brief Function converting Taylor model to GSL array
  void _TM2vec
    ( const TVT*TM, const double*RxL, const double*RxU, const double*A,
      const double*Z, double*vec );

  //! @brief Function converting Taylor model to GSL array
  void _TM2vec
    ( const TVT*TM, const T*D, const double*A, const double*Z, double*vec );

  //! @brief Function converting Taylor model to GSL array
  void _TM2vec
    ( const TVT*TM, const double*A, const double*Z, double*vec );

  //! @brief Function converting GSL array to Taylor model
  bool _vec2var
    ( const double*vec, TVT*TM );

  //! @brief Function converting rotated Taylor remainder terms into original coordinates
  template<typename U> void _d2x
    ( const U*d, U*x, const bool reinit=true ) const;

  //! @brief Function converting rotated Taylor remainder terms into original coordinates
  template<typename U> void _e2x
    ( const double*scaling, const U*d, U*x, const bool reinit=true ) const;

  //! @brief Function converting Taylor remainder ellipsoids into interval remainder bounds
  bool _e2x
    ( const double*Q, E&Ed, T*Id, T*Ix ) const;

  //! @brief Static wrapper to function to calculate the DINEQ-TMs RHS values
  static int MC_GSLRHSTM__
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Function to calculate the DINEQ-TMs RHS values
  int _ODEBND_GSL_RHSTM
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Static wrapper to function to calculate the DINEQ-TMs RHS derivatives
  static int MC_GSLJACTM__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to calculate the DINEQ-TMs RHS derivatives
  int _ODEBND_GSL_JACTM
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm
    ( const double t, const U*x, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm
    ( const double t, const U*x1, const std::string&var1, const U*x2,
      const std::string&var2, std::ostream&os=std::cout ) const;

  //! @brief Position in symmetric matrix stored in lower triangular form
  unsigned int _ndxLT
    ( const unsigned int i, const unsigned int j ) const
    { return( i<j? _ndxLT(j,i): i+j*_nx-j*(j+1)/2 ); }

  //! @brief Private methods to block default compiler methods
  ODEBND_GSL(const ODEBND_GSL&);
  ODEBND_GSL& operator=(const ODEBND_GSL&);
};

template <typename T, typename IVP>
 ODEBND_GSL<T,IVP>* ODEBND_GSL<T,IVP>::pODEBND_GSL = 0;

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_GSL_init
( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver )
{
  // Set GSL driver
  if( driver ) gsl_odeiv2_driver_free( driver );
  switch( options.INTMETH ){
  case Options::RK8PD:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rk8pd,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::MSADAMS:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_msadams,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::MSBDF:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_msbdf,
      options.H0, options.ATOL, options.RTOL );
    break;
  case Options::RKF45: default:
    driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45,
      options.H0, options.ATOL, options.RTOL );
    break;
  }
  gsl_odeiv2_driver_set_hmin( driver, options.HMIN );  
  gsl_odeiv2_driver_set_nmax( driver, options.NMAX );  

  // Set state storage
  delete [] _vec_state;
  _vec_state = new double[ sys.dimension ];
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_GSL_term()
{
  // Get final CPU time
  _final_stats();
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_dp2x
( const U*d, const U*p, U*x ) const
{
  for( unsigned int ix=0; ix<_nx; ix++ ){
    x[ix] = _xref[ix];
    for( unsigned int jx=0; jx<_nx; jx++ )
      x[ix] += d[jx] * _A[jx*_nx+ix];
    for( unsigned int jp=0; jp<_np; jp++ )
      x[ix] += ( p[jp] - _pref[jp] ) * _B[jp*_nx+ix];
  }
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_ep2x
( const double*scaling, const U*d, const U*p, U*x ) const
{
  for( unsigned int ix=0; ix<_nx; ix++ ){   
    x[ix] = _xref[ix];
    if( scaling )
      for( unsigned int jx=0; jx<_nx; jx++ )
        x[ix] += scaling[ix+jx*_nx] * d[jx];
    else
      x[ix] += d[ix];
    for( unsigned int jp=0; jp<_np; jp++ )
      x[ix] += ( p[jp] - _pref[jp] ) * _B[jp*_nx+ix];
  }
  return;
}

template <typename T, typename IVP> inline bool
ODEBND_GSL<T,IVP>::_ep2x
( const double*Q, E&Ed, T*Id, const T*Ip, T*Ix ) const
{
  Ed.set( _nx, _Q );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
  std::cout << "Ed =" << Ed << std::endl;
#endif

#ifdef MC__ODEBND_GSL_ELLTR_SCALED
  if( !Ed.sqrtQ(true).n ) return false;
  for( unsigned int ix=0; ix<_nx; ix++ )
    Id[ix] = T( -1., 1. );
  _ep2x( Ed.sqrtQ(true).array, Id, Ip, Ix );
#else
  for( unsigned int ix=0; ix<_nx; ix++ )
    Id[ix] = T( Ed.l(ix), Ed.u(ix) );
  _ep2x( 0, Id, Ip, Ix );
#endif

#ifdef MC__ODEBND_GSL_LINTR_DEBUG
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "Ix[" << ix << "] = " << Ix[ix] << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return true;
}

template <typename T, typename IVP> inline bool
ODEBND_GSL<T,IVP>::_vec2var
( const double*vec, T*Ix )
{
  unsigned int ivec=0;
  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      Ix[ix] = T( vec[ivec], vec[ivec+1] );
      ivec+=2;
    }
    return true;

  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Id[ix] = T( vec[ivec], vec[ivec+1] );
      ivec+=2;
    }
    for( unsigned int ix=0; ix<_nx; ix++ )
      _xref[ix] = vec[ivec++];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      _A[iA] = vec[ivec++];
    for( unsigned int iZ=0; iZ<_nx*_nx; iZ++ )
      _Z[iZ] = vec[ivec++];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      _B[iB] = vec[ivec++];
    _dp2x( _Id, _Ip, Ix );
    return true;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ )
      _xref[ix] = vec[ivec++];
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      _Q[iQ] = vec[ivec++];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      _B[iB] = vec[ivec++];
    return _ep2x( _Q, _Ed, _Id, _Ip, Ix );

  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_IA2vec
( const double*dL, const double*dU, const double*xref, const double*A,
  const double*Z, const double*B, double*vec )
{
  unsigned int ivec=0;

  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = dL[ix];
      vec[ivec++] = dU[ix];
    }
    if( options.WRAPMIT == Options::NONE ||
        options.WRAPMIT == Options::DINEQ ) return;
    for( unsigned int ix=0; ix<_nx; ix++ )
      vec[ivec++] = xref[ix];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = A[iA];
    for( unsigned int iZ=0; iZ<_nx*_nx; iZ++ )
      vec[ivec++] = Z[iZ];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      vec[ivec++] = B[iB];
    return;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ )
      vec[ivec++] = xref[ix];
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      vec[ivec++] = A[iQ];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      vec[ivec++] = B[iB];
    if( options.WRAPMIT == Options::ELLIPS ) return;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = dL[ix];
      vec[ivec++] = dU[ix];
    }
    return;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_IA2vec
( const T*Id, const double*xref, const double*A, const double*Z,
  const double*B, double*vec )
{
  unsigned int ivec=0;

  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = Op<T>::l( Id[ix] );
      vec[ivec++] = Op<T>::u( Id[ix] );
    }
    if( options.WRAPMIT == Options::NONE ||
        options.WRAPMIT == Options::DINEQ ) return;
    for( unsigned int ix=0; ix<_nx; ix++ )
      vec[ivec++] = xref[ix];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = A[iA];
    for( unsigned int iZ=0; iZ<_nx*_nx; iZ++ )
      vec[ivec++] = Z[iZ];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      vec[ivec++] = B[iB];
    return;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ )
      vec[ivec++] = xref[ix];
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      vec[ivec++] = A[iQ];
    for( unsigned int iB=0; iB<_nx*_np; iB++ )
      vec[ivec++] = B[iB];
    if( options.WRAPMIT == Options::ELLIPS ) return;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = Op<T>::l( Id[ix] );
      vec[ivec++] = Op<T>::u( Id[ix] );
    }
    return;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_init
( T*Ix0 )
{
  switch( options.WRAPMIT){

  case Options::NONE:
  case Options::DINEQ:
    for( unsigned int ix=0; ix<_nx; ix++ )
      Ix0[ix] = IVP().IC( ix, _Ip );
    return _IA2vec( Ix0, 0, 0, 0, 0, _vec_state );

  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ip=0; ip<_np; ip++ ){
      _pref[ip] = Op<T>::mid( _Ip[ip] );
      _MVPp[ip].set( _MVPenv, ip, _Ip[ip] );
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      TVT ICx0 = IVP().IC( ix, _MVPp ).center();
      _xref[ix] = ICx0.constant();
      for( unsigned int jx=0; jx<_nx; jx++ )
        _A[jx*_nx+ix] = _Z[jx*_nx+ix] = (ix==jx? 1.: 0.);
      double*ICx0_linear = ICx0.linear();
      for( unsigned int jp=0; jp<_np; jp++ )
        _B[jp*_nx+ix] = ( ICx0_linear? ICx0_linear[jp]: 0. );
      delete[] ICx0_linear;
      _Id[ix] = ICx0.remainder();
      Ix0[ix] = ICx0.bound();
    }
    return _IA2vec( _Id, _xref, _A, _Z, _B, _vec_state );

  case Options::ELLIPS:
  default:
    for( unsigned int ip=0; ip<_np; ip++ ){
      _pref[ip] = Op<T>::mid( _Ip[ip] );
      _MVPp[ip].set( _MVPenv, ip, _Ip[ip] );
    }
    double trR2 = 0.;
    for( unsigned int ix=0, iQ=0; ix<_nx; ix++ ){
      TVT ICx0 = IVP().IC( ix, _MVPp ).center();
      _xref[ix] = ICx0.constant();
      trR2 += sqr( Op<T>::diam( ICx0.remainder() ) / 2. );
      _Q[iQ++] = Op<T>::diam( ICx0.remainder() ) / 2.;
      for( unsigned int jx=ix+1; jx<_nx; jx++ ) _Q[iQ++] = 0.;
      double*ICx0_linear = ICx0.linear();
      for( unsigned int jp=0; jp<_np; jp++ )
        _B[jp*_nx+ix] = ( ICx0_linear? ICx0_linear[jp]: 0. );
      delete[] ICx0_linear;
      Ix0[ix] = ICx0.bound();
    }
    for( unsigned int ix=0, iQ=0; ix<_nx; ix++, iQ+=_nx-ix )
      _Q[iQ] *= std::sqrt( trR2 );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    E Ex0( _nx, _Q, _xref );
    std::cout << "Ex0 =" << Ex0 << std::endl;
#endif
    return _IA2vec( Ix0, _xref, _Q, 0, _B, _vec_state );
  }
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::_ODEBND_GSL_RHSI
( double t, const double* x, double* xdot, void* user_data )
{
  if( !_vec2var( x, _Ix ) ) return GSL_EBADFUNC;
  stats.numRHS++;

  switch( options.WRAPMIT){
  case Options::NONE:{
    unsigned int ivec=0;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      T Ixdot = IVP().RHS( ix, _Ip, _Ix, t, _istg );
      xdot[ivec++] = Op<T>::l( Ixdot );
      xdot[ivec++] = Op<T>::u( Ixdot );
    }
    return GSL_SUCCESS;  
   }
   
  case Options::DINEQ:{
    unsigned int ivec=0;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      T Ixi = _Ix[ix];
      _Ix[ix] = Op<T>::l( Ixi );
      xdot[ivec++] = Op<T>::l( IVP().RHS( ix, _Ip, _Ix, t, _istg ) );
      _Ix[ix] = Op<T>::u( Ixi );
      xdot[ivec++] = Op<T>::u( IVP().RHS( ix, _Ip, _Ix, t, _istg ) );
      _Ix[ix] = Ixi;
    }
    return GSL_SUCCESS;  
   }
   
  case Options::LINTRANS:
  case Options::DITRANS:{
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    std::cout << "@t=" << t << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "d[" << ix << "] = " << _Id[ix]
                << " : " << Op<T>::mid(_Id[ix]) << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "A[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _A[jx*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "Z[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _Z[jx*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    { int dum; std::cin >> dum; }
#endif

    // Compute differential inequality bound derivatives
    _MVXPenv->options.PROPAGATE_BNDT = false;
    _MVXPenv->options.SCALE_VARIABLES = false;
    _MVXPenv->options.CENTER_REMAINDER = false;
    _MVXPenv->options.REF_MIDPOINT = false;
    for( unsigned int jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Id[jx], 0. );
    for( unsigned int jp=0; jp<_np; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
    _dp2x( _MVXPd, _MVXPp, _MVXPx );
    for( unsigned int jx=0; jx<_nx; jx++ )
      _Iddot[jx] = 0.;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      TVT fi = IVP().RHS( ix, _MVXPp, _MVXPx, t, _istg );
      _xrefdot[ix] = fi.constant();
      for( unsigned int jx=0; jx<_nx; jx++ )
        _Adot[jx*_nx+ix] = fi.linear(jx);
      for( unsigned int jp=0; jp<_np; jp++ )
        _Bdot[jp*_nx+ix] = fi.linear(_nx+jp);
      T Rfi = fi.remainder();
      for( unsigned int iord=2; iord<=_MVXPenv->nord(); iord++ ){
        Rfi += fi.bound( iord );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "Rfi =" << fi.remainder() << std::endl;
        std::cout << "Rhi =" << Rfi << std::endl;
	{ int dum; std::cin >> dum; }
#endif
      }
      for( unsigned int jx=0; jx<_nx; jx++ )
        _Iddot[jx] += _Z[ix*_nx+jx] * Rfi;
      
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
      std::cout << "Adot[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _Adot[jx*_nx+ix] << "  ";
      std::cout << std::endl;
      std::cout << "Bdot[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _Bdot[jx*_nx+ix] << "  ";
      std::cout << std::endl;
#endif
    }

    // Compute Zdot = -Z*Adot*Z
    double* ZAdot = new double[_nx*_nx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    double*AZ = new double[_nx*_nx];
#endif
    for( unsigned int ix=0; ix<_nx; ix++ ){
      for( unsigned int jx=0; jx<_nx; jx++ ){
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        AZ[jx*_nx+ix] = 0.;
#endif
        ZAdot[jx*_nx+ix] = 0.;
        for( unsigned int kx=0; kx<_nx; kx++ ){
	  ZAdot[jx*_nx+ix] += _Z[kx*_nx+ix] * _Adot[jx*_nx+kx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
          AZ[jx*_nx+ix] += _Z[kx*_nx+ix] * _A[jx*_nx+kx];
#endif
        }
      }
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      for( unsigned int jx=0; jx<_nx; jx++ ){
        _Zdot[jx*_nx+ix] = 0.;
        for( unsigned int kx=0; kx<_nx; kx++ )
	  _Zdot[jx*_nx+ix] -= ZAdot[kx*_nx+ix] * _Z[jx*_nx+kx];
      }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "A*Z[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << AZ[jx*_nx+ix] << "  ";
      std::cout << std::endl;
#endif
    }
    delete[] ZAdot;
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    delete[] AZ;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "Adot[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _Adot[jx*_nx+ix] << "  ";
      std::cout << std::endl;
    }
#endif

    // Compute differential inequality bounds
    if( options.WRAPMIT == Options::LINTRANS ){
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _dLdot[ix] = Op<T>::l( _Iddot[ix] );
        _dUdot[ix] = Op<T>::u( _Iddot[ix] );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "dLdot[" << ix << "] = " << _dLdot[ix] << std::endl;
        std::cout << "dUdot[" << ix << "] = " << _dLdot[ix] << std::endl;
#endif
      }
    }

    else{
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _MVXPd[ix].set( _MVXPenv, ix, Op<T>::l( _Id[ix] ), 0. );
        _dp2x( _MVXPd, _MVXPp, _MVXPx );
        _Iddot[ix] = 0.;
        for( unsigned int jx=0; jx<_nx; jx++ ){
          TVT fiL = IVP().RHS( jx, _MVXPp, _MVXPx, t, _istg );
          T RfiL = fiL.remainder();
          for( unsigned int iord=2; iord<=_MVXPenv->nord(); iord++ )
            RfiL += fiL.bound( iord );
          _Iddot[ix] += _Z[jx*_nx+ix] * RfiL;
        }
        _dLdot[ix] = Op<T>::l( _Iddot[ix] );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "dLdot[" << ix << "] = " << _dLdot[ix] << std::endl;
#endif

        _MVXPd[ix].set( _MVXPenv, ix, Op<T>::u( _Id[ix] ), 0. );
        _dp2x( _MVXPd, _MVXPp, _MVXPx );
        _Iddot[ix] = 0.;
        for( unsigned int jx=0; jx<_nx; jx++ ){
          TVT fiU = IVP().RHS( jx, _MVXPp, _MVXPx, t, _istg );
          T RfiU = fiU.remainder();
          for( unsigned int iord=2; iord<=_MVXPenv->nord(); iord++ )
            RfiU += fiU.bound( iord );
          _Iddot[ix] += _Z[jx*_nx+ix] * RfiU;
        }
        _dUdot[ix] = Op<T>::u( _Iddot[ix] );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "dUdot[" << ix << "] = " << _dUdot[ix] << std::endl;
#endif

        _MVXPd[ix].set( _MVXPenv, ix, _Id[ix], 0. );
        _dp2x( _MVXPd, _MVXPp, _MVXPx );
      }
    }

#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    //_dp2x( _MVXPd, _MVXPp, _MVXPx );
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "X[" << ix << "] = " << _MVXPx[ix].B()
                << " : " << Op<T>::mid(_MVXPx[ix].B()) << std::endl;
    {int dum; std::cin >> dum;}
#endif

    // Copy time derivative in GSL array
    _IA2vec( _dLdot, _dUdot, _xrefdot, _Adot, _Zdot, _Bdot, xdot );
    
    return GSL_SUCCESS;
   }
   
  case Options::ELLIPS:{
  default:
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    std::cout << "@t=" << t << std::endl;
    //E::options.PSDCHK = true;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xref[" << ix << "] = " << _xref[ix] << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "B[" << ix << ",#] = ";
      for( unsigned int jp=0; jp<_np; jp++ )
        std::cout << _B[jp*_np+ix] << "  ";
      std::cout << std::endl;
    }
    std::cout << "Ed = " << _Ed << std::endl;
    { int dum; std::cin >> dum; }
#endif

    // Compute differential inequality bound derivatives
    _MVXPenv->options.PROPAGATE_BNDT = false;
    _MVXPenv->options.SCALE_VARIABLES = true;
    _MVXPenv->options.CENTER_REMAINDER = true;
    _MVXPenv->options.REF_MIDPOINT = true;
    for( unsigned int jx=0; jx<_nx; jx++ )
      _MVXPd[jx].set( _MVXPenv, jx, _Id[jx] );
    for( unsigned int jp=0; jp<_np; jp++ )
      _MVXPp[jp].set( _MVXPenv, _nx+jp, _Ip[jp] );
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
    _ep2x( _Ed.sqrtQ(true).array, _MVXPd, _MVXPp, _MVXPx );
#else
    _ep2x( 0, _MVXPd, _MVXPp, _MVXPx );
#endif
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
    { int dum; std::cin >> dum; }
#endif

    double trQ = 0., sumkappa = 0.;
    for( unsigned int ix=0; ix<_nx; ix++ )
      trQ += ( _Q[_ndxLT(ix,ix)]>0? _Q[_ndxLT(ix,ix)]: machprec() );
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
    const double* sqrtQ = _Ed.sqrtQ(true).array;
#endif
     
    for( unsigned int ix=0; ix<_nx; ix++ ){
      TVT TMfi = IVP().RHS( ix, _MVXPp, _MVXPx, t, _istg );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "MVXPf[" << ix << "] = " << TMfi << std::endl;
      { int dum; std::cin >> dum; }
#endif
      _xrefdot[ix] = TMfi.constant();
      _Iddot[ix] = TMfi.remainder();
      for( unsigned int iord=2; iord<=_MVXPenv->nord(); iord++ ){
        _xrefdot[ix] += Op<T>::mid( TMfi.bound( iord ) );
        _Iddot[ix] += TMfi.bound( iord ) - Op<T>::mid( TMfi.bound( iord ) );
      }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
      std::cout << "Iddot[" << ix << "] = " << _Iddot[ix]
                << " : " << Op<T>::mid(_Iddot[ix]) << std::endl;
      std::cout << "kappa[" << ix << "] = "
                << Op<T>::diam( _Iddot[ix] ) / 2.
                 / ( std::sqrt( trQ ) + options.QTOL ) << std::endl;
#endif
      sumkappa += ( Op<T>::diam( _Iddot[ix] ) / 2. )
                / ( std::sqrt( trQ ) + options.QTOL );
      for( unsigned int jx=0; jx<_nx; jx++ )
        _A[ix+jx*_nx] = TMfi.linear(jx);
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "A[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _A[ix+jx*_nx] << "  ";
      std::cout << std::endl;
#endif
      for( unsigned int jp=0; jp<_np; jp++ )
        _Bdot[ix+jp*_nx] = TMfi.linear(_nx+jp);
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      std::cout << "Bdot[" << ix << ",#] = ";
      for( unsigned int jp=0; jp<_np; jp++ )
        std::cout << _Bdot[ix+jp*_nx] << "  ";
      std::cout << std::endl;
#endif
    }

    for( unsigned int jx=0; jx<_nx; jx++ ){
      for( unsigned int ix=jx; ix<_nx; ix++ ){
        _Qdot[_ndxLT(ix,jx)] = sumkappa * _Q[_ndxLT(ix,jx)];
        for( unsigned int kx=0; kx<_nx; kx++ )
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
          _Qdot[_ndxLT(ix,jx)] += sqrtQ[kx+ix*_nx] * _A[jx+kx*_nx]
	                        + _A[ix+kx*_nx] * sqrtQ[kx+jx*_nx];
#else
          _Qdot[_ndxLT(ix,jx)] += _Q[_ndxLT(ix,kx)] * _A[jx+kx*_nx]
	                        + _A[ix+kx*_nx] * _Q[_ndxLT(kx,jx)];
#endif
      }
      _Qdot[_ndxLT(jx,jx)] += ( Op<T>::diam( _Iddot[jx] ) / 2. )
                            * ( std::sqrt( trQ ) + options.QTOL );
    }
     
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xrefdot[" << ix << "] = " << _xrefdot[ix] << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "Qdot[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<=ix; jx++ )
        std::cout << _Qdot[_ndxLT(ix,jx)] << "  ";
      std::cout << std::endl;
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "Bdot[" << ix << ",#] = ";
      for( unsigned int jp=0; jp<_np; jp++ )
        std::cout << _Bdot[jp*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    E Exdot( _nx, _Qdot, _xrefdot );
    std::cout << "Exdot =" << Exdot << std::endl;
#endif

    // Copy time derivative in GSL array
    _IA2vec( _dLdot, _dUdot, _xrefdot, _Qdot, 0, _Bdot, xdot );
    
    return GSL_SUCCESS;
   }
  }
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::MC_GSLRHSI__
( double t, const double* x, double* xdot, void* user_data )
{
  ODEBND_GSL<T,IVP> *pODEBND_GSL = ODEBND_GSL<T,IVP>::pODEBND_GSL;
  int flag = pODEBND_GSL->_ODEBND_GSL_RHSI( t, x, xdot, user_data );
  ODEBND_GSL<T,IVP>::pODEBND_GSL = pODEBND_GSL;
  return flag;
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::_ODEBND_GSL_JACI
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  stats.numJAC++;
  switch( options.WRAPMIT){
#ifdef MC__ODEBND_GSL_FI_USE
  case Options::NONE:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _FIx[ix] = T( x[ix], x[ix+_nx] );
      _FIx[ix].sen(_nx,ix);
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      FIT RHS = IVP().RHS( ix, _Ip, _FIx, t, _istg );
      xdot[ix]     = RHS.l();
      xdot[ix+_nx] = RHS.u();
      for( unsigned int jx=0; jx<_nx; jx++ ){
        jac[ix*2*_nx+jx]           = RHS.lbsenlb()? RHS.lbsenlb(jx): 0.;
        jac[ix*2*_nx+_nx+jx]       = RHS.lbsenub()? RHS.lbsenub(jx): 0.;
        jac[(_nx+ix)*2*_nx+jx]     = RHS.ubsenlb()? RHS.ubsenlb(jx): 0.;
        jac[(_nx+ix)*2*_nx+_nx+jx] = RHS.ubsenub()? RHS.ubsenub(jx): 0.;
      }
    }
    return GSL_SUCCESS;  
#endif
  
  case Options::DINEQ:
#ifdef MC__ODEBND_GSL_FI_USE
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _FIx[ix] = T( x[ix], x[ix+_nx] );
      _FIx[ix].sen(_nx,ix);
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _FIx[ix] = x[ix];
      _FIx[ix].sen(_nx,ix);
      FIT RHSL = IVP().RHS( ix, _Ip, _FIx, t, _istg );
      xdot[ix] = RHSL.l();
      for( unsigned int jx=0; jx<_nx; jx++ ){
        jac[ix*2*_nx+jx]     = RHSL.lbsenlb()? RHSL.lbsenlb(jx): 0.;
        jac[ix*2*_nx+_nx+jx] = RHSL.lbsenub()? RHSL.lbsenub(jx): 0.;
      }
      _FIx[ix] = x[_nx+ix];
      _FIx[ix].sen(_nx,ix);
      FIT RHSU = IVP().RHS( ix, _Ip, _FIx, t, _istg );
      xdot[ix+_nx] = RHSU.u();
      for( unsigned int jx=0; jx<_nx; jx++ ){
        jac[(_nx+ix)*2*_nx+jx]     = RHSU.ubsenlb()? RHSU.ubsenlb(jx): 0.;
        jac[(_nx+ix)*2*_nx+_nx+jx] = RHSU.ubsenub()? RHSU.ubsenub(jx): 0.;
      }
      _FIx[ix] = T( x[ix], x[ix+_nx] );
      _FIx[ix].sen(_nx,ix);
    }
    return GSL_SUCCESS;
#endif

  case Options::LINTRANS:
  case Options::DITRANS:
  case Options::ELLIPS:
  default:
    // Jacobian not yet implemented for certain bounders
    return GSL_EBADFUNC;
  }
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::MC_GSLJACI__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODEBND_GSL<T,IVP> *pODEBND_GSL = ODEBND_GSL<T,IVP>::pODEBND_GSL;
  int flag = pODEBND_GSL->_ODEBND_GSL_JACI( t, x, jac, xdot, user_data );
  ODEBND_GSL<T,IVP>::pODEBND_GSL = pODEBND_GSL;
  return flag;
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_GSL_init
( const T* Ip )
{
  // Define ODE system in GSL format
  _sys_dineqI.function = MC_GSLRHSI__;
  _sys_dineqI.jacobian = MC_GSLJACI__;
  _sys_dineqI.dimension = 0;
  switch( options.WRAPMIT){
  case Options::LINTRANS:
  case Options::DITRANS:
    _sys_dineqI.dimension += _nx*(1+2*_nx+_np);
  case Options::NONE:
  case Options::DINEQ:
    _sys_dineqI.dimension += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    _sys_dineqI.dimension += _nx*(1+_np)+_nx*(_nx+1)/2;
    break;
  }
  _sys_dineqI.params = 0;

  // Set GSL driver
  _GSL_init( _sys_dineqI, _driver_dineqI );

  // Internal Taylor model environment reset
  if( _MVXPenv && ( _MVXPenv->nord() != options.ORDMIT 
                 || _MVXPenv->nvar() != _nx+_np ) ){
    delete[] _MVXPx; _MVXPx = 0;
    delete[] _MVXPd; _MVXPd = 0;
    delete[] _MVXPp; _MVXPp = 0;
    delete _MVXPenv; _MVXPenv = 0;
  }

  // State/parameter arrays
  _Ip = Ip;
  if( !_Ix )  _Ix = new T[_nx];
#ifdef MC__ODEBND_GSL_FI_USE
  if( !_FIx ) _FIx = new FIT[_nx];
#endif

  switch( options.WRAPMIT){
  case Options::ELLIPS:
  default:
    if( !_pref )    _pref    = new double[_np];
    if( !_xref )    _xref    = new double[_nx];
    if( !_A )	    _A       = new double[_nx*_nx];
    if( !_B )	    _B       = new double[_nx*_np];
    if( !_Q )	    _Q       = new double[_nx*(_nx+1)/2];
    if( !_Id )      _Id      = new T[_nx];
    if( !_xrefdot ) _xrefdot = new double[_nx];
    if( !_Bdot )    _Bdot    = new double[_nx*_np];
    if( !_Qdot )    _Qdot    = new double[_nx*(_nx+1)/2];
    if( !_Iddot )   _Iddot   = new T[_nx];
    if( !_dLdot )   _dLdot   = new double[_nx];
    if( !_dUdot )   _dUdot   = new double[_nx];
    if( !_MVXPenv ) _MVXPenv = new TMT( _nx+_np, options.ORDMIT );
    if( !_MVXPx )   _MVXPx   = new TVT[_nx];
    if( !_MVXPd )   _MVXPd   = new TVT[_nx];
    if( !_MVXPp )   _MVXPp   = new TVT[_np];
    _MVXPenv->options.BOUNDER_TYPE = TMT::Options::HYBRID;
    if( !_MVPp )    _MVPp    = new TVT[_np];
    break;

  case Options::LINTRANS:
  case Options::DITRANS:
    if( !_pref )    _pref    = new double[_np];
    if( !_xref )    _xref    = new double[_nx];
    if( !_A )	    _A       = new double[_nx*_nx];
    if( !_Z )	    _Z       = new double[_nx*_nx];
    if( !_B )	    _B       = new double[_nx*_np];
    if( !_Id )      _Id      = new T[_nx];
    if( !_xrefdot ) _xrefdot = new double[_nx];
    if( !_Adot )    _Adot    = new double[_nx*_nx];
    if( !_Zdot )    _Zdot    = new double[_nx*_nx];
    if( !_Bdot )    _Bdot    = new double[_nx*_np];
    if( !_Iddot )   _Iddot   = new T[_nx];
    if( !_dLdot )   _dLdot   = new double[_nx];
    if( !_dUdot )   _dUdot   = new double[_nx];
    if( !_MVXPenv ) _MVXPenv = new TMT( _nx+_np, options.ORDMIT );
    if( !_MVXPx )   _MVXPx   = new TVT[_nx];
    if( !_MVXPd )   _MVXPd   = new TVT[_nx];
    if( !_MVXPp )   _MVXPp   = new TVT[_np];
    _MVXPenv->options.BOUNDER_TYPE = TMT::Options::HYBRID;
    if( !_MVPp )    _MVPp    = new TVT[_np];
    break;

  case Options::NONE:
  case Options::DINEQ:
    break;
  }

  // Reset result record and statistics
  _results_IA.clear();  
  _init_stats();
}

//! @fn template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS ODEBND_GSL<T,IVP>::bounds(
//! const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
//! std::ostream&os )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] interval state enclosures at stage times
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP>
inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::bounds
( const unsigned int ns, const double*tk, const T*Ip, T**Ixk, std::ostream&os )
{
  // Initialize trajectory integration with GSL
  _GSL_init( Ip );
  _t = tk[0];
  _init( Ixk[0] );
  if( options.DISPLAY >= 1 ) _print_interm( tk[0], Ixk[0], "x", os );
  // Record initial results
  if( options.RESRECORD )
    _results_IA.push_back( Results( tk[0], _nx, Ixk[0] ) );

  try{
    // Integrate ODEs through each stage using GSL
    pODEBND_GSL = this;
    for( _istg=0; _istg<ns; _istg++ ){
      if( gsl_odeiv2_driver_apply( _driver_dineqI, &_t, tk[_istg+1], _vec_state )
        != GSL_SUCCESS ){
	_GSL_term();
        if( options.DISPLAY >= 1 ) _print_stats( stats, os );
	return FAILURE;
      }
      _vec2var( _vec_state, Ixk[_istg+1] );
      if( options.DISPLAY >= 1 ) _print_interm( tk[_istg+1], Ixk[_istg+1], "x", os );
      // Record intermediate results
      if( options.RESRECORD )
        _results_IA.push_back( Results( tk[_istg+1], _nx, Ixk[_istg+1] ) );
    }
  }
  catch(...){
    _GSL_term();
    if( options.DISPLAY >= 1 ) _print_stats( stats, os );
    return FAILURE;
  }

  _GSL_term();
  if( options.DISPLAY >= 1 ) _print_stats( stats, os );
  return NORMAL;
}

//! @fn template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS ODEBND_GSL<T,IVP>::hausdorff(
//! const unsigned int ns, const double*tk, const T*Ip, double**Hxk,
//! const unsigned int nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the interval enclosure
//! and the exact reachable set projected onto each variable:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Hxk</a> [output] Hausdorff distance between the interval enclosure
//!     and the exact reachable set projected onto each variable, at stage times
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::hausdorff
( const unsigned int ns, const double*tk, const T*Ip, double**Hxk,
  const unsigned int nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  _ODESLV_GSL->options = options;
  _ODESLV_GSL->options.RESRECORD = false;

  // Compute exact bounds 
  T** Ixk0 = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Ixk0[is] = new T[_nx];
  if( _ODESLV_GSL->bounds( ns, tk, Ip, Ixk0, nsamp, os ) != NORMAL ){
    for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
    delete[] Ixk0;
    return FAILURE;
  }

  // Compute approximate bounds
  T** Ixk = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Ixk[is] = new T[_nx];
  try{ bounds( ns, tk, Ip, Ixk, os ); }
  catch(...){;}
  unsigned int nsf = _istg;
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Ixk[_istg+1][ix] = T( -1e20, 1e20 );

  // Compute Hausdorff metric at each time step
  struct loc{ static double hausdorff
    ( const T&Ix, const T&Ix0 )
    { return std::max( std::fabs(Op<T>::l(Ix)-Op<T>::l(Ix0)),
                       std::fabs(Op<T>::u(Ix)-Op<T>::u(Ix0)) ); }
  };

  options.DISPLAY = DISPLAY_SAVE;
  for( unsigned int is=0; is<ns+1; is++ ){
    for( unsigned int ix=0; ix<_nx; ix++ )
      Hxk[is][ix] = loc::hausdorff( Ixk[is][ix], Ixk0[is][ix] );
    if( options.DISPLAY >= 1 ) _print_interm( tk[is], Hxk[is], "d", os );
  }

  for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
  delete[] Ixk0;
  for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk[is];
  delete[] Ixk;

  return NORMAL;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_d2x
( const U*d, U*x, const bool reinit ) const
{
  for( unsigned int ix=0; ix<_nx; ix++ ){
   if( reinit ) x[ix] = 0.;
    for( unsigned int jx=0; jx<_nx; jx++ )
      x[ix] += d[jx] * _A[jx*_nx+ix];
  }
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_e2x
( const double*scaling, const U*d, U*x, const bool reinit ) const
{
  for( unsigned int ix=0; ix<_nx; ix++ ){   
    if( reinit ) x[ix] = 0.;
    if( scaling )
      for( unsigned int jx=0; jx<_nx; jx++ )
        x[ix] += scaling[ix+jx*_nx] * d[jx];
    else
      x[ix] += d[ix];
  }
  return;
}

template <typename T, typename IVP> inline bool
ODEBND_GSL<T,IVP>::_e2x
( const double*Q, E&Ed, T*Id, T*Ix ) const
{
  Ed.set( _nx, _Q );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
  std::cout << "Ed =" << Ed << std::endl;
#endif

#ifdef MC__ODEBND_GSL_ELLTR_SCALED
  if( !Ed.sqrtQ(true).n ) return false;
  for( unsigned int ix=0; ix<_nx; ix++ )
    Id[ix] = T( -1., 1. );
  _e2x( Ed.sqrtQ(true).array, Id, Ix );
#else
  for( unsigned int ix=0; ix<_nx; ix++ )
    Id[ix] = T( Ed.l(ix), Ed.u(ix) );
  _e2x( 0, Id, Ix );
#endif

#ifdef MC__ODEBND_GSL_LINTR_DEBUG
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "Ix[" << ix << "] = " << Ix[ix] << std::endl;
  { int dum; std::cin >> dum; }
#endif
  return true;
}

template <typename T, typename IVP> inline bool
ODEBND_GSL<T,IVP>::_vec2var
( const double*vec, TVT*TMx )
{
  unsigned int ivec=0;
  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
    for( unsigned int ix=0; ix<_nx; ix++  ){
      TMx[ix].set( _TMenv );
      TMx[ix].set( vec+ivec );
      ivec+=_TMenv->nmon();
      TMx[ix].set( T( vec[ivec], vec[ivec+1] ) );
      ivec+=2;
    }
    return true;

  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ivec+=_TMenv->nmon()+2, ix++  ){
      TMx[ix].set( _TMenv );
      TMx[ix].set( vec+ivec );
      _Id[ix] = T( vec[ivec+_TMenv->nmon()], vec[ivec+_TMenv->nmon()+1] );
    }
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      _A[iA] = vec[ivec++];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      _Z[iA] = vec[ivec++];
    _d2x( _Id, _Ix );
    for( unsigned int ix=0; ix<_nx; ix++  )
      TMx[ix].set( _Ix[ix] );
    return true;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++  ){
      TMx[ix].set( _TMenv );
      TMx[ix].set( vec+ivec );
      ivec+=_TMenv->nmon();
    }
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      _Q[iQ] = vec[ivec++];
    if( !_e2x( _Q, _Ed, _Id, _Ix ) ) return false;
    for( unsigned int ix=0; ix<_nx; ix++  )
      TMx[ix].set( _Ix[ix] );
    return true;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_TM2vec
( const TVT*TM, const double*RxL, const double*RxU, const double*A,
  const double*Z, double*vec )
{
  unsigned int ivec=0;

  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
      vec[ivec++] = RxL[ix];
      vec[ivec++] = RxU[ix];
    }
    if( options.WRAPMIT == Options::NONE ||
        options.WRAPMIT == Options::DINEQ ) return;
  
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = A[iA];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = Z[iA];
    return;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
    }
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      vec[ivec++] = A[iQ];
    if( options.WRAPMIT == Options::ELLIPS ) return;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = RxL[ix];
      vec[ivec++] = RxU[ix];
    }
    return;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_TM2vec
( const TVT*TM, const T*D, const double*A, const double*Z, double*vec )
{
  unsigned int ivec=0;

  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
      vec[ivec++] = Op<T>::l( D[ix] );
      vec[ivec++] = Op<T>::u( D[ix] );
    }
    if( options.WRAPMIT == Options::NONE ||
        options.WRAPMIT == Options::DINEQ ) return;
  
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = A[iA];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = Z[iA];
    return;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
    }
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      vec[ivec++] = A[iQ];
    if( options.WRAPMIT == Options::ELLIPS ) return;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = Op<T>::l( D[ix] );
      vec[ivec++] = Op<T>::u( D[ix] );
    }
    return;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_TM2vec
( const TVT*TM, const double*A, const double*Z, double*vec )
{
  unsigned int ivec=0;

  switch( options.WRAPMIT ){
  case Options::NONE:
  case Options::DINEQ:
  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
      vec[ivec++] = Op<T>::l( TM[ix].remainder() );
      vec[ivec++] = Op<T>::u( TM[ix].remainder() );
    }
    if( options.WRAPMIT == Options::NONE ||
        options.WRAPMIT == Options::DINEQ ) return;
  
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = A[iA];
    for( unsigned int iA=0; iA<_nx*_nx; iA++ )
      vec[ivec++] = Z[iA];
    return;

  case Options::ELLIPS:
  default:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::pair<unsigned int, const double*> TMcoef = TM[ix].coefmon();
      unsigned int imon = 0;
      for( ; imon<TMcoef.first; imon++ )
        vec[ivec++] = TMcoef.second[imon];
      for( ; imon<_TMenv->nmon(); imon++ )
        vec[ivec++] = 0.;
    }
    for( unsigned int iQ=0; iQ<_nx*(_nx+1)/2; iQ++ )
      vec[ivec++] = A[iQ];
    if( options.WRAPMIT == Options::ELLIPS ) return;

    for( unsigned int ix=0; ix<_nx; ix++ ){
      vec[ivec++] = Op<T>::l( TM[ix].remainder() );
      vec[ivec++] = Op<T>::u( TM[ix].remainder() );
    }
    return;
  }
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_init
( TVT*TMx0 )
{
  switch( options.WRAPMIT){

  case Options::NONE:
  case Options::DINEQ:
    for( unsigned int ix=0; ix<_nx; ix++ )
      TMx0[ix] = IVP().IC( ix, _TMp );
    return _TM2vec( TMx0, 0, 0, _vec_state );

  case Options::LINTRANS:
  case Options::DITRANS:
    for( unsigned int ip=0; ip<_np; ip++ )
      _pref[ip] = Op<T>::mid( _TMp[ip].B() );
    for( unsigned int ix=0; ix<_nx; ix++ ){
      TMx0[ix] = IVP().IC( ix, _TMp ).center();
      _Id[ix] = TMx0[ix].remainder();
      for( unsigned int jx=0; jx<_nx; jx++ )
        _A[jx*_nx+ix] = _Z[jx*_nx+ix] = (ix==jx? 1.: 0.);
    }
    return _TM2vec( TMx0, _Id, _A, _Z, _vec_state );

  case Options::ELLIPS:
  default:{
    for( unsigned int ip=0; ip<_np; ip++ )
      _pref[ip] = Op<T>::mid( _TMp[ip].B() );
    double trR2 = 0.;
    for( unsigned int ix=0, iQ=0; ix<_nx; ix++ ){
      TMx0[ix] = IVP().IC( ix, _TMp ).center();
      trR2 += sqr( Op<T>::diam( TMx0[ix].remainder() ) / 2. );
      _Q[iQ++] = Op<T>::diam( TMx0[ix].remainder() ) / 2.;
      for( unsigned int jx=ix+1; jx<_nx; jx++ ) _Q[iQ++] = 0.;
      _Id[ix] = TMx0[ix].remainder();
    }
    for( unsigned int ix=0, iQ=0; ix<_nx; ix++, iQ+=_nx-ix )
      _Q[iQ] *= std::sqrt( trR2 );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    E ERx0( _nx, _Q );
    std::cout << "ERx0 =" << ERx0 << std::endl;
#endif
    return _TM2vec( TMx0, _Q, 0, _vec_state ); }
  }
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::_ODEBND_GSL_RHSTM
( double t, const double* x, double* xdot, void* user_data )
{
  if( !_vec2var( x, _TMx ) ) return GSL_EBADFUNC;
#ifdef MC__ODEBND_GSL_TM_DEBUG
  std::cerr << "TMx Intermediate" << std::endl;
  _print_interm( t, _TMx, "x", std::cerr );
#endif
  stats.numRHS++;

  switch( options.WRAPMIT){
  case Options::NONE:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg );
      _RxLdot[ix] = Op<T>::l( _TMxdot[ix].R() );
      _RxUdot[ix] = Op<T>::u( _TMxdot[ix].R() );
#ifdef MC__ODEBND_GSL_TM_DEBUG
      std::cout << "TMf[" << ix << "] = " << _TMxdot[ix] << std::endl;
#endif
    }
    _TM2vec( _TMxdot,_RxLdot, _RxUdot, 0, 0, xdot );
#ifdef MC__ODEBND_GSL_TM_DEBUG
    {int dum; std::cin >> dum;}
#endif
    return GSL_SUCCESS;  

  case Options::DINEQ:
    for( unsigned int ix=0; ix<_nx; ix++ ){
      T Rxi = _TMx[ix].R();
      _TMx[ix].set( Op<T>::l( Rxi ) );
      _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg );
      _RxLdot[ix] = Op<T>::l( _TMxdot[ix].R() );
      _TMx[ix].set( Op<T>::u( Rxi ) );
      _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg );
      _RxUdot[ix] = Op<T>::u( _TMxdot[ix].R() );
      _TMx[ix].set( Rxi );
    }
    _TM2vec( _TMxdot, _RxLdot, _RxUdot, 0, 0, xdot );
#ifdef MC__ODEBND_GSL_TM_DEBUG
    std::cerr << "TMxdot Intermediate" << std::endl;
    _print_interm( t, _TMxdot, "x", std::cerr );
    {int dum; std::cin >> dum;}
#endif
    return GSL_SUCCESS;  

  case Options::LINTRANS:
  case Options::DITRANS:{
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    std::cout << "@t=" << t << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "d[" << ix << "] = " << _Id[ix]
                << " : " << Op<T>::mid(_Id[ix]) << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "A[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _A[jx*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      std::cout << "Z[" << ix << ",#] = ";
      for( unsigned int jx=0; jx<_nx; jx++ )
        std::cout << _Z[jx*_nx+ix] << "  ";
      std::cout << std::endl;
    }
    { int dum; std::cin >> dum; }
#endif

    // Taylor Expansion of order q in (X,P)
    if( options.ORDMIT ){

      for( unsigned int jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _np+jx, _Id[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _TMx[jx].center().set( T(0.) ), true );
      }
      _d2x( _MVXPd, _MVXPx, false );
      for( unsigned int ix=0; ix<_nx; ix++ )
        _Iddot[ix] = 0.;

      for( unsigned int ix=0; ix<_nx; ix++ ){
        TVT MVXPfi = IVP().RHS( ix, _MVXPp, _MVXPx, t, _istg ).center();
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "MVXPf[" << ix << "] = " << MVXPfi << std::endl;
#endif
        // Extract Taylor model in P and set to 0
        MVXPfi.get( _TMxdot[ix].set(_TMenv), true );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "MVXPf[" << ix << "] = " << MVXPfi << std::endl;
        std::cout << "TMxdot[" << ix << "] =" << _TMxdot[ix] << std::endl;
        { int dum; std::cin >> dum; }
#endif
      // Extract linear part in X and set to 0
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "_Adot[" << ix << ",#] = ";
#endif
        for( unsigned int jx=0; jx<_nx; jx++ ){
          _Adot[ix+jx*_nx] = MVXPfi.linear( _np+jx, true );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
          std::cout << _Adot[ix+jx*_nx] << "  ";
#endif
        }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << std::endl;
#endif
        // Bound remaining terms
        T MVXPfiR = MVXPfi.bound();
        for( unsigned int jx=0; jx<_nx; jx++ )
          _Iddot[jx] += _Z[ix*_nx+jx] * MVXPfiR;
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "Iddot[" << ix << "] = " << _Iddot[ix]
                  << " : " << Op<T>::mid(_Iddot[ix]) << std::endl;
#endif
      }

      // Compute Zdot = -Z*Adot*Z
      double* ZAdot = new double[_nx*_nx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      double*AZ = new double[_nx*_nx];
#endif
      for( unsigned int ix=0; ix<_nx; ix++ ){
        for( unsigned int jx=0; jx<_nx; jx++ ){
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
          AZ[jx*_nx+ix] = 0.;
#endif
          ZAdot[jx*_nx+ix] = 0.;
          for( unsigned int kx=0; kx<_nx; kx++ ){
	    ZAdot[jx*_nx+ix] += _Z[kx*_nx+ix] * _Adot[jx*_nx+kx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
            AZ[jx*_nx+ix] += _Z[kx*_nx+ix] * _A[jx*_nx+kx];
#endif
          }
        }
      }
      for( unsigned int ix=0; ix<_nx; ix++ ){
        for( unsigned int jx=0; jx<_nx; jx++ ){
          _Zdot[jx*_nx+ix] = 0.;
          for( unsigned int kx=0; kx<_nx; kx++ )
	    _Zdot[jx*_nx+ix] -= ZAdot[kx*_nx+ix] * _Z[jx*_nx+kx];
        }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "A*Z[" << ix << ",#] = ";
        for( unsigned int jx=0; jx<_nx; jx++ )
          std::cout << AZ[jx*_nx+ix] << "  ";
        std::cout << std::endl;
#endif
      }
      delete[] ZAdot;
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      delete[] AZ;
      for( unsigned int ix=0; ix<_nx; ix++ ){
        std::cout << "Adot[" << ix << ",#] = ";
        for( unsigned int jx=0; jx<_nx; jx++ )
          std::cout << _Adot[jx*_nx+ix] << "  ";
        std::cout << std::endl;
      }
#endif
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _dLdot[ix] = Op<T>::l( _Iddot[ix] );
        _dUdot[ix] = Op<T>::u( _Iddot[ix] );
      }
    }
    
    // Taylor Expansion of order q in P and mean-value theorem in X
    else{
      T*Ip = new T[_np];
      for( unsigned int ip=0; ip<_np; ip++ )
        Ip[ip] = _TMp[ip].B();

      FT*dIx = new FT[_nx];
      for( unsigned int ix=0; ix<_nx; ix++ ){
        dIx[ix] = _TMx[ix].B(); dIx[ix].diff(ix,_nx);
        for( unsigned int jx=0; jx<_nx; jx++ )
          _Adot[jx*_nx+ix] = _Zdot[jx*_nx+ix] = 0.;
      }

      FT*dIf = new FT[_nx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      double*AZ = new double[_nx*_nx];
#endif
      for( unsigned int ix=0; ix<_nx; ix++ ){
        dIf[ix] = IVP().RHS( ix, Ip, dIx, t, _istg );
        for( unsigned int jx=0; jx<_nx; jx++ ){
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
          AZ[jx*_nx+ix] = 0.;
#endif
          for( unsigned int kx=0; kx<_nx; kx++ ){
            _Adot[jx*_nx+ix] += Op<T>::mid( dIf[ix].d(kx) ) * _A[jx*_nx+kx];
            _Zdot[kx*_nx+jx] -= _Z[ix*_nx+jx] * Op<T>::mid( dIf[ix].d(kx) );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
            AZ[jx*_nx+ix] += _Z[kx*_nx+ix] * _A[jx*_nx+kx];
#endif
          }
        }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "A*Z[" << ix << ",#] = ";
        for( unsigned int jx=0; jx<_nx; jx++ )
          std::cout << AZ[jx*_nx+ix] << "  ";
        std::cout << std::endl;
#endif
    }
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      delete[] AZ;
#endif

      // Compute Taylor model polynomial part
      for( unsigned int ix=0; ix<_nx; ix++ )
        ( _TMx[ix].C() ).set( T(0.) );
      for( unsigned int ix=0; ix<_nx; ix++ )
        _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg );

      // Compute Taylor rotated remainder part
      T*JA = new T[_nx*_nx];
      for( unsigned int ix=0; ix<_nx; ix++ ){
        for( unsigned int jx=0; jx<_nx; jx++ ){
          JA[jx*_nx+ix] = 0.;
          for( unsigned int kx=0; kx<_nx; kx++ )
            JA[jx*_nx+ix] += ( dIf[ix].d(kx) - Op<T>::mid( dIf[ix].d(kx) ) )
	                   * _A[jx*_nx+kx];
        }
      }

      if( options.WRAPMIT == Options::LINTRANS ){
        for( unsigned int ix=0; ix<_nx; ix++ ){
          _Iddot[ix] = 0.;
          for( unsigned int jx=0; jx<_nx; jx++ ){
            T JjAD = 0.;
            for( unsigned int kx=0; kx<_nx; kx++ )
              JjAD += JA[kx*_nx+jx] * _Id[kx];
            _Iddot[ix] += _Z[jx*_nx+ix] * ( _TMxdot[jx].R() + JjAD );
          }
          _dLdot[ix] = Op<T>::l( _Iddot[ix] );
          _dUdot[ix] = Op<T>::u( _Iddot[ix] );
        }
      }

      else{
        for( unsigned int ix=0; ix<_nx; ix++ ){
#ifdef MC__ODEBND_GSL_LINTR_TMJAC
          T _Idi = _Id[ix];
          _Id[ix] = Op<T>::l( _Idi ) * Op<T>::zeroone();
          _d2x( _Id, _Ix );
          for( unsigned int ix2=0; ix2<_nx; ix2++ ){
            dIx[ix2] = _Ix[ix2]; dIx[ix2].diff(ix2,_nx);
          }
          for( unsigned int ix2=0; ix2<_nx; ix2++ ){
            dIf[ix2] = IVP().RHS( ix2, Ip, dIx, t, _istg );
            for( unsigned int jx=0; jx<_nx; jx++ ){
              JA[jx*_nx+ix2] = 0.;
              for( unsigned int kx=0; kx<_nx; kx++ )
                JA[jx*_nx+ix2] += ( dIf[ix2].d(kx) - Op<T>::mid( dIf[ix2].d(kx) ) )
                                * _A[jx*_nx+kx];
            }
          }
#endif
          _Iddot[ix] = 0.;
          for( unsigned int jx=0; jx<_nx; jx++ ){
            T JjADL = 0.;
            for( unsigned int kx=0; kx<_nx; kx++ )
              JjADL += JA[kx*_nx+jx] * ( kx==ix? Op<T>::l(_Id[kx]): _Id[kx] );
            _Iddot[ix] += _Z[jx*_nx+ix] * ( _TMxdot[jx].R() + JjADL );
          }
          _dLdot[ix] = Op<T>::l( _Iddot[ix] );

#ifdef MC__ODEBND_GSL_LINTR_TMJAC
          _Id[ix] = Op<T>::u( _Idi ) * Op<T>::zeroone();
          _d2x( _Id, _Ix );
          for( unsigned int ix2=0; ix2<_nx; ix2++ ){
            dIx[ix2] = _Ix[ix2]; dIx[ix2].diff(ix2,_nx);
          }
          for( unsigned int ix2=0; ix2<_nx; ix2++ ){
            dIf[ix2] = IVP().RHS( ix2, Ip, dIx, t, _istg );
            for( unsigned int jx=0; jx<_nx; jx++ ){
              JA[jx*_nx+ix2] = 0.;
              for( unsigned int kx=0; kx<_nx; kx++ )
                JA[jx*_nx+ix2] += ( dIf[ix2].d(kx) - Op<T>::mid( dIf[ix2].d(kx) ) )
                                * _A[jx*_nx+kx];
            }
          }
#endif
          _Iddot[ix] = 0.;
          for( unsigned int jx=0; jx<_nx; jx++ ){
            T JjADU = 0.;
            for( unsigned int kx=0; kx<_nx; kx++ )
              JjADU += JA[kx*_nx+jx] * ( kx==ix? Op<T>::u(_Id[kx]): _Id[kx] );
            _Iddot[ix] += _Z[jx*_nx+ix] * ( _TMxdot[jx].R() + JjADU );
          }
          _dUdot[ix] = Op<T>::u( _Iddot[ix] );

#ifdef MC__ODEBND_GSL_LINTR_TMJAC
          _Id[ix] = _Idi;
#endif
        }
      }

      delete[] Ip;
      delete[] dIx;
      delete[] dIf;
      delete[] JA;
    }
    
    // Copy time derivative in GSL array
    _TM2vec( _TMxdot, _dLdot, _dUdot, _Adot, _Zdot, xdot );
    
    return GSL_SUCCESS;
   }
   
  case Options::ELLIPS:{
  default:
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    std::cout << "@t=" << t << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "TMx[ " << ix << "] = " << _TMx[ix] << std::endl;
    std::cout << "Ed = " << _Ed << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "Id[" << ix << "] = " << _Id[ix]
                << " : " << Op<T>::mid(_Id[ix]) << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "Ix[" << ix << "] = " << _Ix[ix]
                << " : " << Op<T>::mid(_Ix[ix]) << std::endl;
#endif

    // Taylor Expansion of order q in P and mean-value theorem in X
    if( !options.ORDMIT ){
    
      T*Ip = new T[_np];
      for( unsigned int ip=0; ip<_np; ip++ )
        Ip[ip] = _TMp[ip].B();
      FT*dIx = new FT[_nx];
      F*dx = new F[_nx];
      for( unsigned int ix=0; ix<_nx; ix++ ){
        dIx[ix] = _TMx[ix].center().B();
        dIx[ix].diff(ix,_nx);
        dx[ix] = _TMx[ix].constant();
        dx[ix].diff(ix,_nx);
      }

      // Compute Taylor model polynomial part
      for( unsigned int ix=0; ix<_nx; ix++ )
        ( _TMx[ix].center() ).set( T(0.) );
      for( unsigned int ix=0; ix<_nx; ix++ )
        _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg ).center();

      // Compute Taylor remainder ellipsoidal bound
      double trQ = 0., sumkappa = 0.;
      for( unsigned int ix=0; ix<_nx; ix++ )
        trQ += ( _Q[_ndxLT(ix,ix)]>0? _Q[_ndxLT(ix,ix)]: machprec() );
      const double srqt_trQ = (trQ>0? std::sqrt( trQ ): 0.) + options.QTOL;
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
      const double* sqrtQ = _Ed.sqrtQ(true).array;
#endif

      for( unsigned int ix=0; ix<_nx; ix++ ){
        FT dIfi = IVP().RHS( ix, Ip, dIx, t, _istg );
        F  dfi  = IVP().RHS( ix, _pref, dx, t, _istg );
        for( unsigned int jx=0; jx<_nx; jx++ )
          //_A[ix+jx*_nx] = Op<T>::mid( dIfi.d(jx) );
          _A[ix+jx*_nx] = dfi.d(jx);
        _Iddot[ix] = _TMxdot[ix].remainder();
        for( unsigned int jx=0; jx<_nx; jx++ )
          _Iddot[ix] += ( dIfi.d(jx) - _A[ix+jx*_nx] ) * _Ix[jx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "Iddot[" << ix << "] = " << _Iddot[ix]
                  << " : " << Op<T>::mid(_Iddot[ix]) << std::endl;
        std::cout << "kappa[" << ix << "] = "
                  << Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ ) << std::endl;
#endif
        sumkappa += Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ );
      }

      for( unsigned int jx=0; jx<_nx; jx++ ){
        for( unsigned int ix=jx; ix<_nx; ix++ ){
          _Qdot[_ndxLT(ix,jx)] = sumkappa * _Q[_ndxLT(ix,jx)];
          for( unsigned int kx=0; kx<_nx; kx++ )
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
            _Qdot[_ndxLT(ix,jx)] += sqrtQ[kx+ix*_nx] * _A[jx+kx*_nx]
	                          + _A[ix+kx*_nx] * sqrtQ[kx+jx*_nx];
#else
            _Qdot[_ndxLT(ix,jx)] += _Q[_ndxLT(ix,kx)] * _A[jx+kx*_nx]
	                          + _A[ix+kx*_nx] * _Q[_ndxLT(kx,jx)];
#endif
        }
        _Qdot[_ndxLT(jx,jx)] += Op<T>::diam( _Iddot[jx] ) / 2. * srqt_trQ;
      }

      delete[] Ip;
      delete[] dIx;
      delete[] dx;
    }

    // Taylor Expansion of order q in P and order options.ORDMIT in (X,P)
    else if( options.ORDMIT < _TMenv->nord() ){
    
      T*Ip = new T[_np];
      for( unsigned int ip=0; ip<_np; ip++ )
        Ip[ip] = _TMp[ip].B();
      FT*dIx = new FT[_nx];
      F*dx = new F[_nx];
      for( unsigned int ix=0; ix<_nx; ix++ ){
        dIx[ix] = _TMx[ix].center().B();
        dIx[ix].diff(ix,_nx);
        dx[ix] = _TMx[ix].constant();
        dx[ix].diff(ix,_nx);
      }

      // Compute Taylor model polynomial part
      for( unsigned int ix=0; ix<_nx; ix++ )
        ( _TMx[ix].center() ).set( T(0.) );
      for( unsigned int ix=0; ix<_nx; ix++ )
        _TMxdot[ix] = IVP().RHS( ix, _TMp, _TMx, t, _istg ).center();

      // Compute Taylor remainder ellipsoidal bound
      double trQ = 0., sumkappa = 0.;
      for( unsigned int ix=0; ix<_nx; ix++ )
        trQ += ( _Q[_ndxLT(ix,ix)]>0? _Q[_ndxLT(ix,ix)]: machprec() );
      const double srqt_trQ = (trQ>0? std::sqrt( trQ ): 0.) + options.QTOL;
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
      const double* sqrtQ = _Ed.sqrtQ(true).array;
#endif

      for( unsigned int ix=0; ix<_nx; ix++ ){
        FT dIfi = IVP().RHS( ix, Ip, dIx, t, _istg );
        F  dfi  = IVP().RHS( ix, _pref, dx, t, _istg );
        for( unsigned int jx=0; jx<_nx; jx++ )
          //_A[ix+jx*_nx] = Op<T>::mid( dIfi.d(jx) );
          _A[ix+jx*_nx] = dfi.d(jx);
        _Iddot[ix] = _TMxdot[ix].remainder();
        for( unsigned int jx=0; jx<_nx; jx++ )
          _Iddot[ix] += ( dIfi.d(jx) - _A[ix+jx*_nx] ) * _Ix[jx];
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "Iddot[" << ix << "] = " << _Iddot[ix]
                  << " : " << Op<T>::mid(_Iddot[ix]) << std::endl;
        std::cout << "kappa[" << ix << "] = "
                  << Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ ) << std::endl;
#endif
        sumkappa += Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ );
      }

      for( unsigned int jx=0; jx<_nx; jx++ ){
        for( unsigned int ix=jx; ix<_nx; ix++ ){
          _Qdot[_ndxLT(ix,jx)] = sumkappa * _Q[_ndxLT(ix,jx)];
          for( unsigned int kx=0; kx<_nx; kx++ )
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
            _Qdot[_ndxLT(ix,jx)] += sqrtQ[kx+ix*_nx] * _A[jx+kx*_nx]
	                          + _A[ix+kx*_nx] * sqrtQ[kx+jx*_nx];
#else
            _Qdot[_ndxLT(ix,jx)] += _Q[_ndxLT(ix,kx)] * _A[jx+kx*_nx]
	                          + _A[ix+kx*_nx] * _Q[_ndxLT(kx,jx)];
#endif
        }
        _Qdot[_ndxLT(jx,jx)] += Op<T>::diam( _Iddot[jx] ) / 2. * srqt_trQ;
      }

      delete[] Ip;
      delete[] dIx;
      delete[] dx;
    }

    // Taylor Expansion of order q in (X,P)
    else{

      for( unsigned int jx=0; jx<_nx; jx++ ){
        _MVXPd[jx].set( _MVXPenv, _np+jx, _Id[jx] );
        _MVXPx[jx].set( _MVXPenv ).set( _TMx[jx].center().set( T(0.) ), true );
      }
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
      _e2x( _Ed.sqrtQ(true).array, _MVXPd, _MVXPx, false );
#else
      _e2x( 0, _MVXPd, _MVXPx, false );
#endif
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
      for( unsigned int ix=0; ix<_nx; ix++ )
        std::cout << "MVXPx[ " << ix << "] = " << _MVXPx[ix] << std::endl;
      { int dum; std::cin >> dum; }
#endif

      // Compute Taylor model with remainder ellipsoidal bound
      double trQ = 0., sumkappa = 0.;
      for( unsigned int ix=0; ix<_nx; ix++ )
        trQ += _Q[_ndxLT(ix,ix)];
      const double srqt_trQ = (trQ>0? std::sqrt( trQ ): 0.) + options.QTOL;
#ifdef MC__ODEBND_GSL_ELLTR_SCALED
      const double* sqrtQ = _Ed.sqrtQ(true).array;
#endif

      for( unsigned int ix=0; ix<_nx; ix++ ){
        TVT MVXPfi = IVP().RHS( ix, _MVXPp, _MVXPx, t, _istg ).center();
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "MVXPf[" << ix << "] = " << MVXPfi << std::endl;
#endif
        // Extract Taylor model in P and set to 0
        MVXPfi.get( _TMxdot[ix].set(_TMenv), true );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "TMxdot[" << ix << "] =" << _TMxdot[ix] << std::endl;
        std::cout << "MVXPf[" << ix << "] = " << MVXPfi << std::endl;
        { int dum; std::cin >> dum; }
#endif
        // Extract linear part in X and set to 0
        for( unsigned int jx=0; jx<_nx; jx++ )
          _A[ix+jx*_nx] = MVXPfi.linear( _np+jx, true );
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "A[" << ix << ",#] = ";
        for( unsigned int jx=0; jx<_nx; jx++ )
          std::cout << _A[ix+jx*_nx] << "  ";
        std::cout << std::endl;
        std::cout << "MVXPf[" << ix << "] = " << MVXPfi << std::endl;
#endif
        // Bound remaining terms
        _Iddot[ix] = MVXPfi.bound();
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
        std::cout << "Iddot[" << ix << "] = " << _Iddot[ix]
                  << " : " << Op<T>::mid(_Iddot[ix]) << std::endl;
        std::cout << "kappa[" << ix << "] = "
                  << Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ ) << std::endl;
#endif
        sumkappa += Op<T>::diam( _Iddot[ix] ) / ( 2. * srqt_trQ );
      }

      for( unsigned int jx=0; jx<_nx; jx++ ){
        for( unsigned int ix=jx; ix<_nx; ix++ ){
          _Qdot[_ndxLT(ix,jx)] = sumkappa * _Q[_ndxLT(ix,jx)];
          for( unsigned int kx=0; kx<_nx; kx++ )
            _Qdot[_ndxLT(ix,jx)] += _Q[_ndxLT(ix,kx)] * _A[jx+kx*_nx]
	                          + _A[ix+kx*_nx] * _Q[_ndxLT(kx,jx)];
        }
        _Qdot[_ndxLT(jx,jx)] += Op<T>::diam( _Iddot[jx] ) / 2. * srqt_trQ;
      }
    }    
     
#ifdef MC__ODEBND_GSL_LINTR_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "TMxdot[" << ix << "] =" << _TMxdot[ix] << std::endl;
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "Iddot[" << ix << "] =" << _Iddot[ix] << std::endl;
    E ERxdot( _nx, _Qdot );
    std::cout << "ERxdot =" << ERxdot << std::endl;
    { int dum; std::cin >> dum; }
#endif
   
    // Copy time derivative in GSL array
    _TM2vec( _TMxdot, _RxLdot, _RxUdot, _Qdot, 0, xdot );
    
    return GSL_SUCCESS;
   }
  }
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::MC_GSLRHSTM__
( double t, const double* x, double* xdot, void* user_data )
{
  ODEBND_GSL<T,IVP> *pODEBND_GSL = ODEBND_GSL<T,IVP>::pODEBND_GSL;
  int flag = pODEBND_GSL->_ODEBND_GSL_RHSTM( t, x, xdot, user_data );
  ODEBND_GSL<T,IVP>::pODEBND_GSL = pODEBND_GSL;
  return flag;
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::_ODEBND_GSL_JACTM
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  stats.numJAC++;
  return GSL_EBADFUNC;
}

template <typename T, typename IVP> inline int
ODEBND_GSL<T,IVP>::MC_GSLJACTM__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODEBND_GSL<T,IVP> *pODEBND_GSL = ODEBND_GSL<T,IVP>::pODEBND_GSL;
  int flag = pODEBND_GSL->_ODEBND_GSL_JACTM( t, x, jac, xdot, user_data );
  ODEBND_GSL<T,IVP>::pODEBND_GSL = pODEBND_GSL;
  return flag;
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::_GSL_init
( const TVT* TMp )
{
  // Define ODE system in GSL format
  _sys_dineqTM.function = MC_GSLRHSTM__;
  _sys_dineqTM.jacobian = MC_GSLJACTM__;
  _sys_dineqTM.dimension = _TMenv->nmon()*_nx;
  switch( options.WRAPMIT){
  case Options::LINTRANS:
  case Options::DITRANS:
    _sys_dineqTM.dimension += 2*_nx*_nx;
    // no break here!
  case Options::NONE:
  case Options::DINEQ:
    _sys_dineqTM.dimension += 2*_nx;
    break;
  case Options::ELLIPS:
  default:
    _sys_dineqTM.dimension += _nx*(_nx+1)/2;
    break;
  }
  _sys_dineqTM.params = 0;

  // Set GSL driver
  _GSL_init( _sys_dineqTM, _driver_dineqTM );

  // Internal Taylor model environment reset
  if( _MVXPenv && ( _MVXPenv->nord() != std::min(_TMenv->nord(),options.ORDMIT) 
                 || _MVXPenv->nvar() != _nx+_np ) ){
    delete[] _MVXPx; _MVXPx = 0;
    delete[] _MVXPd; _MVXPd = 0;
    delete[] _MVXPp; _MVXPp = 0;
    delete _MVXPenv; _MVXPenv = 0;
  }

  // Size state/parameter arrays
  _TMp = TMp;
  if( !_TMx )	 _TMx	 = new TVT[_nx];
  if( !_TMxdot ) _TMxdot = new TVT[_nx];
  if( !_RxLdot ) _RxLdot = new double[_nx];
  if( !_RxUdot ) _RxUdot = new double[_nx];

  switch( options.WRAPMIT){
  case Options::ELLIPS:
    if( !_pref )    _pref    = new double[_np];
    if( !_xref )    _xref    = new double[_nx];
    if( !_A )	    _A       = new double[_nx*_nx];
    if( !_Q )	    _Q       = new double[_nx*(_nx+1)/2];
    if( !_Id )      _Id      = new T[_nx];
    if( !_Ix )      _Ix      = new T[_nx];
    if( !_Qdot )    _Qdot    = new double[_nx*(_nx+1)/2];
    if( !_Iddot )   _Iddot   = new T[_nx];
    if( !_MVXPenv ) _MVXPenv = new TMT( _nx+_np, std::min(_TMenv->nord(),options.ORDMIT) );
    _MVXPenv->options = _TMenv->options;
    if( !_MVXPx )   _MVXPx   = new TVT[_nx];
    if( !_MVXPd )   _MVXPd   = new TVT[_nx];
    if( !_MVXPp )   _MVXPp   = new TVT[_np];
    for( unsigned int ip=0; ip<_np; ip++ )
      _MVXPp[ip].set( _MVXPenv, ip, _TMp[ip].B() );
    if( !_MVPp )    _MVPp    = new TVT[_np];
    break;

  case Options::LINTRANS:
  case Options::DITRANS:
    // Size reference and transformation arrays
    if( !_pref )     _pref    = new double[_np];
    if( !_xref )     _xref    = new double[_nx];
    if( !_A )	     _A       = new double[_nx*_nx];
    if( !_Z )	     _Z       = new double[_nx*_nx];
    if( !_Id )       _Id      = new T[_nx];
    if( !_Ix )       _Ix      = new T[_nx];
    if( !_Adot )     _Adot    = new double[_nx*_nx];
    if( !_Zdot )     _Zdot    = new double[_nx*_nx];
    if( !_Iddot )    _Iddot   = new T[_nx];
    if( !_dLdot )    _dLdot   = new double[_nx];
    if( !_dUdot )    _dUdot   = new double[_nx];
    if( !_MVXPenv )  _MVXPenv = new TMT( _nx+_np, std::min(_TMenv->nord(),options.ORDMIT) );
    _MVXPenv->options = _TMenv->options;
    if( !_MVXPx )    _MVXPx   = new TVT[_nx];
    if( !_MVXPd )    _MVXPd   = new TVT[_nx];
    if( !_MVXPp )    _MVXPp   = new TVT[_np];
    for( unsigned int ip=0; ip<_np; ip++ )
      _MVXPp[ip].set( _MVXPenv, ip, _TMp[ip].B() );
    if( !_MVPp )     _MVPp    = new TVT[_np];
    break;

  case Options::NONE:
  case Options::DINEQ:
    break;
  }

  // Reset result record and statistics
  _results_IA.clear();  
  _init_stats();
}

//! @fn template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS ODEBND_GSL<T,IVP>::bounds(
//! const unsigned int ns, const double*tk, const TVT*TMp, TVT**TMxk,
//! std::ostream&os )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>TMp</a> [input] Taylor model of parameter set
//!   - <a>TMxk</a> [output] Taylor model of state enclosures at stage times
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP>
inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::bounds
( const unsigned int ns, const double*tk, const TVT*TMp, TVT**TMxk, std::ostream&os )
{
  // Check Taylor model compatibility and size
  unsigned int kp=_np;
  for( unsigned int ip=0; ip<_np && kp==_np; ip++ )
    if( TMp[ip].env() ) kp = ip;
  if( kp==_np || TMp[kp].env()->nvar()!=_np ) return FATAL;
  _TMenv = TMp[kp].env();  

  try{
    // Initialize trajectory integration with GSL
    _GSL_init( TMp );
    _t = tk[0];
    _init( TMxk[0] );
    if( options.DISPLAY >= 1 ) _print_interm( tk[0], TMxk[0], "x", os );
    // Record initial results
    if( options.RESRECORD )
      _results_IA.push_back( Results( tk[0], _nx, TMxk[0] ) );

    // Integrate ODEs through each stage using GSL
    pODEBND_GSL = this;
    for( _istg=0; _istg<ns; _istg++ ){
      if( gsl_odeiv2_driver_apply( _driver_dineqTM, &_t, tk[_istg+1], _vec_state )
        != GSL_SUCCESS ){
	_GSL_term();
        if( options.DISPLAY >= 1 ) _print_stats( stats, os );
	return FAILURE;
      }
      _vec2var( _vec_state, TMxk[_istg+1] );
      if( options.DISPLAY >= 1 ) _print_interm( tk[_istg+1], TMxk[_istg+1], "x", os );
      // Record intermediate results
      if( options.RESRECORD )
        _results_IA.push_back( Results( tk[_istg+1], _nx, TMxk[_istg+1] ) );
    }
  }
  catch(...){
    _GSL_term();
    if( options.DISPLAY >= 1 ) _print_stats( stats, os );
    return FAILURE;
  }

  _GSL_term(); 
  if( options.DISPLAY >= 1 ) _print_stats( stats, os );
  return NORMAL;
}

template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::_remainders
( const unsigned int ns, const double*tk, const T*Ip, const TVT*const*TMxk,
  T**Rxk, const unsigned int nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  _ODESLV_GSL->options = options; 
  _ODESLV_GSL->options.RESRECORD = false;
  STATUS flag = NORMAL;
  
  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ )
    xk[is] = new double[_nx];
  flag = _ODESLV_GSL->states( ns, tk, p, xk, os );
  if( flag != NORMAL || nsamp <= 1 ){
    delete[] p;
    for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return flag;
  }   
  for( unsigned int is=0; is<=ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Rxk[is][ix] = xk[is][ix] - TMxk[is][ix].polynomial( p );
  
  // Start sampling process
  unsigned int* vsamp = new unsigned int[_np];
  flag = _remainders( ns, tk, Ip, TMxk, Rxk, nsamp, vsamp, 0, p, xk, os );

  options.DISPLAY = DISPLAY_SAVE;
  for( unsigned int is=0; is<=ns; is++ ){
    if( options.DISPLAY >= 1 ) _print_interm( tk[is], Rxk[is], "d", os );
  // Record intermediate results
  if( options.RESRECORD )
    _results_IA.push_back( Results( tk[is], _nx, Rxk[is] ) );
  }
  
  // Clean-up
  delete[] p;
  for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag;
}

template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::_remainders
( const unsigned int ns, const double*tk, const T*Ip, const TVT*const*TMxk,
  T**Rxk, const unsigned int nsamp, unsigned int* vsamp, const unsigned int ip,
  double*p, double**xk, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Update bounds for all sampling points
  for( unsigned int isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < _np ){
      flag = _remainders( ns, tk, Ip, TMxk, Rxk, nsamp, vsamp, ip+1, p, xk, os );
      if( flag != NORMAL ) return flag;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODEBND_GSL_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    flag = _ODESLV_GSL->states( ns, tk, p, xk, os );
    if( flag != NORMAL ) return flag;
    for( unsigned int is=0; is<=ns; is++ )
      for( unsigned int ix=0; ix<_nx; ix++ )
        Rxk[is][ix] = Op<T>::hull( xk[is][ix]-TMxk[is][ix].polynomial(p),
                                   Rxk[is][ix] );
  }

  return flag;
}  

//! @fn template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS ODEBND_GSL<T,IVP>::hausdorff(
//! const unsigned int ns, const double*tk, const TVT*TMp, double**Hxk,
//! double**Rxk, const unsigned int nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the Taylor model bound 
//! and the actual reachable set projected onto each variable as well as the
//! Hausdorff distance between the Taylor model remainder bound and the actual
//! range of the remainder function, using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>TMp</a> [input] Taylor model of parameter set
//!   - <a>Hxk</a> [output] Hausdorff distance between the Taylor model bound 
//!     and the actual reachable set projected onto each variable, at stage times
//!   - <a>Rxk</a> [output] Hausdorff distance between  the Taylor model remainder
//!     bound and the actual range of the remainder function, at stage times
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP> inline typename ODEBND_GSL<T,IVP>::STATUS
ODEBND_GSL<T,IVP>::hausdorff
( const unsigned int ns, const double*tk, const TVT*TMp, double**Hxk,
  double**Rxk, const unsigned int nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  _ODESLV_GSL->options = options;
  _ODESLV_GSL->options.RESRECORD = false;

  // Compute exact bounds 
  T* Ip = new T[_np];
  for( unsigned int ip=0; ip<_np; ip++ ) Ip[ip] = TMp[ip].B();
  T** Ixk0 = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Ixk0[is] = new T[_nx];
  if( _ODESLV_GSL->bounds( ns, tk, Ip, Ixk0, nsamp, os ) != NORMAL ){
    for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
    delete[] Ixk0;
    delete[] Ip;
    return FAILURE;
  }

  // Compute approximate bounds
  TVT** TMxk = new TVT*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) TMxk[is] = new TVT[_nx];
  try{ bounds( ns, tk, TMp, TMxk, os ); }
  catch(...){;}
  unsigned int nsf = _istg;
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      TMxk[is+1][ix] = T( -1e20, 1e20 );

  // Compute remainder bounds 
  T** Rxk0 = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Rxk0[is] = new T[_nx];
  _remainders( nsf, tk, Ip, TMxk, Rxk0, nsamp, os );
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Rxk0[is+1][ix] = T( 0. );

  // Compute Hausdorff metric at each time step
  struct loc{ static double hausdorff
    ( const T&Ix, const T&Ix0 )
    { return std::max( std::fabs(Op<T>::l(Ix)-Op<T>::l(Ix0)),
                       std::fabs(Op<T>::u(Ix)-Op<T>::u(Ix0)) ); }
  };

  options.DISPLAY = DISPLAY_SAVE;
  for( unsigned int is=0; is<ns+1; is++ ){
    for( unsigned int ix=0; ix<_nx; ix++ ){
      Hxk[is][ix] = loc::hausdorff( TMxk[is][ix].B(), Ixk0[is][ix] );
      Rxk[is][ix] = loc::hausdorff( TMxk[is][ix].R(), Rxk0[is][ix] );
    }
    if( options.DISPLAY >= 1 ){
      _print_interm( tk[is], Hxk[is], "d", Rxk[is], "r", os );
    }
  }

  for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
  delete[] Ixk0;
  for( unsigned int is=0; is<ns+1; is++ ) delete[] Rxk0[is];
  delete[] Rxk0;
  for( unsigned int is=0; is<ns+1; is++ ) delete[] TMxk[is];
  delete[] TMxk;
  delete[] Ip;

  return NORMAL;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_print_interm
( const double t, const U*x, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_GSL<T,IVP>::_print_interm
( const double t, const U*x1, const std::string&var1, const U*x2,
  const std::string&var2, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var1.c_str() << "[" << ix << "] = " << x1[ix]
       << "   " << var2.c_str() << "[" << ix << "] = " << x2[ix] << std::endl;
  return;
}

template <typename T, typename IVP> inline void
ODEBND_GSL<T,IVP>::record
( std::ofstream&bndrec, const unsigned int iprec ) const
{
  if( !bndrec ) return;

  // Specify format
  bndrec << std::right << std::scientific << std::setprecision(iprec);

  // Record computed interval bounds at stage times
  typename std::vector< Results >::const_iterator it = _results_IA.begin();
  for( ; it != _results_IA.end(); ++it ){
    bndrec << std::setw(iprec+9) << (*it).t;
    for( unsigned int ix=0; ix<_nx; ix++ )
      bndrec << std::setw(iprec+9) << mc::Op<T>::l( (*it).X[ix] )
             << std::setw(iprec+9) << mc::Op<T>::u( (*it).X[ix] );
    bndrec << std::endl;
  }
}

} // end namescape mc

#endif

