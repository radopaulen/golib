// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODESLV_GSL_HPP
#define MC__ODESLV_GSL_HPP

#undef  MC__ODESLV_GSL_SAMPLE_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "odestruct.hpp"
#include "base_gsl.hpp"

#include "mcop.hpp"

#include "fadiff.h"
#include "gsl/gsl_errno.h"
#include "gsl/gsl_odeiv2.h"

namespace mc
{
//! @brief C++ class computing solutions of parametric ODEs using GSL.
////////////////////////////////////////////////////////////////////////
//! mc::ODESLV_GSL is a C++ class that computes solutions of parametric
//! ordinary differential equations (ODEs) using GSL.
////////////////////////////////////////////////////////////////////////
template <typename T, typename IVP>
class ODESLV_GSL: public BASE_GSL, public ODESTRUCT
{
private:

  typedef fadbad::F< double > F;

  //! @brief GSL data type for ODE system (original ODEs)
  gsl_odeiv2_system _sys_traj;

  //! @brief GSL driver for ODE system (original ODEs)
  gsl_odeiv2_driver *_driver_traj;

  //! @brief states in forward AD
  F *_Fx;

  //! @brief parameters -- do *NOT* free!
  const double *_p;

public:
  /** @ingroup ODESLV
   *  @{
   */

  //! @brief Default constructor
  ODESLV_GSL
    ()
    : BASE_GSL(), ODESTRUCT(IVP())
    {
      // Initalize state/parameter
      _p = 0;
      _Fx = 0;

      // Initalize GSL
      _driver_traj = 0;
    }

  //! @brief Default destructor
  virtual ~ODESLV_GSL()
    {
      // Free state arrays -- Do *NOT* delete _p _Ip _TMp _TMenv
      delete[] _Fx;
      
      // Free GSL arrays
      if( _driver_traj ) gsl_odeiv2_driver_free( _driver_traj );
    }

  //! @brief Integrator options
  struct Options: public BASE_GSL::Options
  {
    //! @brief Constructor
    Options():
      BASE_GSL::Options(), DISPLAY(0), RESRECORD(false)
      {}
    //! @brief Assignment operator
    template <typename U> Options& operator=
      ( U&options ){
        BASE_GSL::Options::operator=(options);
        DISPLAY   = options.DISPLAY;
        RESRECORD = options.RESRECORD;
        return *this;
      }
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
        return "ODESLV_GSL::Exceptions  Error due to calling a not yet implemented feature";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Store integration bounds at a given time instant
  struct Results
  {
    //! @brief Constructors
    Results
      ( const double tk, const unsigned int nxk, const T*Ixk ):
      t( tk ), nx( nxk )
      { X = new T[nx];
        for( unsigned int ix=0; ix<nx; ix++ ) X[ix] = Ixk[ix]; }
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

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned int iprec=5 ) const;

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm
    ( const double t, const U*x, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Computes numerical solution of parametric ODEs
  STATUS states
    ( const unsigned int ns, const double*tk, const double*p, double**xk,
      std::ostream&os=std::cout );

  //! @brief Computes approximate interval enclosure of reachable set of parametric ODEs using parameter sampling
  STATUS bounds
    ( const unsigned int ns, const double*tk, const T*Ip, T**Txk,
      const unsigned int nsamp, std::ostream&os=std::cout );
  /** @} */

  //! @brief static pointer to class
  static ODESLV_GSL<T,IVP> *pODESLV_GSL;

private:

  //! @brief Vector storing interval bound results (see Options::RESRECORD)
  std::vector< Results > _results;

  //! @brief Function to finalize statistics for GSL
  void _GSL_term();

  //! @brief Function to initialize GSL driver
  void _GSL_init
    ( gsl_odeiv2_system &sys, gsl_odeiv2_driver *&driver );

  //! @brief Function to initialize GSL for state trajectories
  void _GSL_init
    ( const double* p );

  //! @brief Static wrapper to function to calculate the ODEs RHS values
  static int MC_GSLRHS__
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Function to calculate the ODEs RHS values
  int _ODESLV_GSL_RHS
    ( double t, const double* x, double* xdot, void* user_data );

  //! @brief Static wrapper to function to calculate the ODEs RHS derivatives
  static int MC_GSLJAC__
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Function to calculate the ODEs RHS derivatives
  int _ODESLV_GSL_JAC
    ( double t, const double* x, double* jac, double* xdot, void* user_data );

  //! @brief Recursive function computing bounds on solutions of IVP in ODEs using sampling
  STATUS _states
    ( const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
      const unsigned int nsamp, unsigned int* vsamp, const unsigned int ip,
      double*p, double**xk, std::ostream&os );

  //! @brief Private methods to block default compiler methods
  ODESLV_GSL(const ODESLV_GSL&);
  ODESLV_GSL& operator=(const ODESLV_GSL&);
};

template <typename T, typename IVP>
ODESLV_GSL<T,IVP>* ODESLV_GSL<T,IVP>::pODESLV_GSL = 0;

template <typename T, typename IVP> inline void
ODESLV_GSL<T,IVP>::_GSL_term()
{
  // Get final CPU time
  _final_stats();
}

template <typename T, typename IVP> inline int
ODESLV_GSL<T,IVP>::_ODESLV_GSL_RHS
( double t, const double* x, double* xdot, void* user_data )
{
  stats.numRHS++;
  for( unsigned int ix=0; ix<_nx; ix++ )
    xdot[ix] = IVP().RHS( ix, _p, x, t, _istg );
  return GSL_SUCCESS;  
}

template <typename T, typename IVP> inline int
ODESLV_GSL<T,IVP>::MC_GSLRHS__
( double t, const double* x, double* xdot, void* user_data )
{
  ODESLV_GSL<T,IVP> *pODESLV_GSL = ODESLV_GSL<T,IVP>::pODESLV_GSL;
  int flag = pODESLV_GSL->_ODESLV_GSL_RHS( t, x, xdot, user_data );
  ODESLV_GSL<T,IVP>::pODESLV_GSL = pODESLV_GSL;
  return flag;
}

template <typename T, typename IVP> inline int
ODESLV_GSL<T,IVP>::_ODESLV_GSL_JAC
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  stats.numJAC++;
  for( unsigned int ix=0; ix<_nx; ix++ ){
    _Fx[ix] = x[ix];
    _Fx[ix].diff(ix,_nx);
  }
  for( unsigned int ix=0; ix<_nx; ix++ ){
    F RHS = IVP().RHS( ix, _p, _Fx, t, _istg );
    xdot[ix] = RHS.x();
    for( unsigned int jx=0; jx<_nx; jx++ )
      jac[ix*_nx+jx] = RHS.d(jx);
  }
  return GSL_SUCCESS;
}

template <typename T, typename IVP> inline int
ODESLV_GSL<T,IVP>::MC_GSLJAC__
( double t, const double* x, double* jac, double* xdot, void* user_data )
{
  ODESLV_GSL<T,IVP> *pODESLV_GSL = ODESLV_GSL<T,IVP>::pODESLV_GSL;
  int flag = pODESLV_GSL->_ODESLV_GSL_JAC( t, x, jac, xdot, user_data );
  ODESLV_GSL<T,IVP>::pODESLV_GSL = pODESLV_GSL;
  return flag;
}

template <typename T, typename IVP> inline void
ODESLV_GSL<T,IVP>::_GSL_init
( const double* p )
{
  // Define ODE system in GSL format
  _sys_traj.function = MC_GSLRHS__;
  _sys_traj.jacobian = MC_GSLJAC__;
  _sys_traj.dimension = _nx;
  _sys_traj.params = 0;
  
  // Set GSL driver
  _GSL_init( _sys_traj, _driver_traj );

  // Size state/parameter arrays
  _p = p;
  if( !_Fx ) _Fx = new F[ _nx ];

  // Initialize statistics
  _init_stats();

  return;
}

template <typename T, typename IVP> inline void
ODESLV_GSL<T,IVP>::_GSL_init
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

  return;
}

//! @fn template <typename T, typename IVP> inline typename ODESLV_GSL<T,IVP>::STATUS ODESLV_GSL<T,IVP>::states(
//! const unsigned int ns, const double*tk, const double*p, double**xk,
//! std::ostream&os )
//!
//! This function computes the solution of the parametric ODEs defined in IVP:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>p</a> [input] parameter values
//!   - <a>xk</a> [output] state values at stage times
//!   - <a>os</a> [input] output stream
//!   .
//! The return value is the status.
template <typename T, typename IVP> inline typename ODESLV_GSL<T,IVP>::STATUS
ODESLV_GSL<T,IVP>::states
( const unsigned int ns, const double*tk, const double*p, double**xk,
  std::ostream&os )
{
  // Initialize trajectory integration with GSL
  _GSL_init( p );
  _t = tk[0];
  for( unsigned int ix=0; ix<_nx; ix++ )
    xk[0][ix] = IVP().IC( ix, _p );
  if( options.DISPLAY >= 1 ) _print_interm( tk[0], xk[0], "x", os );

  // Integrate ODEs through each stage using GSL
  pODESLV_GSL = this;
  for( _istg=0; _istg<ns; _istg++ ){
    //gsl_odeiv2_driver_set_hmax( _driver_traj, tk[_istg+1]-_t+1e-8 );
    for( unsigned int ix=0; ix<_nx; ix++ )
      xk[_istg+1][ix] = xk[_istg][ix];
    if( gsl_odeiv2_driver_apply( _driver_traj, &_t, tk[_istg+1], xk[_istg+1] )
      != GSL_SUCCESS ){
      _GSL_term();
      return FAILURE;
    }
    if( options.DISPLAY >= 1 ) _print_interm( tk[_istg+1], xk[_istg+1], "x", os );
  }

  _GSL_term();
  if( options.DISPLAY >= 1 ) _print_stats( stats, os );
  return NORMAL;
}

//! @fn template <typename T, typename IVP> inline typename ODESLV_GSL<T,IVP>::STATUS ODESLV_GSL<T,IVP>::bounds(
//! const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
//! const unsigned int nsamp, std::ostream&os )
//!
//! This function computes an approximate interval enclosure of the
//! reachable set of the parametric ODEs defined in IVP using equally
//! spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Ip</a> [input] interval parameter set
//!   - <a>Ixk</a> [output] approximate interval state enclosures at stage times
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP> inline typename ODESLV_GSL<T,IVP>::STATUS
ODESLV_GSL<T,IVP>::bounds
( const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
  const unsigned int nsamp, std::ostream&os )
{
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = 0;
  STATUS flag = NORMAL;
  
  // Initialization of sampled bounds at parameter lower bound
  double *p = new double[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    p[ip] = Op<T>::l(Ip[ip]);
  double **xk = new double*[ns+1];
  for( unsigned int is=0; is<=ns; is++ )
    xk[is] = new double[_nx];
  flag = states( ns, tk, p, xk, os );
  if( flag != NORMAL || nsamp <= 1 ){
    delete[] p;
    for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
    return flag;
  }   
  for( unsigned int is=0; is<=ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Ixk[is][ix] = xk[is][ix];
  
  // Start sampling process
  unsigned int* vsamp = new unsigned int[_np];
  flag = _states( ns, tk, Ip, Ixk, nsamp, vsamp, 0, p, xk, os );

  options.DISPLAY = DISPLAY_SAVE;
  for( unsigned int is=0; is<=ns; is++ ){
    if( options.DISPLAY >= 1 ) _print_interm( tk[is], Ixk[is], "x", os );
  // Record intermediate results
  if( options.RESRECORD )
    _results.push_back( Results( tk[is], _nx, Ixk[is] ) );
  }
  
  // Clean-up
  delete[] p;
  for( unsigned int is=0; is<=ns; is++ ) delete[] xk[is]; delete[] xk;
  delete[] vsamp;
  
  return flag;
}

template <typename T, typename IVP> inline typename ODESLV_GSL<T,IVP>::STATUS
ODESLV_GSL<T,IVP>::_states
( const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
  const unsigned int nsamp, unsigned int* vsamp, const unsigned int ip,
  double*p, double**xk, std::ostream&os )
{
  STATUS flag = NORMAL;

  // Update bounds for all sampling points
  for( unsigned int isamp=0; isamp<nsamp; isamp++ ){
    vsamp[ip] = isamp;

    // Continue recursive call
    if( ip+1 < _np ){
      flag = _states( ns, tk, Ip, Ixk, nsamp, vsamp, ip+1, p, xk, os );
      if( flag != NORMAL ) return flag;
      continue;
    }

    // Update bounds for current point
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
    std::cout << "Sample: ";
#endif
    for( unsigned int ip=0; ip<_np; ip++ ){
      p[ip] = Op<T>::l( Ip[ip] ) + vsamp[ip]/(nsamp-1.) * Op<T>::diam( Ip[ip] );
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
      std::cout << p[ip] << "  ";
#endif
    }
#ifdef MC__ODESLV_GSL_SAMPLE_DEBUG
    std::cout << std::endl;
#endif
    flag = states( ns, tk, p, xk, os );
    if( flag != NORMAL ) return flag;
    for( unsigned int is=0; is<=ns; is++ )
      for( unsigned int ix=0; ix<_nx; ix++ )
        Ixk[is][ix] = Op<T>::hull( xk[is][ix], Ixk[is][ix] );
  }

  return flag;
}  

template <typename T, typename IVP> inline void
ODESLV_GSL<T,IVP>::record
( std::ofstream&bndrec, const unsigned int iprec ) const
{
  if( !bndrec ) return;

  // Specify format
  bndrec << std::right << std::scientific << std::setprecision(iprec);

  // Record computed interval bounds at stage times
  typename std::vector< Results >::const_iterator it = _results.begin();
  for( ; it != _results.end(); ++it ){
    bndrec << std::setw(iprec+9) << (*it).t;
    for( unsigned int ix=0; ix<_nx; ix++ )
      bndrec << std::setw(iprec+9) << mc::Op<T>::l( (*it).X[ix] )
             << std::setw(iprec+9) << mc::Op<T>::u( (*it).X[ix] );
    bndrec << std::endl;
  }
}

template <typename T, typename IVP> template<typename U> inline void
ODESLV_GSL<T,IVP>::_print_interm
( const double t, const U*x, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

} // end namescape mc

#endif

