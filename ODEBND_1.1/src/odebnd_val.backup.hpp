// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ODEBND_VAL_HPP
#define MC__ODEBND_VAL_HPP

#undef  MC__ODEBND_VAL_DEBUG

#include <stdexcept>
#include <cassert>
#include <cmath>
#include <fstream>
#include <vector>
#include <sys/time.h>

#include "odestruct.hpp"
#include "odetaylor.hpp"
#include "odeslv_gsl.hpp"
#include "ellipsoid_cpplapack.hpp"
#include "tmodel.hpp"
#include "mcop.hpp"
#include "mcfadbad.hpp"

namespace mc
{
//! @brief C++ class computing enclosures of the reachable set of parametric ODEs using MC++ and (validated) Taylor expansion
////////////////////////////////////////////////////////////////////////
//! mc::ODEBND_VAL is a C++ class that computes enclosures of the
//! reachable set of parametric ordinary differential equations
//! (ODEs) using MC++. It implements a validated method based on Taylor
//! expansion where either intervals or ellispoids are propagated, as
//! well as their combination with Taylor models to enable high-order
//! convergence. The use of ellipsoids enables stability of the
//! enclosures when the parameter host is sufficiently small.
////////////////////////////////////////////////////////////////////////
template <typename T, typename IVP>
class ODEBND_VAL: public virtual ODESTRUCT, public ODETAYLOR<T,IVP>,
                  public ODETAYLOR<TVar<T>,IVP>, public BASE_ODE
                   
{
private:

  typedef fadbad::F< double > F;
  typedef fadbad::F< T > FT;
  typedef Ellipsoid E;
  typedef TModel<T> TMT;
  typedef TVar<T> TVT;
  typedef fadbad::F< TVT > FTVT;

  //! @brief ODE solver based on GSL
  ODESLV_GSL<T,IVP>* _ODESLV_GSL;

  //! @brief parameter interval bounds
  T *_Ip;

  //! @brief parameter ellipsoidal bounds
  E _Ep;

  //! @brief state interval bounds
  T *_Ix;

  //! @brief joint state-parameter ellipsoidal bounds
  E _Exp;

  //! @brief joint state-parameter remainder bounds in Taylor model
  T *_Rxp;

  //! @brief state interval bounds for stepsize validation
  T *_Ix_val;

  //! @brief RHS derivatives w.r.t. state-parameter
  CPPL::dgematrix _fxp;

  //! @brief stepsize
  double _h;

  //! @brief maximal valid stepsize
  double _hmax;

  //! @brief Maximal propagation time
  double _tmax;

  //! @brief Scaling coefficients
  double* _scaling;

  //! @brief Step count
  unsigned long _nsteps;

  //! @brief mc::TModel environment
  TMT *_TMenv;

  //! @brief parameter Taylor model
  TVT *_TMp;

  //! @brief state Taylor model
  TVT *_TMx;

  //! @brief state interval remainder bounds in Taylor model
  T *_Rx;

  //! @brief state interval remainder bounds in Taylor model at expansion point
  T *_Rx_TE;

  //! @brief state ellipsoidal remainder bounds in Taylor model
  E _Ex;

  //! @brief state ellipsoidal remainder bounds in Taylor model at expansion point
  E _Ex_TE;

  //! @brief RHS derivatives w.r.t. state remainders
  CPPL::dgematrix _fx;

  //! @brief mc::TModel environment for ODE right-hand side
  TMT *_TMRHSenv;

  //! @brief State in right-hand side Taylor expansion
  TVT *_TMRHSx;

  //! @brief Parameter in right-hand side Taylor expansion
  TVT *_TMRHSp;

  //! @brief mc::TModel environment for ODE intial value
  TMT *_TMIVenv;

  //! @brief Parameter in initial value Taylor expansion
  TVT *_TMIVp;

  //! @brief Structure holding Taylor expansion arrays
  template <typename U>
  struct TExpansion
  {
    //! @brief Taylor expansion order
    unsigned int ordmax;
    //! @brief Taylor expansion coefficient values
    U** f;
    //! @brief Current time
    double t;
    //! @brief Default constructor
    TExpansion(): ordmax(0), f(0), t(0./0.) {}
    //! @brief Destructor
    virtual ~TExpansion()
    {
      reset();
    }
    //! @brief Reset arrays
    void reset()
    {
      for( unsigned int q=0; q<=ordmax; q++ )
        if( f ) delete[] f[q];
      delete[] f; f = 0;
      t = 0./0.;
    }
    //! @brief Reset arrays
    void set
    ( const unsigned int order )
    {
      if( order <= ordmax ) return;
      reset();
      ordmax = order;
      f = new U*[ordmax+1];
      const unsigned int nx = IVP().nx();
      for( unsigned int q=0; q<=order; q++ )
        f[q] = new U[nx];
    }
  };

  //! @brief Structure holding Taylor expansion arrays
  template <typename U>
  struct TFExpansion: public TExpansion<U>
  {
    //! @brief Taylor expansion coefficient derivatives w.r.t. states
    U** dfdx;
    //! @brief Taylor expansion coefficient derivatives w.r.t. parameters
    U** dfdp;
    //! @brief Default constructor
    TFExpansion(): TExpansion<U>(), dfdx(0), dfdp(0) {}
    //! @brief Destructor
    virtual ~TFExpansion()
    {
      reset();
    }
    //! @brief Reset arrays
    void reset()
    {
      for( unsigned int q=0; q<=TExpansion<U>::ordmax; q++ ){
        if( dfdx ) delete[] dfdx[q];
        if( dfdp ) delete[] dfdp[q];
      }
      delete[] dfdx; delete[] dfdp;
      dfdx = dfdp = 0;
      TExpansion<U>::reset();
    }
    //! @brief Reset arrays
    void set
    ( const unsigned int order )
    {
      if( order <= TExpansion<U>::ordmax ) return;
      reset();
      TExpansion<U>::set(order);
      dfdx = new U*[TExpansion<U>::ordmax+1];
      dfdp = new U*[TExpansion<U>::ordmax+1];
      const unsigned int nx = IVP().nx(), np = IVP().np();
      for( unsigned int q=0; q<=order; q++ ){
        dfdx[q] = new U[nx*nx];
        dfdp[q] = ( np? new U[nx*np]: 0 );
      }
    }
  };

  //! @brief Structure holding Taylor expansion arrays
  template <typename U>
  struct TFFdExpansion: public TExpansion<U>
  {
    //! @brief Taylor expansion coefficient directional second derivatives w.r.t. states
    U** d2fdx2dd;
    //! @brief Default constructor
    TFFdExpansion(): TExpansion<U>(), d2fdx2dd(0) {}
    //! @brief Destructor
    virtual ~TFFdExpansion()
    {
      reset();
    }
    //! @brief Reset arrays
    void reset()
    {
      for( unsigned int q=0; q<=TExpansion<U>::ordmax; q++ )
        if( d2fdx2dd ) delete[] d2fdx2dd[q];
      delete[] d2fdx2dd; d2fdx2dd = 0;
      TExpansion<U>::reset();
    }
    //! @brief Reset arrays
    void set
    ( const unsigned int order )
    {
      if( order <= TExpansion<U>::ordmax ) return;
      reset();
      TExpansion<U>::set(order);
      d2fdx2dd = new U*[TExpansion<U>::ordmax+1];
      const unsigned int nx = IVP().nx();
      for( unsigned int q=0; q<=order; q++ )
        d2fdx2dd[q] = new U[nx*nx];
    }
  };

  //! @brief Taylor expansion arrays for T arithmetic
  TExpansion<T> _T_TE;

  //! @brief Taylor expansion arrays for T arithmetic for stepsize validation
  TExpansion<T> _T_TE_val;

  //! @brief Taylor expansion arrays for TVT arithmetic
  TExpansion<TVT> _TVT_TE;

  //! @brief Taylor expansion arrays with first-order derivatives for T arithmetic
  TFExpansion<T> _T_TFE;

  //! @brief Taylor expansion arrays with directional second-order derivatives for T arithmetic
  TFFdExpansion<T> _T_TFFdE;

  //! @brief Structure holding invariant arrays
  template <typename U>
  struct Invariant
  {
    //! @brief Current time
    double t;
    //! @brief Invariant values
    U* inv;
    //! @brief Default constructor
    Invariant(): t(0./0.)
    {
      const unsigned int ni = IVP().ni();
      inv = ( ni? new U[ni]: 0 );
    }
    //! @brief Destructor
    virtual ~Invariant()
    {
      delete[] inv;
    }
  };

  //! @brief Structure holding invariant arrays
  template <typename U>
  struct FInvariant: public Invariant<U>
  {
    //! @brief Invariant derivatives w.r.t. states
    U* dinvdx;
    //! @brief Invariant derivatives w.r.t. parameters
    U* dinvdp;
    //! @brief Default constructor
    FInvariant(): Invariant<U>()
    {
      const unsigned int ni = IVP().ni();
      const unsigned int nx = IVP().nx();
      const unsigned int np = IVP().np();
      dinvdx = ( ni? new U[ni*nx]: 0 );
      dinvdp = ( ni&&np? new U[ni*np]: 0 );
    }
    //! @brief Destructor
    virtual ~FInvariant()
    {
      delete[] dinvdx;
      delete[] dinvdp;
    }
  };

  //! @brief Structure holding invariant arrays
  template <typename U>
  struct FFdInvariant: public Invariant<U>
  {
    //! @brief Invariant directional second derivatives w.r.t. states
    U* d2invdx2dd;
    //! @brief Default constructor
    FFdInvariant(): Invariant<U>()
    {
      const unsigned int ni = IVP().ni();
      d2invdx2dd = ( ni? new U[ni*ni]: 0 );
    }
    //! @brief Destructor
    virtual ~FFdInvariant()
    {
      delete[] d2invdx2dd;
    }
  };

  //! @brief Invariant arrays for TVT arithmetic
  Invariant<T> _TVT_I;

  //! @brief Invariant arrays with first-order derivatives for T arithmetic
  FInvariant<T> _T_FI;

  //! @brief Invariant arrays with directional second-order derivatives for T arithmetic
  FFdInvariant<T> _T_FFdI;

public:
  /** @ingroup ODEBND
   *  @{
   */
  //! @brief Default constructor
  ODEBND_VAL
    ()
    : ODESTRUCT(IVP()), ODETAYLOR<T,IVP>(), ODETAYLOR<TVT,IVP>(), BASE_ODE()
    {
      // Initalize state/parameter
      _Ip = _Ix = _Ix_val = _Rxp = _Rx = _Rx_TE = 0;
      _TMenv = 0;
      _TMp = _TMx = 0;
      _scaling = new double[_nx];

      // Initialize Taylor expansion in states/parameters
      _TMRHSenv = new TMT( _nx+_np, options.ORDMIT );
      _TMRHSenv->options.BOUNDER_TYPE = options.BNDMIT;
      _TMRHSx = _TMRHSp = 0;
      _TMIVenv = new TMT( _np, options.ORDMIT );
      _TMIVenv->options.BOUNDER_TYPE = options.BNDMIT;
      _TMIVp = 0;

      // Initialize ellipsoidal calculus
#ifdef MC__ODEBND_VAL_DEBUG
      E::options.PSDCHK = true;
#else
      E::options.PSDCHK = false;
#endif

      // Initalize GSL
      _ODESLV_GSL = new ODESLV_GSL<T,IVP>();
    }

  //! @brief Default destructor
  virtual ~ODEBND_VAL()
    {
      // Free arrays -- Do *NOT* delete _TMenv
      delete[] _Ip;
      delete[] _Ix;
      delete[] _Ix_val;
      delete[] _Rxp;
      delete[] _TMRHSx;
      delete[] _TMRHSp;
      delete _TMRHSenv;
      delete[] _TMIVp;
      delete _TMIVenv;
      delete[] _TMp;
      delete[] _TMx;
      delete[] _Rx;
      delete[] _Rx_TE;
      delete[] _scaling;

      // Free GSL arrays
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
      BASE_GSL::Options(), TSORDER(10), WRAPMIT(ELLIPS), ORDMIT(1),
      BNDMIT(TMT::Options::HYBRID), HMIN(1e-8), HMAX(1e-8), HREDUC(0.8),
      TOL(1e-6), SCALING(true), ATOL(1e-8), QTOL(1e-8), TSTOP(true),
      USEINV(true), DISPLAY(0), RESRECORD(false)
      {}
    //! @brief Enumeration of wrapping mitigation strategies
    enum WRAPPING_STRATEGY{
      NONE=0,		//!< No wrapping mitigation
      ELLIPS		//!< Ellipsoidal contractor with linear preconditioning [Default]
    };
    //! @brief Order of Taylor series expansion (Default: 5)
    unsigned int TSORDER;
    //! @brief Wrapping mitigation strategy
    WRAPPING_STRATEGY WRAPMIT;
    //! @brief Order of wrapping mitigation strategy (Default: 2)
    unsigned int ORDMIT;
    //! @brief Bounder type for Taylor model in wrapping mitigation strategy (Default: TMT::Options::HYBRID)
    typename TMT::Options::BOUNDER BNDMIT;
    //! @brief Minimum step-size, \f$h_{\rm min}\f$ (Default: 1e-8)
    double HMIN;
    //! @brief Maximum step-size, \f$h_{\rm max}\f$ (Default: 1e8)
    double HMAX;
    //! @brief Reduction factor for step-size between 0-1 (default: 0.8)
    double HREDUC;
    //! @brief Tolerance on truncation error for step size selection (Default: 1e-6)
    double TOL;
    //! @brief Whether or not to adjust scaling for state components (default: true)
    bool SCALING;
    //! @brief Absolute tolerance for scaling (Default: 1e-8)
    double ATOL;
    //! @brief Tolerance when dividing by trace of shape matrix in ellipsoidal bounds (Default: 1e-8)
    double QTOL;
    //! @brief Whether to stop and reinitialize the integrator at time steps (default: true)
    bool TSTOP;
    //! @brief Whether or not to use the specified invariants for bounds contraction (default: true)
    bool USEINV;
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
      SIZE=1,	//!< Error due to inconsistent size
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
        return "ODEBND_VAL::Exceptions  Error due to calling a not yet implemented feature";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Computes interval and ellipsoidal enclosures of reachable set of parametric ODEs
  template <typename U> STATUS bounds
    ( const unsigned int ns, const double*tk, const U&Up, T**Ixk, E*Exk,
      std::ostream&os=std::cout );

  //! @brief Computes Taylor model enclosure of reachable set of parametric ODEs
  STATUS bounds
    ( const unsigned int ns, const double*tk, const TVT*TMp, const E&Ep,
      TVT**TMxk, E*Exk, std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distance between interval enclosure and actual reachable set of parametric ODEs, using parameter sampling
  template <typename U> STATUS hausdorff
    ( const unsigned int ns, const double*tk, const U&Up, double**Hxk,
      const unsigned int nsamp, const typename ODESLV_GSL<T,IVP>::Options&optGSL,
      std::ostream&os=std::cout );

  //! @brief Computes Hausdorff distances between Taylor model enclosure and actual reachable set and between Taylor model remainder bound and actual remainder function range, using parameter sampling
  STATUS hausdorff
    ( const unsigned int ns, const double*tk, const TVT*TMp, double**HBxk,
      double**HRxk, const unsigned int nsamp,
      const typename ODESLV_GSL<T,IVP>::Options&optGSL, std::ostream&os=std::cout );

  //! @brief Record results in file <a>bndrec</a>, with accuracy of <a>iprec</a> digits
  void record
    ( std::ofstream&bndrec, const unsigned int iprec=5 ) const;
  /** @} */

private:

  //! @brief Vector storing interval bound results (upon request only)
  std::vector< Results > _results_IA;

  //! @brief Function to convert T* array in center-radius form
  std::pair<CPPL::dcovector,CPPL::dcovector> _I2rc
    ( const unsigned int n, const T*I ) const;

  //! @brief Function to compute the Taylor series expansion for all states in T arithmetic
  bool _expansion
    ( TExpansion<T>&T_TE, const T*Ix, bool force=false );

  //! @brief Function to compute the Taylor series expansion for all states in T arithmetic with first-order derivatives
  bool _expansion
    ( TFExpansion<T>&T_TFE, const T*Ix, bool force=false );

  //! @brief Function to compute the Taylor series expansion for all states in T arithmetic with (directional) second-order derivatives
  bool _expansion
    ( TFFdExpansion<T>&T_TFFdE, const T*Ix, const T*Id, bool force=false );

  //! @brief Function to compute the Taylor series expansion for all states in TVT arithmetic
  bool _expansion
    ( TExpansion<TVT>&TVT_TE, const T*Ix, bool force=false );

  //! @brief Function to compute the Taylor series expansion for all states in TVT arithmetic
  bool _expansion
    ( TExpansion<TVT>&TVT_TE, const TVT*TMx, bool force=false );

  //! @brief Function to compute the Taylor series expansion for all states in TVT arithmetic
  bool _expansion
    ( TExpansion<TVT>&TVT_TE, TVT*Px, const T*Rx, bool force=false );

  //! @brief Function to modify Taylor expansion coefficients
  template <typename U> void _expansion_add
    ( U**f, const T*const*dfdx, const T*Rx, const T*const*d2fdx2dd=0 );

  //! @brief Function to compute all the invariants in T arithmetic with first-order derivatives
  bool _invariant
    ( FInvariant<T>&T_FI, const T*Ix, bool force=false );

  //! @brief Function to compute all the invariants in T arithmetic with (directional) second-order derivatives
  bool _invariant
    ( FFdInvariant<T>&T_FFdI, const T*Ix, const T*Id, bool force=false );

  //! @brief Function to compute all the invariants in TVT arithmetic
  bool _invariant
    ( Invariant<TVT>&TVT_I, TVT*Px, const T*Rx, bool force=false );

  //! @brief Function to modify invariants
  template <typename U> void _invariant_add
    ( U*inv, const T*dinvdx, const T*Rx, const T*d2invdx2dd=0 );

  //! @brief Function to update scaling
  template <typename U> void _scaling_update
    ( const U*Ux );

  //! @brief Function to guess stepsize
  template <typename U> double _stepsize_init
    ( TExpansion<U>&U_TE );

  //! @brief Function to validate stepsize
  template <typename U> bool _stepsize_valid
    ( TExpansion<U>&U_TE, const double&h );

  //! @brief Function to select and validate stepsize in T arithmetic
  template <typename U> double _stepsize
    ( TExpansion<T>&T_TE, TExpansion<U>&U_TE_val );

  //! @brief Function to select and validate stepsize in TVT arithmetic
  template <typename U> double _stepsize
    ( TExpansion<TVT>&TVT_TE, TExpansion<U>&U_TE_val );

  //! @brief Function to compute the Taylor-based predictor for state component <a>ix</a> within a (possibly interval valued) step-size <a>h</a> in T arithmetic
  template <typename U, typename V> U _predictor
    ( const TExpansion<U>&U_TE, const unsigned int ix, const V&h );

  //! @brief Function to compute the Taylor-based predictor for derivative component <a>jx</a> of state component <a>ix</a> within a (possibly interval valued) step-size <a>h</a> in T arithmetic
  template <typename U, typename V> std::pair<U,double> _predictor
    ( const TFExpansion<U>&U_TFE, const unsigned int ix, const unsigned jx,
      const V&h );

  //! @brief Function to compute the Taylor-based predictor for directional second derivative of state component <a>ix</a> within a (possibly interval valued) step-size <a>h</a> in T arithmetic
  template <typename U, typename V> U _predictor
    ( const TFFdExpansion<U>&U_TFFdE, const unsigned int ix, const V&h );

  //! @brief Function to initialize parameter interval bounds
  void _T_prepare
    ( const T*Ip );

  //! @brief Function to initialize parameter ellipsoidal bounds
  void _T_prepare
    ( const E&Ep );

  //! @brief Function to prepare for state bounding
  void _T_prepare();

  //! @brief Function to initialize state interval bounds
  void _T_init();

  //! @brief Function to propagate enclosures in T arithmetic
  STATUS _T_propagate
    ( const double tnext, std::ostream&os );

  //! @brief Function to display final statistic in T arithmetic
  void _T_print_stats
    ( const Stats&stats, std::ostream&os=std::cout ) const;

  //! @brief Function to initialize parameter Taylor models
  void _TVT_prepare
    ( const TVT*TMp, const E&Ep );

  //! @brief Function to prepare for Taylor model-based state bounding
  void _TVT_prepare();

  //! @brief Function to initialize state Taylor models
  void _TVT_init();

  //! @brief Function to propagate enclosures in TMT arithmetic
  STATUS _TVT_propagate
    ( const double tnext, std::ostream&os );

  //! @brief Function to display final statistic in TVT arithmetic
  void _TVT_print_stats
    ( const Stats&stats, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm_ptr
    ( const double t, const U*x, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U, typename V> void _print_interm_ptr_ref
    ( const double t, const U*x, const V&r, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm_ref
    ( const double t, const U&x, const std::string&var, std::ostream&os=std::cout ) const;

  //! @brief Function to display intermediate results
  template<typename U> void _print_interm_ptr
    ( const double t, const U*x1, const std::string&var1, const U*x2,
      const std::string&var2, std::ostream&os=std::cout ) const;

  //! @brief Private methods to block default compiler methods
  ODEBND_VAL(const ODEBND_VAL&);
  ODEBND_VAL& operator=(const ODEBND_VAL&);
};

template <typename T, typename IVP> inline std::pair<CPPL::dcovector,CPPL::dcovector>
ODEBND_VAL<T,IVP>::_I2rc
( const unsigned int n, const T*I ) const
{
  assert( n && I );
  CPPL::dcovector r(n), c(n);  
  for( unsigned int i=0; i<n; i++ ){
    r(i) = Op<T>::diam(I[i])/2.;
    c(i) = Op<T>::mid(I[i]);
  }
  return std::make_pair( r, c );
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TExpansion<T>&T_TE, const T*Ix, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, T_TE.t ) ) return false;
  ODETAYLOR<T,IVP>::T_expand( _istg, _t, Ix, _Ip,
    options.TSORDER+1, T_TE.f );
  T_TE.t = _t;
  return true;
}

template <typename T, typename IVP>
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TExpansion<TVT>&TVT_TE, const T*Ix, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, TVT_TE.t ) ) return false;
  for( unsigned int ix=0; ix<_nx; ix++ )
    _TMRHSx[ix].set( _TMRHSenv, ix, Ix[ix] );
  for( unsigned int ip=0; ip<_np; ip++ ){
    _TMRHSp[ip].set( _TMRHSenv, _nx+ip, _Ip[ip] );
    _Rxp[_nx+ip] = _TMRHSp[ip].constant();
  }
  ODETAYLOR<TVT,IVP>::T_expand( _istg, _t, _TMRHSx, _TMRHSp,
    options.TSORDER+1, TVT_TE.f );
  TVT_TE.t = _t;
  return true;
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TExpansion<TVT>&TVT_TE, const TVT*TMx, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, TVT_TE.t ) ) return false;
  ODETAYLOR<TVT,IVP>::T_expand( _istg, _t, TMx, _TMp,
    options.TSORDER+1, TVT_TE.f );
  TVT_TE.t = _t;
  return true;
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TFExpansion<T>&T_TFE, const T*Ix, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, T_TFE.t ) ) return false;
  ODETAYLOR<T,IVP>::TF_expand( _istg, _t, Ix, _Ip,
    options.TSORDER+1, T_TFE.f, T_TFE.dfdx );
  T_TFE.t = _t;
  return true;
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TFFdExpansion<T>&T_TFFdE, const T*Ix, const T*Id, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, T_TFFdE.t ) ) return false;
  ODETAYLOR<T,IVP>::TFFd_expand( _istg, _t, Ix, _Ip,
    options.TSORDER+1, T_TFFdE.f, Id, T_TFFdE.d2fdx2dd );
  T_TFFdE.t = _t;
  return true;
}

template <typename T, typename IVP>
inline bool
ODEBND_VAL<T,IVP>::_expansion
( TExpansion<TVT>&TVT_TE, TVT*Px, const T*Rx, bool force )
{
  // Update Taylor expansion results if needed
  if( !force && isequal( _t, TVT_TE.t ) ) return false;
  for( unsigned int ip=0; ip<_np; ip++ ){
    _TMRHSp[ip].set( _TMRHSenv, ip, _Ip[ip] );
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_TMRHSp[" << ip << "] = " << _TMRHSp[ip] << std::endl;
#endif
  }
  for( unsigned int jx=0; jx<_nx; jx++ ){
    // Dummy variables to initialize expansion w.r.t. remainder term
    TVT TMRHSrj( _TMRHSenv, _np+jx, Rx[jx] );
    _TMRHSx[jx].set( _TMRHSenv ).set( Px[jx], true );
    _TMRHSx[jx] += TMRHSrj;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Px[" << jx << "] = " << Px[jx] << std::endl;
    std::cout << "_TMRHSx[" << jx << "] = " << _TMRHSx[jx] << std::endl;
#endif
  }
  ODETAYLOR<TVT,IVP>::T_expand( _istg, _t, _TMRHSx, _TMRHSp,
    options.TSORDER+1, TVT_TE.f );
  TVT_TE.t = _t;
  return true;
}

template <typename T, typename IVP> template <typename U>
inline void
ODEBND_VAL<T,IVP>::_expansion_add
( U**f, const T*const*dfdx, const T*Rx, const T*const*d2fdx2dd )
{
  // Update Taylor expansion coefficients
  for( unsigned int q=0; q<=options.TSORDER; q++ ){
    for( unsigned int ix=0; ix<_nx; ix++ ){
      for( unsigned int jx=0; jx<_nx; jx++ )
        f[q][ix] += (dfdx[q][ix+jx*_nx]-Op<T>::mid(dfdx[q][ix+jx*_nx]))*Rx[jx];
      if( d2fdx2dd && d2fdx2dd[q] )
        f[q][ix] += 0.5*d2fdx2dd[q][ix];
    }
  }
  return;
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_invariant
( FInvariant<T>&T_FI, const T*Ix, bool force )
{
  // Update invariant results if needed
  if( !force && isequal( _t, T_FI.t ) ) return false;
  ODETAYLOR<T,IVP>::F_invariant( _istg, _t, Ix, _Ip, T_FI.inv, T_FI.dinvdx );
  T_FI.t = _t;
  return true;
}

template <typename T, typename IVP> 
inline bool
ODEBND_VAL<T,IVP>::_invariant
( FFdInvariant<T>&T_FFdI, const T*Ix, const T*Id, bool force )
{
  // Update invariant results if needed
  if( !force && isequal( _t, T_FFdI.t ) ) return false;
  ODETAYLOR<T,IVP>::FFd_invariant( _istg, _t, Ix, _Ip, T_FFdI.inv, Id,
    T_FFdI.d2invdx2dd );
  T_FFdI.t = _t;
  return true;
}

template <typename T, typename IVP>
inline bool
ODEBND_VAL<T,IVP>::_invariant
( Invariant<TVT>&TVT_I, TVT*Px, const T*Rx, bool force )
{
  // Update invariant results if needed
  if( !force && isequal( _t, TVT_I.t ) ) return false;
  for( unsigned int ip=0; ip<_np; ip++ ){
    _TMRHSp[ip].set( _TMRHSenv, ip, _Ip[ip] );
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "_TMRHSp[" << ip << "] = " << _TMRHSp[ip] << std::endl;
#endif
  }
  for( unsigned int jx=0; jx<_nx; jx++ ){
    // Dummy variables to initialize expansion w.r.t. remainder term
    TVT TMRHSrj( _TMRHSenv, _np+jx, Rx[jx] );
    _TMRHSx[jx].set( _TMRHSenv ).set( Px[jx], true );
    _TMRHSx[jx] += TMRHSrj;
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Px[" << jx << "] = " << Px[jx] << std::endl;
    std::cout << "_TMRHSx[" << jx << "] = " << _TMRHSx[jx] << std::endl;
#endif
  }
  ODETAYLOR<TVT,IVP>::invariant( _istg, _t, _TMRHSx, _TMRHSp, TVT_I.inv );
  TVT_I.t = _t;
  return true;
}

template <typename T, typename IVP> template <typename U>
inline void
ODEBND_VAL<T,IVP>::_invariant_add
( U*inv, const T*dinvdx, const T*Rx, const T*d2invdx2dd )
{
  // Update Taylor expansion coefficients
  for( unsigned int ix=0; ix<_nx; ix++ ){
    for( unsigned int jx=0; jx<_nx; jx++ )
      inv[ix] += (dinvdx[ix+jx*_nx]-Op<T>::mid(dinvdx[ix+jx*_nx]))*Rx[jx];
    if( d2invdx2dd )
      inv[ix] += 0.5*d2invdx2dd[ix];
  }
  return;
}

template <typename T, typename IVP> template <typename U>
inline void
ODEBND_VAL<T,IVP>::_scaling_update
( const U*Ux )
{
  // Adjust scaling heuristic
  for( unsigned int ix=0; ix<_nx; ix++ )
    _scaling[ix] = ( options.SCALING? Op<U>::diam(Ux[ix])/2. + options.ATOL/options.TOL: 1. );
}

template <typename T, typename IVP> template <typename U>
inline double
ODEBND_VAL<T,IVP>::_stepsize_init
( TExpansion<U>&U_TE )
{
//#ifdef MC__ODEBND_VAL_DEBUG_STEPSIZE
  std::cout << "U_TE.f[options.TSORDER+1]:";
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "  " << U_TE.f[options.TSORDER+1][ix];
  std::cout << std::endl;
//#endif
  double maxnorm = 0.;
  // Compute maximum norm of abs(sigma^{-1}*R(0))
  for( unsigned int ix=0; ix<_nx; ix++ )
    maxnorm = std::max( maxnorm, Op<U>::abs(U_TE.f[options.TSORDER+1][ix])/_scaling[ix] );
  return options.HREDUC * std::pow( options.TOL / maxnorm, 1./(double)options.TSORDER );
}

template <typename T, typename IVP> template <typename U>
inline bool
ODEBND_VAL<T,IVP>::_stepsize_valid
( TExpansion<U>&U_TE, const double&h )
{
//#ifdef MC__ODEBND_VAL_DEBUG_STEPSIZE
  std::cout << "U_TE.f[options.TSORDER+1]:";
  for( unsigned int ix=0; ix<_nx; ix++ )
    std::cout << "  " << U_TE.f[options.TSORDER+1][ix];
  std::cout << std::endl;
  std::cout << "current step: h=" << h << " <=? ";
//#endif
  for( unsigned int ix=0; ix<_nx; ix++ ){
    // Compute components of abs(sigma^{-1}*R(0))
    double normi = Op<U>::abs(U_TE.f[options.TSORDER+1][ix])/_scaling[ix];
    // Check stepsize validity for component ix
//#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "  " << std::pow( options.TOL / normi, 1./(double)options.TSORDER );
//#endif
    if( h > std::pow( options.TOL / normi, 1./(double)options.TSORDER ) )
      return false;
  }
  return true;
}

template <typename T, typename IVP> template <typename U>
inline double
ODEBND_VAL<T,IVP>::_stepsize
( TExpansion<T>&T_TE, TExpansion<U>&U_TE_val)
{
  // Stepsize validation
  for( double h = _stepsize_init( T_TE ); true; h *= options.HREDUC ){
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "\nStepsize attempt: h=" << h << std::endl;
#endif
    // Check minimum stepsize criterion
    if( h < options.HMIN ) return h;
    // Compute interval enclosure on full step [0,h]
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Ix_val[ix] = _predictor( T_TE, ix, h*T(0.,1.) );
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_Ix_val[ix] =" << _Ix_val[ix] << std::endl;
#endif
    }
    // Bound Taylor expansion remainder for interval enclosure on [0,h]
    try{ _expansion( U_TE_val, _Ix_val, true ); }
    catch(...){
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "Exception thrown during expansion for stepsize selection"
                << "- h=" << h << std::endl;
      int dum; std::cin>>dum;
#endif
      continue;
    }
    if( _stepsize_valid( U_TE_val, h ) ) return h;
  }
}

template <typename T, typename IVP> template <typename U>
inline double
ODEBND_VAL<T,IVP>::_stepsize
( TExpansion<TVT>&TVT_TE, TExpansion<U>&U_TE_val)
{
  // Stepsize validation
  for( double h = _stepsize_init( TVT_TE ); true; h *= options.HREDUC ){
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "\nStepsize attempt: h=" << h << std::endl;
#endif
    // Check minimum stepsize criterion
    if( h < options.HMIN ) return h;
    // Compute interval enclosure on full step [0,h]
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _Ix_val[ix] = _predictor( TVT_TE, ix, h*T(0.,1.) ).bound();
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_Ix_val[ix] =" << _Ix_val[ix] << std::endl;
#endif
    }
    // Bound Taylor expansion remainder for interval enclosure on [0,h]
    try{ _expansion( U_TE_val, _Ix_val, true ); }
    catch(...){
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "Exception thrown during expansion for stepsize selection"
                << "- h=" << h << std::endl;
      int dum; std::cin>>dum;
#endif
      continue;
    }
    if( _stepsize_valid( U_TE_val, h ) ) return h;
  }
}

template <typename T, typename IVP> template <typename U, typename V>
inline U
ODEBND_VAL<T,IVP>::_predictor
( const TExpansion<U>&U_TE, const unsigned int ix, const V&h )
{
#ifdef MC__ODETAYLOR_DEBUG
  for( unsigned int q=0; q<=options.TSORDER; q++ )
    std::cout << "ODEBND_VAL::_predictor *** f(" << ix << ")[" << q << "] = "
              << U_TE.f[q][ix] << std::endl;
#endif
  // Compute predictor of xi(_t+h)
  unsigned int q = options.TSORDER;
  U xi = U_TE.f[q][ix];
  for( ; q>0; q-- ){ xi *= h; xi += U_TE.f[q-1][ix]; }
  return xi;
}

template <typename T, typename IVP> template <typename U, typename V>
inline std::pair<U,double>
ODEBND_VAL<T,IVP>::_predictor
( const TFExpansion<U>&U_TFE, const unsigned int ix, const unsigned jx,
  const V&h )
{
  // Compute predictor of dxi(_t+h)/dx
  unsigned int q = options.TSORDER;
  U dxi = U_TFE.dfdx[q][ix+jx*_nx]-Op<U>::mid(U_TFE.dfdx[q][ix+jx*_nx]);
  double dxci = Op<U>::mid(U_TFE.dfdx[q][ix+jx*_nx]);
  for( ; q>0; q-- ){
    dxi *= h; dxi += U_TFE.dfdx[q-1][ix+jx*_nx]-Op<U>::mid(U_TFE.dfdx[q-1][ix+jx*_nx]);
    dxci *= h; dxci += Op<U>::mid(U_TFE.dfdx[q-1][ix+jx*_nx]);
  }
  return std::make_pair(dxi,dxci);
}

template <typename T, typename IVP> template <typename U, typename V>
inline U
ODEBND_VAL<T,IVP>::_predictor
( const TFFdExpansion<U>&U_TFFdE, const unsigned int ix, const V&h )
{
  // Compute predictor of d^T.d2xi(_t+h)/dx2.d
  unsigned int q = options.TSORDER;
  U d2xid = U_TFFdE.d2fdx2dd[q][ix];
  for( ; q>0; q-- ){ d2xid *= h; d2xid += U_TFFdE.d2fdx2dd[q-1][ix]; }
  return d2xid;
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_T_prepare
( const T* Ip )
{
  _T_prepare();

  // Parameters
  if( !_Ip ) _Ip = new T[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    _Ip[ip] = Ip[ip];
  _Ep.set();
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_T_prepare
( const E& Ep )
{
  _T_prepare();

  // Parameters
  if( !_Ip ) _Ip = new T[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    _Ip[ip] = T( Ep.l(ip), Ep.u(ip) );
  _Ep = Ep;
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_T_prepare()
{
  // Internal Taylor model environment reset
  if( _TMIVenv && ( _TMIVenv->nord() != options.ORDMIT 
                 || _TMIVenv->nvar() != _np ) ){
    delete[] _TMIVp; _TMIVp = 0;
    delete _TMIVenv; _TMIVenv = 0;
  }
  if( _TMRHSenv && ( _TMRHSenv->nord() != options.ORDMIT 
                  || _TMRHSenv->nvar() != _nx+_np ) ){
    delete[] _TMRHSx; _TMRHSx = 0;
    delete[] _TMRHSp; _TMRHSp = 0;
    delete _TMRHSenv; _TMRHSenv = 0;
  }

  // Define environments and size arrays as appropriate
  switch( options.WRAPMIT){
  case Options::NONE:
    _T_TE.set( options.TSORDER+1 );
    _T_TE_val.set( options.TSORDER+1 );
    break;

  case Options::ELLIPS:
  default:
    if( !_TMIVenv )  _TMIVenv  = new TMT( _np, options.ORDMIT );
    _TMIVenv->options.BOUNDER_TYPE = options.BNDMIT;
    _TMIVenv->options.SCALE_VARIABLES = false;
    _TMIVenv->options.CENTER_REMAINDER = true;
    if( !_TMIVp )    _TMIVp    = new TVT[_np];
    if( !_TMRHSenv ) _TMRHSenv = new TMT( _nx+_np, options.ORDMIT );
    _TMRHSenv->options.BOUNDER_TYPE = options.BNDMIT;
    _TMRHSenv->options.SCALE_VARIABLES = false;
    _TMRHSenv->options.CENTER_REMAINDER = true;
    if( !_TMRHSx )   _TMRHSx   = new TVT[_nx];
    if( !_TMRHSp )   _TMRHSp   = new TVT[_np];
    if( !_fxp.n )    _fxp.resize(_nx+_np,_nx+_np); _fxp.identity();
    if( !_Rxp )      _Rxp = new T[_nx+_np];
    _TVT_TE.set( options.TSORDER+1 );
    _T_TE_val.set( options.TSORDER+1 );
    break;
  }

  // States
  if( !_Ix ) _Ix = new T[_nx];
  if( !_Ix_val ) _Ix_val = new T[_nx];

  // Reset result record and statistics
  _results_IA.clear();  
  _init_stats();
  _nsteps = 0;
  ODETAYLOR<T,IVP>::reset_calls();
  ODETAYLOR<TVT,IVP>::reset_calls();
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_T_init()
{
  switch( options.WRAPMIT){

  case Options::NONE:{
    for( unsigned int ix=0; ix<_nx; ix++ )
      _Ix[ix] = IVP().IC( ix, _Ip );
    _Exp.set();
    return; }

  case Options::ELLIPS:
  default:{
    CPPL::dgematrix x0p(_nx+_np,_np); x0p.zero();
    for( unsigned int ip=0; ip<_np; ip++ ){
      // COULD BE IMPROVED BY PRECONDITIONNING WITH Q^{1/2}
      _TMIVp[ip].set( _TMIVenv, ip, _Ip[ip] );
      for( unsigned int jp=0; jp<_np; jp++ )
        x0p(_nx+jp,ip) = ( ip==jp? 1.: 0. );
      _Rxp[_nx+ip] = _TMIVp[ip].constant() + _TMIVp[ip].remainder();
    }
    for( unsigned int ix=0; ix<_nx; ix++ ){
      TVT TMIVx0i = IVP().IC( ix, _TMIVp ).center();
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "TMIVx0[" << ix << "] =" << TMIVx0i << std::endl;
#endif
      for( unsigned int ip=0; ip<_np; ip++ )
        x0p(ix,ip) = TMIVx0i.linear(ip,true); // resets to 0
      _Rxp[ix] = TMIVx0i.bound();
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_Rxp[" << ix << "] =" << _Rxp[ix] << std::endl;
#endif
    }
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Ep =" << _Ep << std::endl;
#endif
    // Case parameter set given as an interval box
    if( !_Ep.Q().n ){
      std::pair<CPPL::dcovector,CPPL::dcovector> rc = _I2rc( _np, _Ip );
      _Exp = minksum_ea( mtimes(E(rc.first),x0p), _I2rc(_nx+_np,_Rxp), options.QTOL );
    }
    // Case parameter set given as an ellipsoid
    else
      _Exp = minksum_ea( mtimes(_Ep.O(),x0p), _I2rc(_nx+_np,_Rxp), options.QTOL );
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Exp0 =" << _Exp << std::endl;
#endif
    for( unsigned int ix=0; ix<_nx; ix++ )
      _Ix[ix] = T( _Exp.l(ix), _Exp.u(ix) );
    return; }
  }
}

template <typename T, typename IVP>
inline typename ODEBND_VAL<T,IVP>::STATUS
ODEBND_VAL<T,IVP>::_T_propagate
( const double tnext, std::ostream&os )
{
  // Propagate validated enclosure until tnext
  for( _h = 0; tnext > _t; _t += _h, ++_nsteps ){
    switch( options.WRAPMIT){

    case Options::NONE:{
      // Taylor series expansion at current point
      _expansion( _T_TE, _Ix );

      // Stepsize selection and validation
      _scaling_update( _Ix );
      _h = _stepsize( _T_TE, _T_TE_val );
      if( _h < options.HMIN ) return FAILURE;
      // Modify stepsize if stage time is exceeded
      if( _t+_h > tnext ) _h = tnext-_t;

      // Enclosure propagation
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _Ix[ix] = _predictor( _T_TE, ix, _h );
        // Remainder is 2nd- and higher-order terms *plus* truncation error
        _Ix[ix] += _h*options.TOL*_scaling[ix]*T(-1.,1.);
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "Ix(" << _t+_h << ")[" << ix << "] = " << _Ix[ix] << std::endl;
#endif
      }
      _Exp.set();

      break;
     }

    case Options::ELLIPS:
    default:{
      // Taylor series expansion at current point
      _expansion( _TVT_TE, _Ix );

      // Stepsize selection and validation
      _scaling_update( _Ix );
      _h = _stepsize( _TVT_TE, _T_TE_val );
      if( _h < options.HMIN ) return FAILURE;
      // Modify stepsize if stage time is exceeded
      if( _t+_h > tnext ) _h = tnext-_t;

      // Enclosure propagation
      for( unsigned int ix=0; ix<_nx; ix++ ){
        TVT TMRHSxhi = _predictor( _TVT_TE, ix, _h ).center();
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "TMRHSxh[" << ix << "]) =" << TMRHSxhi;
#endif
        for( unsigned int jxp=0; jxp<_nx+_np; jxp++ )
          _fxp(ix,jxp) = TMRHSxhi.linear(jxp,true); // resets to 0
        // Remainder is 2nd- and higher-order terms *plus* truncation error
        _Rxp[ix] = TMRHSxhi.bound() + _h*options.TOL*_scaling[ix]*T(-1.,1.);
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "diam(_Rxp[" << ix << "]) =" << Op<T>::diam(_Rxp[ix]) << std::endl;
#endif
      }
      
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_Exp(" << _t << ") =" << _Exp << std::endl;
      std::cout << "_fxp(" << _t << ") =\n" << _fxp << std::endl;
#endif
     //_Exp = mtimes(_Exp,_fxp);
     _Exp = minksum_ea( mtimes(_Exp.O(),_fxp), _I2rc(_nx+_np,_Rxp), options.QTOL );
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "_Exp(" << _t+_h << ") =" << _Exp << std::endl;
      {int dum; std::cin >> dum;}
#endif
      if( _Exp.info() ) return FAILURE;

      // Other updates
      for( unsigned int ix=0; ix<_nx; ix++ )
        _Ix[ix] = T( _Exp.l(ix), _Exp.u(ix) );

      break;
     }
    }
  }
  return NORMAL;
}

//! @fn template <typename T, typename IVP> inline typename ODEBND_VAL<T,IVP>::STATUS ODEBND_VAL<T,IVP>::bounds(
//! const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
//! std::ostream&os )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Up</a> [input] parameter set, either niterval box or ellipsoid
//!   - <a>Ixk</a> [output] interval state enclosures at stage times
//!   - <a>Exk</a> [output] ellipsoidal state enclosures at stage times (only
//!     if ODEBND_VAL::Options::ELLIPS selected for wrapping mitigation)
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename IVP> template <typename U>
inline typename ODEBND_VAL<T,IVP>::STATUS
ODEBND_VAL<T,IVP>::bounds
( const unsigned int ns, const double*tk, const U&Up, T**Ixk, E*Exk, std::ostream&os )
{
  // Initialize trajectory integration
  _T_prepare( Up );
  _t = tk[0];
  _T_init();
  for( unsigned int ix=0; ix<_nx; ix++ ) Ixk[0][ix] = _Ix[ix];
  Exk[0] = _Exp;
  if( options.DISPLAY >= 1 )
    switch( options.WRAPMIT ){
      case Options::NONE:   _print_interm_ptr( tk[0], _Ix, "x", os ); break;
      case Options::ELLIPS: _print_interm_ref( tk[0], _Exp, "x", os ); break;
    }
  // Record initial results
  if( options.RESRECORD )
    _results_IA.push_back( Results( tk[0], _nx, _Ix ) );

  try{
    // Integrate ODEs through each stage
    for( _istg=0; _istg<ns; _istg++ ){
      if( _T_propagate( tk[_istg+1], os ) != NORMAL ){
        _final_stats();
        if( options.DISPLAY >= 1 ) _T_print_stats( stats, os );
	return FAILURE;
      }
      for( unsigned int ix=0; ix<_nx; ix++ ) Ixk[_istg+1][ix] = _Ix[ix];
      Exk[_istg+1] = _Exp;
      if( options.DISPLAY >= 1 )
        switch( options.WRAPMIT ){
          case Options::NONE:   _print_interm_ptr( tk[_istg+1], _Ix, "x", os ); break;
          case Options::ELLIPS: _print_interm_ref( tk[_istg+1], _Exp, "x", os ); break;
        }
      // Record intermediate results
      if( options.RESRECORD ) _results_IA.push_back( Results( tk[_istg+1], _nx, _Ix ) );
    }
  }
  catch(...){
    _final_stats();
    if( options.DISPLAY >= 1 ) _T_print_stats( stats, os );
    return FAILURE;
  }

  _final_stats();
  if( options.DISPLAY >= 1 ) _T_print_stats( stats, os );
  return NORMAL;
}

//! @fn template <typename T, typename IVP> template <typename U> inline typename ODEBND_VAL<T,IVP>::STATUS ODEBND_VAL<T,IVP>::hausdorff(
//! const unsigned int ns, const double*tk, const U&Up, double**Hxk,
//! const unsigned int nsamp, std::ostream&os )
//!
//! This function computes the Hausdorff distance between the interval enclosure
//! and the exact reachable set projected onto each variable:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>Up</a> [input] parameter set
//!   - <a>Hxk</a> [output] Hausdorff distance between the interval enclosure
//!     and the exact reachable set projected onto each variable, at stage times
//!   - <a>nsamp</a> [input] number of samples for each parameter
//!   - <a>os</a> [input] output stream
//! .
//! The return value is the status.
template <typename T, typename IVP> template <typename U> inline typename ODEBND_VAL<T,IVP>::STATUS
ODEBND_VAL<T,IVP>::hausdorff
( const unsigned int ns, const double*tk, const U&Up, double**Hxk,
  const unsigned int nsamp, const typename ODESLV_GSL<T,IVP>::Options&optGSL,
  std::ostream&os )
{
  _ODESLV_GSL->options = optGSL;
  int DISPLAY_SAVE = options.DISPLAY;
  options.DISPLAY = _ODESLV_GSL->options.DISPLAY = 0;
  _ODESLV_GSL->options.RESRECORD = false;

  // Compute approximate bounds
  T** Ixk = new T*[ns+1];
  E*  Exk = new E[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Ixk[is] = new T[_nx];
  try{ bounds( ns, tk, Up, Ixk, Exk, os ); }
  catch(...){;}
  unsigned int nsf = _istg;
  for( unsigned int is=nsf; is<ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++ )
      Ixk[_istg+1][ix] = T( -1e20, 1e20 );

  // Compute exact bounds 
  T** Ixk0 = new T*[ns+1];
  for( unsigned int is=0; is<ns+1; is++ ) Ixk0[is] = new T[_nx];
  if( _ODESLV_GSL->bounds( ns, tk, _Ip, Ixk0, nsamp, os ) != ODESLV_GSL<T,IVP>::NORMAL ){
    for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
    delete[] Ixk0;
    return FAILURE;
  }

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
    if( options.DISPLAY >= 1 ) _print_interm_ptr( tk[is], Hxk[is], "d", os );
  }

  for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk0[is];
  delete[] Ixk0;
  for( unsigned int is=0; is<ns+1; is++ ) delete[] Ixk[is];
  delete[] Ixk;
  delete[] Exk;

  return NORMAL;
}


template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_TVT_prepare
( const TVT* TMp, const E& Ep )
{
  _TVT_prepare();

  // Parameters
  if( !_Ip ) _Ip = new T[_np];
  if( !_TMp ) _TMp = new TVT[_np];
  for( unsigned int ip=0; ip<_np; ip++ ){
    _Ip[ip] = TMp[ip].bound();
    _TMp[ip] = TMp[ip];
  }
  if( !Ep.n() )
    _Ep.set();
  else if ( _TMenv->nvar() != Ep.n() )
    throw Exceptions( Exceptions::SIZE );
  else
    _Ep = Ep;
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_TVT_prepare()
{
  // Internal Taylor model environment reset
  if( options.ORDMIT <= 2 || ( _TMIVenv
   && ( _TMIVenv->nord() != _TMenv->nord() || _TMIVenv->nvar() != _np ) ) ){
    delete[] _TMIVp; _TMIVp = 0;
    delete _TMIVenv; _TMIVenv = 0;
  }
  if( options.ORDMIT <= 2 || ( _TMRHSenv
   && ( _TMRHSenv->nord() != _TMenv->nord() || _TMRHSenv->nvar() != _nx+_np ) ) ){
    delete[] _TMRHSx; _TMRHSx = 0;
    delete[] _TMRHSp; _TMRHSp = 0;
    delete _TMRHSenv; _TMRHSenv = 0;
  }

  // Define environments and size arrays as appropriate
  switch( options.WRAPMIT){
  case Options::NONE:
    _TVT_TE.set( options.TSORDER+1 );
    _T_TE_val.set( options.TSORDER+1 );
    break;

  case Options::ELLIPS:
  default:
    _TVT_TE.set( options.TSORDER+1 );
    _T_TE_val.set( options.TSORDER+1 );
    if( options.ORDMIT <= 2 ){
      _T_TFE.set( options.TSORDER+1 );
      _T_TFFdE.set( options.TSORDER+1 );
      if( !_Rx_TE ) _Rx_TE = new T[_nx];
    }
    else{
      if( !_TMIVenv )  _TMIVenv  = new TMT( _np, _TMenv->nord() );
      _TMIVenv->options = _TMenv->options;
      //_TMIVenv->options.SCALE_VARIABLES = false;
      //_TMIVenv->options.CENTER_REMAINDER = true;
      if( !_TMIVp )    _TMIVp    = new TVT[_np];
      if( !_TMRHSenv ) _TMRHSenv = new TMT( _nx+_np, _TMenv->nord() );
      _TMRHSenv->options = _TMenv->options;
      //_TMRHSenv->options.SCALE_VARIABLES = false;
      //_TMRHSenv->options.CENTER_REMAINDER = true;
      if( !_TMRHSx )   _TMRHSx   = new TVT[_nx];
      if( !_TMRHSp )   _TMRHSp   = new TVT[_np];
    }
    if( !_fx.n ) _fx.resize(_nx,_nx);
    if( !_Rx )   _Rx = new T[_nx];
    break;
  }

  // States
  if( !_Ix ) _Ix = new T[_nx];
  if( !_TMx ) _TMx = new TVT[_nx];
  if( !_Ix_val ) _Ix_val = new T[_nx];

  // Reset result record and statistics
  _results_IA.clear();  
  _init_stats();
  _nsteps = 0;
  ODETAYLOR<T,IVP>::reset_calls();
  ODETAYLOR<TVT,IVP>::reset_calls();
}

template <typename T, typename IVP> inline void
ODEBND_VAL<T,IVP>::_TVT_init()
{
  switch( options.WRAPMIT){

  case Options::NONE:
    for( unsigned int ix=0; ix<_nx; ix++ )
      _TMx[ix] = IVP().IC( ix, _TMp ).center();
    _Ex.set();
    return;

  case Options::ELLIPS:
  default:{
    for( unsigned int ix=0; ix<_nx; ix++ ){
      _TMx[ix] = IVP().IC( ix, _TMp ).center();
      _Rx[ix] = _TMx[ix].remainder();
      _Ix[ix] = _TMx[ix].bound();
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "TMx0[" << ix << "] =" << _TMx[ix] << std::endl;
#endif
    }
    std::pair<CPPL::dcovector,CPPL::dcovector> rc = _I2rc( _nx, _Rx );
    _Ex.set( rc.first, rc.second ); // rc.second should be all zero as Taylor model is centered
#ifdef MC__ODEBND_VAL_DEBUG
    std::cout << "Ex0 =" << _Ex << std::endl;
#endif
    // THERE MIGHT BE A BETTER WAY OF INITIALIZING IF AN ELLIPSOIDAL
    // REMAINDER BOUND IS KNOWN FOR THE PARAMETERS
    return; }
  }
}

template <typename T, typename IVP>
inline typename ODEBND_VAL<T,IVP>::STATUS
ODEBND_VAL<T,IVP>::_TVT_propagate
( const double tnext, std::ostream&os )
{
  // Propagate validated enclosure until tnext
  for( _h = 0; tnext > _t; _t += _h ){
    switch( options.WRAPMIT){

    case Options::NONE:
      // Taylor series expansion at current point
      _expansion( _TVT_TE, _TMx );

      // Stepsize selection and validation
      _scaling_update( _TMx );
      _h = _stepsize( _TVT_TE, _T_TE_val );
      if( _h < options.HMIN ) return FAILURE;
      if( _h > options.HMAX ) _h = options.HMAX;
      // Modify stepsize if stage time is exceeded
      if( _t+_h > tnext ) _h = tnext-_t;
      ++_nsteps;

      // Enclosure propagation
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _TMx[ix] = _predictor( _TVT_TE, ix, _h );
        // Remainder is 2nd- and higher-order terms *plus* truncation error
        _TMx[ix] += _h*options.TOL*_scaling[ix]*T(-1.,1.);
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "TMx(" << _t+_h << ")[" << ix << "] = " << _TMx[ix] << std::endl;
#endif
      }
      _Ex.set();

      break;

    case Options::ELLIPS:
    default:
     // In this variant a bound on the Jacobian matrix is computed and the
     // linear part is taken as the mid-point of this matrix
     if( options.ORDMIT <= 2 ){

      if( options.TSTOP || isequal( _t, _tmax ) ){

        // Start by updating scaling (if applicable)
        _scaling_update( _TMx );

        // Taylor series expansion at current point for stepsize selection,
        // first- and/or second-order derivative of coefficients, and polynomial part
        _Ex_TE = _Ex; // Copy of _Ex at expansion point for updates later on
        for( unsigned int ix=0; ix<_nx; ix++ ){
          _Rx_TE[ix] = _TMx[ix].remainder(); // Copy of Taylor remainder at expansion point for updates later on
          _TMx[ix].center().set( T(0.) ); // get centered polynomial part
        }
        _expansion( _TVT_TE, _TMx );
        if( options.ORDMIT <= 1 ){
          _expansion( _T_TFE, _Ix ); // 1st-order terms in T arithmetic to reduce CPU time
          _expansion_add( _TVT_TE.f, _T_TFE.dfdx, _Rx_TE );
          // Stepsize selection and validation
          _hmax = _h = _stepsize( _T_TFE, _T_TE_val );
        }
        else{
          _expansion( _T_TFFdE, _Ix, _Rx_TE ); // 2nd-order terms in T arithmetic to reduce CPU time
          for( unsigned int ix=0; ix<_nx; ix++ )
            _Ix[ix] = _TMx[ix].bound(); // Bound on polynomial part
          _expansion( _T_TFE, _Ix ); // 1st-order terms in T arithmetic to reduce CPU time
          _expansion_add( _TVT_TE.f, _T_TFE.dfdx, _Rx_TE, _T_TFFdE.d2fdx2dd );
          // Stepsize selection and validation
          _hmax = _h = _stepsize( _T_TFFdE, _T_TE_val );
        }

        // Return error if stepsize too small
        if( _hmax < options.HMIN ) return FAILURE;
        if( _hmax > options.HMAX ) _hmax = _h = options.HMAX;
        // Modify stepsize if stage time is exceeded
        _tmax = _t+_hmax;
        ++_nsteps;
//#ifdef MC__ODEBND_VAL_DEBUG
       std::cout << "\nCurrent step: hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
       { int dum; std::cin >> dum; }
//#endif
      }

      // By-pass stepsize validation if previous integration step can be
      // carried out further and ODE RHS is continuous
      else{
        _h = _hmax;
      }

      // Limit stepsize so as to stop at stage time tnext
      // Do not use _t since might be different from Taylor expansion time
      if( _tmax > tnext ) _h = tnext-_TVT_TE.t;
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "\nCurrent step: t=" << _t << "  h=" << _h << "  tnext=" << tnext
                << "  hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
#endif

      // Enclosure propagation
      for( unsigned int ix=0; ix<_nx; ix++ ){
        _TMx[ix] = _predictor( _TVT_TE, ix, _h ).center();
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "_TMx[" << ix << "] =" << _TMx[ix];
#endif
        // Remainder := Remainder of (modified) Taylor model propagation
        _Rx[ix] = _TMx[ix].remainder();
        // Linear part for ellipsoidal remainder update
        for( unsigned int jx=0; jx<_nx; jx++ )
          _fx(ix,jx) = _predictor( _T_TFE, ix, jx, _h ).second;
        // Remainder += truncation error
        _Rx[ix] += _h*options.TOL*_scaling[ix]*T(-1.,1.);
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "_Rx[" << ix << "] =" << _Rx[ix]
                  << "  diam(_Rx[" << ix << "]) =" << Op<T>::diam(_Rx[ix])
                  << "  midp(_Rx[" << ix << "]) =" << Op<T>::mid(_Rx[ix])
                  << std::endl;
#endif
      }
     }

     // In this variant a Taylor model in the joint state-parameter and
     // of the same order as the parameter Taylor model is computed
     else{

      if( options.TSTOP || isequal( _t, _tmax ) ){

        // Start by updating scaling (if applicable)
        _scaling_update( _TMx );

        // Taylor series expansion at current point
        for( unsigned int ix=0; ix<_nx; ix++ ){
          _Rx[ix] = _TMx[ix].center().remainder(); // get centered remainder
          _TMx[ix].set( T(0.) ); // cancel remainder term to get polynomial part
#ifdef MC__ODEBND_VAL_DEBUG
          std::cout << "_TMx[" << ix << "] =" << _TMx[ix];
          std::cout << "_Rx[" << ix << "] =" << _Rx[ix] << std::endl;
#endif
        }
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "Expansion performed: t=" << _t << std::endl;
#endif
        _expansion( _TVT_TE, _TMx, _Rx );
        _Ex_TE = _Ex; // Copy of _Ex at expansion point since reused for update

        // Stepsize selection and validation
        _hmax = _h = _stepsize( _TVT_TE, _T_TE_val );
        if( _hmax < options.HMIN ) return FAILURE;
        if( _hmax > options.HMAX ) _hmax = _h = options.HMAX;
        // Modify stepsize if stage time is exceeded
        _tmax = _t+_hmax;
        ++_nsteps;
      }

      // By-pass stepsize validation if previous integration step can be
      // carried out further and ODE RHS is continuous
      else{
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "Expansion by-passed: t=" << _t << "  tmax=" << _tmax << "  hmax=" << _hmax << std::endl;
#endif
        _h = _hmax;
      }

      // Limit stepsize so as to stop at stage time tnext
      // Do not use _t since might be different from Taylor expansion time
      if( _tmax > tnext ) _h = tnext-_TVT_TE.t;
#ifdef MC__ODEBND_VAL_DEBUG
      std::cout << "\nCurrent step: t=" << _t << "  h=" << _h << "  tnext=" << tnext
                << "  hmax=" << _hmax << "  tmax=" << _tmax << std::endl;
#endif

      // Enclosure propagation
      for( unsigned int ix=0; ix<_nx; ix++ ){
        TVT TMRHSxhi = _predictor( _TVT_TE, ix, _h ).center();
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "TMRHSxh[" << ix << "] =" << TMRHSxhi;
#endif
        // Extract Taylor model polynomial part in p and set to 0
        TMRHSxhi.get( _TMx[ix].set(_TMenv), true );
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "_TMx[" << ix << "] =" << _TMx[ix];
        std::cout << "TMRHSxh[" << ix << "] =" << TMRHSxhi;
#endif
        // Extract first-order terms and set to 0
        for( unsigned int jx=0; jx<_nx; jx++ )
          _fx(ix,jx) = TMRHSxhi.linear(_np+jx,true);
        // Remainder is higher-order terms *plus* truncation error
        _Rx[ix] = TMRHSxhi.bound() + _h*options.TOL*_scaling[ix]*T(-1.,1.);
#ifdef MC__ODEBND_VAL_DEBUG
        std::cout << "_Rx[" << ix << "] =" << _Rx[ix]
                  << "  diam(_Rx[" << ix << "]) =" << Op<T>::diam(_Rx[ix])
                  << "  midp(_Rx[" << ix << "]) =" << Op<T>::mid(_Rx[ix])
                  << std::endl;
#endif
      }
    }

#ifdef MC__ODEBND_VAL_DEBUG
     std::cout << "_fx(" << _TVT_TE.t << ") =\n" << _fx << std::endl;
     std::cout << "_Ex(" << _TVT_TE.t << ") =" << _Ex_TE << std::endl;
#endif
     // Correct stepsize in case stepsize selection/validation was by-passed
     _h -= (_t-_TVT_TE.t); 
#ifdef MC__ODEBND_VAL_DEBUG
     std::cout << "\nAdjusted step: h=" << _h << "  texp=" << _TVT_TE.t << std::endl;
#endif
     //_Ex = mtimes(_Ex_TE,_fx);
     _Ex = minksum_ea( mtimes(_Ex_TE,_fx), _I2rc(_nx,_Rx), options.QTOL );
     for( unsigned int ix=0; ix<_nx; ix++ )
       _TMx[ix] += _Ex_TE.c(ix);
     _Ex = _Ex.O();
#ifdef MC__ODEBND_VAL_DEBUG
     std::cout << "_Ex(" << _t+_h << ") =" << _Ex << std::endl;
     {int dum; std::cin >> dum;}
#endif
     if( _Ex.info() ) return FAILURE;

/*
     // Contract ellipsoidal remainder using ODE invariants
     if( options.USEINV ){
       for( unsigned int ix=0; ix<_nx; ix++ ) _TMx[ix].set( T(0.) );
       for( unsigned int ii=0; ii<_ni; ii++ ){
         TVT TMinv = INV( ii, _TMp, _TMx, _t, _istg );
	 if( Op<T>::diam(TMinv.T()) < machprec() ){
	   
	 }
       }
     }
*/     
     // Update interval remainder
     for( unsigned int ix=0; ix<_nx; ix++ ){
       _TMx[ix].set( T( _Ex.l(ix), _Ex.u(ix) ) );
       _Ix[ix] = _TMx[ix].bound();
     }
 
     break;
    }
  }

  return NORMAL;
}

//! @fn template <typename T, typename IVP> inline typename ODEBND_VAL<T,IVP>::STATUS ODEBND_VAL<T,IVP>::bounds(
//! const unsigned int ns, const double*tk, const T*Ip, T**Ixk,
//! std::ostream&os )
//!
//! This function computes an interval enclosure of the reachable set of 
//! the parametric ODEs defined in IVP using equally spaced samples:
//!   - <a>ns</a> [input] number of time stages
//!   - <a>tk</a> [input] stage times, including the initial time
//!   - <a>TMp</a> [input] Taylor model parameter enclosure
//!   - <a>Ep</a> [input] ellipsoidal remainders in Taylor model parameter enclosure
//!   - <a>Txk</a> [output] Taylor model state enclosures at stage times
//!   - <a>Exk</a> [output] ellipsoidal remainders in Taylor model state enclosures at stage times
//!   - <a>os</a> [input] output stream (default: std::cout)
//! .
//! The return value is the status.
template <typename T, typename IVP>
inline typename ODEBND_VAL<T,IVP>::STATUS
ODEBND_VAL<T,IVP>::bounds
( const unsigned int ns, const double*tk, const TVT*TMp, const E&Ep,
  TVT**TMxk, E*Exk, std::ostream&os )
{
  // Check Taylor model compatibility and size
  unsigned int kp=_np;
  for( unsigned int ip=0; ip<_np && kp==_np; ip++ )
    if( TMp[ip].env() ) kp = ip;
  if( kp==_np || TMp[kp].env()->nvar()!=_np ) return FATAL;
  _TMenv = TMp[kp].env();  

  // Initialize trajectory integration
  _TVT_prepare( TMp, Ep );
  _t = _tmax = tk[0];
  _TVT_init();
  for( unsigned int ix=0; ix<_nx; ix++ ) TMxk[0][ix] = _TMx[ix];
  Exk[0] = _Ex;
  if( options.DISPLAY >= 1 )
    switch( options.WRAPMIT ){
      case Options::NONE:   _print_interm_ptr( tk[0], _TMx, "x", os ); break;
      case Options::ELLIPS: _print_interm_ptr_ref( tk[0], _TMx, _Ex, "x", os ); break;
    }
  // Record initial results
  if( options.RESRECORD )
    _results_IA.push_back( Results( tk[0], _nx, _TMx ) );

  try{
    // Integrate ODEs through each stage
    for( _istg=0; _istg<ns; _istg++ ){
      if( _TVT_propagate( tk[_istg+1], os ) != NORMAL ){
        _final_stats();
        if( options.DISPLAY >= 1 ) _TVT_print_stats( stats, os );
	return FAILURE;
      }
      for( unsigned int ix=0; ix<_nx; ix++ ) TMxk[_istg+1][ix] = _TMx[ix];
      Exk[_istg+1] = _Ex;
      if( options.DISPLAY >= 1 )
        switch( options.WRAPMIT ){
          case Options::NONE:   _print_interm_ptr( tk[_istg+1], _TMx, "x", os ); break;
          case Options::ELLIPS: _print_interm_ptr_ref( tk[_istg+1], _TMx, _Ex, "x", os ); break;
        }
      // Record intermediate results
      if( options.RESRECORD ) _results_IA.push_back( Results( tk[_istg+1], _nx, _TMx ) );
    }
  }
  catch(...){
    _final_stats();
    if( options.DISPLAY >= 1 ) _TVT_print_stats( stats, os );
    return FAILURE;
  }

  _final_stats();
  if( options.DISPLAY >= 1 ) _TVT_print_stats( stats, os );
  return NORMAL;
}

template <typename T, typename IVP>
inline void
ODEBND_VAL<T,IVP>::_T_print_stats
( const Stats&stats, std::ostream&os ) const
{
  // Final statistics
  os << " No TIME STEPS      " << _nsteps << std::endl;
  os << " No RHS EXPANSIONS  ";
  switch( options.WRAPMIT){
  case Options::NONE:
    os << "IA:" << ODETAYLOR<T,IVP>::T_calls() << std::endl;
    break;
  case Options::ELLIPS:
    os << "TM:" << ODETAYLOR<TVT,IVP>::T_calls() << std::endl;
    break;
  }
  os << " CPU TIME (SEC)     " << std::fixed << std::left
     << std::setprecision(5) << stats.cputime << std::endl;
  return;
}

template <typename T, typename IVP>
inline void
ODEBND_VAL<T,IVP>::_TVT_print_stats
( const Stats&stats, std::ostream&os ) const
{
  // Final statistics
  os << " No TIME STEPS      " << _nsteps << std::endl;
  os << " No RHS EXPANSIONS  ";
  switch( options.WRAPMIT){
  case Options::NONE:
    os << "IA:" << ODETAYLOR<T,IVP>::T_calls()
       << "  TM:" << ODETAYLOR<TVT,IVP>::T_calls() << std::endl;
    break;
  case Options::ELLIPS:
    os << "IA:" << ODETAYLOR<T,IVP>::T_calls()
       << "  TM:" << ODETAYLOR<TVT,IVP>::T_calls() << std::endl;
    break;
  }
  os << " CPU TIME (SEC)     " << std::fixed << std::left
     << std::setprecision(5) << stats.cputime << std::endl;
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_VAL<T,IVP>::_print_interm_ptr
( const double t, const U*x, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  return;
}

template <typename T, typename IVP> template<typename U, typename V> inline void
ODEBND_VAL<T,IVP>::_print_interm_ptr_ref
( const double t, const U*x, const V&r, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  for( unsigned int ix=0; ix<_nx; ix++ )
    os << " " << var.c_str() << "[" << ix << "] = " << x[ix] << std::endl;
  os << " " << "R" << var.c_str() << " =" << r << std::endl;
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_VAL<T,IVP>::_print_interm_ref
( const double t, const U&x, const std::string&var, std::ostream&os ) const
{
  os << " @t = " << std::scientific << std::setprecision(4)
                 << std::left << t << " :" << std::endl;
  os << " " << var.c_str() << " =" << x << std::endl;
  return;
}

template <typename T, typename IVP> template<typename U> inline void
ODEBND_VAL<T,IVP>::_print_interm_ptr
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
ODEBND_VAL<T,IVP>::record
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
