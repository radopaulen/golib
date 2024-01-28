// Copyright (C) 2012, 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_ELLIPSOID Ellipsoidal Calculus
\author Benoit Chachuat

The mc::Ellipsoid class provides a minimal set of operators and functions for ellipsoids together with constructors and data accessing functions. An ellipsoid with center \f$c \in \mathbb R^n\f$ and shape matrix \f$Q \in \mathbb S_{+}^n\f$ is defined as 
\f{align*}
  \mathcal E(c,Q) := & \left\{ \left. c + Q^\frac{1}{2} v \ \right| \ \exists v \in \mathbb R^{n}: \ v^T v \leq 1 \right\} \subseteq \mathbb R^{n} .
\f}

In order to define the ellipsoid \f$\mathcal E(c,Q)\f$ with
\f{align*}
  c = & \left(\begin{array}{c} 1.\\-1.\end{array}\right),\ \text{and} & 
  Q = & \left(\begin{array}{cc} 0.5 & -0.4 \\ -0.4 & 0.5 \end{array}\right).
\f}

we proceed as follows:

\code
    const unsigned int n = 2;
    CPPL::dcovector c(n); CPPL::csymatrix Q(n);
    c(0) =  1.;  Q(0,0) =  0.5;
    c(1) = -1.;  Q(1,0) = -0.4;  Q(1,1) = 0.5 ;
    mc::Ellipsoid E( Q, c ); 
\endcode

This ellipsoid is simply displayed as:

\code
    std::cout << E << std::endl;
\endcode

In the present case, the following information is displayed:

\verbatim
center:
 1.00000e+00
 5.00000e-01

shape:
 5.00000e-01 {-4.00000e-01}
 -4.00000e-01  5.00000e-01 
\endverbatim


\section sec_ELL_opt How Are the Options Set for Ellipsoidal Calculus?

The options are defined in the structure mc::Ellipsoid::Options. Current options are as follows:

<TABLE border="1">
<CAPTION><EM>Options in mc::Ellipsoid: name, type and
description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>PSDCHK</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether or not to check positive semi-definiteness of shape matrices
     <TR><TH><tt>PSDTOL</tt> <TD><tt>double</tt> <TD>1e2*MACHPREC
         <TD>Absolute tolerance for positive semi-definiteness check of shape matrices
     <TR><TH><tt>RKTOLA</tt> <TD><tt>double</tt> <TD>MACHPREC
         <TD>Absolute tolerance for rank and regularization of shape matrices
     <TR><TH><tt>RKTOLR</tt> <TD><tt>double</tt> <TD>MACHPREC*1e6
         <TD>Relative tolerance for rank and regularization of shape matrices
</TABLE>

The mc::Ellipsoid class has a public static member called mc::Ellipsoid::options that can be used to set/modify the options; for example:

\code
      mc::Ellipsoid::options.PSDCHK = true;
\endcode

\section sec_ELL_err Errors Encountered during Ellipsoidal Calculus?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::Ellipsoid::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a natural interval extension, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will stop.

Possible errors encountered during ellipsoidal calculus:

<TABLE border="1">
<CAPTION><EM>Errors with mc::Ellipsoid</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Non positive-semi definite shape matrix
</TABLE>

\section sec_ELL_refs References

Kurzhanskiy, A., A.. and P. Varaiya, <A href="http://code.google.com/p/ellipsoids">"Ellipsoidal Toolbox"</A>, Technical Report UCB/EECS-2006-46, EECS Department, University of California, Berkeley, May 2006.

*/

#ifndef MC__ELLIPSOID_HPP
#define MC__ELLIPSOID_HPP

#include <iostream>
#include <iomanip>
#include <stdarg.h>
#include <cassert>

#include "cpplapack_dsysv.hpp"
#include "mcfunc.hpp"

#undef  MC__ELLIPSOID_DEBUG_SQRT
#undef  MC__ELLIPSOID_DEBUG_HPINTERSECTION
#undef  MC__ELLIPSOID_DEBUG_INTERSECTION_EA
#undef  MC__ELLIPSOID_DEBUG


namespace mc
{
//! @brief C++ class for ellipsoidal calculus
////////////////////////////////////////////////////////////////////////
//! mc::Ellipsoid is a C++ class for ellipsoidal calculus. Round-off
//! errors are not accounted for in the computations (non-verified
//! implementation).
////////////////////////////////////////////////////////////////////////
class Ellipsoid
////////////////////////////////////////////////////////////////////////
{
  friend std::ostream& operator<<
    ( std::ostream&, const Ellipsoid& );
  friend Ellipsoid ell_unitball
    ( const unsigned int n );
  friend Ellipsoid mtimes
    ( const Ellipsoid&E, const CPPL::dgematrix&A, const CPPL::dcovector&b );
  friend Ellipsoid minksum_ea
    ( const std::vector<Ellipsoid>&E, const CPPL::dcovector&D );
  friend std::vector<Ellipsoid> minksum_ea
    ( const std::vector<Ellipsoid>&E, const std::vector<CPPL::dcovector>&D );
  friend Ellipsoid minksum_ea
    ( const Ellipsoid&E, const std::pair<CPPL::dcovector,CPPL::dcovector>&I,
      const double&TOL, const double&EPS );
  friend Ellipsoid inv
    ( const Ellipsoid&E );
  friend std::vector<Ellipsoid> inv
    ( const std::vector<Ellipsoid>&E );
  friend double dist
    ( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP );
  friend Ellipsoid hpintersection
    ( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP );
  friend Ellipsoid hpintersection
    ( const Ellipsoid&E, const std::vector< std::pair<CPPL::dcovector,double> >&HP );
  friend Ellipsoid intersection_ea
    ( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP,
  const double&TOL );
  friend Ellipsoid intersection_ea
    ( const Ellipsoid&E, const std::vector< std::pair<CPPL::dcovector,double> >&HP,
  const double&TOL );

private:

  //! @brief Ellipsoid shape matrix (dense symmetric format)
  CPPL::dsymatrix _Q;
  //! @brief Ellipsoid center
  CPPL::dcovector _c;
  //! @brief Pointer to eigenvalues
  std::pair< CPPL::dcovector, CPPL::dgematrix > _eigQ;
  //! @brief Pointer to square root of shape matrix
  CPPL::dsymatrix _sqrtQ;
  //! @brief Pointer to singular values
  std::pair< CPPL::dcovector, std::pair<CPPL::dgematrix,CPPL::dgematrix> > _svdQ;
  //! @brief Pointer to inverse of shape matrix
  CPPL::dsymatrix _invQ;

  //! @brief Whether shape matrix was checked for postive semi-definiteness
  bool _PSDchecked;

  //! @brief Delete pointers in class
  void _reset_auxiliary()
    {
      _eigQ.first.clear(); _eigQ.second.clear();
      _sqrtQ.clear();
      _svdQ.first.clear(); _svdQ.second.first.clear(); _svdQ.second.second.clear();
      _invQ.clear();
      _PSDchecked = false;
    }

public:

  /** @defgroup ELLIPSOID Ellipsoidal Calculus
   *  @{
   */
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in \mathbb R^n\f$, returns \f$(\mathcal E + d)\f$
  Ellipsoid& operator+=
    ( const CPPL::dcovector& d )
    {
      _c += d;
      return *this;
    }
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in \mathbb R^n\f$, returns \f$(\mathcal E - d)\f$
  Ellipsoid& operator-=
    ( const CPPL::dcovector& d )
    {
      _c -= d;
      return *this;
    }
  //! @brief Assignment operator
  Ellipsoid& operator=
    ( const Ellipsoid&E )
    {
      _Q = E._Q;
      _c = E._c;
      _reset_auxiliary();
      return *this;
    }

  //! @brief Ellipsoid options
  static struct Options
  {
    //! @brief Constructor
    Options():
      PSDCHK(true), PSDTOL(machprec()*1e2), RKTOLA(machprec()), RKTOLR(1e6*mc::machprec()),
      ROOTTOL(1e-10), ROOTSECANT(false), ROOTMAXIT(0)
      {}
    //! @brief Whether or not to check positive semi-definiteness of shape matrices (default=false)
    bool PSDCHK;
    //! @brief Absolute tolerance for positive semi-definiteness check of shape matrix (default=1e2*mc::machprec())
    double PSDTOL;
    //! @brief Absolute tolerance for rank and regularization of shape matrix (default=mc::machprec())
    double RKTOLA;
    //! @brief Relative tolerance for rank and regularization of shape matrix (default=1e6*mc::machprec())
    double RKTOLR;
    //! @brief Absolute stopping tolerance for root-finding method (objective function value less than ROOTTOL; default=1e-8)
    double ROOTTOL;
    //! @brief Whether to use the secant method for root finding (default=true)
    bool ROOTSECANT;
    //! @brief Maximum number of iteration for root-finding method (default=0 - no maximum)
    unsigned int ROOTMAXIT;
  } options;

  //! @brief Ellipsoid exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for Interval exception handling
    enum TYPE{
      NONPSD=1,	//!< Non positive-semi definite shape matrix
      LAPACK,	//!< Linear algebra routine in LAPACK failed
      ROOT,	//!< Root-finding routine failed
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Return error flag
    int ierr(){ return _ierr; }
    //! @brief Return error description
    std::string what(){
      switch( _ierr ){
      case NONPSD:
        return "Not positive semi-definite shape matrix";
      case LAPACK: default:
        return "Failure in LAPACK components";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Default constructor (needed to declare arrays of Ellipsoid class)
  Ellipsoid
    ():
    _Q(), _c(), _PSDchecked(false)
    {}

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$
  Ellipsoid
    ( const CPPL::dsymatrix& Q, const CPPL::dcovector& c=CPPL::dcovector() ):
    _Q(Q), _c(c), _PSDchecked(false)
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( !c.l || Q.n == c.l );
#endif
      if( options.PSDCHK ){
        bool PSD; _isPSD( Q, _eigQ.first, PSD );
        if( !PSD ){
	  std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
      if( _c.l != _Q.n ){ _c.resize(_Q.n); _c.zero(); }
    }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$ (lower triangular part stored contiguously and columnwise)
  Ellipsoid
    ( const unsigned int n, const double*Q, const double*c=0 ):
    _PSDchecked(false)
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( n && Q );
#endif
      _Q.resize(n);
      for( unsigned int j=0,k=0; j<n; j++ )
        for( unsigned int i=j; i<n; i++ )
          _Q(i,j) = Q[k++];
      _c.resize(_Q.n);
      for( unsigned int i=0; i<n; i++ )
        _c(i) = ( c? c[i]: 0. );

      if( options.PSDCHK ){
        bool PSD; _isPSD( _Q, _eigQ.first, PSD );
        if( !PSD ){
	  std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
    }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ enclosing interval vector of radius \f$r\f$ centered at \f$c\f$
  Ellipsoid
    ( const CPPL::dcovector& r, const CPPL::dcovector& c=CPPL::dcovector() ):
    _Q(r.l), _c(c), _PSDchecked(false)
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( !c.l || r.l == c.l );
#endif
      if( options.PSDCHK ){
        for( unsigned int i=0; i<r.l; i++ ) if( r(i) < 0 ){
          std::cout << "Interval radius: " << r(i) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
        }
        _PSDchecked = true;
      }
      _Q.zero();
      double nrm2_r = CPPL::nrm2( r );
      for( unsigned int i=0; i<r.l; i++ ) _Q(i,i) = r(i)*nrm2_r;
      if( _c.l != _Q.n ){ _c.resize(_Q.n); _c.zero(); }
    }

  //! @brief Copy constructor
  Ellipsoid
    ( const Ellipsoid&E ):
    _Q(E._Q), _c(E._c), _PSDchecked(false)
    {}

  //! @brief Destructor
  ~Ellipsoid()
    {}

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$
  Ellipsoid& set
    ( const CPPL::dsymatrix& Q=CPPL::dsymatrix(), const CPPL::dcovector& c=CPPL::dcovector() )
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( !c.l || Q.n == c.l );
#endif
      _reset_auxiliary();
      if( Q.n && options.PSDCHK ){
        bool PSD; _isPSD( Q, _eigQ.first, PSD );
        if( !PSD ){
	  std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
      _Q = Q;
      _c = c;
      if( _c.l != _Q.n ){ _c.resize(_Q.n); _c.zero(); }
      return *this;
    }

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$ (lower triangular part stored contiguously and columnwise)
  Ellipsoid& set
    ( const unsigned int n, const double*Q, const double*c=0 )
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( n && Q );
#endif
      _Q.resize(n);
      for( unsigned int j=0,k=0; j<n; j++ )
        for( unsigned int i=j; i<n; i++ )
          _Q(i,j) = Q[k++];
      _c.resize(_Q.n);
      for( unsigned int i=0; i<n; i++ )
        _c(i) = ( c? c[i]: 0. );
 
      _reset_auxiliary();
      if( options.PSDCHK ){
        bool PSD; _isPSD( _Q, _eigQ.first, PSD );
        if( !PSD ){
	  std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
      return *this;
    }

  //! @brief Define ellipsoid of dimension \f$n\f$ enclosing interval vector of radius \f$r\f$ centered at \f$c\f$
  Ellipsoid& set
    ( const CPPL::dcovector& r, const CPPL::dcovector& c=CPPL::dcovector() )
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( !c.l || r.l == c.l );
#endif
      _reset_auxiliary();
      if( options.PSDCHK ){
        for( unsigned int i=0; i<r.l; i++ ) if( r(i) < 0 ){
          std::cout << "Interval radius: " << r(i) << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
        }
        _PSDchecked = true;
      }
      _Q.resize( r.l );
      _Q.zero();
      double nrm2_r = nrm2( r );
      for( unsigned int i=0; i<r.l; i++ ) _Q(i,i) = r(i)*nrm2_r;
      _c = c;
      if( _c.l != _Q.n ){ _c.resize(_Q.n); _c.zero(); }
      return *this;
    }

  //! @brief Extend dimension by one, by appending row <a>Qi</a> to shape matrix and entry <a>ci</a> to center
  Ellipsoid& extend
    ( const CPPL::drovector& Qi, const double& ci=0. )
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( Qi.l == _Q.n+1 );
#endif
      _reset_auxiliary();
      CPPL::dsymatrix Qext(_Q.n+1);
      CPPL::dcovector cext(_Q.n+1);
      for( unsigned int i=0; i<_Q.n; i++ ){
        cext(i) = _c(i);
        for( unsigned int j=0; j<=i; j++ )
          Qext(i,j) = _Q(i,j);
      }
      cext(_Q.n) = ci;
      for( unsigned int i=0; i<=_Q.n; i++ )
        Qext(_Q.n,i) = Qi(i);
      _c = cext;
      _Q = Qext;
      return *this;
    }

  //! @brief Recenter ellipsoid at the origin by canceling out the centre
  Ellipsoid O() const
    {
      return Ellipsoid(_Q);
    }
  //! @brief Return dimension of ellipsoid
  unsigned int n() const
    {
      return _Q.n;
    }
  //! @brief Return center of ellipsoid
  const CPPL::dcovector& c() const
    {
      return _c;
    }
  //! @brief Return center of ellipsoid
  CPPL::dcovector& c()
    {
      return _c;
    }
  //! @brief Return center coefficient
  double c
    ( unsigned int i ) const
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i<_Q.n );
#endif
      return _c(i);
    }
  //! @brief Return shape matrix of ellipsoid
  const CPPL::dsymatrix& Q() const
    {
      return _Q;
    }
  //! @brief Return shape matrix coefficient
  double Q
    ( unsigned int i, unsigned int j ) const
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i<_Q.n && j<_Q.n );
#endif
      return _Q(i,j);
    }
  //! @brief Return/set shape matrix coefficient
  double& Q
    ( unsigned int i, unsigned int j )
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i<_Q.n && j<_Q.n );
#endif
      return _Q(i,j);
    }

  //! @brief Return trace of shape matrix
  double trQ() const
    {
      if( !_Q.n ) return 0.;
      double tr(_Q(0,0));
      for( unsigned int i=1; i<_Q.n; i++ ) tr += _Q(i,i);
      return tr;
    }

  //! @brief Return eigenvalues and eigenvectors of shape matrix
  const std::pair< CPPL::dcovector, CPPL::dgematrix >& eigQ()
    {
      if( !_eigQ.first.l || !_eigQ.second.n )   
        _eigen( _Q, _eigQ.first, _eigQ.second );
      return _eigQ;
    }

  //! @brief Return square root of shape matrix
  bool psdQ()
    {
      bool PSD; _isPSD( _Q, _eigQ.first, PSD );
      return PSD;
    }

  //! @brief Return square root of shape matrix
  const CPPL::dsymatrix& sqrtQ
    ( const bool complete=false )
    {
      if( !_sqrtQ.n ){
        if( !_eigQ.first.l || !_eigQ.second.n )   
          _eigen( _Q, _eigQ.first, _eigQ.second );
        _sqrt( _Q, _eigQ.first, _eigQ.second, _sqrtQ, _PSDchecked );
      }
      if( complete ) _sqrtQ.complete();
      return _sqrtQ;
    }

  //! @brief Return rank of shape matrix
  unsigned int rankQ()
    {
      svdQ();
      unsigned int rank = _Q.n;
      for( unsigned int i=0; i<_Q.n; i++, rank-- )
        if( _svdQ.first(_Q.n-i-1) > options.RKTOLA
         && _svdQ.first(_Q.n-i-1) > options.RKTOLR*_svdQ.first(0) ) break;
      return rank;
    }

  //! @brief Return singular value decomposition of shape matrix
  const std::pair< CPPL::dcovector, std::pair<CPPL::dgematrix,CPPL::dgematrix> >& svdQ()
    {
      if( !_svdQ.first.l || !_svdQ.second.first.n || !_svdQ.second.second.n )
        _svd(  _Q, _svdQ.first, _svdQ.second.first, _svdQ.second.second );
      return _svdQ;
    }

  //! @brief Return pointer to regularized shape matrix
  const CPPL::dsymatrix& regQ()
    {
      const unsigned int r = rankQ();
#ifdef MC__ELLIPSOID_DEBUG_REGQ
      std::cout << "***Ellipsoid::regQ -- r = " << r << std::endl;
#endif
      if( r == _Q.n ) return _Q;

      CPPL::dgbmatrix E(_Q.n,_Q.n,0,0); E.zero();
      const double eps = std::max( options.RKTOLA, options.RKTOLR*svdQ().first(0) );
      for( unsigned int i=0; i<_Q.n-r; i++ ) E(r+i,r+i) = eps;
      CPPL::dgematrix UEUT = svdQ().second.first * E * t(svdQ().second.first);
#ifdef MC__ELLIPSOID_DEBUG_REGQ
      std::cout << "***Ellipsoid::regQ -- UEUT = " << UEUT << std::endl;
#endif
      for( unsigned int i=r; i<_Q.n; i++ )
        _svdQ.first(i) += eps;
      for( unsigned int i=0; i<_Q.n; i++ ){
        for( unsigned int j=0; j<=i; j++ )
          _Q(i,j) += 0.5*(UEUT(i,j)+UEUT(j,i));
      }
#ifdef MC__ELLIPSOID_DEBUG_REGQ
      _reset_auxiliary();
#endif
      return _Q;
    }

  //! @brief Return pointer to regularized shape matrix
  const CPPL::dsymatrix& invQ()
    {
      if( !_invQ.n ) _inv( _Q, _invQ );
      return _invQ;
    }

  //! @brief Computes an orthogonal matrix rotating the vector <a>x</a> so that it is parallel to the vector <a>v</a>
  CPPL::dgematrix align
    ( const CPPL::dcovector&v, const CPPL::dcovector&x ) const
    {
       if( v.l != x.l )
        return CPPL::dgematrix();

       CPPL::dgematrix vmat( v.l, 1 ), xmat( x.l, 1 );
       for( unsigned int i=0; i<v.l; i++ ){
         vmat( i, 0 ) = v( i );
         xmat( i, 0 ) = x( i );
       }
       CPPL::dcovector Sv, Sx;
       CPPL::dgematrix Uv, Ux, VTv, VTx;
       _svd(  vmat, Sv, Uv, VTv );
       _svd(  xmat, Sx, Ux, VTx );
       return (Uv * VTv(0,0)) * t( Ux * VTx(0,0) );
    }

  //! @brief Return lower bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double l
    ( const unsigned int i ) const
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i>=0 && i<_Q.n );
#endif
      return _c(i) - std::sqrt(std::max(0.,_Q(i,i)));
    }
  //! @brief Return upper bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double u
    ( const unsigned int i ) const
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i>=0 && i<_Q.n );
#endif
      return _c(i) + std::sqrt(std::max(0.,_Q(i,i)));
    }
  //! @brief Return maximum radius for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double r
    ( const unsigned int i ) const
    {
#ifdef MC__ELLIPSOID_DEBUG
      assert( i>=0 && i<_Q.n );
#endif
      return std::sqrt(std::max(0.,_Q(i,i)));
    }
  /** @} */
  
private:

  //! @brief Wrapper to LAPACK function <TT>dsyev_</TT> doing eigenvalue decomposition of a symmetric matrix
  static void _eigen
    ( const CPPL::dsymatrix&Q, CPPL::dcovector&D, CPPL::dgematrix&U );

  //! @brief Wrapper to LAPACK function <TT>dsyev_</TT> doing eigenvalue decomposition of a symmetric matrix (eigenvalues only)
  static void _eigen
    ( const CPPL::dsymatrix&Q, CPPL::dcovector&D );

  //! @brief Check if a symmetric matrix is positive definite
  static void _isPSD
    ( const CPPL::dsymatrix&Q, CPPL::dcovector&D, bool&PSD );

  //! @brief Compute the square-root of a symmetric matrix
  static void _sqrt
    ( const CPPL::dsymatrix&Q, CPPL::dcovector&D, CPPL::dgematrix&U,
      CPPL::dsymatrix&R, bool&PSDchecked );

  //! @brief Wrapper to LAPACK function <TT>_dgesvd</TT> doing singular value decomposition of a general matrix
  static void _svd
    ( const CPPL::dgematrix&Q, CPPL::dcovector&S, CPPL::dgematrix&U,
      CPPL::dgematrix&VT );

  //! @brief Wrapper to LAPACK function <TT>_dgesvd</TT> doing singular value decomposition of a general matrix
  static void _svd
    ( const CPPL::dsymatrix&Q, CPPL::dcovector&S, CPPL::dgematrix&U,
      CPPL::dgematrix&VT );

  //! @brief Wrapper to LAPACK functions <TT>_dsytrf</TT> and <TT>_dsytri</TT> doing factorization and inversion of a symmetric indefinite matrix
  static void _inv
    ( const CPPL::dsymatrix&Q, CPPL::dsymatrix&Qinv );

  //! @brief Prototype function for finding junction points in convex/concave envelopes of univariate terms
  typedef double (puniv)
    ( const double x, const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
      const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const unsigned int n );
  //! @brief Secant method for root finding 
  static double _secant
    ( const double x0, const double x1, const double xL, const double xU,
      puniv f, const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
      const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const unsigned int n );
  //! @brief Golden section search method for root finding 
  static double _goldsect
    ( const double xL, const double xU, puniv f, const CPPL::dcovector&c1,
      const CPPL::dsymatrix&W1, const CPPL::dcovector&c2, const CPPL::dsymatrix&W2,
      const unsigned int n );
  //! @brief Golden section search iterations 
  static double _goldsect_iter
    ( const bool init, const double a, const double fa, const double b,
      const double fb, const double c, const double fc, puniv f,
      const CPPL::dcovector&c1, const CPPL::dsymatrix&W1, const CPPL::dcovector&c2,
      const CPPL::dsymatrix&W2, const unsigned int n );
  //! @brief Function whose root in the interval (0, 1) determines the minimal volume ellipsoid overapproximating the intersection of two ellipsoids
  static double _ell_fusionlambda
    ( const double a, const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
      const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const unsigned int n );
  //! @brief Function to find parameter value for minimal volume ellipsoid in ellispoid intersection
  static double _ell_get_lambda
    ( const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
      const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const bool isHP );

  //! @brief Function to calculate the trace of matrix <a>W</a>
  template <typename M> static double _trace
    ( const M&W )
    {
      double trace = 0.;
      for( unsigned int i=0; i<W.n; i++ ) trace += W(i,i);
      return trace;
    }

  //! @brief Function to calculate the determinant of symmetric matrix given its eignevalues <a>eigW</a>
  static double _det
    ( const CPPL::dcovector&eigW )
    {
      double det = 1.;
      for( unsigned int i=0; i<eigW.l; i++ ) det *= eigW(i);
      return det;
    }

  //! @brief Pause the program execution and prompt the user
  static void _pause()
    { double tmp;
      std::cout << "ENTER <1> TO CONTINUE" << std::endl;
      std::cin  >> tmp; }
};

////////////////////////////////////////////////////////////////////////

Ellipsoid::Options Ellipsoid::options;

inline void
Ellipsoid::_eigen
( const CPPL::dsymatrix&Q, CPPL::dcovector&D )
{
  if( !Q.n ){
    D.clear();
    return;
  }
  CPPL::dsymatrix S = Q;
  std::vector<double> vD;
  if( S.dsyev( vD ) ) throw Exceptions( Exceptions::LAPACK );
  D.resize( Q.n );
  typename std::vector<double>::const_iterator Di = vD.begin();
  for( unsigned int i=0; Di!=vD.end(); ++Di, i++ ) D(i) = *Di;
  return;
}

inline void
Ellipsoid::_eigen
( const CPPL::dsymatrix&Q, CPPL::dcovector&D, CPPL::dgematrix&U )
{
  if( !Q.n ){
    D.clear(); U.clear();
    return;
  }
  CPPL::dsymatrix S = Q;
  std::vector<double> vD;
  std::vector<CPPL::dcovector> vU;
  if( S.dsyev( vD, vU ) ) throw Exceptions( Exceptions::LAPACK );
  D.resize( Q.n );
  typename std::vector<double>::const_iterator Di = vD.begin();
  for( unsigned int i=0; Di!=vD.end(); ++Di, i++ ) D(i) = *Di;
  U.resize( Q.n, Q.m );
  typename std::vector<CPPL::dcovector>::const_iterator Uj = vU.begin();
  for( unsigned int j=0; Uj!=vU.end(); ++Uj, j++ )
    for( unsigned int i=0; i<(*Uj).l; i++ ) U(i,j) = (*Uj)(i);
  return;
}

inline void
Ellipsoid::_isPSD
( const CPPL::dsymatrix&Q, CPPL::dcovector&D, bool&PSD )
{
  if( !Q.n ){
    PSD = false;
    D.clear();
    return;
  }
  D.resize( Q.n );
  _eigen( Q, D );

  PSD = true;
  if( D(0) < -std::fabs( options.PSDTOL ) ){
#ifdef MC__ELLIPSOID_DEBUG
    std::cout << "Min. eigenvalue: " << D(0) << " < 0 !" << std::endl;
    _pause();
#endif
    PSD = false;
  }
  for( unsigned int i=0; i<D.l; i++ ){
    if( D(i) < 0. ) D(i) = 0.;
    else break; // b/c eigenvalues are returned in increasing order
  }
  return;
}

inline void
Ellipsoid::_sqrt
( const CPPL::dsymatrix&Q, CPPL::dcovector&D, CPPL::dgematrix&U,
  CPPL::dsymatrix&sqrtQ, bool&PSDchecked )
{
  bool PSD; _isPSD( Q, D, PSD );
  if( options.PSDCHK && !PSDchecked && !PSD ){
    PSDchecked = true;
    throw Exceptions( Exceptions::NONPSD );
  }
  if( !D.l || !U.n || !U.m ){ sqrtQ.clear(); return; }

  CPPL::dgbmatrix sqrtD(D.l,D.l,0,0);
  for( unsigned int i=0; i<sqrtD.n; i++ )
    sqrtD(i,i) = std::sqrt( D(i) );
  CPPL::dgematrix sqrtQ_ge = U*sqrtD*t(U);
  sqrtQ.resize( sqrtQ_ge.n );
  for( unsigned int i=0; i<sqrtQ_ge.n; i++ )
    for( unsigned int j=0; j<=i; j++ )
      sqrtQ(i,j) = 0.5*(sqrtQ_ge(i,j)+sqrtQ_ge(j,i));

#ifdef MC__ELLIPSOID_DEBUG_SQRT
  std::cout << "sqrtQ*sqrtQ:\n" << sqrtQ*sqrtQ;
#endif
  return;
}

inline void
Ellipsoid::_svd
( const CPPL::dgematrix&Q, CPPL::dcovector&S, CPPL::dgematrix&U,
  CPPL::dgematrix&VT )
{
  if( !Q.n || !Q.m ){
    S.clear(); U.clear(); VT.clear();
    return;
  }
  CPPL::dgematrix R( Q ); // local copy; otherwise Q entries are overwritten...
  if( R.dgesvd( S, U, VT ) ) throw Exceptions( Exceptions::LAPACK );
  return;
}

inline void
Ellipsoid::_svd
( const CPPL::dsymatrix&Q, CPPL::dcovector&S, CPPL::dgematrix&U,
  CPPL::dgematrix&VT )
{
  return _svd( Q.to_dgematrix(), S, U, VT );
}

inline void
Ellipsoid::_inv
( const CPPL::dsymatrix&Q, CPPL::dsymatrix&Qinv )
{
  if( !Q.n ){
    Qinv.clear();
    return;
  }
  //CPPL::dsymatrix R = i( Q );
  CPPL::dsymatrix R;
  if( CPPL::dsysv( Q, R ) ) throw Exceptions( Exceptions::LAPACK );
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "\n***Ellispoid: Q\n" << Q;
  std::cout << "***Ellispoid: inv(Q) (ill-conditioned)\n" << R;
#endif
  CPPL::dgematrix invRQ_R = i( R * Q ) * R;
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "***Ellispoid: inv(R*Q)*R\n" << invRQ_R;
#endif
  //Qinv = R;
  Qinv.resize( Q.n );
  for( unsigned int i=0; i<Q.n; i++ )
    for( unsigned int j=0; j<=i; j++ )
      Qinv(i,j) = 0.5 * ( invRQ_R(i,j) + invRQ_R(j,i) ); 
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "***Ellispoid: inv(Q) (proper)\n" << invRQ_R;
#endif

  return;
}

inline Ellipsoid mtimes
( const Ellipsoid&E, const CPPL::dgematrix&A, const CPPL::dcovector&b=CPPL::dcovector() )
{
#ifdef MC__ELLIPSOID_DEBUG
  assert( E._Q.n == A.n );
#endif  
  Ellipsoid AE;
  if( !A.m ) return AE;

  // Transformed center
  AE._c = A * E._c;
  if( b.l == A.m ) AE._c += b;

  // Transformed shape
  CPPL::dgematrix AQAT = A * E._Q * t(A);
  AE._Q.resize( A.m );
  for( unsigned int i=0; i<A.m; i++ )
    for( unsigned int j=0; j<=i; j++ )
      AE._Q(i,j) = 0.5*(AQAT(i,j)+AQAT(j,i));
#ifdef MC__ELLIPSOID_DEBUG
  bool PSD; AE._isPSD( AE._Q, AE._eigQ.first, PSD ); assert( PSD );
#endif
  return AE;
}

inline Ellipsoid minksum_ea
( const std::vector<Ellipsoid>&E, const CPPL::dcovector&D )
{
  Ellipsoid EA;
  if( !E.size() ) return EA;

  // Check size
  const unsigned int n = (*E.begin())._Q.n;
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
#ifdef MC__ELLIPSOID_DEBUG
  for( ; Ei!=E.end(); ++Ei ) assert( (*Ei)._Q.n == n );
  assert( D.l == n );
#endif

  // Construct center and shape matrix
  Ei = E.begin(); CPPL::dcovector c_EA( (*Ei)._c );
  for( ++Ei; Ei!=E.end(); ++Ei ) c_EA += (*Ei)._c;
  CPPL::dsymatrix Q_EA( n );
  double sum_sqrt_lTQil = 0.;
  for( Ei=E.begin(); Ei!=E.end(); ++Ei ){
    double sqrt_lTQil = std::sqrt(D%((*Ei)._Q*D));
    sum_sqrt_lTQil += sqrt_lTQil;
    if( Ei==E.begin() ) Q_EA  = (*Ei)._Q/sqrt_lTQil;
    else                Q_EA += (*Ei)._Q/sqrt_lTQil;
  }
  Q_EA *= sum_sqrt_lTQil;
  EA.set(Q_EA,c_EA);

  return EA;  
}

inline std::vector<Ellipsoid> minksum_ea
( const std::vector<Ellipsoid>&E, const std::vector<CPPL::dcovector>&D )
{
  std::vector<Ellipsoid> EA;
  if( !E.size() || !D.size() ) return EA;

  // Check size
  const unsigned int n = (*E.begin())._Q.n;
  typename std::vector<CPPL::dcovector>::const_iterator Di = D.begin();
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
#ifdef MC__ELLIPSOID_DEBUG
  for( ; Ei!=E.end(); ++Ei ) assert( (*Ei)._Q.n == n );
  for( ; Di!=D.end(); ++Di ) assert( (*Di).l == n );
#endif

  // Construct centers and shape matrices
  Ei = E.begin(); CPPL::dcovector c_EA( (*Ei)._c );
  for( ++Ei; Ei!=E.end(); ++Ei ) c_EA += (*Ei)._c;
  CPPL::dsymatrix Q_EA( n );
  for( Di=D.begin(); Di!=D.end(); ++Di ){
    double sum_sqrt_lTQil = 0.;
    for( Ei=E.begin(); Ei!=E.end(); ++Ei ){
      double sqrt_lTQil = std::sqrt((*Di)%((*Ei)._Q*(*Di)));
      sum_sqrt_lTQil += sqrt_lTQil;
      if( Ei==E.begin() ) Q_EA  = (*Ei)._Q/sqrt_lTQil;
      else                Q_EA += (*Ei)._Q/sqrt_lTQil;
    }
    Q_EA *= sum_sqrt_lTQil;
    EA.push_back( Ellipsoid(Q_EA,c_EA) );
  }

  return EA;  
}

inline Ellipsoid minksum_ea
( const Ellipsoid&E, const std::pair<CPPL::dcovector,CPPL::dcovector>&I,
  const double&TOL, const double&EPS=machprec() )
{
#ifdef MC__ELLIPSOID_DEBUG
  assert( E._Q.n == I.first.l && ( !I.second.l || E._Q.n == I.second.l ) );
#endif
  Ellipsoid EA;
  if( !E._Q.n ) return EA;

  // Construct center and shape matrix
  CPPL::dcovector c_EA( E._c );
  if( I.second.l ) c_EA += I.second;

  double trQ = 0.0;
  CPPL::dcovector sqrR(E._Q.n);
  for( int i=0; i<E._Q.n; i++ ){
    trQ += E._Q(i,i)/(E._Q(i,i)+TOL);
    sqrR(i) = I.first(i)/std::sqrt(E._Q(i,i)+TOL);
  }
  trQ = std::sqrt(trQ);
  double kappa = trQ;
  for( int i=0; i<E._Q.n; i++ ) kappa += sqrR(i);

  CPPL::dsymatrix Q_EA( E._Q*kappa/(trQ+EPS) );
  for( int i=0; i<E._Q.n; i++ )
    Q_EA(i,i) += I.first(i)*I.first(i)*kappa/(sqrR(i)+EPS)+EPS;

  return EA.set(Q_EA,c_EA);  
}

inline Ellipsoid inv
( const Ellipsoid&E )
{
  Ellipsoid Einv( E );
  if( !E._Q.n ) return Einv;
  Einv.regQ();
  Einv._Q = Einv.invQ();
  Einv._reset_auxiliary();
  return Einv;
}  

inline std::vector<Ellipsoid> inv
( const std::vector<Ellipsoid>&E )
{
  std::vector<Ellipsoid> Einv;
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
  for( ; Ei!=E.end(); ++Ei ) Einv.push_back( inv(*Ei) );
  return Einv;
}  

inline double dist
( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP )
{
  if( E._Q.n != HP.first.l ){
    std::cerr << "mc::dist: ellipsoid and hyperplane must be of the same dimension.\n";
    return 0./0.; // return NaN
  }

  return ( std::fabs( HP.second - E._c % HP.first )
         - std::sqrt( HP.first % ( E._Q * HP.first ) ) )
         / ( std::sqrt( HP.first % HP.first ) );
}

inline Ellipsoid intersection_ea
( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP,
  const double&TOL=machprec() )
{
  if( E._Q.n != HP.first.l ){
    std::cerr << "mc::intersection_ea: ellipsoid and hyperplane must be of the same dimension.\n";
    return Ellipsoid();
  }
  const unsigned int n = E._Q.n;
  if( !n ) return Ellipsoid();

  CPPL::dcovector v = -HP.first / std::sqrt( HP.first % HP.first );
  double c = -HP.second / std::sqrt( HP.first % HP.first );
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "dist(E,HP): " << dist( E, HP ) << std::endl;
  std::cout << "HP.vT E.c > 0 ? " << ( v % E._c > c? 'T': 'F' ) << std::endl;
#endif
  if( dist( E, HP ) > 0 && v % E._c > c )
    return E;
  if( dist( E, HP ) > 0 && v % E._c < c )
    return Ellipsoid();
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "E:\n" << E << std::endl;
  std::cout << "HP:\n" << HP.first << "; " << HP.second << std::endl;
#endif

  Ellipsoid EM( E );
  CPPL::dsymatrix W2(n);
  CPPL::dgematrix vmat(n,1);
  for( unsigned int i=0; i<n; i++ ) vmat(i,0) = v(i);
  CPPL::dgematrix vmatvmatT = vmat * t(vmat);
  for( unsigned int i=0; i<n; i++ )
    for( unsigned int j=0; j<=i; j++ )
      W2(i,j) = 0.5 * ( vmatvmatT(i,j) + vmatvmatT(j,i) )
              / ( 4.*EM.eigQ().first(n-1) ); 
  Ellipsoid Einter = hpintersection( EM, HP );
  double h = 2.*std::sqrt( EM.eigQ().first(n-1) );
  CPPL::dcovector c2 = Einter._c + h*v;
  Ellipsoid EM2( W2, c2 ); EM2.regQ(); W2 = EM2._Q;
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "W2:\n" << W2 << std::endl;
  std::cout << "c2:\n" << c2 << std::endl;
#endif

  const CPPL::dcovector& c1 = EM._c;
  EM.regQ(); const CPPL::dsymatrix& W1 = EM.invQ();
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "E (reg):\n" << EM << std::endl;
  std::cout << "rank(E._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "inv(E._Q):\n" << EM.invQ() << std::endl;
#endif

  double a = Ellipsoid::_ell_get_lambda( c1, W1, c2, W2, true );
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "a (ell_get_lambda):" << a << std::endl;
#endif
  CPPL::dsymatrix W = a*W1 + (1-a)*W2;
  CPPL::dsymatrix Winv( n ); Ellipsoid::_inv( W, Winv );
  double k = 1. - a*(1-a)*(c2-c1)%(W2*Winv*W1*(c2-c1));
  CPPL::dcovector q = Winv*(a*W1*c1 + (1-a)*W2*c2);
  CPPL::dsymatrix Q = (1+TOL)*k*Winv;
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "W:\n" << W << std::endl;
  std::cout << "Winv:\n" << Winv << std::endl;
  std::cout << "k:" << k << std::endl;
  std::cout << "q:" << q << std::endl;
  std::cout << "Q:\n" << Q << std::endl;
#endif
  return Ellipsoid( Q, q );
}

inline double
Ellipsoid::_ell_get_lambda
( const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
  const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const bool isHP )
{
  assert( c1.l == W1.n && c2.l == W2.n && c1.l == c2.l );
  const unsigned int n = c1.l;

  double a;
  if( !options.ROOTSECANT ){
    a= _goldsect( -.1, 1., _ell_fusionlambda, c1, W1, c2, W2, n );
    return( a<=0.? ( isHP? 1.: 0. ): a );
  }

  try{
    a = _secant( .9, 1., -.1, 1., _ell_fusionlambda, c1, W1, c2, W2, n );
  }
  catch( Ellipsoid::Exceptions ){
    a = _goldsect( -.1, 1., _ell_fusionlambda, c1, W1, c2, W2, n );
  }
    return( a<=0.? ( isHP? 1.: 0. ): a );
}

inline double
Ellipsoid::_secant
( const double x0, const double x1, const double xL, const double xU,
  puniv f, const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
  const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const unsigned int n )
{
  double xkm = std::max(xL,std::min(xU,x0));
  double fkm = f(xkm,c1,W1,c2,W2,n);
  double xk = std::max(xL,std::min(xU,x1));
  
  for( unsigned int it=0; !options.ROOTMAXIT || it<options.ROOTMAXIT; it++ ){
    double fk = f(xk,c1,W1,c2,W2,n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
    std::cout << "xk (fk):" << xk << "  (" << fk << ")" << std::endl;
#endif
    if( std::fabs(fk) < options.ROOTTOL ) return xk;
    double Bk = (fk-fkm)/(xk-xkm);
    if( Bk == 0 ) throw Exceptions( Exceptions::ROOT );
    if( isequal(xk,xL) && fk/Bk>0 ) return xk;
    if( isequal(xk,xU) && fk/Bk<0 ) return xk;
    xkm = xk;
    fkm = fk;
    xk = std::max(xL,std::min(xU,xk-fk/Bk));
  }

  throw Exceptions( Exceptions::ROOT );
}

inline double
Ellipsoid::_goldsect
( const double xL, const double xU, puniv f, const CPPL::dcovector&c1,
  const CPPL::dsymatrix&W1, const CPPL::dcovector&c2, const CPPL::dsymatrix&W2,
  const unsigned int n )
{
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  const double fL = f(xL,c1,W1,c2,W2,n), fU = f(xU,c1,W1,c2,W2,n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "x (f):" << xL << "  (" << fL << ")" << std::endl;
  std::cout << "x (f):" << xU << "  (" << fU << ")" << std::endl;
#endif
  if( fL*fU > 0 ) throw Exceptions( Exceptions::ROOT );
  const double xm = xU-phi*(xU-xL), fm = f(xm,c1,W1,c2,W2,n);
  return _goldsect_iter( true, xL, fL, xm, fm, xU, fU, f, c1, W1, c2, W2, n );
}

inline double
Ellipsoid::_goldsect_iter
( const bool init, const double a, const double fa, const double b,
  const double fb, const double c, const double fc, puniv f,
  const CPPL::dcovector&c1, const CPPL::dsymatrix&W1, const CPPL::dcovector&c2,
  const CPPL::dsymatrix&W2, const unsigned int n )
// a and c are the current bounds; the minimum is between them.
// b is a center point
{
  static unsigned int iter;
  iter = ( init? 1: iter+1 );
  const double phi = 2.-(1.+std::sqrt(5.))/2.;
  bool b_then_x = ( c-b > b-a );
  double x = ( b_then_x? b+phi*(c-b): b-phi*(b-a) );
  //if( std::fabs(c-a) < options.ROOTTOL*(std::fabs(b)+std::fabs(x)) 
  // || iter > options.ROOTMAXIT ) return (c+a)/2.;
  double fx = f(x,c1,W1,c2,W2,n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "x (f):" << x << "  (" << fx << ")" << std::endl;
#endif
  if( std::fabs(fx) < options.ROOTTOL || ( options.ROOTMAXIT && iter > options.ROOTMAXIT ) )
    return x;
  if( b_then_x )
    return( fa*fx<0? _goldsect_iter( false, a, fa, b, fb, x, fx, f, c1, W1, c2, W2, n ):
                     _goldsect_iter( false, b, fb, x, fx, c, fc, f, c1, W1, c2, W2, n ) );
  return( fa*fb<0? _goldsect_iter( false, a, fa, x, fx, b, fb, f, c1, W1, c2, W2, n ):
                   _goldsect_iter( false, x, fx, b, fb, c, fc, f, c1, W1, c2, W2, n ) );
}

inline double
Ellipsoid::_ell_fusionlambda
( const double a, const CPPL::dcovector&c1, const CPPL::dsymatrix&W1,
  const CPPL::dcovector&c2, const CPPL::dsymatrix&W2, const unsigned int n )
{
  CPPL::dsymatrix W = a*W1 + (1-a)*W2;
  CPPL::dsymatrix Winv( n ); Ellipsoid::_inv( W, Winv );
  double k = 1. - a*(1-a)*(c2-c1)%(W2*Winv*W1*(c2-c1));
  CPPL::dcovector c = Winv*(a*W1*c1 + (1-a)*W2*c2);

  CPPL::dcovector eigW; _eigen( W, eigW );
  double detW = Ellipsoid::_det(eigW);
  double traW = Ellipsoid::_trace(detW*Winv*(W1-W2));
  return k*detW*traW - n*(detW*detW) * ( c%(2*W1*c1-2*W2*c2+(W2-W1)*c)
         - c1%(W1*c1) + c2%(W2*c2) );
}

inline Ellipsoid hpintersection
( const Ellipsoid&E, const std::pair<CPPL::dcovector,double>&HP )
{
  if( E._Q.n != HP.first.l ){
    std::cerr << "mc::intersection_ea: ellipsoid and hyperplane must be of the same dimension.\n";
    return Ellipsoid();
  }
  const unsigned int n = E._Q.n;
  if( !n ) return Ellipsoid();

#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "dist(E,HP): " << dist( E, HP ) << std::endl;
#endif
  if( dist( E, HP ) > 0 ) return Ellipsoid();
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "E:\n" << E << std::endl;
  std::cout << "HP:\n" << HP.first << "; " << HP.second << std::endl;
#endif

  CPPL::dcovector e1( n ); e1.zero(); e1(0) = 1.;
  CPPL::dgematrix T = E.align( e1, HP.first );
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "T:\n" << T << std::endl;
  std::cout << "T*d:\n" << T*HP.first << std::endl;  
#endif
  CPPL::dcovector f = ( HP.second / ( HP.first % HP.first ) ) * T * HP.first;
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "f: " << f << std::endl;
#endif
  Ellipsoid EM = mtimes( E, T, -f );
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "EM:\n" << EM << std::endl;
  std::cout << "rank(EM._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "svd(EM._Q):\n" << EM.svdQ().first << std::endl;
#endif
  EM.regQ();
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "EM (reg):\n" << EM << std::endl;
  std::cout << "rank(EM._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "svd(EM._Q):\n" << EM.svdQ().first << std::endl;
  std::cout << "inv(EM._Q):\n" << EM.invQ() << std::endl;
#endif

  CPPL::dsymatrix W( n-1 );
  CPPL::dcovector w( n-1 );
  for( unsigned int i=1; i<n; i++ ){
    for( unsigned int j=i; j<n; j++ )
      W(i-1,j-1) = EM.invQ()(i,j);
    w(i-1) = EM.invQ()(i,0);
  }
  const double w11 = EM.invQ()(0,0);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "W:\n" << W << std::endl;
  std::cout << "w:\n" << w << std::endl;
  std::cout << "w11:\n" << w11 << std::endl;
#endif
  CPPL::dsymatrix Winv( n-1 );
  Ellipsoid::_inv( W, Winv );
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "Winv:\n" << Winv << std::endl;
#endif

  CPPL::dcovector Winv_w( Winv * w );
  const double h = EM.c(0)*EM.c(0) * ( w11 - w % Winv_w );
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "h: " << h << std::endl;
#endif

  CPPL::dsymatrix Z( n ); Z.zero();
  CPPL::dcovector z( EM._c ); z(0) = 0.; 
  for( unsigned int i=1; i<n; i++ ){
    for( unsigned int j=i; j<n; j++ )
      Z(i,j) = (1.-h) * Winv(i-1,j-1);
    z(i) += EM._c(0) * Winv_w(i-1);
  }
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "Z:\n" << Z << std::endl;
  std::cout << "z:\n" << z << std::endl;
#endif
  Ellipsoid I( Z, z );
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "I: " << I;
  std::cout << "eig(I): " << I.eigQ().first;
#endif
  I += f;
  return mtimes( I, t(T) );
}

inline Ellipsoid hpintersection
( const Ellipsoid&E, const std::vector< std::pair<CPPL::dcovector,double> >&HP )
{
  Ellipsoid Ei( E );
  typename std::vector< std::pair<CPPL::dcovector,double> >::const_iterator HPi = HP.begin();
  for( ; HPi!=HP.end(); ++HPi )
    Ei = hpintersection( Ei, *HPi );
  return Ei; 
}

inline Ellipsoid ell_unitball
( const unsigned int n )
{
  Ellipsoid Eunit;
  Eunit._c.resize(n); Eunit._c.zero();
  Eunit._Q.resize(n); Eunit._Q.identity();
  return Eunit;
}

inline std::ostream&
operator<<
( std::ostream&out, const Ellipsoid&E)
{
  if( !E.n() ) return out;
  const int iprec = 5;
  out << std::scientific << std::setprecision(iprec);
  out << "\ncenter:\n" << E._c;
  out << "shape:\n" << E._Q;
  return out;
}

} // namespace mc

#endif
