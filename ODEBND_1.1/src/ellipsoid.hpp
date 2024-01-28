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

Suppose we want to define an ellipsoid \f$\mathcal E(c,Q)\f$ with
\f{align*}
  c = & \left(\begin{array}{c} 1.\\-1.\end{array}\right),\ \text{and} & 
  Q = & \left(\begin{array}{cc} 0.5 & -0.4 \\ -0.4 & 0.5 \end{array}\right).
\f}

First, we specify whether we want to define the shape matrix in upper- or lower-triangular form:

\code
    mc::Ellipsoid::options.STORAGE = mc::Ellipsoid::Options::UPPER;
\endcode

Then, we define the ellipsoid as follows:

\code
    const unsigned int n = 2;
    double Q[n*(n-1)] = { 0.5, -0.4, 0.5 };
    double c[n] = { 1., -1. };
    mc::Ellipsoid E( n, Q, c ); 
    std::cout << E << std::endl;
\endcode

An ellipsoid is simply displayed as:

\code
    std::cout << E << std::endl;
\endcode

In the present case, the following information is displayed:

\verbatim
 n: 2     c:  1.00000e+00   Q:  5.00000e-01 -4.00000e-01
             -1.00000e+00                    5.00000e-01
\endverbatim


\section sec_ELL_opt How Are the Options Set for Ellipsoidal Calculus?

The options are defined in the structure mc::Ellipsoid::Options. Current options are as follows:

<TABLE border="1">
<CAPTION><EM>Options in mc::Ellipsoid: name, type and
description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>STORAGE</tt> <TD><tt>mc::Ellipsoid::Options::SHAPE_MATRIX</tt> <TD>mc::Ellipsoid::Options::LOWER
         <TD>Storage format for shape matrix
     <TR><TH><tt>PSDCHK</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether or not to check positive semi-definiteness of shape matrices
     <TR><TH><tt>PSDTOL</tt> <TD><tt>double</tt> <TD>1e2*MACHPREC
         <TD>Tolerance for positive semi-definiteness check of shape matrices
</TABLE>

The mc::Ellipsoid class has a public static member called mc::Ellipsoid::options that can be used to set/modify the options; for example:

\code
      mc::Ellipsoid::options.PSDCHK = true;
      mc::Ellipsoid::options.STORAGE = mc::Ellipsoid::Options::UPPER;
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

#include "mcfunc.hpp"

#undef  MC__ELLIPSOID_DEBUG_INV
#undef  MC__ELLIPSOID_DEBUG_SVD
#undef  MC__ELLIPSOID_DEBUG_EIGEN
#undef  MC__ELLIPSOID_DEBUG_SQRT

extern "C" void dsyev_
( const char*jobz, const char*uplo, const unsigned int*n, double*a,
  const unsigned int*lda, double*w, double*work, const int*lwork, int*info );

extern "C" void dgesvd_
( const char*jobu, const char*jobvt, const unsigned int*m, const unsigned int*n,
  double*a, const unsigned int*lda, double*s, double*u, const unsigned int*lu,
  double*vt, const unsigned int*lvt, double*work, const int*lwork, int*info );

extern "C" void dsytrf_
( const char*uplo, const unsigned int*n, double*a, const unsigned int*lda,
  int*ipiv, double*work, const int*lwork, int*info );

extern "C" void dsytri_
( const char*uplo, const unsigned int*n, double*a, const unsigned int*lda,
  const int*ipiv, double*work, int*info );

extern "C" double dnrm2_
( const unsigned int*n, const double*a, const unsigned int*incx );

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
  friend Ellipsoid inv
    ( const Ellipsoid&E );
  friend Ellipsoid mtimes
    ( const Ellipsoid&E, const unsigned int m, const double*A,
      const double*b );
  friend Ellipsoid* minksum_ea
    ( const unsigned int nE, const Ellipsoid*E, const unsigned int nD,
      const double*D );
  friend Ellipsoid intersection_ea
    ( const Ellipsoid&E1, const Ellipsoid&E2 );
  friend Ellipsoid ell_unitball
    ( const unsigned int n );

private:

  //! @brief Ellipsoid dimension
  unsigned int _n;
  //! @brief Ellipsoid storage
  double *_p;
  //! @brief Pointer to ellipsoid shape matrix
  double *_Q;
  //! @brief Pointer to ellipsoid center
  double *_c;
  //! @brief Pointer to eigenvalues
  std::pair< double*, double* > _eigQ;
  //! @brief Pointer to singular values
  std::pair< double*, std::pair<double*,double*> > _svdQ;
  //! @brief Pointer to square root of shape matrix
  double *_sqrtQ;
  //! @brief Pointer to inverse of shape matrix
  double *_invQ;
  //! @brief Whether shape matrix was checked for postive semi-definiteness
  bool _PSDchecked;

  //! @brief Return the ith diagonal element of shape matrix
  double _diagQ
    ( const unsigned int i ) const;
  //! @brief Return/set the ith diagonal element of shape matrix
  double& _diagQ
    ( const unsigned int i );

  //! @brief Delete pointers in class
  void _reset_auxiliary()
    {
      delete[] _eigQ.first;
      delete[] _eigQ.second;
      delete[] _sqrtQ;
      delete[] _svdQ.first;
      delete[] _svdQ.second.first;
      delete[] _svdQ.second.second;
      delete[] _invQ;
      _eigQ.first = _eigQ.second = _sqrtQ = _svdQ.first = _svdQ.second.first =
      _svdQ.second.second = _invQ = 0;
      _PSDchecked = false;
    }

public:

  /** @defgroup ELLIPSOID Ellipsoidal Calculus
   *  @{
   */
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in \mathbb R^n\f$, returns \f$(\mathcal E + d)\f$
  Ellipsoid& operator+=
    ( const double *d )
    {
      for( unsigned int i=0; i<_n; i++ )
        _c[i] += d[i];
      return *this;
    }
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in \mathbb R^n\f$, returns \f$(\mathcal E - d)\f$
  Ellipsoid& operator-=
    ( const double *d )
    {
      for( unsigned int i=0; i<_n; i++ )
        _c[i] -= d[i];
      return *this;
    }

  //! @brief Ellipsoid options
  static struct Options
  {
    //! @brief Constructor
    Options():
      STORAGE(LOWER), PSDCHK(false), PSDTOL(machprec()*1e2), RANKTOL(machprec()*1e2)
      {}
    //! @brief Enumeration type for shape matrix storage
    enum SHAPE_MATRIX{
      LOWER=0,		//!< Lower triangular part of shape matrix stored row-wise
      UPPER		//!< Upper triangular part of shape matrix stored row-wise
    };
    //! @brief Storage format for shape matrix (default=LOWER)
    SHAPE_MATRIX STORAGE;
    //! @brief Whether or not to check positive semi-definiteness of shape matrices (default=false)
    bool PSDCHK;
    //! @brief Tolerance for positive semi-definiteness check of shape matrices (default=1e2*mc::machprec())
    double PSDTOL;
    //! @brief Absolute tolerance for rank of shape matrices (default=1e2*mc::machprec())
    double RANKTOL;
  } options;

  //! @brief Ellipsoid exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for Interval exception handling
    enum TYPE{
      NONPSD=1,	//!< Non positive-semi definite shape matrix
      SIZE,	//!< Operation on ellipsoids of different dimension not permitted
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Return error flag
    int ierr(){ return _ierr; }
    //! @brief Return error description
    std::string what(){
      switch( _ierr ){
      case NONPSD: default:
        return "Nonpositive semi-definite shape matrix";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Default constructor (needed to declare arrays of Ellipsoid class)
  Ellipsoid
    ():
    _n(0), _p(0), _Q(0), _c(0), _eigQ(NULL,NULL), _svdQ(NULL,_eigQ), _sqrtQ(0),
    _invQ(0), _PSDchecked(false)
    {}

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$
  Ellipsoid
    ( const unsigned int n, const double *Q=0, const double *c=0 ):
    _n(n), _eigQ(NULL,NULL), _svdQ(NULL,_eigQ), _sqrtQ(0), _invQ(0),
    _PSDchecked(false)
    {
      if( options.PSDCHK ){
        if( !_isPSD( n, Q, _eigQ, options.STORAGE ) ){
	  std::cout << "Min. eigenvalue: " << *_eigQ.first << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
      
      _p = new double[_n+_n*(_n+1)/2];
      _c = _p;
      _Q = _p+_n;
      for( unsigned int i=0, k=0; i<n; i++ ){
        _c[i] = ( c? c[i]: 0.);
	for( unsigned int j=i; j<n; j++, k++ )
	  _Q[k] = ( Q? Q[k]: 0. );
      }
    }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ enclosing interval vector of radius \f$r\f$ centered at \f$c\f$
  Ellipsoid
    ( const unsigned int n, const double *r=0, const double *c=0 ):
    _n(n), _eigQ(NULL,NULL), _svdQ(NULL,_eigQ), _sqrtQ(0), _invQ(0),
    _PSDchecked(false)
    {
      if( options.PSDCHK ){
        for( unsigned int i=0; i<n; i++ ) if( r[i] < 0 ){
          std::cout << "Interval radius: " << r[i] << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
        }
        _PSDchecked = true;
      }
      
      _p = new double[_n+_n*(_n+1)/2];
      _c = _p;
      _Q = _p+_n;
      double tr2 = 0.;
      for( unsigned int i=0, k=0; i<n; i++ ){
        tr2 += ( r? r[i]*r[i]: 0. );
	for( unsigned int j=i; j<n; j++, k++ )
	  _Q[k] = 0.;
      }
      for( unsigned int i=0, k=0; i<n; i++ ){
        _c[i] = ( c? c[i]: 0. );
        _diagQ[i] = ( r? r[i]*std::sqrt(tr2): 0. );
      }
    }

  //! @brief Copy constructor
  Ellipsoid
    ( const Ellipsoid&E ):
    _n(E._n), _eigQ(NULL,NULL), _svdQ(NULL,_eigQ), _sqrtQ(0), _invQ(0),
    _PSDchecked(false)
    {
      if( !_n ) return;
      _p = new double[_n+_n*(_n+1)/2];
      _c = _p;
      _Q = _p+_n;
      for( unsigned int i=0, k=0; i<_n; i++ ){
        _c[i] = E._c[i];
	for( unsigned int j=i; j<_n; j++, k++ ) _Q[k] = E._Q[k];
      }
    }

  //! @brief Destructor
  ~Ellipsoid()
    {
      delete[] _eigQ.first;
      delete[] _eigQ.second;
      delete[] _sqrtQ;
      delete[] _svdQ.first;
      delete[] _svdQ.second.first;
      delete[] _svdQ.second.second;
      delete[] _invQ;
      delete[] _p;
    }

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$
  Ellipsoid& set
    ( const unsigned int n, const double *Q=0, const double *c=0 )
    {
      _reset_auxiliary();
      if( options.PSDCHK ){
        if( !_PSDchecked && !_isPSD( n, Q, _eigQ, options.STORAGE ) ){
	  std::cout << "Min. eigenvalue: " << *_eigQ.first << std::endl;
          throw Exceptions( Exceptions::NONPSD );
	}
        _PSDchecked = true;
      }
      
      if( _n != n ){
        _n = n;
        delete[] _p;
        _p = new double[_n+_n*(_n+1)/2];
        _c = _p;
        _Q = _p+_n;
      }
      for( unsigned int i=0, k=0; i<n; i++ ){
        _c[i] = ( c? c[i]: 0.);
	for( unsigned int j=i; j<n; j++, k++ )
	  _Q[k] = ( Q? Q[k]: 0. );
      }
      return *this;
    }

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape matrix \f$Q\f$
  Ellipsoid& set
    ( const unsigned int n, const double *r=0, const double *c=0 )
    {
      _reset_auxiliary();
      if( options.PSDCHK ){
        for( unsigned int i=0; i<n; i++ ) if( r[i] < 0 ){
          std::cout << "Interval radius: " << r[i] << " < 0 !" << std::endl;
          throw Exceptions( Exceptions::NONPSD );
        }
        _PSDchecked = true;
      }
      
      if( _n != n ){
        delete[] _p;
        _n = n;
        _p = new double[_n+_n*(_n+1)/2];
        _c = _p;
        _Q = _p+_n;
      }
      double tr2 = 0.;
      for( unsigned int i=0, k=0; i<n; i++ ){
        tr2 += ( r? r[i]*r[i]: 0. );
	for( unsigned int j=i; j<n; j++, k++ )
	  _Q[k] = 0.;
      }
      for( unsigned int i=0, k=0; i<n; i++ ){
        _c[i] = ( c? c[i]: 0. );
        _diagQ[i] = ( r? r[i]*std::sqrt(tr2): 0. );
      }
      return *this;
    }

  //! @brief Return dimension of ellipsoid
  unsigned int n() const
    {
      return _n;
    }
  //! @brief Return center of ellipsoid
  const double* c() const
    {
      return _c;
    }
  //! @brief Return shape matrix of ellipsoid
  const double* Q() const
    {
      return _Q;
    }
  //! @brief Return shape matrix coefficient
  double Q
    ( unsigned int i, unsigned int j ) const
    {
      assert( i<_n && j<_n );
      switch( options.STORAGE ){
      case Ellipsoid::Options::LOWER:
        return( j<=i? _Q[_n*j-j*(j-1)/2+i-j]: _Q[_n*i-i*(i-1)/2+j-i] );
      case Ellipsoid::Options::UPPER: default:
        return( j>=i? _Q[j*(j+1)/2+i]: _Q[i*(i+1)/2+j] );
      }
    }
  //! @brief Return/set shape matrix coefficient
  double& Q
    ( unsigned int i, unsigned int j )
    {
      assert( i<_n && j<_n );
      switch( options.STORAGE ){
      case Ellipsoid::Options::LOWER:
        return( j<=i? _Q[_n*j-j*(j-1)/2+i-j]: _Q[_n*i-i*(i-1)/2+j-i] );
      case Ellipsoid::Options::UPPER: default:
        return( j>=i? _Q[j*(j+1)/2+i]: _Q[i*(i+1)/2+j] );
      }
    }
  //! @brief Return square root of shape matrix
  const double* sqrtQ()
    {
      if( !_sqrtQ ) _sqrtQ = _sqrt( _n, _Q, _eigQ, options.STORAGE, _PSDchecked );
      return _sqrtQ;
    }
  //! @brief Return eigenvalues and eigenvectors of shape matrix
  std::pair< const double*, const double*> eigQ()
    {
      if( !_eigQ.first || !_eigQ.second ){
        delete[] _eigQ.first;
        delete[] _eigQ.second;
        _eigQ = _eigen( _n, _Q, options.STORAGE, true );
      }
      return _eigQ;
    }
  //! @brief Return rank of shape matrix
  unsigned int rankQ()
    {
      if( !_svdQ.first ){
        delete[] _svdQ.first;
        delete[] _svdQ.second.first;
        delete[] _svdQ.second.second;
        _svdQ = _svd( _n, _Q, options.STORAGE, false, false );
      }
      unsigned int rank = _n;
      for( unsigned int i=0; i<_n && _svdQ.first[_n-i-1]<options.RANKTOL; i++, rank-- ){}
      return rank;
    }
  //! @brief Return singular value decomposition of shape matrix
  std::pair< const double*, std::pair<const double*,const double*> > svdQ()
    {
      if( !_svdQ.first || !_svdQ.second.first || !_svdQ.second.second ){
        delete[] _svdQ.first;
        delete[] _svdQ.second.first;
        delete[] _svdQ.second.second;
        _svdQ = _svd( _n, _Q, options.STORAGE, true, true );
      }
      return _svdQ;
    }
  //! @brief Return pointer to regularized shape matrix
  const double* regQ()
    {
      const unsigned int r = rankQ();
      if( r == _n ) return _Q;

      _reset_auxiliary();
      const double*U = svdQ().second.first;
      const double eps = 2.*options.RANKTOL;
      for( unsigned int i=0; i<r; i++ ){
        _svdQ.first[i] += eps;
	for( unsigned int j=0; j<=i; j++ )
	  for( unsigned int k=0; k<_n-r; k++ )
            Q(i,j) += eps * ( U[i+(r+k)*_n] * U[j+(r+k)*_n] );
      }
      for( unsigned int i=0; i<_n-r; i++ ){
        _svdQ.first[r+i] += eps;
        for( unsigned int j=0; j<=i; j++ )
          for( unsigned int k=0; k<_n-r; k++ )
            Q(r+i,r+j) += eps * ( U[r+i+(r+k)*_n] * U[r+j+(r+k)*_n] );
      }
      for( unsigned int i=0; i<_n-r; i++ )
        for( unsigned int j=0; j<r; j++ )
          for( unsigned int k=0; k<_n-r; k++ )
            Q(r+i,j) += eps * ( U[r+i+(r+k)*_n] * U[j+(r+k)*_n] );

      return _Q;
    }
  //! @brief Return pointer to regularized shape matrix
  const double* invQ()
    {
      if( !_invQ )
        _invQ = _inv( _n, _Q, options.STORAGE );
      return _invQ;
    }

  //! @brief Return lower bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double l
    ( const unsigned int i )
    {
      return( i<0 || i>=_n ? 0.: -std::sqrt(std::max(0.,_diagQ(i))) );
    }
  //! @brief Return upper bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double u
    ( const unsigned int i )
    {
      return( i<0 || i>=_n ? 0.: std::sqrt(std::max(0.,_diagQ(i))) );
    }
  /** @} */
  
private:

  //! @brief Wrapper to LAPACK function <TT>_dsyev</TT> doing eigenvalue decomposition of a symmetric matrix
  static std::pair<double*,double*> _eigen
    ( const unsigned int n, const double*Q, Options::SHAPE_MATRIX STORAGE,
      bool eigv=true );

  //! @brief Compute the square-root of a symmetric matrix
  static double* _sqrt
    ( const unsigned int n, const double*Q, std::pair<double*,double*>&eigQ,
      Options::SHAPE_MATRIX STORAGE, bool&PSDchecked );

  //! @brief Check if a symmetric matrix is positive definite
  static bool _isPSD
    ( const unsigned int n, const double*Q, std::pair<double*,double*>&eigQ,
      Options::SHAPE_MATRIX STORAGE );

  //! @brief Wrapper to LAPACK function <TT>_dgesvd</TT> doing singular value decomposition of a general matrix
  static std::pair< double*,std::pair<double*,double*> > _svd
    ( const unsigned int n, const double*Q, Options::SHAPE_MATRIX STORAGE,
      const bool getu=true, const bool getvt=true );

  //! @brief Wrapper to LAPACK functions <TT>_dsytrf</TT> and <TT>_dsytri</TT> doing factorization and inversion of a symmetric indefinite matrix
  static double* _inv
    ( const unsigned int n, const double*Q,
      Ellipsoid::Options::SHAPE_MATRIX STORAGE );
  
  //! @brief Display the elements of a 2D array
  static void _display
    ( const unsigned int m, const unsigned int n, const double*a,
      const unsigned int lda, const std::string&stra,
      std::ostream&os=std::cout );
  
  //! @brief Pause the program execution and prompt the user
  static void _pause()
    { double tmp;
      std::cout << "ENTER <1> TO CONTINUE" << std::endl;
      std::cin  >> tmp; }
};

////////////////////////////////////////////////////////////////////////

Ellipsoid::Options Ellipsoid::options;

inline bool
Ellipsoid::_isPSD
( const unsigned int n, const double*Q, std::pair<double*,double*>&eigQ,
  Options::SHAPE_MATRIX STORAGE )
{
  if( !eigQ.first || !eigQ.second ){
    delete[] eigQ.first;
    delete[] eigQ.second;
    eigQ = _eigen( n, Q, STORAGE );
  }
  if( !eigQ.first ) return false;
  bool isPSD = true;
  if( *eigQ.first < -std::fabs( options.PSDTOL ) ){
    std::cout << "Min. eigenvalue: " << *eigQ.first << std::endl;
    int dum; std::cin >> dum;
    isPSD = false;
  }
  for( unsigned int i=0; i<n; i++ ){
    if( eigQ.first[i] < 0. ) eigQ.first[i] = 0.;
    else break;
  }
  return isPSD;
}

inline double*
Ellipsoid::_sqrt
( const unsigned int n, const double*Q, std::pair<double*,double*>&eigQ,
  Ellipsoid::Options::SHAPE_MATRIX STORAGE, bool&PSDchecked )
{
  bool isPSD = _isPSD( n, Q, eigQ, STORAGE );
  if( options.PSDCHK ){
    if( !PSDchecked && !isPSD )
      throw Exceptions( Exceptions::NONPSD );
    PSDchecked = true;
  }
  if( !eigQ.first ) return 0;

  double *sqrtQ = new double[n*n];
  for( unsigned int i=0; i<n; i++ )
    for( unsigned int j=0; j<n; j++ ){
      sqrtQ[i+j*n] = 0.;
      for( unsigned int k=0; k<n; k++ )
        sqrtQ[i+j*n] += ( eigQ.second[i+k*n] * std::sqrt( eigQ.first[k] ) )
	                * eigQ.second[j+k*n];
    }

#ifdef MC__ELLIPSOID_DEBUG_SQRT
  double* Qchk = new double[n*n];
  for( unsigned int i=0; i<n; i++ )
    for( unsigned int j=0; j<n; j++ ){
      Qchk[i+j*n] = 0.;
      for( unsigned int k=0; k<n; k++ )
        Qchk[i+j*n] += sqrtQ[i+k*n] * sqrtQ[k+j*n];
    }
  _display( n, n, Qchk, n, "Matrix sqrtQ*sqrtQ", std::cout );
  delete[] Qchk;
#endif
  return sqrtQ;
}

inline std::pair<double*,double*>
Ellipsoid::_eigen
( const unsigned int n, const double*Q, Ellipsoid::Options::SHAPE_MATRIX STORAGE,
  const bool eigv )
{
  if( !n ) return std::make_pair((double*)0,(double*)0);

  int info;
  double*D = new double[n];
  double*U = new double[n*n];
  char JOBZ = (eigv?'V':'N'), UPLO = 'U';
  for( unsigned int j=0, k=0; j<n; j++ )
    switch( STORAGE ){
    case Options::LOWER:
      for( unsigned int i=j; i<n; i++ )  U[j*n+i] = U[i*n+j] = Q[k++]; break;
    case Options::UPPER:
      for( unsigned int i=0; i<=j; i++ ) U[j*n+i] = U[i*n+j] = Q[k++]; break;
    }
#ifdef MC__ELLIPSOID_DEBUG_EIGEN
  _display( n, n, U, n, "Matrix Q", std::cout );
#endif

  // get optimal size
  double worktmp;
  int lwork = -1;
  dsyev_( &JOBZ, &UPLO, &n, U, &n, D, &worktmp, &lwork, &info );

  // perform eigenvalue decomposition
  lwork = (int)worktmp;
  double*work = new double[lwork];
  dsyev_( &JOBZ, &UPLO, &n, U, &n, D, work, &lwork, &info );
#ifdef MC__ELLIPSOID_DEBUG_EIGEN
  _display( n, n, U, n, "Matrix U", std::cout );
  _display( 1, n, D, 1, "Matrix D", std::cout );
#endif
  delete[] work;

#ifdef MC__ELLIPSOID_DEBUG_EIGEN
  std::cout << "INFO: " << info << std::endl;
  _pause();
#endif
  //double* pnull = 0;
  if( info  ){ delete[] D; delete[] U;
               return std::make_pair((double*)0,(double*)0); }
  if( !eigv ){ delete[] U; U = 0; }
  return std::make_pair( D, U );
}

inline std::pair< double*,std::pair<double*,double*> >
Ellipsoid::_svd
( const unsigned int n, const double*Q, Ellipsoid::Options::SHAPE_MATRIX STORAGE,
  const bool getu, const bool getvt )
{
  if( !n ) return std::make_pair((double*)0,std::make_pair((double*)0,(double*)0));

  int info;
  double*A = new double[n*n];
  for( unsigned int j=0, k=0; j<n; j++ )
    switch( STORAGE ){
    case Options::LOWER:
      for( unsigned int i=j; i<n; i++ )  A[j*n+i] = A[i*n+j] = Q[k++]; break;
    case Options::UPPER:
      for( unsigned int i=0; i<=j; i++ ) A[j*n+i] = A[i*n+j] = Q[k++]; break;
    }
#ifdef MC__ELLIPSOID_DEBUG_SVD
  _display( n, n, A, n, "Matrix Q", std::cout );
#endif
  double*S = new double[n];
  char JOBU = ( getu? 'A': 'N' ), JOBVT = ( getvt? 'A': 'N' );
  double*U  = ( getu?  new double[n*n]: 0 );
  double*VT = ( getvt? new double[n*n]: 0 );

  // get optimal size
  double worktmp;
  int lwork = -1;
  dgesvd_( &JOBU, &JOBVT, &n, &n, A, &n, S, U, &n, VT, &n, &worktmp, &lwork, &info );

  // perform singular value decomposition
  lwork = (int)worktmp;
  double*work = new double[lwork];
  dgesvd_( &JOBU, &JOBVT, &n, &n, A, &n, S, U, &n, VT, &n, work, &lwork, &info );
#ifdef MC__ELLIPSOID_DEBUG_SVD
  _display( n, n, U, n, "Matrix U", std::cout );
  _display( 1, n, S, 1, "Matrix S", std::cout );
  _display( n, n, VT, n, "Matrix VT", std::cout );
#endif
  delete[] work;
  delete[] A;

#ifdef MC__ELLIPSOID_DEBUG_SVD
  std::cout << "INFO: " << info << std::endl;
  _pause();
#endif
  if( info ){ delete[] S; delete[] U; delete[] VT;
              return std::make_pair((double*)0,std::make_pair((double*)0,(double*)0)); }
  return std::make_pair( S, std::make_pair( U, VT ) );
}

inline double*
Ellipsoid::_inv
( const unsigned int n, const double*Q,
  Ellipsoid::Options::SHAPE_MATRIX STORAGE )
{
  if( !n ) return (double*)0;

  int info;
  int *ipiv = new int[n];
  double*A = new double[n*n];
  char UPLO = 'U';
  for( unsigned int j=0, k=0; j<n; j++ )
    switch( STORAGE ){
    case Options::LOWER:
      for( unsigned int i=j; i<n; i++ )  A[j*n+i] = A[i*n+j] = Q[k++]; break;
    case Options::UPPER:
      for( unsigned int i=0; i<=j; i++ ) A[j*n+i] = A[i*n+j] = Q[k++]; break;
    }
#ifdef MC__ELLIPSOID_DEBUG_INV
  _display( n, n, A, n, "Matrix Q", std::cout );
#endif

  // get optimal size
  double worktmp;
  int lwork = -1;
  dsytrf_( &UPLO, &n, A, &n, ipiv, &worktmp, &lwork, &info );

  // perform pivoting and inversion
  lwork = (int)worktmp;
  double*work = new double[lwork];
  dsytrf_( &UPLO, &n, A, &n, ipiv, work, &lwork, &info );
  if( !info ) dsytri_( &UPLO, &n, A, &n, ipiv, work, &info );
#ifdef MC__ELLIPSOID_DEBUG_INV
  _display( n, n, A, n, "Matrix Q^-1", std::cout );
#endif
  delete[] work;
  delete[] ipiv;

  if( info ){ delete[] A; return (double*)0; }
  double *Qinv = new double[n*n];
  for( unsigned int j=0, k=0; j<n; j++ )
    switch( STORAGE ){
    case Options::LOWER:
      for( unsigned int i=j; i<n; i++ )  Qinv[k++] = A[i*n+j]; break;
    case Options::UPPER:
      for( unsigned int i=0; i<=j; i++ ) Qinv[k++] = A[i*n+j]; break;
    }
  delete[] A;
  return Qinv;
}

inline void
Ellipsoid::_display
( const unsigned int m, const unsigned int n, const double*a,
  const unsigned int lda, const std::string&stra, std::ostream&os )
{
  os << stra << " =" << std::endl << std::scientific
     << std::setprecision(5);
  for( unsigned int im=0; a && im<m; im++ ){
    for( unsigned int in=0; in<n; in++ ){
      os << a[in*lda+im] << "  ";
    }
    os << std::endl;
  }
  os << std::endl;

  if( os == std::cout || os == std::cerr ) _pause();
}

inline double
Ellipsoid::_diagQ
( const unsigned int i ) const
{
  switch( options.STORAGE ){
  case Options::LOWER:
    return _Q[_n*i-i*(i-1)/2];
  case Options::UPPER: default:
    return _Q[i*(i+1)/2+i];
  }
}

inline double&
Ellipsoid::_diagQ
( const unsigned int i )
{
  switch( options.STORAGE ){
  case Options::LOWER:
    return _Q[_n*i-i*(i-1)/2];
  case Options::UPPER: default:
    return _Q[i*(i+1)/2+i];
  }
}

inline Ellipsoid* minksum_ea
( const unsigned int nE, const Ellipsoid*E, const unsigned int nD,
  const double*D )
{
  if( !nE || !E || !nD || !D ) return 0;

  const unsigned int n = E[0]._n;
  for( unsigned int iE=1; iE< nE; iE++ )
    if( E[iE]._n != n )
      throw Ellipsoid::Exceptions( Ellipsoid::Exceptions::SIZE );

  Ellipsoid *EA = new Ellipsoid[nD];
  for( unsigned int iD=0; iD< nD; iD++ ) EA[iD].set( n );
    
  // Centers
  for( unsigned int i=0; i<E[0]._n; i++ ){
    double ci = E[0]._c[0];
    for( unsigned int iE=1; iE< nE; iE++ ) ci += E[iE]._c[i];
    for( unsigned int iD=0; iD< nD; iD++ ) EA[iD]._c[i] = ci;
  }

  // Shape matrices
  for( unsigned int iD=0; iD< nD; iD++ ){
    double sum_sqrt_lTQl = 0.;
    for( unsigned int iE=0; iE< nE; iE++ ){
      double lTQl = 0.;
      for( unsigned int i=0; i<n; i++ ){
        double Ql_i = 0.;
        for( unsigned int j=0; j<n; j++ )
	  Ql_i += E[iE].Q(i,j) * D[iD*n+j];
        lTQl += D[iD*n+i] * Ql_i;
      }
      sum_sqrt_lTQl += std::sqrt(lTQl);
      for( unsigned int i=0; i<n; i++ )
        for( unsigned int j=0; j<=i; j++ ){
          if( iE ) EA[iD].Q(i,j) += E[iE].Q(i,j) / std::sqrt(lTQl);
	  else     EA[iD].Q(i,j) = E[0].Q(i,j) / std::sqrt(lTQl);
	}
    }
    for( unsigned int i=0; i<n*(n+1)/2; i++ )
      EA[iD]._Q[i] *= sum_sqrt_lTQl;
  }

  return EA;  
}

inline Ellipsoid mtimes
( const Ellipsoid&E, const unsigned int m, const double*A, const double*b=0 )
{
  Ellipsoid mE( m );
  if( !m || !A ) return mE;

  // Transformed center
  for( unsigned int i=0; i<m; i++ ){
    mE._c[i] = ( b? b[i]: 0. );
    for( unsigned int j=0; j<E._n; j++ )
      mE._c[i] += A[i+m*j] * E._c[j];
  }
  
  // Transformed shape
  double*QAT = new double[E._n*m];
  for( unsigned int i=0; i<E._n; i++ )
    for( unsigned int j=0; j<m; j++ ){
      QAT[i+E._n*j] = 0.;
      for( unsigned int k=0; k<E._n; k++ )
        QAT[i+E._n*j] += E.Q(i,k) * A[j+m*k];
    }
  for( unsigned int i=0; i<m; i++ )
    for( unsigned int j=0; j<=i; j++ ){
      mE.Q(i,j) = 0.;
      for( unsigned int k=0; k<E._n; k++ )
        mE.Q(i,j) += A[i+m*k] * QAT[k+E._n*j];
    }
  delete[] QAT;
  
  return mE;
}

inline Ellipsoid inv
( const Ellipsoid&E )
{
  Ellipsoid Einv( E );
  if( !E._n ) return Einv;
  if( Einv.rankQ() < Einv._n ) Einv.regQ();
  const double*Qinv = Einv.invQ();
  for( unsigned int k=0; k<Einv._n*(Einv._n+1)/2; k++ )
    Einv._Q[k] = Qinv[k];
  Einv._reset_auxiliary();
  return Einv;
}  

inline Ellipsoid intersection_ea
( const Ellipsoid&E1, const Ellipsoid&E2 )
{
  if( E1._n != E2._n )
    throw Ellipsoid::Exceptions( Ellipsoid::Exceptions::SIZE );
  const unsigned int n = E1._n;
  
  Ellipsoid Einter( n );
  
  return Einter;
}

inline Ellipsoid ell_unitball
( const unsigned int n )
{
  Ellipsoid Eunit( n );
  for( unsigned int i=0; i<n; i++ ) Eunit._diagQ(i) = 1.;
  return Eunit;
}

inline std::ostream&
operator<<
( std::ostream&out, const Ellipsoid&E)
{
  const int iprec = 5;
  out << std::scientific << std::setprecision(iprec);
  out << "\n n: "  << std::left << std::setw(3) << E._n << std::right;
  if( !E._n ) return out;
  switch( E.options.STORAGE ){
  case Ellipsoid::Options::LOWER:
    for( unsigned int i=0; i<E._n; i++ ){
      if( !i ){ out << "   c:" << std::setw(iprec+8) << E._c[0]
                    << "   Q:" << std::setw(iprec+8) << E._Q[0]; continue; }
      out << std::endl
          << std::setw(12+iprec+8) << E._c[i]
	  << std::setw(5+iprec+8)  << E._Q[i];
      for( unsigned int j=1; j<=i; j++ )
        out << std::setw(iprec+8) << E._Q[E._n*j-j*(j-1)/2+i-j];
    }
    out << std::endl;
    break;
  case Ellipsoid::Options::UPPER: default:
    for( unsigned int i=0; i<E._n; i++ ){
      if( !i ) out << "   c:" << std::setw(iprec+8) << E._c[0]
                   << "   Q:";
      else out << std::setw(12+iprec+8) << E._c[i]
               << std::setw(5+(iprec+8)*i) << " ";      
      for( unsigned int j=i; j<E._n; j++ )
        out << std::setw(iprec+8) << E._Q[j*(j+1)/2+i];
      out << std::endl;
    }
    break;
  }
  return out;
}

} // namespace mc

#endif
