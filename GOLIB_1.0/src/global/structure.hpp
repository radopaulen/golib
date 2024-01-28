// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_STRUCTURE Structure Detection for Factorable Functions
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

mc::Structure is a C++ class that determines the structure of mathematical expressions, namely their sparsity pattern and linearity, for a given set of participating variables. It relies on the operator overloading and function overloading mechanisms of C++. The overloaded operators are: `+', `-', `*', and `/'; the overloaded functions are: `exp', `log', `sqr', `pow', `sqrt', `fabs', `xlog', `min', `max', `cos', `sin', `tan', `acos', `asin' and `atan'.


\section sec_StructureEval How Do I Determine the Structure of a Factorable Function?

Suppose you are given 4 variables \f$x_1,\ldots,x_4\f$ and want to determine the sparsity pattern and linearity of the vector following function
\f{eqnarray*}
  {\bf f}({\bf x}) = \left(\begin{array}{c} f_1({\bf x})\\ f_2({\bf x})\end{array}\right) = \left(\begin{array}{c} x_3 x_4+x_1\\x_1(\exp(x_3-x_4))^2+x_2 \end{array}\right)
\f}

First, define the variables \f$x_1,\ldots,x_4\f$ as

\code
      const int NX = 4;
      Structure X[NX];
      for( int i=0; i<NX; i++ ) X[i].indep(i);
\endcode

Essentially, the first line means that <tt>X</tt> is an array of mc::Structure class objects, and the second line defines X[0],...X[NX-1] as independent variables with indices 0,...,NX-1, respectively.

Once the independent variables \f${\bf x}\f$ have been defined, determine the structure of \f${\bf f}({\bf x})\f$ simply as

\code
      const int NF = 2;
      Structure F[NF] = { X[2]*X[3]+X[0],
                          X[0]*pow(exp(X[2]-X[3]),2)+X[1] };
\endcode

Retrieve the structure - both the sparsity pattern and the linearity - of \f$f_1\f$ and \f$f_2\f$ as

\code
      std::map<int,bool> F0_dep = F[0].dep();
      std::map<int,bool> F1_dep = F[1].dep();
\endcode

You can also display the structure as

\code
      std::cout << "Variable dependence of F[0]: " << F[0] << std::endl;
      std::cout << "Variable dependence of F[1]: " << F[1] << std::endl;
\endcode

The corresponding output is

\verbatim
      Variable dependence of F[0]: { 0L 2NL 3NL }
      Variable dependence of F[1]: { 0NL 1L 2NL 3NL }
\endverbatim

which indicates that X[0], X[2] and X[3] participate in F[0], but not X[1]. Moreover, X[0] participates linearly, unlike X[2] and X[3].

\section sec_StructureErr Errors Encountered in Determining the Structure of a Factorable Function?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type Structure::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during a calculation, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will stop.

Possible errors encountered in determining the structure of a factorable function are:

<TABLE border="1">
<CAPTION><EM>Errors during Structure Determination</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>-33</tt> <TD>Error due to calling a feature not yet implemented in mc::Structure
</TABLE>

*/

#ifndef MC__SPARSE_HPP
#define MC__SPARSE_HPP

#include <iostream>
#include <map>

namespace mc
{

//! @brief C++ class for evaluation of the sparsity pattern of a factorable function
////////////////////////////////////////////////////////////////////////
//! Interval is a C++ class for evaluating the sparsity pattern of a
//! factorable function
////////////////////////////////////////////////////////////////////////
class Structure
////////////////////////////////////////////////////////////////////////
{
  // friends of class Structure for operator overloading
  friend Structure operator+
    ( const Structure& );
  friend Structure operator+
    ( const Structure&, const Structure& );
  friend Structure operator+
    ( const double, const Structure& );
  friend Structure operator+
    ( const Structure&, const double );
  friend Structure operator-
    ( const Structure& );
  friend Structure operator-
    ( const Structure&, const Structure& );
  friend Structure operator-
    ( const double, const Structure& );
  friend Structure operator-
    ( const Structure&, const double );
  friend Structure operator*
    ( const Structure&, const Structure& );
  friend Structure operator*
    ( const Structure&, const double );
  friend Structure operator*
    ( const double, const Structure& );
  friend Structure operator/
    ( const Structure&, const Structure& );
  friend Structure operator/
    ( const Structure&, const double );
  friend Structure operator/
    ( const double, const Structure& );
  friend std::ostream& operator<<
    ( std::ostream&, const Structure& );
  friend bool operator==
    ( const Structure&, const Structure& );
  friend bool operator!=
    ( const Structure&, const Structure& );
  friend bool operator<=
    ( const Structure&, const Structure& );
  friend bool operator>=
    ( const Structure&, const Structure& );
  friend bool operator<
    ( const Structure&, const Structure& );
  friend bool operator>
    ( const Structure&, const Structure& );

  // friends of class Structure for function overloading
  friend Structure inv
    ( const Structure& );
  friend Structure sqr
    ( const Structure& );
  friend Structure exp
    ( const Structure& );
  friend Structure log
    ( const Structure& );
  friend Structure cos
    ( const Structure& );
  friend Structure sin
    ( const Structure& );
  friend Structure tan
    ( const Structure& );
  friend Structure acos
    ( const Structure& );
  friend Structure asin
    ( const Structure& );
  friend Structure atan
    ( const Structure& );
  friend Structure fabs
    ( const Structure& );
  friend Structure sqrt
    ( const Structure& );
  friend Structure xlog
    ( const Structure& );
  friend Structure erf
    ( const Structure& );
  friend Structure erfc
    ( const Structure& );
  friend Structure arh
    ( const Structure&, const double );
  friend Structure pow
    ( const Structure&, const int );
  friend Structure pow
    ( const Structure&, const double );
  friend Structure pow
    ( const Structure&, const Structure& );
  friend Structure min
    ( const Structure&, const Structure& );
  friend Structure max
    ( const Structure&, const Structure& );
  friend Structure min
    ( const unsigned int, const Structure* );
  friend Structure max
    ( const unsigned int, const Structure* );
  friend Structure sum
    ( const unsigned int, const Structure* );
  friend Structure prod
    ( const unsigned int, const Structure* );

public:

  //! @brief Structure exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for Structure exception handling
    enum TYPE{
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented in mc::Structure
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}

    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
  private:
    TYPE _ierr;
  };

  typedef std::map<int,bool> t_Structure;
  typedef t_Structure::iterator it_Structure;
  typedef t_Structure::const_iterator cit_Structure;

  // other operator overloadings (inlined)
  Structure& operator=
    ( const double c )
    { _dep.clear(); return *this; }
  Structure& operator=
    ( const Structure&S )
    { if( this != &S ) _dep = S._dep; return *this; }
  Structure& operator+=
    ( const double c )
    { return *this; }
  Structure& operator+=
    ( const Structure&S )
    { return combine( S ); }
  Structure& operator-=
    ( const double c )
    { return *this; }
  Structure& operator-=
    ( const Structure&S )
    { return combine( S ); }
  Structure& operator*=
    ( const double c )
    { return *this; }
  Structure& operator*=
    ( const Structure&S )
    { return combine( S, false ); }
  Structure& operator/=
    ( const double c )
    { return *this; }
  Structure& operator/=
    ( const Structure&S )
    { return combine( S, false ); }

  /** @defgroup Structure Structure Detection for Factorable Functions
   *  @{
   */
  //! @brief Default constructor (needed to declare arrays of Structure class)
  Structure
    ( const double c=0. )
    {}
  //! @brief Copy constructor
  Structure
    ( const Structure&S ):
    _dep(S._dep)
    {}
  //! @brief Destructor
  ~Structure()
    {}

  //! @brief Sets as independent with index <a>ind</a>
  Structure& indep
    ( const int ind )
    { _dep.clear(); _dep.insert( std::make_pair(ind,true) ); return *this; }

  //! @brief Determines if the current object is dependent on the variable of index <a>ind</a>
  std::pair<bool,bool> dep
    ( const int ind )
    { cit_Structure it = _dep.find(ind);
      return( it==_dep.end()? std::make_pair(false,true):
        std::make_pair(true,it->second) ); }
//   //! @brief Determines if the current object is dependent on the variable <a>S</a>
//   bool dep
//     ( const Structure&S )
//     { return( _dep.size()==combine(*this,S)._dep.size()? true: false); }

  //! @brief Returns the dependency set
  const t_Structure& dep() const
    { return _dep; }
  t_Structure& dep()
    { return _dep; }

  //! @brief Combines with the dependency sets of another variable
  Structure& combine
    ( const Structure&S, const bool linear=true );
  //! @brief Combines the dependency sets of two variables
  static Structure combine
    ( const Structure&S1, const Structure&S2, const bool linear=true );

  //! @brief Turns current dependent variables into nonlinear
  Structure& nonlinear();
  //! @brief Turns current dependent variables into nonlinear
  static Structure nonlinear
    ( const Structure&S );
  /** @} */
  
private:

  //! @brief Dependency set
  t_Structure _dep;
};

////////////////////////////////////////////////////////////////////////

inline Structure&
Structure::nonlinear()
{
  it_Structure it = _dep.begin();
  for( ; it != _dep.end(); ++it ) it->second = false;
  return *this;
}

inline Structure
Structure::nonlinear
( const Structure&S )
{
  Structure S2( S );
  return S2.nonlinear(); 
}

inline Structure&
Structure::combine
( const Structure&S, const bool linear )
{
  cit_Structure cit = S._dep.begin();
  for( ; cit != S._dep.end(); ++cit ){
    std::pair<it_Structure,bool> ins = _dep.insert( *cit );
    if( !ins.second ) ins.first->second = ( ins.first->second && cit->second );
  }
  return( linear? *this: nonlinear() );
}

inline Structure
Structure::combine
( const Structure&S1, const Structure&S2, const bool linear )
{
  Structure S3( S1 );
  return S3.combine( S2, linear );
}

inline std::ostream&
operator<<
( std::ostream&out, const Structure&S)
{
  out << "{ ";
  Structure::cit_Structure iS = S._dep.begin();
  for( ; iS != S._dep.end(); ++iS )
    out << iS->first << (iS->second?"L":"NL") << " ";
  out << "}";
  return out;
}

inline Structure
operator+
( const Structure&S )
{
  return S;
}

inline Structure
operator-
( const Structure&S )
{
  return S;
}

inline Structure
operator+
( const double c, const Structure&S )
{
  return S;
}

inline Structure
operator+
( const Structure&S, const double c )
{
  return S;
}

inline Structure
operator+
( const Structure&S1, const Structure&S2 )
{ 
  return Structure::combine( S1, S2 );
}

inline Structure
sum
( const unsigned int n, const Structure*S )
{
  if( n==2 ) return S[0] + S[1];
  return S[0] + sum( n-1, S+1 );
}

inline Structure
operator-
( const double c, const Structure&S )
{
  return S;
}

inline Structure
operator-
( const Structure&S, const double c )
{
  return S;
}

inline Structure
operator-
( const Structure&S1, const Structure&S2 )
{
  return Structure::combine( S1, S2 );
}

inline Structure
operator*
( const double c, const Structure&S )
{
  return S;
}

inline Structure
operator*
( const Structure&S, const double c )
{
  return S;
}

inline Structure
operator*
( const Structure&S1, const Structure&S2 )
{
  if( S1._dep.empty() ) return S2;
  if( S2._dep.empty() ) return S1;
  return Structure::combine( S1, S2, false );
}

inline Structure
prod
( const unsigned int n, const Structure*S )
{
  if( n==2 ) return S[0] * S[1];
  return S[0] * prod( n-1, S+1 );
}

inline Structure
operator/
( const Structure&S, const double c )
{
  return S;
}

inline Structure
operator/
( const double c, const Structure&S )
{
  return S;
}

inline Structure
operator/
( const Structure&S1, const Structure&S2 )
{
  if( S1._dep.empty() ) return inv( S2 );
  if( S2._dep.empty() ) return S1;
  return Structure::combine( S1, S2, false );
}

inline Structure
inv
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
sqr
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
exp
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
arh
( const Structure &S, const double a )
{
  return Structure::nonlinear( S );
}

inline Structure
log
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
xlog
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
erf
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
erfc
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
sqrt
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
fabs
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
pow
( const Structure&S, const int n )
{
  if( n == 0 ){ Structure S2; return S2; }
  if( n == 1 ) return S;
  return Structure::nonlinear( S );
}

inline Structure
pow
( const Structure&S, const double a )
{
  return Structure::nonlinear( S );
}

inline Structure
pow
( const Structure&S1, const Structure&S2 )
{
  if( S1._dep.empty() ) return Structure::nonlinear( S2 );
  if( S2._dep.empty() ) return Structure::nonlinear( S1 );
  return Structure::combine( S1, S2, false );
}

inline Structure
min
( const Structure&S1, const Structure&S2 )
{
  if( S1._dep.empty() ) return S2;
  if( S2._dep.empty() ) return S1;
  return Structure::combine( S1, S2, false );
}

inline Structure
max
( const Structure&S1, const Structure&S2 )
{
  if( S1._dep.empty() ) return S2;
  if( S2._dep.empty() ) return S1;
  return Structure::combine( S1, S2, false );
}

inline Structure
min
( const unsigned int n, const Structure*S )
{
  if( n==2 ) return min( S[0], S[1] );
  return min( S[0], min( n-1, S+1 ) );
}

inline Structure
max
( const unsigned int n, const Structure*S )
{
  if( n==2 ) return max( S[0], S[1] );
  return max( S[0], max( n-1, S+1 ) );
}

inline Structure
cos
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
sin
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
tan
( const Structure&S )
{
  return Structure::nonlinear( S );
}

inline Structure
acos
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
asin
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline Structure
atan
( const Structure &S )
{
  return Structure::nonlinear( S );
}

inline bool
operator==
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep == S2._dep );
}

inline bool
operator!=
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep != S2._dep );
}

inline bool
operator<=
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep <= S2._dep );
}

inline bool
operator>=
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep >= S2._dep );
}

inline bool
operator<
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep < S2._dep );
}

inline bool
operator>
( const Structure&S1, const Structure&S2 )
{
  return( S1._dep > S2._dep );
}

} // namespace mc

#endif
