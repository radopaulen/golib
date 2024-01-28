// Copyright (C) 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.

/*!
\page page_CHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
\author Jai Rajyaguru, Mario E. Villanueva, Beno&icirc;t Chachuat

A \f$q\f$th-order Chebyshev model of a multivariate function \f$f:\mathbb{R}^n\to\mathbb{R}\f$ that is (at least) \f$(q+1)\f$-times continuously differentiable on the domain \f$D\f$, consists of the \f$q^{\rm th}\f$-order multivariate polynomial \f$\mathcal P\f$, matching either the truncated Chebyshev expansion of \f$\mathcal P\f$ up to order \f$q^{\rm th}\f$ or the \f$q^{\rm th}\f$-order Chebyshev interpolated polynomial, plus a remainder term \f$\mathcal R\f$, so that
\f{align*}
  f({x}) \in \mathcal P({x}-\hat{x}) \oplus \mathcal R, \quad \forall {x}\in D.
\f}
The polynomial part \f$\mathcal P\f$ is propagated symbolically and accounts for functional dependencies. The remainder term \f$\mathcal R\f$, on the other hand, is traditionally computed using interval analysis [Brisebarre & Joldes, 2010]; see figure below. More generally, convex/concave bounds or an ellipsoidal enclosure can be computed for the remainder term of vector-valued functions too. In particular, it can be established that the remainder term has convergence order (no less than) \f$q+1\f$ with respect to the diameter of the domain set \f$D\f$ under mild conditions [Bompadre <I>et al.</I>, 2012].

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html Chebyshev_model.png</TD>
</TR>
</TABLE></CENTER>

The classes mc::CModel and mc::CVar provide an implementation of Chebyshev model arithmetic. We note that mc::CModel / mc::CVar is <b>not a verified implementation</b> in the sense that rounding errors are not accounted for in propagating the coefficients in the multivariate polynomial part, which are treated as floating-point numbers.

The implementation of mc::CModel and mc::CVar relies on the operator/function overloading mechanism of C++. This makes the computation of Chebyshev models both simple and intuitive, similar to computing function values in real arithmetics or function bounds in interval arithmetic (see \ref page_INTERVAL). Moreover, mc::CVar can be used as the template parameter of other available types in MC++; for instance, mc::CVar can be used in order to propagate the underlying interval bounds in mc::McCormick. Likewise, mc::CVar can be used as the template parameter of the types fadbad::F, fadbad::B and fadbad::T of <A href="http://www.fadbad.com/fadbad.html">FADBAD++</A> for computing Chebyshev models of either the partial derivatives or the Chebyshev coefficients of a factorable function (see \ref sec_CHEBYSHEV_fadbad).

mc::CModel and mc::CVar themselves are templated in the type used to propagate bounds on the remainder term. By default, mc::CModel and mc::CVar can be used with the non-verified interval type mc::Interval of MC++. For reliability, however, it is strongly recommended to use verified interval arithmetic such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.hpp</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.hpp</tt>). As already noted, convex/concave bounds on the remainder term can also be propagated by using the type mc::McCormick of MC++, thereby enabling McCormick-Chebyshev models.

As well as propagating Chebyshev models for factorable functions, mc::CModel and mc::CVar provide support for computing bounds on the Chebyshev model range (multivariate polynomial part). We note that computing exact bounds for multivariate polynomials is a hard problem in general. Instead, a number of computationally tractable, yet typically conservative, bounding approaches are implemented in mc::CModel and mc::CVar, which include:
- Bounding every monomial term independently and adding these bounds;
- Bounding the first- and diagonal second-order terms exactly and adding bounds for the second-order off-diagonal and higher-order terms computed independently [Lin & Stadtherr, 2007];
- Bounding the terms up to order 2 based on an eigenvalue decomposition of the corresponding Hessian matrix and adding bounds for the higher-order terms computed independently;
- Expressing the multivariate polynomial in Bernstein basis, thereby providing bounds as the minimum/maximum among all Bernstein coefficients [Lin & Rokne, 1995; 1996].
.

Examples of Chebyshev models (blue lines) constructed with mc::CModel and mc::CVar are shown on the figure below for the factorable function \f$f(x)=x \exp(-x^2)\f$ (red line) for \f$x\in [-0.5,1]\f$. Also shown on these plots are the interval bounds computed from the Chebyshev models.

<CENTER><TABLE BORDER=0>
<TR>
<TD>\image html CM-1D.png</TD>
</TR>
</TABLE></CENTER>


\section sec_CHEBYSHEV_I How do I compute a Chebyshev model with interval remainder bound of a factorable function?

Suppose we want to compute a 4th-order Chebyshev model for the real-valued function \f$f(x,y)=x\exp(x+y^2)-y^2\f$ with \f$(x,y)\in [1,2]\times[0,1]\f$. For simplicity, bounds on the remainder terms are computed using the default interval type mc::Interval here:

\code
      #include "interval.hpp"
      #include "cmodel.hpp"
      typedef mc::Interval I;
      typedef mc::CModel<I> CM;
      typedef mc::CVar<I> CV;
\endcode

First, the number of independent variables in the factorable function (\f$x\f$ and \f$y\f$ here) as well as the order of the Chebyshev model (4th order here) are specified by defining an mc::CModel object as:

\code
      CM mod( 2, 4 );
\endcode

Next, the variables \f$x\f$ and \f$y\f$ are defined as follows:

\code
      CV X( &mod, 0, I(1.,2.) );
      CV Y( &mod, 1, I(0.,1.) );
\endcode

Essentially, the first line means that <tt>X</tt> is a variable of class mc::CVar, participating in the Chebyshev model <tt>mod</tt>, belonging to the interval \f$[1,2]\f$, and having index 0 (indexing in C/C++ start at 0 by convention!). The same holds for the Chebyshev variable <tt>Y</tt>, participating in the model <tt>mod</tt>, belonging to the interval \f$[0,1]\f$, and having index 1.

Having defined the variables, a Chebyshev model of \f$f(x,y)=x\exp(x+y^2)-y^2\f$ on \f$[1,2]\times[0,1]\f$ at the mid-point \f$(\frac{3}{2},\frac{1}{2})\f$ is simply computed as:

\code
      CV F = X*exp(X+pow(Y,2))-pow(Y,2);
\endcode

This model can be displayed to the standard output as:

\code
      std::cout << "f Chebyshev model: " << F << std::endl;
\endcode

which produces the following output:

\verbatim
f Chebyshev model: 
   a0    =  8.38199e+00     0  0
   a1    =  1.90755e+00     0  1
   a2    =  3.59621e+00     1  0
   a3    =  7.47482e-01     0  2
   a4    =  9.00782e-01     1  1
   a5    =  6.30186e-01     2  0
   a6    =  1.56945e-01     0  3
   a7    =  3.35238e-01     1  2
   a8    =  1.55141e-01     2  1
   a9    =  6.67468e-02     3  0
   a10   =  3.49519e-02     0  4
   a11   =  6.58449e-02     1  3
   a12   =  6.04330e-02     2  2
   a13   =  1.80397e-02     3  1
   a14   =  5.41191e-03     4  0
   R     =  [ -2.09182e+00 :  2.22652e+00 ]
   B     =  [ -1.02564e+01 :  3.93973e+01 ]
\endverbatim

<tt>a0</tt>,...,<tt>a14</tt> refer to the coefficients of the monomial terms in the Chebyshev model, with the corresponding variable orders given in the subsequent columns. The remainder term as well as the Chebyshev model range estimator are reported next.

Other operations involve retreiving the remainder bound, centering the remainder term in a Chebyshev model, or computing the value of its polynomial part at a given point:

\code
      I B = F.B();
      F.C();
      double x[2] = { 0.5, 1.5 };
      double Pval = F.P( x );
\endcode

See the documentations of mc::CModel and mc::CVar for a complete list of member functions. 


\section sec_CHEBYSHEV_fct Which functions are overloaded for Chebyshev model arithmetic?

mc::CVar overloads the usual functions <tt>exp</tt>, <tt>log</tt>, <tt>sqr</tt>, <tt>sqrt</tt>, <tt>pow</tt>, <tt>inv</tt>, <tt>cos</tt>, <tt>sin</tt>, <tt>tan</tt>, <tt>acos</tt>, <tt>asin</tt>, <tt>atan</tt>. Unlike mc::Interval and mc::McCormick, the functions <tt>min</tt>, <tt>max</tt> and <tt>fabs</tt> are not overloaded in mc::CVar as they are nonsmooth. Moreover, mc::CVar defines the following functions:
- <tt>inter(x,y,z)</tt>, computing a Chebyshev model of the intersection \f$x = y\cap z\f$ of two Chebyshev models and returning true/false if the intersection is nonempty/empty. With Chebyshev models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this intersection is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & [\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C})] \cap [\mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})]\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ the Chebyshev model range bounder, and \f$\eta\f$ a real scalar in \f$[0,1]\f$. Choosing \f$\eta=1\f$ amounts to setting the polynomial part \f$\mathcal P_{x}\f$ as \f$\mathcal P_y\f$, whereas \f$\eta=0\f$ sets \f$\mathcal P_{x}\f$ as \f$\mathcal P_z\f$. The parameter \f$\eta\f$ can be defined in mc::CModel::Options::REF_POLY.
- <tt>hull(x,y)</tt>, computing a Chebyshev model of the union \f$x = y\cup z\f$ of two Chebyshev models. With Chebyshev models \f$\mathcal P_y\oplus\mathcal R_y\f$ and \f$\mathcal P_z\oplus\mathcal R_z\f$, this union is computed as follows:
\f{align*}
  \mathcal P_{x} =\ & (1-\eta) \mathcal P_y^{\rm C} + \eta \mathcal P_z^{\rm C}\\
  \mathcal R_{x} =\ & {\rm hull}\{\mathcal R_y^{\rm C}\oplus\eta\mathcal{B}(\mathcal P_y^{\rm C}-\mathcal P_z^{\rm C}), \mathcal R_z^{\rm C}\oplus (1-\eta)\mathcal{B}(\mathcal P_z^{\rm C}-\mathcal P_y^{\rm C})\}\,.
\f}
with \f$\mathcal{B}(\cdot)\f$ and \f$\eta\f$ as previously.


\section sec_CHEBYSHEV_opt How are the options set for the computation of a Chebyshev model?

The class mc::CModel has a public member called mc::CModel::options that can be used to set/modify the options; e.g.,

\code
      model.options.BOUNDER_TYPE = CM::Options::EIGEN;
      model.options.SCALE_VARIABLES = true;
\endcode

The available options are the following:

<TABLE border="1">
<CAPTION><EM>Options in mc::CModel::Options: name, type and description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>BOUNDER_TYPE</tt> <TD><tt>mc::CModel::Options::BOUNDER</tt> <TD>mc::CModel::Options::LSB
         <TD>Chebyshev model range bounder.
     <TR><TH><tt>BOUNDER_ORDER</tt> <TD><tt>unsigned int</tt> <TD>0
         <TD>Order of Bernstein polynomial for Chebyshev model range bounding, when mc::CModel::options::BOUNDER_TYPE = mc::CModel::options::BERNSTEIN is selected. Only values greater than the actual Chebyshev model order are accounted for; see [Lin & Rokne, 1996].
     <TR><TH><tt>CENTER_REMAINDER</tt> <TD><tt>bool</tt> <TD>false
         <TD>Whether to center the remainder term during Chebyshev model propagation.
     <TR><TH><tt>REF_POLY</tt> <TD><tt>double</tt> <TD>0.
         <TD>Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_CHEBYSHEV_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
     <TR><TH><tt>DISPLAY_DIGITS</tt> <TD><tt>unsigned int</tt> <TD>5
         <TD>Number of digits in output stream for Chebyshev model coefficients.
</TABLE>


\section sec_CM_err Errors What errors can I encounter during computation of a Chebyshev model?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, a class object of type mc::CModel::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during the computation of a Chebyshev model, and then make the appropriate changes. Should an exception be thrown and not caught by the calling program, the execution will abort.

Possible errors encountered during the computation of a Chebyshev model are:

<TABLE border="1">
<CAPTION><EM>Errors during the Computation of a Chebyshev Model</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Division by zero
     <TR><TH><tt>2</tt> <TD>Failed to compute eigenvalue decomposition in range bounder CModel::Options::EIGEN
     <TR><TH><tt>3</tt> <TD>Failed to compute the maximum gap between a univariate term and its Bernstein model
     <TR><TH><tt>-1</tt> <TD>Number of variable in Chebyshev model must be nonzero
     <TR><TH><tt>-2</tt> <TD>Failed to construct Chebyshev variable
     <TR><TH><tt>-3</tt> <TD>Chebyshev model bound does not intersect with bound in template parameter arithmetic
     <TR><TH><tt>-4</tt> <TD>Operation between Chebyshev variables linked to different Chebyshev models
     <TR><TH><tt>-5</tt> <TD>Maximum size of Chebyshev model reached (monomials indexed as unsigned int)
     <TR><TH><tt>-33</tt> <TD>Feature not yet implemented in mc::CModel
</TABLE>

Moreover, exceptions may be thrown by the template parameter class itself.


\section sec_CM_refs References

- Brisebarre, N., and M. Joldes, <A href="http://hal.archives-ouvertes.fr/docs/00/48/17/37/PDF/RRLIP2010-13.pdf">Chebyshev Interpolation Polynomial-based Tools for Rigorous Computing</A>, <i>Research Report No RR2010-13</i>, Ecole Normale Sup&eaccute;rieure de Lyon, Unit&eaccute; Mixte de Recherche CNRS-INRIA-ENS LYON-UCBL No 5668, 2010
.
*/

#ifndef MC__CMODEL_H
#define MC__CMODEL_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <stdarg.h>
#include <cassert>
#include <climits>
#include <limits>
#include <stdlib.h>
#include <complex>

#include "mcfunc.hpp"
#include "mcop.hpp"
#include "mclapack.hpp"

#undef  MC__CMODEL_DEBUG
#undef  MC__CMODEL_DEBUG_SCALE
#define MC__CMODEL_CHECK
#undef  MC__CVAR_DEBUG_EXP
#undef  MC__CVAR_DEBUG_BERSTEIN

namespace mc
{

template <typename T> class CVar;
typedef unsigned long long CM_size;

//! @brief C++ class for Chebyshev model computation of factorable function - Chebyshev model environment
////////////////////////////////////////////////////////////////////////
//! mc::CModel is a C++ class for definition of Chebyshev model
//! environment. Propagation of Chebyshev models for factorable functions
//! is via the C++ class mc::CVar. The template parameter corresponds
//! to the type used to propagate the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class CModel
////////////////////////////////////////////////////////////////////////
{
  friend class CVar<T>;
  template <typename U> friend class CModel;

  template <typename U> friend CVar<U> pow
    ( const CVar<U>&, const int );

public:

  /** @addtogroup CHEBYSHEV Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Chebyshev model environment for <tt>nvar</tt> variables and order <tt>nord</tt>
  CModel
    ( const unsigned int nvar, const unsigned int nord )
    { _size( nvar, nord ); }
  //! @brief Destructor of Chebyshev model environment

  ~CModel()
    { _cleanup(); }

  //! @brief Number of variables in Chebyshev model environment
  unsigned int nvar() const
    { return _nvar; };

  //! @brief Order of Chebyshev model environment
  unsigned int nord() const
    { return _nord; };

  //! @brief Total number of monomial terms in Chebyshev variable
  unsigned int nmon() const
    { return _nmon; };

  //! @brief Const pointer to array of size <tt>nmon()*nvar()</tt> with variable exponents on every term. The exponent for variable <tt>ivar</tt> in monomial term <tt>imon</tt> is at position <tt>imon*nvar()+ivar</tt>.
  const unsigned int* expmon() const
    { return _expmon; };

  //! @brief Index of term whose variable exponents match those in array <tt>iexp</tt> (of size <tt>nvar()</tt>)
  unsigned int loc_expmon
    ( const unsigned int *iexp ) const
    { return _loc_expmon( iexp ); }

  //! @brief Get const pointer to array of size <tt>nvar</tt> with original variable bounds
  const T* bndvar() const
    { return _bndvar; }

  //! @brief Exceptions of mc::CModel
  class Exceptions
  {
  public:
    //! @brief Enumeration type for CModel exception handling
    enum TYPE{
      DIV=1,	//!< Division by zero scalar
      EIGEN,	//!< Failed to compute eigenvalue decomposition in range bounder CModel::Options::EIGEN
      BERNSTEIN,//!< Failed to compute the maximum gap between a univariate term and its Bernstein model
      SIZE=-1,	//!< Number of variable in Chebyshev model must be nonzero
      INIT=-2,	//!< Failed to construct Chebyshev variable
      INCON=-3, //!< Chebyshev model bound does not intersect with bound in template parameter arithmetic
      CMODEL=-4,//!< Operation between Chebyshev variables linked to different Chebyshev models
      MAXSIZE=-5,//!< Maximum size of Chebyshev model reached (monomials indexed as unsigned int)
      UNDEF=-33 //!< Feature not yet implemented in mc::CModel
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Error flag
    int ierr(){ return _ierr; }
    //! @brief Error description
    std::string what(){
      switch( _ierr ){
      case DIV:
        return "mc::CModel\t Division by zero scalar";
      case EIGEN:
        return "mc::CModel\t Range bounder with eigenvalue decomposition failed";
      case BERNSTEIN:
        return "mc::CModel\t Bernstein remainder evalution failed";
      case SIZE:
        return "mc::CModel\t Inconsistent Chebyshev model dimension";
      case INIT:
        return "mc::CModel\t Chebyshev variable initialization failed";
      case INCON:
        return "mc::CModel\t Inconsistent bounds with template parameter arithmetic";
      case CMODEL:
        return "mc::CModel\t Operation between Chebyshev variables in different Chebyshev model environment not allowed";
      case MAXSIZE:
        return "mc::CModel\t Maximum size in Chebyshev model reached";
      case UNDEF:
        return "mc::CModel\t Feature not yet implemented in mc::CModel class";
      default:
        return "mc::CModel\t Undocumented error";
      }
    }

  private:
    TYPE _ierr;
  };

  //! @brief Options of mc::CModel
  struct Options
  {
    //! @brief Constructor of mc::CModel::Options
    Options():
      BOUNDER_TYPE(LSB), BOUNDER_ORDER(0), 
      CENTER_REMAINDER(false), REF_POLY(0.), DISPLAY_DIGITS(5)
      {}
    //! @brief Copy constructor of mc::CModel::Options
    template <typename U> Options
      ( U&options )
      : BOUNDER_TYPE( options.BOUNDER_TYPE ),
        BOUNDER_ORDER( options.BOUNDER_ORDER ),
        CENTER_REMAINDER( options.CENTER_REMAINDER ),
        REF_POLY(options.REF_POLY),
	DISPLAY_DIGITS(options.DISPLAY_DIGITS)
      {}
    //! @brief Assignment of mc::CModel::Options
    template <typename U> Options& operator =
      ( U&options ){
        BOUNDER_TYPE     = options.BOUNDER_TYPE;
        BOUNDER_ORDER    = options.BOUNDER_ORDER;
        CENTER_REMAINDER = options.CENTER_REMAINDER;
        REF_POLY         = options.REF_POLY,
	DISPLAY_DIGITS   = options.DISPLAY_DIGITS;
        return *this;
      }
    //! @brief Chebyshev model range bounder option
    enum BOUNDER{
      NAIVE=0,	//!< Naive polynomial range bounder
      LSB,	//!< Lin & Stadtherr range bounder
      EIGEN,	//!< Eigenvalue decomposition-based bounder
      BERNSTEIN,//!< Bernstein range bounder
      HYBRID	//!< Hybrid LSB + EIGEN range bounder
    };
    //! @brief Chebyshev model range bounder - See \ref sec_CHEBYSHEV_opt
    BOUNDER BOUNDER_TYPE;
    //! @brief Order of Bernstein polynomial for Chebyshev model range bounding (no less than Chebyshev model order!). Only if mc::CModel::options::BOUNDER_TYPE is set to mc::CModel::options::BERNSTEIN.
    unsigned int BOUNDER_ORDER;
    //! @brief Array of Chebyshev model range bounder names (for display)
    static const std::string BOUNDER_NAME[5];
    //! @brief Whether to center the remainder term during Chebyshev model propagation
    bool CENTER_REMAINDER;
    //! @brief Scalar in \f$[0,1]\f$ related to the choice of the polynomial part in the overloaded functions mc::inter and mc::hull (see \ref sec_TAYLOR_fct). A value of 0. amounts to selecting the polynomial part of the left operand, whereas a value of 1. selects the right operand.
    double REF_POLY;
    //! @brief Number of digits in output stream for Chebyshev model coefficients.
    unsigned int DISPLAY_DIGITS;
  } options;
  /** @} */

private:  
  //! @brief Order of Chebyshev model
  unsigned int _nord;
  //! @brief Number of variables in Chebyshev model
  unsigned int _nvar;
  //! @brief Total number of monomial terms in Chebyshev model
  unsigned int _nmon;
  //! @brief Array of size <tt>_nord</tt> with indices of first monomial term of order <tt>iord=1,...,_nord</tt> in Chebyshev model
  unsigned int *_posord;
  //! @brief Order used to size _posord
  unsigned _posord_size;
  //! @brief Array of size <tt>_nmon*_nvar</tt> with variable exponents in all terms. The exponent for variable <tt>ivar</tt> in term <tt>imon</tt> is at location <tt>imon*nvar()+ivar</tt>.
  unsigned int *_expmon;
  //! @brief Number of monomial coefficients used to size _expmon
  unsigned _expmon_size;
  //! @brief Triple array of size <tt>(_nmon+1,<=_nmon,2^_nvar)</tt> with indices of terms from product of two Chebyshev terms <tt>imon=1,...,_nmon</tt> and <tt>jmon=1,...,_nmon</tt>.
  unsigned int ***_prodmon;
  //! @brief Array of size <tt>_nvar</tt> with bounds on original variables <tt>ivar=1,...,_nvar</tt>
  T *_bndvar;
  //! @brief Have any of the model variables been modified?
  bool _modvar;
  //! @brief Array of <tt>(_nvar+_nord-1)*(_nord+1)</tt> contining binomial coefficients
  CM_size *_binom;
  //! @brief Maximum binomial coefficients in array _binom
  std::pair<unsigned int, unsigned int> _binom_size;

  //! @brief Internal Chebyshev variable to speed-up computations and reduce dynamic allocation
  CVar<T>* _CV;

  //! @brief Set Chebyshev model order <tt>nord</tt> and number of variables <tt>nvar</tt>
  void _size
    ( const unsigned int nvar, const unsigned int nord );

  //! @brief Populate array <tt>_bndmon</tt>
  void _set_bndmon();

  //! @brief Populate array <tt>_posord</tt> up to order <tt>nord</tt>
  void _set_posord
    ( const unsigned int nord );

  //! @brief Extend array <tt>_posord</tt> up to maximum order <tt>maxord</tt>
  void _ext_posord
    ( const unsigned int maxord );

  //! @brief Populate array <tt>_expmon</tt> up to order <tt>nord</tt>
  void _set_expmon
    ( const unsigned int nord );

  //! @brief Extend array <tt>_expmon</tt> up to order <tt>maxord</tt> to accomodate <tt>maxmon</tt> coefficients
  void _ext_expmon
    ( const unsigned int maxord, const bool full=false );
  
  //! @brief Generate variable exponents <tt>iexp</tt> for subsequent monomial order <tt>iord</tt>
  void _next_expmon
    ( unsigned int *iexp, const unsigned int iord ) const;

  //! @brief Populate array _prodmon w/ exponents resulting from the product of two monomial terms 1,...,nmon
  void _set_prodmon();

  //! @brief Generate indices of basis functions in product monics recursively
  void _set_prodmon_exp
    ( unsigned int&ivar, const unsigned int*iexp, const unsigned int iord,
      unsigned int*iprodexp, unsigned int&nprodmon );

  //! @brief Get index of monomial term with variable exponents <tt>iexp</tt> in <tt>1,...,_nmon</tt>
  unsigned int _loc_expmon
    ( const unsigned int *iexp ) const;
    
  //! @brief Populate array <tt>_binom</tt> with binomial coefficients up to order <tt>nord</tt>
  void _set_binom
    ( const unsigned int nord );
    
  //! @brief Extend array <tt>_binom</tt> with binomial coefficients up to order <tt>nord</tt>
  void _ext_binom
    ( const unsigned int nord );

  //! @brief Get binomial coefficient \f$\left(\stackrel{n}{k}\right)\f$
  CM_size _get_binom
    ( const unsigned int n, const unsigned int k ) const;

  //! @brief Populate array <tt>_bndvar</tt> for variable <tt>ivar</tt> with range <tt>X</tt>
  void _set_bndvar
    ( const unsigned int ivar, const T&X );

  //! @brief Reset the variable bound arrays
  void _reset();

  //! @brief Clean up the arrays
  void _cleanup();

  //! @brief Recursive calculation of nonnegative integer powers
  CVar<T> _intpow
    ( const CVar<T>&CV, const int n );

} ;

template <typename T> const std::string CModel<T>::Options::BOUNDER_NAME[5]
  = { "NAIVE", "LSB", "EIGEN", "BERNSTEIN", "HYBRID" };

//! @brief C++ class for Chebyshev model computation of factorable function - Chebyshev model propagation
////////////////////////////////////////////////////////////////////////
//! mc::CVar is a C++ class for propagation of Chebyshev models through
//! factorable functions. The template parameter corresponds to the
//! type used in computing the remainder bound.
////////////////////////////////////////////////////////////////////////
template <typename T>
class CVar
////////////////////////////////////////////////////////////////////////
{
  template <typename U> friend class CVar;
  template <typename U> friend class CModel;

  template <typename U> friend CVar<U> operator+
    ( const CVar<U>& );
  template <typename U, typename V> friend CVar<U> operator+
    ( const CVar<U>&, const CVar<V>& );
  template <typename U, typename V> friend CVar<U> operator+
    ( const CVar<U>&, const V& );
  template <typename U, typename V> friend CVar<U> operator+
    ( const V&, const CVar<U>& );
  template <typename U> friend CVar<U> operator+
    ( const double, const CVar<U>& );
  template <typename U> friend CVar<U> operator+
    ( const CVar<U>&, const double );
  template <typename U> friend CVar<U> operator-
    ( const CVar<U>& );
  template <typename U, typename V> friend CVar<U> operator-
    ( const CVar<U>&, const CVar<V>& );
  template <typename U, typename V> friend CVar<U> operator-
    ( const CVar<U>&, const V& );
  template <typename U, typename V> friend CVar<U> operator-
    ( const V&, const CVar<U>& );
  template <typename U> friend CVar<U> operator-
    ( const double, const CVar<U>& );
  template <typename U> friend CVar<U> operator-
    ( const CVar<U>&, const double );
  template <typename U> friend CVar<U> operator*
    ( const CVar<U>&, const CVar<U>& );
  template <typename U> friend CVar<U> operator*
    ( const double, const CVar<U>& );
  template <typename U> friend CVar<U> operator*
    ( const CVar<U>&, const double );
  template <typename U> friend CVar<U> operator*
    ( const U&, const CVar<U>& );
  template <typename U> friend CVar<U> operator*
    ( const CVar<U>&, const U& );
  template <typename U> friend CVar<U> operator/
    ( const CVar<U>&, const CVar<U>& );
  template <typename U> friend CVar<U> operator/
    ( const double, const CVar<U>& );
  template <typename U> friend CVar<U> operator/
    ( const CVar<U>&, const double );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const CVar<U>& );

  template <typename U> friend CVar<U> inv
    ( const CVar<U>& );
  template <typename U> friend CVar<U> sqr
    ( const CVar<U>& );
  template <typename U> friend CVar<U> sqrt
    ( const CVar<U>& );
  template <typename U> friend CVar<U> exp
    ( const CVar<U>& );
  template <typename U> friend CVar<U> log
    ( const CVar<U>& );
  template <typename U> friend CVar<U> xlog
    ( const CVar<U>& );
  template <typename U> friend CVar<U> pow
    ( const CVar<U>&, const int );
  template <typename U> friend CVar<U> pow
    ( const CVar<U>&, const double );
  template <typename U> friend CVar<U> pow
    ( const double, const CVar<U>& );
  template <typename U> friend CVar<U> pow
    ( const CVar<U>&, const CVar<U>& );
  template <typename U> friend CVar<U> monomial
    ( const unsigned int, const CVar<U>*, const int* );
  template <typename U> friend CVar<U> cos
    ( const CVar<U>& );
  template <typename U> friend CVar<U> sin
    ( const CVar<U>& );
  template <typename U> friend CVar<U> tan
    ( const CVar<U>& );
  template <typename U> friend CVar<U> acos
    ( const CVar<U>& );
  template <typename U> friend CVar<U> asin
    ( const CVar<U>& );
  template <typename U> friend CVar<U> atan
    ( const CVar<U>& );
  template <typename U> friend CVar<U> hull
    ( const CVar<U>&, const CVar<U>& );
  template <typename U> friend bool inter
    ( CVar<U>&, const CVar<U>&, const CVar<U>& );

private:
  //! @brief Pointer to Chebyshev model environment
  CModel<T> *_CM;

  //! @brief Pointer to internal Chebyshev variable in Chebyshev model environment
  CVar<T>* _CV() const
    { return _CM->_CV; };
  //! @brief Order of Chebyshev model environment
  unsigned int _nord() const
    { return _CM->_nord; };
  //! @brief Number of variables in Chebyshev model environment
  unsigned int _nvar() const
    { return _CM->_nvar; };
  //! @brief Total number of monomial terms in Chebyshev model
  unsigned int _nmon() const
    { return _CM->_nmon; };
  //! @brief Index of first term of order <tt>iord</tt> in Chebyshev model
  unsigned int _posord
    ( const unsigned int iord ) const
    { return _CM->_posord[iord]; };
  //! @brief Const pointer to array of size <tt>nvar()</tt> of variable exponents in term <tt>imon</tt>
  const unsigned int* _expmon
    ( const unsigned int imon ) const
    { return _CM->_expmon+imon*_CM->_nvar; };
  //! @brief Indices of terms from product of two terms <tt>imon</tt> and <tt>jmon</tt> (size 2^nvar())
  unsigned int* _prodmon
    ( const unsigned int imon, const unsigned int jmon ) const
    { return _CM->_prodmon[imon][jmon]; };
  //! @brief Original bounds on variable <tt>ivar</tt>
  const T& _bndvar
    ( const unsigned int ivar ) const
    { return _CM->_bndvar[ivar]; };

public:
  /** @addtogroup TAYLOR Chebyshev Model Arithmetic for Factorable Functions
   *  @{
   */
  //! @brief Constructor of Chebyshev variable for a real scalar
  CVar
    ( const double d=0. );
  //! @brief Constructor of Chebyshev variable for a remainder bound
  CVar
    ( const T&B );
  //! @brief Constructor of Chebyshev variable with index <a>ix</a> (starting from 0) and bounded by <a>X</a>
  CVar
    ( CModel<T>*CM, const unsigned int ix, const T&X );
  //! @brief Copy constructor of Chebyshev variable
  CVar
    ( const CVar<T>&CV );
  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with implicit type conversion)
  template <typename U> CVar
    ( CModel<T>*&CM, const CVar<U>&CV );
  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with explicit type conversion as given by class member function <a>method</a>)
  template <typename U> CVar
    ( CModel<T>*&CM, const CVar<U>&CV, const T& (U::*method)() const );
  //! @brief Copy constructor of Chebyshev variable in different Chebyshev model environment (with explicit type conversion as given by non-class member function <a>method</a>)
  template <typename U> CVar
    ( CModel<T>*&CM, const CVar<U>&CV, T (*method)( const U& ) );

  //! @brief Destructor of Chebyshev variable
  ~CVar()
    { delete [] _coefmon; delete [] _bndord; }

  //! @brief Set Chebyshev variable with index <a>ix</a> (starting from 0),  bounded by <a>X</a>, and with reference point at mid-point <a>Op<T>::mid(X)</a>
  CVar<T>& set
    ( CModel<T>*CM, const unsigned int ix, const T&X )
    { *this = CVar( CM, ix, X ); return *this; }

  //! @brief Set Chebyshev model environment in Chebyshev variable to <tt>env</tt>
  CVar<T>& set
    ( CModel<T>*env )
    { *this = CVar( env ); return *this; }

  //! @brief Set multivariate polynomial coefficients in Chebyshev variable to <tt>coefmon</tt>
  CVar<T>& set
    ( const double*coefmon )
    { for( unsigned int imon=0; imon<(_CM?_nmon():1); imon++ )
        _coefmon[imon] = coefmon[imon];
      _update_bndord(); return *this; }

  //! @brief Set multivariate polynomial coefficients in Chebyshev variable to <tt>coefmon</tt> - only first <tt>coefmon.first</tt> coefficients are set
  CVar<T>& set
    ( std::pair<unsigned int, const double*> coefmon )
    { for( unsigned int imon=0; imon<(_CM?_nmon():1); imon++ )
        _coefmon[imon] = ( imon<coefmon.first && coefmon.second?
	                   coefmon.second[imon]: 0. );
      _update_bndord(); return *this; }

  //! @brief Set remainder term in Chebyshev variable to <tt>bndrem</tt>
  CVar<T>& set
    ( const T&bndrem )
    { *_bndrem = bndrem; _update_bndord(); return *this; }

  //! @brief Set multivariate polynomial coefficients and remainder term to those of Chebyshev variable <tt>CV</tt>, possibly defined in another Chebyshev model environment with less variables and a lower expansion order. Coefficients involving other variables or higher order are initialized to 0 if <tt>reset=true</tt> (default), otherwise they are left unmodified.
  CVar<T>& set
    ( CVar<T>& CV, const bool reset=true )
    { if( !_CM || ( CV._CM && ( _nvar() < CV._nvar() || _nord() < CV._nord() ) ) )
        return *this;
      // Reset coefficients to 0
      for( unsigned int imon=0; reset && imon<_nmon(); imon++ )
        _coefmon[imon] = 0.;
      // Copy coefficients from CV --> *this
      _coefmon[0] = CV._coefmon[0];
      if( !CV._CM ){ *_bndrem = *CV._bndrem; _update_bndord(); return *this; }
      unsigned int*iexp = new unsigned int[_nvar()];
      for( unsigned int jmon=1; jmon<CV._nmon(); jmon++ ){
        for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
	  iexp[ivar] = (ivar<CV._nvar()? CV._CM->_expmon[jmon*CV._nvar()+ivar]: 0);
        _coefmon[_CM->_loc_expmon(iexp)] = CV._coefmon[jmon];
      }
      delete[] iexp;
      *_bndrem = *CV._bndrem; _update_bndord(); return *this; }

  //! @brief Copy multivariate polynomial coefficients from current Chebyshev variable into Chebyshev variable <tt>CV</tt>, possibly defined in another Chebyshev model environment with less variables and a lower expansion order. Copied coefficients are reset to 0 in current Chebyshev variable if <tt>reset=true</tt>, otherwise they are left unmodified (default).
  CVar<T>& get
    ( CVar<T>&CV, const bool reset=false )
    { if( !_CM ){
        CV._coefmon[0] = _coefmon[0];
        // Reset coefficients to 0
	if( reset ) _coefmon[0] = 0.;
        return *this;
      }
      if( !CV._CM || _nvar() < CV._nvar() || _nord() < CV._nord() )
        return *this;
      // Copy coefficients from *this --> CV
      unsigned int*iexp = new unsigned int[_nvar()];
      for( unsigned int jmon=0; jmon<(CV._CM?CV._nmon():1); jmon++ ){
        for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
	  iexp[ivar] = (ivar<CV._nvar()? CV._CM->_expmon[jmon*CV._nvar()+ivar]: 0);
        CV._coefmon[jmon] = _coefmon[_CM->_loc_expmon(iexp)];
        // Reset coefficients to 0
	if( reset ) _coefmon[_CM->_loc_expmon(iexp)] = 0.;
      }
      delete[] iexp;
      if( reset ) _update_bndord();
      *CV._bndrem = 0.; CV._update_bndord(); return *this; }

  //! @brief Get pointer to Chebyshev model environment
  CModel<T>* env() const
    { return _CM; }

  //! @brief Compute bound on all terms of (total) order <tt>iord</tt> in Chebyshev variable
  T bound
    ( const unsigned int iord ) const
    { return (!iord || (_CM && iord<=_nord()))? _bndord[iord]: 0.; }

  //! @brief Compute bound on Chebyshev variable using specified bounder <a>type</a>
  T bound
    ( typename CModel<T>::Options::BOUNDER type ) const
    { return _polybnd(type) + *_bndrem; }

  //! @brief Compute bound on Chebyshev variable using bounder type in mc::CModel::Options::BOUNDER_TYPE)
  T bound() const
    { return bound(_CM?_CM->options.BOUNDER_TYPE:CModel<T>::Options::LSB); }

  //! @brief Return remainder term of Chebyshev variable
  T remainder() const
    { return( *_bndrem ); }

  //! @brief Center remainder term of Chebyshev variable
  CVar<T>& center()
    { _center_CM(); return *this; }

  //! @brief Return new Chebyshev variable with same multivariate polynomial part and zero remainder
  CVar<T> polynomial() const
  { CVar<T> CV = *this; *(CV._bndrem) = 0.; return CV; }

  //! @brief Evaluate polynomial part at <tt>x</tt>
  double polynomial
    ( const double*x ) const;

  //! @brief Shortcut to mc::CVar::bound
  T B
    ( const unsigned int iord ) const
    { return bound(iord); }

  //! @brief Shortcut to mc::CVar::bound
  T B() const
    { return bound(); }

  //! @brief Shortcut to mc::CVar::bound
  T B
    ( typename CModel<T>::Options::BOUNDER type ) const
    { return bound(type); }

  //! @brief Shortcut to mc::CVar::remainder
  T R() const
    { return remainder(); }

  //! @brief Shortcut to mc::CVar::center
  CVar<T>& C()
    { return center(); }

  //! @brief Shortcut to mc::CVar::polynomial
  CVar<T> P() const
    { return polynomial(); }

  //! @brief Shortcut to mc::CVar::polynomial
  double P
    ( const double*x ) const
    { return polynomial( x ); }

  //! @brief Get (possibly scaled) coefficient in monomial term with variable exponents as given in <a>iexp</a>
  double coefmon
    ( const unsigned int*iexp ) const;

  //! @brief Get pair of size of, and const pointer to, array of (possibly scaled) monomial coefficients in multivariate polynomial of Chebyshev variable
  std::pair<unsigned int, const double*> coefmon() const;

  //! @brief Get pair of size of, and const pointer to, array of monomial exponents in multivariate polynomial of Chebyshev variable
  std::pair<unsigned int, const unsigned int*> expmon() const;
  /** @} */

  CVar<T>& operator =
    ( const double );
  CVar<T>& operator =
    ( const CVar<T>& );
  CVar<T>& operator =
    ( const T& );
  template <typename U> CVar<T>& operator +=
    ( const CVar<U>& );
  template <typename U> CVar<T>& operator +=
    ( const U& );
  CVar<T>& operator +=
    ( const double );
  template <typename U> CVar<T>& operator -=
    ( const CVar<U>& );
  template <typename U> CVar<T>& operator -=
    ( const U& );
  CVar<T>& operator -=
    ( const double );
  CVar<T>& operator *=
    ( const CVar<T>& );
  CVar<T>& operator *=
    ( const double );
  CVar<T>& operator *=
    ( const T& );
  CVar<T>& operator /=
    ( const CVar<T>& );
  CVar<T>& operator /=
    ( const double );

private:

  //! @brief Private constructor for real scalar in Chebyshev model environment <tt>CM</tt>
  CVar
    ( CModel<T>*CM, const double d=0. );
  //! @brief Private constructor for remainder bound in Chebyshev model environment <tt>CM</tt>
  CVar
    ( CModel<T>*CM, const T&B );

  //! @brief Array of size <tt>_nmon</tt> with (possibly scaled) Chebyshev coefficients
  double *_coefmon;
  //! @brief Array of size <tt>_nord+2</tt> with bounds for all terms of degrees <tt>iord=0,...,_nord</tt> as well as the remainder bound at position <tt>_nord+1</tt>
  mutable T * _bndord; 
  //! @brief Pointer to the remainder bound
  mutable T * _bndrem;

  //! @brief Update bounds for all terms of degrees <tt>iord=0,...,_nord</tt> in <tt>_bndord</tt>
  void _update_bndord();
  //! @brief Center remainder error term <tt>_bndrem</tt>
  void _center_CM();

  //! @brief Range bounder - using specified bounder type
  T _polybnd
    ( typename CModel<T>::Options::BOUNDER type ) const;
  //! @brief Range bounder - using bounder type in CModel<T>::options.BOUNDER_TYPE
  T _polybnd() const;
  //! @brief Range bounder - naive approach
  T _polybnd_naive() const;
  //! @brief Range bounder - Lin & Stadtherr approach
  T _polybnd_LSB() const;
  //! @brief Range bounder - eigenvalue decomposition-based approach
  T _polybnd_eigen() const;
  //! @brief Range bounder - Bernstein approach
  T _polybnd_bernstein() const;
  //! @brief Compute Bernstein coefficient for variable with index <tt>jmon</tt> and exponents <tt>jexp</tt>, given coefficients in Chebyshev form <tt>coefmon</tt> and maximum order <tt>maxord</tt>
  double _coef_bernstein
    ( const double*coefmon, const unsigned int*jexp,
      const unsigned int jmon, const unsigned int maxord ) const;

  //! @brief Initialize private members
  void _init();
  //! @brief Reinitialize private members
  void _reinit();
  //! @brief Clean up private members
  void _clean();

};

////////////////////////////////// CModel //////////////////////////////////////

template <typename T> inline void
CModel<T>::_size
( const unsigned int nvar, const unsigned int nord )
{
  if( !nvar ) throw Exceptions( Exceptions::SIZE );

  //_cleanup();
  _nvar = nvar;
  _nord = nord; 
  _binom = new CM_size[(nvar+nord-1)*(nord+1)];
  _binom_size = std::make_pair( nvar+nord-1, nord+1 );
  _set_binom( nord );
  _posord = new unsigned int[nord+2];
  _posord_size = nord;
  _set_posord( nord );
  _nmon = _posord[_nord+1];
  _expmon = new unsigned int[_nmon*nvar];
  _expmon_size = _nmon;
  _set_expmon( nord );
  _set_prodmon();
  _bndvar = new T[_nvar];
  for( unsigned int i=0; i<_nvar; i++ )
    _bndvar[i] = 0.;
  _modvar = true;

  _CV = new CVar<T>( this );
}

template <typename T> inline void
CModel<T>::_set_posord
( const unsigned int nord )
{
  _posord[0] = 0;
  _posord[1] = 1;
  for( unsigned int i=1; i<=nord; i++ ){
    CM_size _posord_next = _posord[i] + _get_binom( _nvar+i-1, i );
    if( _posord_next > UINT_MAX )
      throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::MAXSIZE );
    _posord[i+1] = _posord_next;
  }

#ifdef MC__CMODEL_DEBUG
  mc::display( 1, nord+2, _posord, 1, "_posord", std::cout );
#endif
}
    
template <typename T> inline void
CModel<T>::_ext_posord
( const unsigned int maxord )
{
  if( maxord < _posord_size ) return;
  delete[] _posord;
  _posord = new unsigned int[maxord+2];
  _posord_size = maxord;
  _set_posord( maxord );
}

template <typename T> inline void
CModel<T>::_set_expmon
( const unsigned int nord )
{
  unsigned int *iexp = new unsigned int[_nvar] ;
  for( unsigned int k=0; k<_nvar; k++ ) _expmon[k] = 0;
  for( unsigned int i=1; i<=nord; i++ ){
    for( unsigned int j=0; j<_nvar; j++ ) iexp[j] = 0;
    for( unsigned int j=_posord[i]; j<_posord[i+1]; j++ ){
      _next_expmon( iexp, i );
      for( unsigned int k=0; k<_nvar; k++ )
        _expmon[j*_nvar+k] = iexp[k];
    }
  }
  delete[] iexp;

#ifdef MC__CMODEL_DEBUG
  mc::display( _nvar, _expmon_size, _expmon, _nvar, "_expmon", std::cout );
#endif
}
  
template <typename T> inline void
CModel<T>::_next_expmon
( unsigned int *iexp, const unsigned int iord ) const
{
  unsigned int curord;
  do{
    iexp[_nvar-1] += iord;
    unsigned int j = _nvar;
    while( j > 0 && iexp[j-1] > iord ){
      iexp[j-1] -= iord + 1;
      j-- ;
      iexp[j-1]++;
    }
    curord = 0;
    for( unsigned int i=0; i<_nvar; i++ ) curord += iexp[i];
  } while( curord != iord );
}
    
template <typename T> inline void
CModel<T>::_ext_expmon
( const unsigned int maxord, const bool full )
{
  _ext_binom( maxord );
  _ext_posord( maxord ); 
  delete[] _expmon;
  if( full ) _expmon_size = std::pow(maxord+1,_nvar);
  else       _expmon_size = _posord[maxord+1];
  _expmon = new unsigned int[ _expmon_size*_nvar];
  _set_expmon( maxord );
  if( !full ) return;

  unsigned int *iexp = new unsigned int[_nvar];
  for( unsigned int iord=maxord+1, jmon=_posord[maxord+1];
       jmon<_expmon_size; iord++ ){
    _ext_binom( iord );
    _ext_posord( iord );
    if( _posord[iord] >= _posord[iord+1] )
      throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::MAXSIZE );
    for( unsigned int ivar=0; ivar<_nvar; ivar++ ) iexp[ivar] = 0;
    for( unsigned int kmon=_posord[iord]; kmon<_posord[iord+1]; kmon++ ){
      _next_expmon( iexp, iord );
      bool lt_maxord = true;
      for( unsigned int ivar=0; ivar<_nvar; ivar++ ){
        if( iexp[ivar] > maxord ){ lt_maxord = false; break; }
        _expmon[jmon*_nvar+ivar] = iexp[ivar];
      }
      if( lt_maxord ){
#ifdef MC__CMODEL_DEBUG
	std::cout << jmon << ":";
        for( unsigned int ivar=0; ivar<_nvar; ivar++ )
	  std::cout << "  " << _expmon[jmon*_nvar+ivar];
	std::cout << std::endl;
#endif
        jmon++;
      }
    }
  }
  delete[] iexp;
  return;
}

template <typename T> inline void
CModel<T>::_set_prodmon_exp
( unsigned int&kvar, const unsigned int*iexp, const unsigned int iord,
  unsigned int*iprodexp, unsigned int&nprodmon )
{
  if( kvar == _nvar || !nprodmon ) return;

  int nprodnew = 0;
  for( unsigned int imon=0; imon<nprodmon; imon++ ){
    // Product with Chebyshev basis function T0
    if( !iexp[kvar] || !iprodexp[imon*(_nvar+1)+kvar+1] ){
      iprodexp[imon*(_nvar+1)] += iprodexp[imon*(_nvar+1)+kvar+1] + iexp[kvar];
      if( iprodexp[imon*(_nvar+1)] <= _nord )
        iprodexp[imon*(_nvar+1)+kvar+1] += iexp[kvar];
    }
    // Product with Chebyshev basis function Ti, i>0
    else{
      // Only add term if resulting monic index remain less than or equal to _nord
      if( iprodexp[imon*(_nvar+1)] + iprodexp[imon*(_nvar+1)+kvar+1] + iexp[kvar] <= _nord ){
        for( unsigned int ivar=0; ivar<=_nvar; ivar++ )
          iprodexp[(nprodmon+nprodnew)*(_nvar+1)+ivar] = iprodexp[imon*(_nvar+1)+ivar];
        iprodexp[(nprodmon+nprodnew)*(_nvar+1)] += iprodexp[imon*(_nvar+1)+kvar+1] + iexp[kvar];
        iprodexp[(nprodmon+nprodnew)*(_nvar+1)+kvar+1] += iexp[kvar];
        nprodnew++;
      }
      if( iexp[kvar] <= iprodexp[imon*(_nvar+1)+kvar+1] ){
        iprodexp[imon*(_nvar+1)] += iprodexp[imon*(_nvar+1)+kvar+1] - iexp[kvar];
        if( iprodexp[imon*(_nvar+1)] <= _nord )
          iprodexp[imon*(_nvar+1)+kvar+1] -= iexp[kvar];
      }
      else{
        iprodexp[imon*(_nvar+1)] += iexp[kvar] - iprodexp[imon*(_nvar+1)+kvar+1];      
        if( iprodexp[imon*(_nvar+1)] <= _nord )
          iprodexp[imon*(_nvar+1)+kvar+1] = iexp[kvar] - iprodexp[imon*(_nvar+1)+kvar+1];
      }
    }
  }
  nprodmon += nprodnew;
  return _set_prodmon_exp( ++kvar, iexp, iord, iprodexp, nprodmon );
}

template <typename T> inline void
CModel<T>::_set_prodmon()
{
  _prodmon = new unsigned int**[_nmon];
  size_t prodmon_max = _nvar+1;
  for( unsigned int i=0; i<_nvar && i<_nord; i++ ){  // <-- SIZE COULD BE TOO RESTRICTIVE!!!
    if( prodmon_max > std::numeric_limits<size_t>::max()/2 )
      throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::MAXSIZE );
    prodmon_max *= 2;
  }
  unsigned int *iexp = new unsigned int[prodmon_max];

  // Loop over all monic terms
  for( unsigned int iord=0; iord<=_nord; iord++ ){    
    for( unsigned int imon=_posord[iord]; imon<_posord[iord+1]; imon++ ){

      // Loop over all monics whose indices are lower than or equal to imon
      _prodmon[imon] = new unsigned int*[imon+1];
      for( unsigned int jord=0, jmon=0; jord<=iord; jord++ ){    
        for( ; jmon<_posord[jord+1] && jmon<=imon; jmon++ ){

          // Current monic operand
          iexp[0] = 0; // used to store partial order
          for( unsigned int ivar=0; ivar<_nvar; ivar++ )
            iexp[ivar+1] = _expmon[imon*_nvar+ivar];
#ifdef MC__CMODEL_DEBUG
          mc::display( 1, _nvar, iexp+1, 1, "iexp:", std::cout );
          mc::display( 1, _nvar, _expmon+jmon*_nvar, 1, "jexp:", std::cout );
#endif
          // Construct product monics recursively and only keep those of order <= _nord
          unsigned int nprodmax = 1, nprodmon = 1, kvar = 0;
          for( unsigned int ivar=0; ivar<_nvar; ivar++ )
            if( iexp[ivar+1] && _expmon[jmon*_nvar+ivar] ) nprodmax*=2;
          _set_prodmon_exp( kvar, _expmon+jmon*_nvar, jord, iexp, nprodmon );
          _prodmon[imon][jmon] = new unsigned int[nprodmon+2];
          _prodmon[imon][jmon][0] = nprodmax;
          _prodmon[imon][jmon][1] = 0;
          for( unsigned int kmon=0; kmon<nprodmon; kmon++ ){
#ifdef MC__CMODEL_DEBUG
            mc::display( 1, _nvar+1, iexp+kmon*(_nvar+1), 1, "ijexp:", std::cout );
#endif
            // Eliminate monic terms with order > _nord
            if( iexp[kmon*(_nvar+1)] > _nord ) continue;
            _prodmon[imon][jmon][++_prodmon[imon][jmon][1]+1] = _loc_expmon( iexp+kmon*(_nvar+1)+1 );
          }
#ifdef MC__CMODEL_DEBUG
          std::ostringstream oheadi;
          oheadi << "_prodmon[" << imon << "," << jmon << "]:  max=" << nprodmax
                 << "  cur=" << _prodmon[imon][jmon][1];
          mc::display( 1, _prodmon[imon][jmon][1], _prodmon[imon][jmon]+2, 1, oheadi.str(), std::cout );
#endif
        }
      }
    }
  }

  delete[] iexp;
}
    
template <typename T> inline unsigned int
CModel<T>::_loc_expmon
( const unsigned int *iexp ) const
{
  unsigned int ord = 0;
  for( unsigned int i=0; i<_nvar; i++ ) ord += iexp[i];
  assert( ord<_nord+2 );
  unsigned int pos = _posord[ord];
  
  unsigned int p = _nvar ; 
  for( unsigned int i=0; i<_nvar-1; i++ ){
    p--;
    for( unsigned int j=0; j<iexp[i]; j++ )
      pos += _get_binom( p-1+ord-j, ord-j );
    ord -= iexp[i];
  }

  return pos;    
}

template <typename T> inline void
CModel<T>::_set_bndvar
( const unsigned int ivar, const T&X )
{
  _bndvar[ivar] = X;
}

template <typename T> inline void
CModel<T>::_set_binom
( const unsigned int nord )
{
  CM_size *p;
  unsigned int k;
  for( unsigned int i=0; i<_nvar+nord-1; i++ ){
    p = &_binom[i*(nord+1)];
    *p = 1;
    p++;
    *p = i+1;
    p++;
    k = ( i+1<nord? i+1: nord );
    for( unsigned int j=2; j<=k; j++, p++ ) *p = *(p-1) * (i+2-j)/j;
    for( unsigned int j=k+1; j<=nord; j++, p++ ) *p = 0.;
  }
#ifdef MC__CMODEL_DEBUG
  mc::display( _binom_size.second, _binom_size.first, _binom,
    _binom_size.second, "_binom", std::cout );
#endif
}
    
template <typename T> inline void
CModel<T>::_ext_binom
( const unsigned int maxord )
{
  if( maxord < _binom_size.second ) return;
  delete[] _binom;
  _binom = new CM_size[(_nvar+maxord-1)*(maxord+1)];
  _binom_size = std::make_pair( _nvar+maxord-1, maxord+1 );
  _set_binom( maxord );
}

template <typename T> inline CM_size
CModel<T>::_get_binom
( const unsigned int n, const unsigned int k ) const
{
#ifdef MC__CMODEL_CHECK
  assert( n<=_binom_size.first );
  assert( k<=_binom_size.second );
  assert( k<=n );
#endif
  return( n? _binom[(n-1)*_binom_size.second+k]: 1. );
}

template <typename T> inline void
CModel<T>::_reset()
{
  for( unsigned int i=0; i<_nvar; i++ ){
    _bndvar[i] = 0;
  }
}

template <typename T> inline void
CModel<T>::_cleanup()
{
  for( unsigned int i=0; i<_nmon; i++ ){
    for( unsigned int j=0; j<=i; j++ )
      delete[] _prodmon[i][j];
    delete[] _prodmon[i];
  }
  delete[] _prodmon;
  delete[] _expmon;
  delete[] _posord;
  delete[] _bndvar;
  delete[] _binom;
  delete _CV;
}

////////////////////////////////// CVar ///////////////////////////////////////

template <typename T> inline void
CVar<T>::_init()
{
  if( !_CM ){
    _coefmon = new double[1];
    _bndord  = new T[1];
    _bndrem  = _bndord;
    return;
  }
  _coefmon = new double[_nmon()];
  _bndord  = new T[_nord()+2];
  _bndrem  = _bndord + _nord()+1;
}

template <typename T> inline void
CVar<T>::_clean()
{
  delete [] _coefmon; delete [] _bndord;
  _coefmon = 0; _bndord = _bndrem = 0;
}

template <typename T> inline void
CVar<T>::_reinit()
{
  _clean(); _init();
}

template <typename T> inline
CVar<T>::CVar
( const double d )
: _CM( 0 )
{
  _init();
  _coefmon[0] = d;
  _bndord[0] = 0.;
}

template <typename T> inline CVar<T>&
CVar<T>::operator =
( const double d )
{
  if( _CM ){ _CM = 0; _reinit(); }
  _coefmon[0] = d;
  _bndord[0] = 0.;
  return *this;
}

template <typename T> inline
CVar<T>::CVar
( CModel<T>*CM, const double d )
: _CM( CM )
{
  if( !_CM ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INIT );
  _init();
  _coefmon[0] = d;
  for( unsigned int i=1; i<_nmon(); i++ ) _coefmon[i] = 0.;
  _bndord[0] = d;
  for( unsigned int i=1; i<_nord()+2; i++) _bndord[i] = 0.;
}

template <typename T> inline
CVar<T>::CVar
( const T&B )
: _CM( 0 )
{
  _init();
  _coefmon[0] = 0.;
  _bndord[0] = B;
}

template <typename T> inline CVar<T>&
CVar<T>::operator =
( const T&B )
{
  if( _CM ){ _CM = 0; _reinit(); }
  _coefmon[0] = 0.;
  _bndord[0] = B;
  return *this;
}

template <typename T> inline
CVar<T>::CVar
( CModel<T>*CM, const T&B )
: _CM( CM )
{
  if( !_CM ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INIT );
  _init();
  for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] = 0.;
  for( unsigned int i=0; i<_nord()+1; i++) _bndord[i] = 0.;
  *_bndrem = B;
  if( _CM->options.CENTER_REMAINDER ) _center_CM();
}

template <typename T> inline
CVar<T>::CVar
( const CVar<T>&CV )
: _CM(0)
{
  _init();
  *this = CV;
}

template <typename T> inline CVar<T>&
CVar<T>::operator =
( const CVar<T>&CV )
{
  // Same CVar
  if( this == &CV ) return *this;

  // Reinitialization needed?
  if( _CM != CV._CM ){ _CM = CV._CM; _reinit(); }

  // Set to CVar not linked to CModel (either scalar or range)
  if( !_CM ){
    _coefmon[0] = CV._coefmon[0];
    _bndord[0] = CV._bndord[0];
    return *this; 
  }
  // Set to CVar linked to CModel
  for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] = CV._coefmon[i];
  for( unsigned int i=0; i<_nord()+2; i++) _bndord[i] = CV._bndord[i];
  return *this;
}

template <typename T> template <typename U> inline
CVar<T>::CVar
( CModel<T>*&CM, const CVar<U>&CV )
: _CM(CM), _coefmon(0), _bndord(0), _bndrem(0)
{
  _init();
  CVar<U> CVtrunc( CV );
  _coefmon[0] = CVtrunc._coefmon[0];
  CVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _CM && i<_nmon(); i++ ){
    if( CVtrunc._CM && i < CVtrunc._nmon() ){
      _coefmon[i] = CVtrunc._coefmon[i];
      CVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  CVtrunc._update_bndord();
  *_bndrem = T( CVtrunc.B() );
  if( _CM ) _update_bndord();
  return;
}

template <typename T> template <typename U> inline
CVar<T>::CVar
( CModel<T>*&CM, const CVar<U>&CV, T (*method)( const U& ) )
: _CM(CM), _coefmon(0), _bndord(0), _bndrem( 0 )
{
  if( !method ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INIT );
  _init();
  CVar<U> CVtrunc( CV );
  _coefmon[0] = CVtrunc._coefmon[0];
  CVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _CM && i<_nmon(); i++ ){
    if( CVtrunc._CM && i < CVtrunc._nmon() ){
      _coefmon[i] = CVtrunc._coefmon[i];
      CVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  CVtrunc._update_bndord();
  *_bndrem = (*method)( CVtrunc.B() );
  if( _CM ) _update_bndord();
  return;
}

template <typename T> template <typename U> inline
CVar<T>::CVar
( CModel<T>*&CM, const CVar<U>&CV, const T& (U::*method)() const )
: _CM(CM), _coefmon(0), _bndord(0), _bndrem( 0 )
{
  if( !method ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INIT );
  _init();
  CVar<U> CVtrunc( CV );
  _coefmon[0] = CVtrunc._coefmon[0];
  CVtrunc._coefmon[0] = 0. ;
  for( unsigned int i=1; _CM && i<_nmon(); i++ ){
    if( CVtrunc._CM && i < CVtrunc._nmon() ){
      _coefmon[i] = CVtrunc._coefmon[i];
      CVtrunc._coefmon[i] = 0.;
    }
    else
      _coefmon[i] = 0.;
  }
  CVtrunc._update_bndord();
  *_bndrem = (CVtrunc.B().*method)();
  if( _CM ) _update_bndord();
  return;
}

template <typename T> inline
CVar<T>::CVar
( CModel<T>*CM, const unsigned int ivar, const T&X )
: _CM( CM )
{
  if( !CM ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INIT );

  // Keep track of variable bounds in CModel
  _CM->_set_bndvar( ivar, X );
  _init();

  // Populate _coefmon w/ CVar coefficients
  _coefmon[0] = Op<T>::mid(X);
  for( unsigned int i=1; i<_nmon(); i++ ) _coefmon[i] = 0.;
  if( _nord() > 0 ) _coefmon[_nvar()-ivar] = Op<T>::diam(X)/2.;

  // Populate _bndord w/ bounds on CVar terms
  _bndord[0] = _coefmon[0];
  _bndord[1] = X-_coefmon[0];
  for( unsigned int i=2; i<_nord()+2; i++) _bndord[i] = 0.;
}

template <typename T> inline void
CVar<T>::_update_bndord()
{
  if( !_CM ) return;
  _bndord[0] = _coefmon[0];
  for( unsigned int i=1; i<=_nord(); i++ ){
    _bndord[i] = 0.; 
    for( unsigned int j=_posord(i); j<_posord(i+1); j++ )
      _bndord[i] += _coefmon[j] * T(-1.,1.);
  }
}

template <typename T> inline void
CVar<T>::_center_CM()
{
  const double remmid = Op<T>::mid(*_bndrem);
  _coefmon[0] += remmid;
  if( _CM ) _bndord[0] = _coefmon[0];
  *_bndrem -= remmid;
}

template <typename T> inline T
CVar<T>::_polybnd_eigen() const
{
  static const double TOL = 1e-8;

  T bndpol = _coefmon[0];
  if( _nord() == 1 ) bndpol += _bndord[1];

  else if( _nord() > 1 ){
    T bndmon = Op<T>::zeroone()*2.-1.;
    double*U = new double[_nvar()*_nvar()];
    for( unsigned int i=_posord(2); i<_posord(3); i++ ){
      unsigned int i1=0, i2=_nvar();
      const unsigned int*iexp=_expmon(i);
      for( ; i1<_nvar(); i1++ ) if( iexp[i1] ) break;
      if( iexp[i1] == 2 ){
        U[_nvar()*i1+i1] = 2.*_coefmon[i];
        bndpol -= _coefmon[i];
        continue;
      }
      for( i2=i1+1; i2<_nvar(); i2++ ) if( iexp[i2] ) break;
      U[_nvar()*i1+i2] = 0.;
      U[_nvar()*i2+i1] = _coefmon[i]/2.;
    }
#ifdef MC__CVAR_DEBUG_EIGEN
    display( _nvar(), _nvar(), U, _nvar(), "Matrix U", std::cout );
#endif
    double*D = mc::dsyev_wrapper( _nvar(), U, true );
    if( !D ){
      delete[] U;
      throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::EIGEN );
    }

    for( unsigned int i=0; i<_nvar(); i++ ){
      double linaux = 0.;
      T bndaux(0.);
      for( unsigned int k=0; k<_nvar(); k++ ){
        linaux += U[i*_nvar()+k] * _coefmon[_nvar()-k];
        bndaux += U[i*_nvar()+k] * bndmon;
      }
#ifdef MC__CVAR_DEBUG_EIGEN
      std::cout << i << ": LINAUX = " << linaux
                << "  BNDAUX = " << bndaux << std::endl;
#endif
      if( std::fabs(D[i]) > TOL )
        bndpol += D[i] * Op<T>::sqr( linaux/D[i]/2. + bndaux )
                - linaux*linaux/D[i]/4.;
      else     
        bndpol += linaux * bndaux + D[i] * Op<T>::sqr( bndaux );
#ifdef MC__CVAR_DEBUG_EIGEN
        std::cout << "BNDPOL: " << bndpol << std::endl;
#endif
    }
    delete[] U;
    delete[] D;
  }
#ifdef MC__CVAR_DEBUG_EIGEN
  int tmp; std::cin >> tmp;
#endif

  for( unsigned int i=3; i<=_nord(); i++ ) bndpol += _bndord[i];
  return bndpol;
}

template <typename T> inline T
CVar<T>::_polybnd_bernstein() const
{
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );
/*
  // Scale Chebyshev model in unit hypercube with reference point at the origin
  double *coeftrans = new double[_nmon()];
  for( unsigned int imon=0; imon<_nmon(); imon++ )
    coeftrans[imon] = _coefmon[imon]; 
  for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
    _CM->_scale( ivar, coeftrans );
  
  // Expand binomial coefficient and exponent arrays if needed
  const unsigned int maxord = (_CM->options.BOUNDER_ORDER>_nord()? 
    _CM->options.BOUNDER_ORDER: _nord() );
  const unsigned int maxmon = std::pow(maxord+1,_nvar());
  _CM->_ext_expmon( maxord, true );

  // Compute min/max amongst all Bernstein coefficients
  bndmod = coeftrans[0];
#ifdef  MC__CVAR_DEBUG_BERSTEIN
  std::cout << "\n0:  " << bndmod << std::endl;
#endif
  for( unsigned int jmon=1; jmon<maxmon; jmon++ ){
    const unsigned int*jexp = _CM->_expmon + jmon*_nvar();
    const double coefbern = _coef_bernstein( coeftrans, jexp, jmon, maxord );
    bndmod = Op<T>::hull( bndmod, coefbern );
#ifdef  MC__CVAR_DEBUG_BERNSTEIN
    std::cout << jmon << " ["; 
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      std::cout << std::setw(3) << jexp[ivar];
    std::cout << "] : " << coefbern << " : " << bndmod << std::endl;
#endif
  }

  delete[] coeftrans;
  bndmod += *_bndrem;
*/
  return T(0.);
}

/*
template <typename T> inline double
CVar<T>::_coef_bernstein
( const double*coefmon, const unsigned int*jexp,
  const unsigned int jmon, const unsigned int maxord ) const
{
  // Compute bernstein coefficient with variables indices <tt>jexp</tt>
  double coefbern = coefmon[0];
  for( unsigned int imon=1; imon<=std::min(jmon,_nmon()-1); imon++ ){
    const unsigned int*iexp = _CM->_expmon + imon*_nvar();
    // Only append term if monomial degrees are lower
    bool inrange = true;
    for( unsigned int ivar=0; ivar<_nvar() && inrange; ivar++ )
      if( iexp[ivar] > jexp[ivar] ) inrange = false;
    if( !inrange ) continue;
    double termbern = coefmon[imon];
#ifdef  MC__CVAR_DEBUG_BERSTEIN
    std::cout << "  ["; 
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      std::cout << std::setw(3) << iexp[ivar];
    std::cout << "]";
#endif
    for( unsigned int ivar=0; ivar<_nvar(); ivar++ )
      termbern *= (double)_CM->_get_binom(jexp[ivar],iexp[ivar])
                / (double)_CM->_get_binom(maxord,iexp[ivar]);
    coefbern += termbern;
  }
#ifdef  MC__CVAR_DEBUG_BERSTEIN
  std::cout << std::endl;
#endif
  return coefbern;
}
*/

template <typename T> inline T
CVar<T>::_polybnd_LSB() const
{
  static const double TOL = 1e-8;

  T bndpol = _coefmon[0];
  if( _nord() == 1 ) bndpol += _bndord[1];

  else if( _nord() > 1 ){
    T bndmon = Op<T>::zeroone()*2.-1.;
    for( unsigned int i=_posord(2); i<_posord(3); i++ ){
      unsigned int k=0;
      const unsigned int*iexp=_expmon(i);
      for( unsigned int j=0; j<_nvar(); j++ ){
        k = ( iexp[j]==1? _nvar(): j );
        if( iexp[j] ) break;
      }
      // off-diagonal quadratic terms
      if( k == _nvar() ){
        bndpol += _coefmon[i] * bndmon;
        continue;
      }
      // linear and diagonal quadratic terms
      const double ai  = _coefmon[_nvar()-k], aii = _coefmon[i];
      if( std::fabs(aii) > TOL )
        bndpol += (2.*aii)*Op<T>::sqr(bndmon+ai/(aii*4.))-aii-ai*ai/8./aii;
      else
        bndpol += ai*bndmon+aii*bndmon;
    }
  }
  // higher-order terms
  for( unsigned int i=3; i<=_nord(); i++ ) bndpol += _bndord[i];
  return bndpol;
}

template <typename T> inline T
CVar<T>::_polybnd_naive() const
{
  T bndpol = _coefmon[0];
  for( unsigned int i=1; i<=_nord()+1; i++ ) bndpol += _bndord[i];
  return bndpol;
}

template <typename T> inline T
CVar<T>::_polybnd
( typename CModel<T>::Options::BOUNDER type ) const
{
  if( !_CM ) return _coefmon[0];

  switch( type ){
  case CModel<T>::Options::NAIVE:
    return _polybnd_naive();
  case CModel<T>::Options::BERNSTEIN:
    return _polybnd_bernstein();
  case CModel<T>::Options::EIGEN:
    return _polybnd_eigen();
  case CModel<T>::Options::HYBRID:{
    T bndpol;
    if( !Op<T>::inter( bndpol, _polybnd_LSB(), _polybnd_eigen() ) )
      throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::INCON );
    return bndpol;
   }
  case CModel<T>::Options::LSB: default:
    return _polybnd_LSB();
  }
}

template <typename T> inline T
CVar<T>::_polybnd() const
{
  return _polybnd( _CM->options.BOUNDER_TYPE );
}

template <typename T> inline double
CVar<T>::polynomial
( const double*x ) const
{
  if( !_CM ) return _coefmon[0];

  static const double TOL = machprec()*1e2;
  double xs[_nvar()];
  for( unsigned k=0; k<_nvar(); k++ ){
    const double dXk = Op<T>::diam(_CM->_bndvar[k]);
    xs[k] = ( dXk>TOL? 2.*(x[k]-Op<T>::l(_CM->_bndvar[k]))/dXk-1.
                     : Op<T>::mid(_CM->_bndvar[k]) );
  }
  double Pval = _coefmon[0];
  for( unsigned int i=1; i<_nmon(); i++ ){
    const unsigned int *iexp = _expmon(i);
    double valmon = 1.;
    for( unsigned int k=0; k<_nvar(); k++ ){
      if( !iexp[k] ) continue;
      if( iexp[k] == 1 ) valmon *= xs[k];
      else if( iexp[k] == 2 ) valmon *= 2*xs[k]*xs[k]-1;
      else valmon *= std::cos(iexp[k]*std::acos(xs[k]));
    }
    Pval += _coefmon[i] * valmon;
  }
  return Pval;
}

template <typename T> inline double
CVar<T>::coefmon
( const unsigned int*iexp ) const
{
  if( !_CM ) return 0;
  const unsigned int imon = _CM->_loc_expmon( iexp );
  return( imon<_nmon()? _coefmon[imon]: 0. );
}

template <typename T> inline std::pair<unsigned int, const double*>
CVar<T>::coefmon() const
{
  return std::make_pair( (_CM?_nmon():1), _coefmon );
}

template <typename T> inline std::pair<unsigned int, const unsigned int*>
CVar<T>::expmon() const
{
  return std::make_pair( (_CM?_nmon()*_nvar():1), _CM->_expmon );
}

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const CVar<T>&CV )
{
  out << std::endl
      << std::scientific << std::setprecision(5)
      << std::right;

  // Constant model
  if( !CV._CM ){
    out << "   a0    = " << std::right << std::setw(12) << CV._coefmon[0]
        << "      0  0"
        << std::endl
        << "   R     = " << *(CV._bndrem) << std::endl;
  }

  // Chebyshev coefficients and corresponding exponents
  else{
    out << std::setprecision(CV._CM->options.DISPLAY_DIGITS);
    for( unsigned int i=0; i<CV._nmon(); i++ ){
      out << "   a" << std::left << std::setw(4) << i << " = "
          << std::right << std::setw(CV._CM->options.DISPLAY_DIGITS+7)
	  << CV._coefmon[i] << "   ";
      for( unsigned int k=0; k<CV._nvar(); k++ )
        out << std::setw(3) << CV._expmon(i)[k];
      out << std::endl;
    }
    // Remainder term
    out << std::right << "   R     =  " << *(CV._bndrem)
        << std::endl;
  }

  // Range bounder
  out << std::right << "   B     =  " << CV.B()
      << std::endl;

  return out;
}

template <typename T> inline CVar<T>
operator +
( const CVar<T>&CV )
{
  return CV;
}

template <typename T> template <typename U> inline CVar<T>&
CVar<T>::operator +=
( const CVar<U>&CV )
{
  if( !CV._CM ){
    _coefmon[0] += CV._coefmon[0];
    *_bndrem += *(CV._bndrem);
  }
  else if( !_CM ){
    CVar<T> CV2(*this);
    *this = CV; *this += CV2;
  }
  else{
    if( _CM != CV._CM ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::CMODEL );
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] += CV._coefmon[i];
    *_bndrem += *(CV._bndrem);
    _update_bndord();
  }
  if( _CM && _CM->options.CENTER_REMAINDER ) _center_CM();
  return *this;
}

template <typename T, typename U> inline CVar<T>
operator +
( const CVar<T>&CV1, const CVar<U>&CV2 )
{
  CVar<T> CV3( CV1 );
  CV3 += CV2;
  return CV3;
}

template <typename T> inline CVar<T>&
CVar<T>::operator +=
( const double c )
{
  _coefmon[0] += c;
  if( _CM ) _bndord[0] = _coefmon[0];
  return *this;
}

template <typename T> inline CVar<T>
operator +
( const CVar<T>&CV1, const double c )
{
  CVar<T> CV3( CV1 );
  CV3 += c;
  return CV3;
}

template <typename T> inline CVar<T>
operator +
( const double c, const CVar<T>&CV2 )
{
  CVar<T> CV3( CV2 );
  CV3 += c;
  return CV3;
}

template <typename T> template <typename U> inline CVar<T>&
CVar<T>::operator +=
( const U&I )
{
  *_bndrem += I;
  if( _CM && _CM->options.CENTER_REMAINDER ) _center_CM();
  return *this;
}

template <typename T, typename U> inline CVar<T>
operator +
( const CVar<T>&CV1, const U&I )
{
  CVar<T> CV3( CV1 );
  CV3 += I;
  return CV3;
}

template <typename T, typename U> inline CVar<T>
operator +
( const U&I, const CVar<T>&CV2 )
{
  CVar<T> CV3( CV2 );
  CV2 += I;
  return CV3;
}

template <typename T> inline CVar<T>
operator -
( const CVar<T>&CV )
{
  if( !CV._CM ){
    CVar<T> CV2;
    CV2._coefmon[0] = -CV._coefmon[0];
    CV2._bndord[0] = -CV._bndord[0];
    return CV2;
  }
  CVar<T>& CV2 = *CV._CV();
  for( unsigned int i=0; i<CV._nmon(); i++ ) CV2._coefmon[i] = -CV._coefmon[i];
  for( unsigned int i=0; i<CV._nord()+2; i++ ) CV2._bndord[i] = -CV._bndord[i];
  if( CV._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> template <typename U> inline CVar<T>&
CVar<T>::operator -=
( const CVar<U>&CV )
{
  if( !CV._CM ){
    _coefmon[0] -= CV._coefmon[0];
    *_bndrem -= *(CV._bndrem);
  }
  else if( !_CM ){
    CVar<T> CV2(*this);
    *this = -CV; *this += CV2;
  }
  else{
    if( _CM != CV._CM ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::CMODEL );
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] -= CV._coefmon[i];
    *_bndrem -= *(CV._bndrem);
    _update_bndord();
  }
  if( _CM && _CM->options.CENTER_REMAINDER ) _center_CM();
  return *this;
}

template <typename T, typename U> inline CVar<T>
operator-
( const CVar<T>&CV1, const CVar<U>&CV2 )
{
  CVar<T> CV3( CV1 );
  CV3 -= CV2;
  return CV3;
}

template <typename T> inline CVar<T>&
CVar<T>::operator -=
( const double c )
{
  _coefmon[0] -= c;
  if( _CM ) _bndord[0] = _coefmon[0];
  return *this;
}

template <typename T> inline CVar<T>
operator -
( const CVar<T>&CV1, const double c )
{
  CVar<T> CV3( CV1 );
  CV3 -= c;
  return CV3;
}

template <typename T> inline CVar<T>
operator -
( const double c, const CVar<T>&CV2 )
{
  CVar<T> CV3( -CV2 );
  CV3 += c;
  return CV3;
}

template <typename T> template <typename U> inline CVar<T>&
CVar<T>::operator -=
( const U&I )
{
  *_bndrem -= I;
  if( _CM && _CM->options.CENTER_REMAINDER ) _center_CM();
  return *this;
}

template <typename T, typename U> inline CVar<T>
operator -
( const CVar<T>&CV1, const U&I )
{
  CVar<T> CV3( CV1 );
  CV3 -= I;
  return CV3;
}

template <typename T, typename U> inline CVar<T>
operator -
( const U&I, const CVar<T>&CV2 )
{
  CVar<T> CV3( -CV2 );
  CV3 += I;
  return CV3;
}

template <typename T> inline CVar<T>&
CVar<T>::operator *=
( const CVar<T>&CV )
{
   CVar<T> CV2( *this );
   *this = CV * CV2;
   return *this;
}

template <typename T> inline CVar<T>
operator *
( const CVar<T>&CV1, const CVar<T>&CV2 )
{
  if( !CV2._CM )      return( CV1 * CV2._coefmon[0] + CV1 * *(CV2._bndrem) );
  else if( !CV1._CM ) return( CV2 * CV1._coefmon[0] + CV2 * *(CV1._bndrem) );
  else if( &CV1 == &CV2 ) return sqr(CV1);

  if( CV1._CM != CV2._CM )
    throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::CMODEL );
  CVar<T>& CV3 = *CV1._CV();
  for( unsigned int i=0; i<CV3._nmon(); i++ ) CV3._coefmon[i] = 0.;
  *(CV3._bndrem) = 0.;
  T bndmon = Op<T>::zeroone()*2.-1.;

  // Populate _coefmon and _bndrem for product of polynomial parts
  for( unsigned int i=0; i<CV3._nmon(); i++ ){
    for( unsigned int j=0; j<i; j++ ){
      const unsigned int*prodij = CV3._prodmon(i,j);
      for( unsigned int k=0; k<prodij[1]; k++ )
        CV3._coefmon[prodij[2+k]] += ( CV1._coefmon[i] * CV2._coefmon[j]
          + CV1._coefmon[j] * CV2._coefmon[i] ) / (double)prodij[0];
      *(CV3._bndrem) += bndmon * ( ( CV1._coefmon[i] * CV2._coefmon[j]
          + CV1._coefmon[j] * CV2._coefmon[i] )
          * ( 1. - (double)prodij[1] / (double)prodij[0] ) );
    }
    const unsigned int*prodii = CV3._prodmon(i,i);
    for( unsigned int k=0; k<prodii[1]; k++ )
      CV3._coefmon[prodii[2+k]] += CV1._coefmon[i] * CV2._coefmon[i]
        / (double)prodii[0];
    *(CV3._bndrem) += bndmon * ( CV1._coefmon[i] * CV2._coefmon[i]
      * ( 1. - (double)prodii[1] / (double)prodii[0] ) );
  }

  // Populate _coefmon and _bndrem for product of polynomial and remainder parts
  T P1 = CV1._polybnd(), P2 = CV2._polybnd();
  //T P1 = 0., P2 = 0.;
  //for( unsigned int i=0; i<=CV3._nord(); i++ ){
  //  P1 += CV1._bndord[i];
  //  P2 += CV2._bndord[i];
  //}
  T R1 = *(CV3._bndrem) + ( P1 + *(CV1._bndrem) ) * *(CV2._bndrem) + P2 * *(CV1._bndrem);
  T R2 = *(CV3._bndrem) + P1 * *(CV2._bndrem) + ( P2 + *(CV2._bndrem) ) * *(CV1._bndrem);
  if( !Op<T>::inter( *(CV3._bndrem), R1, R2) ) *(CV3._bndrem) = R1;

  // Update _bndord for product term (except remainder term)
  CV3._update_bndord();
  if( CV3._CM->options.CENTER_REMAINDER ) CV3._center_CM();
  return CV3;
}

template <typename T> inline CVar<T>
sqr
( const CVar<T>&CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] *= CV2._coefmon[0];
    *(CV2._bndrem) *= 2. + *(CV2._bndrem);
    return CV2;
  }

  CVar<T>& CV2 = *CV._CV();
  for( unsigned int i=0; i<CV2._nmon(); i++ ) CV2._coefmon[i] = 0.;
  *(CV2._bndrem) = 0.;
  T bndmon = Op<T>::zeroone()*2.-1.;

  // Populate _coefmon and _bndrem for product of polynomial parts
  for( unsigned int i=0; i<CV2._nmon(); i++ ){
    for( unsigned int j=0; j<i; j++ ){
      const unsigned int*prodij = CV2._prodmon(i,j);
      for( unsigned int k=0; k<prodij[1]; k++ )
        CV2._coefmon[prodij[2+k]] += 2. * CV._coefmon[i] * CV._coefmon[j]
          / (double)prodij[0];
      *(CV2._bndrem) += bndmon * ( 2. * CV._coefmon[i] * CV._coefmon[j]
          * ( 1. - (double)prodij[1] / (double)prodij[0] ) );
    }
    const unsigned int*prodii = CV2._prodmon(i,i);
    for( unsigned int k=0; k<prodii[1]; k++ )
      CV2._coefmon[prodii[2+k]] += CV._coefmon[i] * CV._coefmon[i]
        / (double)prodii[0];
    *(CV2._bndrem) += bndmon * ( CV._coefmon[i] * CV._coefmon[i]
      * ( 1. - (double)prodii[1] / (double)prodii[0] ) );
  }

  // Populate _coefmon and _bndrem for product of polynomial and remainder parts
  //T PB = 0.;
  //for( unsigned int i=0; i<=CV3._nord(); i++ ) PB += CV._bndord[i];
  T PB = CV._polybnd();
  *(CV2._bndrem) += ( 2. * PB + *(CV._bndrem) ) * *(CV._bndrem);

  // Populate _bndord for product term (except remainder term)
  CV2._update_bndord();
  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>&
CVar<T>::operator *=
( const double c )
{
  if( !_CM ){
    _coefmon[0] *= c;
    *(_bndrem) *= c;
  }
  else{
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] *= c;
    for( unsigned int i=0; i<_nord()+2; i++ ) _bndord[i] *= c;
  }
  return *this;
}

template <typename T> inline CVar<T>
operator *
( const CVar<T>&CV1, const double c )
{
  CVar<T> CV3( CV1 );
  CV3 *= c;
  return CV3;
}

template <typename T> inline CVar<T>
operator *
( const double c, const CVar<T>&CV2 )
{
  CVar<T> CV3( CV2 );
  CV3 *= c;
  return CV3;
}

template <typename T> inline CVar<T>&
CVar<T>::operator *=
( const T&I )
{
  if( !_CM ){
    *(_bndrem) += _coefmon[0];
    _coefmon[0] = 0.;
    *(_bndrem) *= I;
  }
  else{
    const double Imid = Op<T>::mid(I);
    T Icur = bound();
    for( unsigned int i=0; i<_nmon(); i++ ) _coefmon[i] *= Imid;
    for( unsigned int i=0; i<_nord()+2; i++ ) _bndord[i] *= Imid;
    *_bndrem += (I-Imid)*Icur;
  }
  if( _CM && _CM->options.CENTER_REMAINDER ) _center_CM();
  return (*this);
}

template <typename T> inline CVar<T>
operator *
( const CVar<T>&CV1, const T&I )
{
  CVar<T> CV3( CV1 );
  CV3 *= I;
  return CV3;
}

template <typename T> inline CVar<T>
operator *
( const T&I, const CVar<T>&CV2 )
{
  CVar<T> CV3( CV2 );
  CV3 *= I;
  return CV3;
}

template <typename T> inline CVar<T>&
CVar<T>::operator /=
( const CVar<T>&CV )
{
   *this *= inv(CV);
   return *this;
}

template <typename T> inline CVar<T>
operator /
( const CVar<T>&CV1, const CVar<T>&CV2 )
{
  return CV1 * inv(CV2);
}

template <typename T> inline CVar<T>&
CVar<T>::operator /=
( const double c )
{
  if( isequal( c, 0. ) ) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::DIV );
   *this *= (1./c);
   return *this;
}

template <typename T> inline CVar<T>
operator /
( const CVar<T>&CV, const double c )
{
  if ( isequal( c, 0. )) throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::DIV );
  return CV * (1./c);
}

template <typename T> inline CVar<T>
operator /
( const double c, const CVar<T>&CV )
{
  return inv(CV) * c;
}

template <typename T> inline CVar<T>
inv
( const CVar<T>&CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::inv(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
sqrt
( const CVar<T>&CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::sqrt(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
exp
( const CVar<T>&CV )
{ 
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::exp(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
log
( const CVar<T>&CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::log(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
xlog
( const CVar<T>&CV )
{
  return CV * log( CV );
}

template <typename T> inline CVar<T>
pow
( const CVar<T>&CV, const int n )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::pow(CV._coefmon[0] + *(CV._bndrem), n);
    CV2._update_bndord();
    return CV2;
  }

  if( n < 0 ) return pow( inv( CV ), -n );
  CVar<T> CV2( CV._CM->_intpow( CV, n ) );
  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
CModel<T>::_intpow
( const CVar<T>&CV, const int n )
{
  if( n == 0 ) return 1.;
  else if( n == 1 ) return CV;
  return n%2 ? sqr( _intpow( CV, n/2 ) ) * CV : sqr( _intpow( CV, n/2 ) );
}

template <typename T> inline CVar<T>
pow
( const CVar<T> &CV, const double a )
{
  return exp( a * log( CV ) );
}

template <typename T> inline CVar<T>
pow
( const CVar<T> &CV1, const CVar<T> &CV2 )
{
  return exp( CV2 * log( CV1 ) );
}

template <typename T> inline CVar<T>
pow
( const double a, const CVar<T> &CV )
{
  return exp( CV * std::log( a ) );
}

template <typename T> inline CVar<T>
monomial
(const unsigned int n, const CVar<T>*CV, const int*k)
{
  if( n == 0 ){
    return 1.;
  }
  if( n == 1 ){
    return pow( CV[0], k[0] );
  }
  return pow( CV[0], k[0] ) * monomial( n-1, CV+1, k+1 );
}

template <typename T> inline CVar<T>
cos
( const CVar<T> &CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::cos(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
sin
( const CVar<T> &CV )
{
  return cos( CV - PI/2. );
}

template <typename T> inline CVar<T>
asin
( const CVar<T> &CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::asin(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
acos
( const CVar<T> &CV )
{
  return PI/2. - asin( CV );
}

template <typename T> inline CVar<T>
tan
( const CVar<T> &CV )
{
  return sin(CV) / cos(CV);
}

template <typename T> inline CVar<T>
atan
( const CVar<T> &CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::atan(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
fabs
( const CVar<T> &CV )
{
  if( !CV._CM ){
    CVar<T> CV2( CV );
    CV2._coefmon[0] = 0.;
    *(CV2._bndrem) = Op<T>::fabs(CV._coefmon[0] + *(CV._bndrem));
    CV2._update_bndord();
    return CV2;
  }

  CVar<T> CV2( CV._CM, 0. );
  throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF );

  //////////// TO BE IMPLEMENTED \\\\\\\\\\\\

  if( CV2._CM->options.CENTER_REMAINDER ) CV2._center_CM();
  return CV2;
}

template <typename T> inline CVar<T>
hull
( const CVar<T>&CV1, const CVar<T>&CV2 )
{
  // Neither operands associated to CModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1._coefmon[0] + *(CV1._bndrem);
    T R2 = CV2._coefmon[0] + *(CV2._bndrem);
    return Op<T>::hull(R1, R2);
  }

  // First operand not associated to CModel
  else if( !CV1._CM )
    return hull( CV2, CV1 );

  // Second operand not associated to CModel
  else if( !CV2._CM ){
    CVar<T> CVR = CV1.P();
    return CVR + Op<T>::hull( CV1.R(), CV2._coefmon[0]+*(CV2._bndrem)-CVR.B() );
  }

  // CModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::CMODEL );

  // Perform union
  CVar<T> CV1C( CV1 ), CV2C( CV2 );
  const double eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  T BCVD = (CV1C-CV2C).B();
  return (1.-eta)*CV1C + eta*CV2C + Op<T>::hull( R1C+eta*BCVD, R2C+(eta-1.)*BCVD );
}

template <typename T> inline bool
inter
( CVar<T>&CVR, const CVar<T>&CV1, const CVar<T>&CV2 )
{
  // Neither operands associated to CModel -- Make intersection in T type     
  if( !CV1._CM && !CV2._CM ){
    T R1 = CV1._coefmon[0] + CV1._bndord[0];
    T R2 = CV2._coefmon[0] + CV2._bndord[0];
    T RR( 0. );
    bool flag = Op<T>::inter(RR, R1, R2);
    CVR = RR;
    return flag;
  }

  // First operand not associated to CModel
  else if( !CV1._CM )
    return inter( CVR, CV2, CV1 );

  // Second operand not associated to CModel
  else if( !CV2._CM ){
    T R1 = CV1.R(), B2 = CV2.B();
    CVR = CV1.P();
    if( !Op<T>::inter(*(CVR._bndrem), R1, B2-CVR.B()) )
      return false;
    if( CVR._CM->options.CENTER_REMAINDER ) CVR._center_CM();
    return true;
  }

  // CModel for first and second operands are inconsistent
  else if( CV1._CM != CV2._CM )
    throw typename CModel<T>::Exceptions( CModel<T>::Exceptions::CMODEL );

  // Perform intersection
  CVar<T> CV1C( CV1 ), CV2C( CV2 );
  const double eta = CV1._CM->options.REF_POLY;
  T R1C = CV1C.C().R(), R2C = CV2C.C().R(); 
  CV1C.set(T(0.));
  CV2C.set(T(0.));
  CVR = (1.-eta)*CV1C + eta*CV2C;
  CV1C -= CV2C;
  T BCVD = CV1C.B();
  if( !Op<T>::inter( *(CVR._bndrem), R1C+eta*BCVD, R2C+(eta-1.)*BCVD ) )
    return false;
  if( CVR._CM->options.CENTER_REMAINDER ) CVR._center_CM();
  return true;
}

} // namespace mc

#include "mcop.hpp"

namespace mc
{

//! @brief C++ structure for specialization of the mc::Op templated structure to allow usage of the Chebyshev model type mc::CVar inside other MC++ type, e.g. mc::McCormick
template <> template<typename T> struct Op< mc::CVar<T> >
{
  typedef mc::CVar<T> CV;
  static CV point( const double c ) { return CV(c); }
  static CV zeroone() { return CV( mc::Op<T>::zeroone() ); }
  static void I(CV& x, const CV&y) { x = y; }
  static double l(const CV& x) { return mc::Op<T>::l(x.B()); }
  static double u(const CV& x) { return mc::Op<T>::u(x.B()); }
  static double abs (const CV& x) { return mc::Op<T>::abs(x.B());  }
  static double mid (const CV& x) { return mc::Op<T>::mid(x.B());  }
  static double diam(const CV& x) { return mc::Op<T>::diam(x.B()); }
  static CV inv (const CV& x) { return mc::inv(x);  }
  static CV sqr (const CV& x) { return mc::sqr(x);  }
  static CV sqrt(const CV& x) { return mc::sqrt(x); }
  static CV log (const CV& x) { return mc::log(x);  }
  static CV xlog(const CV& x) { return x*mc::log(x); }
  static CV fabs(const CV& x) { return mc::fabs(x); }
  static CV exp (const CV& x) { return mc::exp(x);  }
  static CV sin (const CV& x) { return mc::sin(x);  }
  static CV cos (const CV& x) { return mc::cos(x);  }
  static CV tan (const CV& x) { return mc::tan(x);  }
  static CV asin(const CV& x) { return mc::asin(x); }
  static CV acos(const CV& x) { return mc::acos(x); }
  static CV atan(const CV& x) { return mc::atan(x); }
  static CV erf (const CV& x) { throw typename mc::CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF ); }
  static CV erfc(const CV& x) { throw typename mc::CModel<T>::Exceptions( CModel<T>::Exceptions::UNDEF ); }
  static CV hull(const CV& x, const CV& y) { return mc::hull(x,y); }
  static CV min (const CV& x, const CV& y) { return mc::Op<T>::min(x.B(),y.B());  }
  static CV max (const CV& x, const CV& y) { return mc::Op<T>::max(x.B(),y.B());  }
  static CV arh (const CV& x, const double k) { return mc::exp(-k/x); }
  template <typename X, typename Y> static CV pow(const X& x, const Y& y) { return mc::pow(x,y); }
  static CV monomial (const unsigned int n, const T* x, const int* k) { return mc::monomial(n,x,k); }
  static bool inter(CV& xIy, const CV& x, const CV& y) { return mc::inter(xIy,x,y); }
  static bool eq(const CV& x, const CV& y) { return mc::Op<T>::eq(x.B(),y.B()); }
  static bool ne(const CV& x, const CV& y) { return mc::Op<T>::ne(x.B(),y.B()); }
  static bool lt(const CV& x, const CV& y) { return mc::Op<T>::lt(x.B(),y.B()); }
  static bool le(const CV& x, const CV& y) { return mc::Op<T>::le(x.B(),y.B()); }
  static bool gt(const CV& x, const CV& y) { return mc::Op<T>::gt(x.B(),y.B()); }
  static bool ge(const CV& x, const CV& y) { return mc::Op<T>::ge(x.B(),y.B()); }
};

} // namespace mc

#endif
