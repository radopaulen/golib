// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_ODEGPE Guaranteed Parameter Estimation in Dynamic Systems Described by Parametric IVPs in ODEs via Set Inversion
\author Benoit C. Chachuat, Radoslav Paulen, Mario Villanueva
\version 0.1
\date 2012
\bug No known bugs.


Consider a dynamic systems by parametric ODEs of the form
\f{align*}
 \dot{\bf x}(t) =&\ {\bf f} ({\bf x}(t), {\bf p}), \qquad {\bf x}(t_0) = {\bf h} ({\bf p}),\\
 \hat{\bf y}(t) =&\ {\bf g}({\bf x}(t), {\bf p}),
\f}
where \f${\bf x}\f$ denotes the \f$n_x\f$-dimensional vector of states; \f${\bf p}\f$, the \f$n_p\f$-dimensional vector of (a priori unknown) parameters; and \f$\hat{\bf y}\f$, the \f$n_y\f$-dimensional vector of model outputs (predictions). Suppose that a set of output measurements \f${\bf y}_{\rm m}\f$ at \f$N\f$ time points \f$t_1,\ldots,t_N\f$ are available, which provide estimates of the actual outputs \f${\bf y}_{\rm p}\f$ within some bounded measurement errors \f${\bf e}\in {\bf E} := [{\bf e}^L, {\bf e}^U]\f$, so that
\f{align*}
 {\bf y}_{\rm p}(t_i) \in {\bf y}_{\rm m}(t_i) + [{\bf e}^L, {\bf e}^U] =: {\bf Y}_i.
\f}

The problem addressed in <I>guaranteed (or bounded-error) parameter estimation</I> is to estimate the set \f${\bf P}_{\rm e}\f$ of <I>all</I> possible parameter values \f${\bf p}\f$ such that \f$\hat{\bf y}(t_i) \in {\bf Y}_i\f$ for every \f$i=1,\ldots,N\f$; that is, 
\f{align*}
 {\bf P}_{\rm e} := \left\{{\bf p}\in {\bf P}_0 \left|\begin{array}{l}
 \exists {\bf x}, {\bf y}:\\
 \dot{\bf x}(t) = {\bf f}({\bf x}(t), {\bf p}),\ {\bf x}(t_0) = {\bf h}({\bf p}),\\
 \hat{\bf y}(t_i) = {\bf g}({\bf x}(t_i), {\bf p}),\\
 \hat{\bf y}(t_i) \in {\bf Y}_i,\\
 \forall t\in[t_0,t_N], \forall i \in \{1, \dots, N\}
 \end{array}\right.\right\}.
\f}

\image html GPE_principle.png

Depicted in red on the left plot in the figure above is the set of all output trajectories satisfying \f$\hat{\bf y}(t_i) \in {\bf Y}_i\f$, \f$i=1,\ldots,N\f$, and on the right plot the corresponding set \f${\bf P}_{\rm e}\f$ projected in the \f$(p_1,p_2)\f$ space. Obtaining an exact characterization of the set \f${\bf P}_{\rm e}\f$ is not possible in general, and one has to resort to approximation techniques to make the problem computationally tractable. The class mc::ODEGPE in MC++ is an algorithm providing such approximations based on set inversion techniques. An over-approximation of the solution set of the parametric ODEs is computed using the class mc::ODEBND or mc::ODEBND_GSL, and domain reduction techniques are used to enhance the convergence of the algorithm by constructing LP relaxations using the class mc::FPRelax.

The class mc::ODEGPE is templated in two variable types T and PE:
- T is the class that performs the underlying interval arithmetics. Both verified and non-verified types are supported. For verified computations, the libraries <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> and <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> are supported; the non-verified interval type mc::Interval of MC++ can also be used. More generally, any interval type can be used provided that the templated structure mc::Op is instantiated accordingly -- see the header files <tt>mcprofil.h</tt> and <tt>mcfilib.h</tt> for examples.
- PE is the class that defines the parameter estimation problem, which has to be a class derived from mc::PESTRUCT. For example, suppose we want to find guaranteed parameter estimates for the following dynamic system involving two state variables \f${\bf x}=(x_1, x_2)^{\sf T}\f$ and three uncertain parameters \f${\bf p}=(p_1, p_2, p_3)^{\sf T}\in[0.01,1]^3\f$:
\f{align*}
 \dot x_1(t)=&-(p_1+p_3)x_1(t)+p_2x_2(t), \qquad x_1(0)=1,\\
 \dot x_2(t)=& p_1x_1(t)-p_2x_2(t), \qquad x_2(0)=0,
\f}
and having a single output variable \f$y\f$ corresponding to the state variable \f$x_2\f$,
\f{align*}
 \hat y(t_k)=&x_2(t_k), \qquad t_k=1,\ldots,15.
\f}
The following class is defined:
.
\code
      #include "odegpe.h"

      // Number of time stages
      const unsigned int NS = 15;
      // Number of model parameters
      const unsigned int NP = 3;
      // Number of state variables
      const unsigned int NX = 2;
      // Number of output variables
      const unsigned int NY = 1;
      // Number of constraints
      const unsigned int NC = 0;

      class PE : public virtual mc::PESTRUCT
      {
      public:
        PE()
          : mc::PESTRUCT( NP, NX, NY, NC )
          {}

        // ODEs right-hand side
	template <typename TX, typename TP>
        TX RHS
          ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
          {
            assert( ix < nx() );
            switch( ix ){
              case 0: return -(p[0]+p[2])*x[0]+p[1]*x[1];
              case 1: return p[0]*x[0]-p[1]*x[1];
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
              case 0: return 1.;
              case 1: return 0.;
              default: throw std::runtime_error("invalid size");
            }
          }
      
        // Initial conditions
        template <typename T>
        T OUT
          ( const unsigned int iy, const T*p, const TX*x, const unsigned int is )
          {
            assert( iy < ny() );
            switch( iy ){
              case 0: return x[1];
              default: throw std::runtime_error("invalid size");
            }
          }

        // Constraints
        template <typename T>
        std::pair<T,t_CTR> CTR
          ( const unsigned int ic, const T*p, T* const*x, const unsigned int ns )
          {
            throw std::runtime_error("invalid size");
          }
      };
\endcode

\section sec_ODEGPE_solve How to Solve a GPE Problem using mc::ODEGPE?

Start by defining the initial parameter set \f${\bf P}\f$, and the time stages \f$t_i\f$, \f$i=0,\ldots,15\f$ (including the initial time and defined in such a way that they match the measurement times---see below):

\code
  typedef mc::Interval I;
  I P[NP];
  for( unsigned int ip=0; ip<NP; ip++ )
    P[ip] = I( 1e-2, 1e0 );

  double tk[NS+1];
  for( unsigned int is=0; is<NS+1; is++ )
    tk[is] = is;
\endcode

Also define the set of measurements \f${\bf Y}\f$, where each record consists of the measurement time, output variable and measurement range:

\code
  typename std::list< mc::ODEGPE<I,PE>::Data > Y;
  double yk[NS] = { .36, .47, .49, .48, .46, .44, .42, .4, .38, .36, .35, .33, .31, .3, .28 };
  double e = 5e-3;
  for( unsigned int is=0; is<NS; is++ )
    Y.push_back( typename mc::ODEGPE<I,PE>::Data( is+1, 0, yk[is]-e, yk[is]+e ) );
\endcode

Then, define an mc::ODEGPE object:

\code
  mc::ODEGPE<I,PE> pPE;
\endcode

The guaranteed parameter estimation problem is solved as follows:

\code
  std::pair<double,double> optimum = pPE.solve( NS, tk, P, Y );
\endcode

where the first and second arguments are the number and value of the time stages, the third one is the parameter bounds, and the fourth one is the measurement data. The return value is of a pair, whose elements are the volumes of the computed inner-approximation and boundary of the guaranteed parameter set, respectively. 

The following result is obtained for this problem:

\verbatim
  INDEX      VOLUME    STACK   CUMUL TIME       INNER         OPEN     PARENT        ACTION

      1  9.702990e-01      1  0.000000e+00  0.000000e+00  9.702990e-01      0       BRANCH2
      2  9.170232e-01      2  5.600300e-02  0.000000e+00  9.170232e-01      1       BRANCH1
      3  9.170232e-01      3  1.040060e-01  0.000000e+00  8.705861e-01      1       BRANCH1
      4  4.120745e-01      4  2.080130e-01  0.000000e+00  7.074116e-01      2       BRANCH0
      5  4.120745e-01      5  2.480150e-01  0.000000e+00  6.697194e-01      2       BRANCH0
      6  2.953371e-01      6  3.360210e-01  0.000000e+00  5.926902e-01      3       BRANCH0
      7  2.953371e-01      7  3.800240e-01  0.000000e+00  5.748107e-01      3         OUTER
      8  1.683451e-01      6  5.000310e-01  0.000000e+00  4.271422e-01      4       BRANCH1
      9  1.683451e-01      7  5.840360e-01  0.000000e+00  4.026270e-01      4       BRANCH1
     10  1.297890e-01      8  6.720420e-01  0.000000e+00  3.663817e-01      6       BRANCH2
     11  1.297890e-01      9  7.600470e-01  0.000000e+00  3.405974e-01      6       BRANCH1
     12  1.290081e-01     10  8.040500e-01  0.000000e+00  3.405974e-01      5       BRANCH1
     13  1.290081e-01     11  8.960560e-01  0.000000e+00  3.091632e-01      5       BRANCH1
     14  6.489451e-02     12  9.400590e-01  0.000000e+00  3.091632e-01     11         OUTER
     15  6.489451e-02     11  1.076067e+00  0.000000e+00  2.767160e-01     11         OUTER
     16  6.450404e-02     10  1.136071e+00  0.000000e+00  2.442687e-01     13       BRANCH2
     17  6.450404e-02     11  1.180074e+00  0.000000e+00  2.428053e-01     13         OUTER
     18  5.965739e-02     10  1.320082e+00  0.000000e+00  2.105533e-01      8       BRANCH2
     19  5.965739e-02     11  1.356085e+00  0.000000e+00  2.101915e-01      8       BRANCH2
     20  4.792719e-02     12  1.444090e+00  0.000000e+00  1.988221e-01      9         OUTER
\endverbatim
...
\verbatim
    829  4.098610e-08    783  4.187462e+01  0.000000e+00  1.008018e-05    386       BRANCH0
    830  4.076584e-08    784  4.193862e+01  0.000000e+00  1.006774e-05    398       BRANCH2
    831  4.076584e-08    785  4.197062e+01  0.000000e+00  1.006646e-05    398       BRANCH2
    832  4.073699e-08    786  4.200662e+01  0.000000e+00  1.006644e-05    437       BRANCH1
    833  4.073699e-08    787  4.207063e+01  0.000000e+00  1.006199e-05    437       BRANCH1
    834  4.049112e-08    788  4.210663e+01  0.000000e+00  1.005813e-05    453       BRANCH1
    835  4.049112e-08    789  4.213863e+01  0.000000e+00  1.005494e-05    453       BRANCH1
    836  4.035653e-08    790  4.217064e+01  0.000000e+00  1.005110e-05    400       BRANCH0
    837  4.035653e-08    791  4.220264e+01  0.000000e+00  1.005110e-05    400       BRANCH0
    838  4.026665e-08    792  4.227064e+01  0.000000e+00  1.004094e-05    295       BRANCH0
    839  4.026665e-08    793  4.233465e+01  0.000000e+00  1.002652e-05    295       BRANCH1
    840  4.013921e-08    794  4.237065e+01  0.000000e+00  1.002651e-05    502       BRANCH2
    841  4.013921e-08    795  4.240265e+01  0.000000e+00  1.002650e-05    502       BRANCH2
    842  3.977716e-08    796  4.243465e+01  0.000000e+00  1.002643e-05    464       BRANCH0
    843  3.977716e-08    797  4.246665e+01  0.000000e+00  1.002643e-05    464       BRANCH0
    844  3.975389e-08    798  4.253466e+01  0.000000e+00  1.001954e-05    454       BRANCH2
    845  3.975389e-08    799  4.256666e+01  0.000000e+00  1.001634e-05    454       BRANCH2
    846  3.974725e-08    800  4.260266e+01  0.000000e+00  1.001355e-05    513       BRANCH0
    847  3.974725e-08    801  4.266667e+01  0.000000e+00  1.000899e-05    513       BRANCH0
    848  3.957155e-08    802  4.273467e+01  0.000000e+00  1.000385e-05    492       BRANCH2
    849  3.957155e-08    803  4.277067e+01  0.000000e+00  1.000039e-05    492       BRANCH2

#  NORMAL TERMINATION:      42.770673 CPU SEC
#  INNER APPROXIMATION:     VOLUME = 0.000000e+00,  NODES = 0
#  BOUNDARY APPROXIMATION:  VOLUME = 9.990558e-06,  NODES = 804
#  TOTAL NUMBER OF NODES:   848
#  MAX OPEN NODES IN STACK: 804
\endverbatim

The partitions \f$\mathbb{P}_{\rm int}\f$ and \f$\mathbb{P}_{\rm bnd}\f$ of inner and boundary parameter subsets, respectively -- see Figure above -- can be recorded to files after completion as follows:

\code
  #include <fstream>
  std::ofstream os_open(  "OPEN_NODES", std::ios_base::out );
  std::ofstream os_inner( "INNER_NODES", std::ios_base::out );

  GPE.output_stacks( os_open, os_inner );

  os_open.close();
  os_inner.close();
\endcode

Each line consists of the lower and upper bounds for all the parameters.

\section sec_ODEGPE_opt How to Set the Options for the Solution of a GPE Problem using mc::ODEGPE?

All of the options are defined in the public structure mc::ODEGPE::Options. Current options are:

<TABLE border="1">
<CAPTION><EM>Options in mc::ODEBND::Options: name, type and
description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>IVP_BOUNDING_PROPAGATION</tt>   <TD><tt>mc::ODEGPE::Options::IVPMETH</tt>  <TD>mc::ODEGPE::Options::VERIF <TD>Set the IVP bounding method
     <TR><TH><tt>IVP_BOUNDING_TYPE</tt>          <TD><tt>mc::ODEGPE::Options::IVPBOUND</tt> <TD>mc::ODEGPE::Options::TM    <TD>Set the IVP bounder
     <TR><TH><tt>TAYLOR_MODEL_ORDER</tt>         <TD><tt>unsigned int</tt>                  <TD>2                          <TD>Order of Taylor/McCormick-Taylor model for IVP bounding (only if mc::ODEGPE::Options::TM or mc::ODEGPE::Options::MCTM selected for IVP_BOUNDING_TYPE)
     <TR><TH><tt>USE_CONSTRAINT_PROPAGATION</tt> <TD><tt>bool</tt>                          <TD>false                      <TD>Whether to pre-process nodes using constraint propagation (make sure to use a validated interval type if applied)
     <TR><TH><tt>USE_DOMAIN_REDUCTION</tt>       <TD><tt>bool</tt>                          <TD>true                       <TD>Whether to post-process nodes using optimality-based domain reduction
     <TR><TH><tt>USE_RLT_TMODEL</tt>             <TD><tt>bool</tt>                          <TD>true                       <TD>Whether to tighten Taylor model relaxations using RLT constraints
     <TR><TH><tt>DOMAIN_REDUCTION_MAXLOOPS</tt>  <TD><tt>unsigned int</tt>                  <TD>4                          <TD>Maximum number of domain reduction loops per node
     <TR><TH><tt>DOMAIN_REDUCTION_THRESHOLD</tt> <TD><tt>double</tt>                        <TD>2e-1                       <TD>Threshold for repeating domain reduction (minimum domain reduction ratio for any variable)
</TABLE>

These options can be set/modified using the public (non-static) member mc::ODEGPE::options as follows:

\code
   GPE.options.TAYLOR_MODEL_ORDER          = 3;
   GPE.options.IVP_BOUNDING_TYPE           = mc::ODEGPE<I,PE>::Options::TM;
   GPE.options.IVP_BOUNDING_PROPAGATION    = mc::ODEGPE<I,PE>::Options::DINEQ;
   GPE.options.USE_CONSTRAINT_PROPAGATION  = false;
   GPE.options.USE_DOMAIN_REDUCTION        = true;
   GPE.options.DOMAIN_REDUCTION_MAXLOOPS   = 10;      
   GPE.options.USE_RLT_TMODEL              = true;
\endcode

Options for other components of MC++ can be set/modified by using the public member functions: mc::ODEGPE::options_SetInv(), mc::ODEGPE::options_ODEBND(), mc::ODEGPE::options_ODEBND_GSL(), mc::ODEGPE::options_FPRelax(), and mc::ODEGPE::options_TModel(), which return references to the option structures in, respectively, mc::SetInv, mc::ODEBND, mc::ODEBND_GSL, mc::FPRelax, and mc::TModel.

For instance, the options in the class mc::SetInv can be set/modified as:

\code
   GPE.options_SetInv().DISPLAY            = 2;
   GPE.options_SetInv().MAX_CPU_TIM        = 1e4;
   GPE.options_SetInv().MAX_NODES          = 1000000;
   GPE.options_SetInv().ABSOLUTE_TOLERANCE = 1e-6;
   GPE.options_SetInv().RELATIVE_TOLERANCE = 1e-6;
\endcode

\section sec_ODEGPE_refs References

- Jaulin, L., and E. Walter, <A href="http://dx.doi.org/10.1016/0005-1098(93)90106-4">Set inversion via interval analysis for nonlinear bounded-error estimation</A>, <i>Automatica</i>, <b>29</b>(4):1053--1064, 1993.
- Lin, Y., and M. A. Stadtherr, <A href="http://dx.doi.org/10.1021/ie0707725"> Guaranteed state and parameter estimation for nonlinear continuous-time systems with bounded-error measurements,</A> <i>Industrial & Engineering Chemistry Research</I>, <B>46</B>:7198--7207, 2007.
- Kieffer, M., and E. Walter, <A href="http://dx.doi.org/10.1002/acs.1194">Guaranteed estimation of the parameters of nonlinear continuous-time models: Contributions of interval analysis</A>, <i>Intermation Journal of Adaptive Control & Signal Processing</i>, <b>25</b>:191--207, 2011.
.
*/

#ifndef MC__ODEGPE_HPP
#define MC__ODEGPE_HPP

#undef USE_CPLEX

#include <stdexcept>
#include <cassert>
#include "setinv.hpp"
#include "pestruct.hpp"
#include "odebnd_val.hpp"
#include "odebnd_gsl.hpp"
#ifdef USE_CPLEX
  #include "fprcplex.hpp"
#else
  #include "fprgurobi.hpp"
#endif
#include "structure.hpp"

#undef MC__ODEGPE_DEBUG
#undef MC__ODEGPE_SHOW_FAILURE
#undef MC__ODEGPE_SHOW_REDUCTION
#undef MC__ODEGPE_TRACE

/* TO DO:
- Describe options in documentation

BUGS:
*/

namespace mc
{
//! @brief C++ class for guaranteed parameter estimation in dynamic systems using set inversion techniques
////////////////////////////////////////////////////////////////////////
//! mc::ODEGPE is a C++ class doing guaranteed parameter estimation in
//! dynamic systems using set-inversion techniques
////////////////////////////////////////////////////////////////////////
template <typename T, typename PE>
class ODEGPE: protected SetInv<T>, public PESTRUCT
{
  // Overloading stdout operator
  template <typename U, typename PU> friend std::ostream& operator<<
    ( std::ostream&os, const ODEGPE<U,PU>& );

public:
  typedef Ellipsoid    E;
  typedef TModel<T>    TMT;
  typedef TVar<T>      TVT;
#ifdef USE_CPLEX
  typedef FPRCplex<T>  FPRT;
#else
  typedef FPRGurobi<T> FPRT;
#endif
  typedef FPVar<T>     FPVT;
  typedef TModel<FPVT> TMFPVT;

  /** @defgroup ODEGPE Guaranteed Parameter Estimation in Dynamic Systems using Set Inversion
   *  @{
   */
  //! @brief Constructor
  ODEGPE();

  //! @brief Destructor
  virtual ~ODEGPE();

  //! @brief ODEGPE options
  struct Options
  {
    //! @brief Constructor
    Options():
      IVP_BOUNDING_TYPE(TM), IVP_BOUNDING_PROPAGATION(DINEQ),
      TAYLOR_MODEL_ORDER(2), USE_CONSTRAINT_PROPAGATION(false),
      USE_DOMAIN_REDUCTION(true), USE_RLT_TMODEL(true),
      DOMAIN_REDUCTION_MAXLOOPS(4), DOMAIN_REDUCTION_THRESHOLD(0.2)
      {}
    //! @brief Display
    void display
      ( std::ostream&out ) const;
    //! @brief IVP bounder
    enum IVPBOUND{
      IA=0,	//!< Interval analysis
      TM,	//!< Taylor model
    };
    //! @brief IVP bounding method
    enum IVPMETH{
      DINEQ=0,	//!< Differential methods not accouting for truncation errors
      VERIF	//!< Verified integration accounting for truncation errors
    };
    //! @brief IVP bounder
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
   //! @brief Maximum number of domain reduction loops per node
    unsigned int DOMAIN_REDUCTION_MAXLOOPS;
    //! @brief Threshold for repeating domain reduction (minimum domain reduction ratio)
    double DOMAIN_REDUCTION_THRESHOLD;
  } options;

  //! @brief ODEGPE exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for ODEGPE exception handling
    enum TYPE{
      DATA=0,	//!< Error in measurement data specification
      REDUC,	//!< Error during variable domain reduction
      UNDEF=-33	//!< Error due to calling a function/feature not yet implemented
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }
    std::string what(){
      switch( _ierr ){
      case DATA:
        return "Error in measurement data specification";
      case REDUC:
        return "Error during variable domain reduction";
      case UNDEF: default:
        return "Error due to calling a feature not yet implemented in ODEGPE";
      }
    }
  private:
    TYPE _ierr;
  };

  //! @brief Structure holding current statistics
  struct Stats{
    double cumul_ODEBND;
    double cumul_FPREL;
    double cumul_REDUC;
    double cumul_SETINV;
  } stats;

  //! @brief Structure to hold measurement data
  struct Data{
    Data( unsigned int is, unsigned int iy, double yl, double yu ):
      stage(is), index(iy), lower(yl), upper(yu)
      {}
    unsigned int stage;
    unsigned int index;
    double lower;
    double upper;
  };

  //! @brief Reference to SetInv options
  typename SetInv<T>::Options& options_SetInv()
    { return SetInv<T>::options; }
  //! @brief Reference to ODEBND options
  typename ODEBND_VAL<T,PE>::Options& options_ODEBND_VAL()
    { return _pODEBND_VAL->options; }
  //! @brief Reference to ODEBND_GSL options
  typename ODEBND_GSL<T,PE>::Options& options_ODEBND_GSL()
    { return _pODEBND_GSL->options; }
  //! @brief Reference to FPRelax options
  typename FPRelax<T>::Options& options_FPRelax()
    { return _pFPR->options; }
  //! @brief Reference to Taylor model options
  typename TModel<T>::Options& options_TModel()
    { return *_pTMOpt; }

  //! @brief Solve GPE problem -- return value is pair of inner-approximation and boundary volumes
  std::pair<double,double> solve
    ( const unsigned int ns, const double*tk, const T*P,
      const std::list< Data >*Y, std::ostream&os=std::cout );
  //! @brief Append all open nodes and inner nodes to, respectively, <a>os_open</a> and <a>os_inner</a> - Number of significant digits is set via <a>DPREC</a> (default=6)
  void output_stacks
    ( std::ostream&os_open, std::ostream&os_inner,
      const unsigned int DPREC=6 ) const;
  /** @} */

private:
  //! @brief ODE Bounder
  ODEBND_VAL<T,PE>* _pODEBND_VAL;
  //! @brief ODE Bounder
  ODEBND_GSL<T,PE>* _pODEBND_GSL;
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
  //! @brief TMFP environment for IVP polyhedral relaxation
  TModel<FPVT>* _pTMFP;
  //! @brief Array holding variable bounds
  std::pair<double,double>* _pP;
  //! @brief TM options
  typename TModel<T>::Options* _pTMOpt;

  //! @brief Variable holding PE problem structure
  Structure _PEstruct[2];
  //! @brief Set of variables participating in IVP
  std::set<unsigned int> _IVPdepend;
  //! @brief Array holding IVP parameter indices
  unsigned int *_IVPparam;

  //! @brief Structure holding problem data for internal use
  struct Internal{
    unsigned int ns;
    const double*tk;  
    const T*P;
    const std::list< Data >*Y;
  } _data;

  //! @brief User-function to subproblem assessment
  typename SetInv<T>::STATUS assess
    ( const unsigned int np, T*P );

  //! @brief Set PE problem structure
  void _set_PE_structure();
  //! @brief Initialize IVP relaxation
  void _init_IVP_relax();

  //! @brief Create output and constraint polyhedral relaxations
  void _relaxation_setup
    ( const T*P );
  //! @brief Create state polyhedral relaxations via interval analysis
  typename SetInv<T>::STATUS _bound_IVP_IA
    ( const T*P, FPVT**FPVxk, typename Options::IVPMETH meth );
  //! @brief Bound states using Taylor models and polyhedral relaxations
  typename SetInv<T>::STATUS _bound_IVP_TM
    ( const T*P, FPVT**FPVxk, typename Options::IVPMETH meth );

  //! @brief Bound given output and test enclosure
  typename SetInv<T>::STATUS _bound_output
    ( const T&Iyk, const double&lo, const double&up );
  //! @brief Bound outputs and test enclosure (interval analyis wrapper)
  typename SetInv<T>::STATUS _bound_output
    ( const T*Ip, T* const*Ixk );
  //! @brief Bound outputs and test enclosure (Taylor model wrapper)
  typename SetInv<T>::STATUS _bound_output
    ( const TVT*TVp, TVT* const*TVxk );

  //! @brief Create output and constraint polyhedral relaxations
  typename SetInv<T>::STATUS _relax_outputs_constraints
    ( FPVT*const*FPVxk );
  //! @brief Perform optimization-based domain reduction
  typename SetInv<T>::STATUS _domain_reduction
    ( FPVT*const*FPVxk, T*P, double &maxred );

  //! @brief Private methods to block default compiler methods
  ODEGPE(const ODEGPE&);
  ODEGPE& operator=(const ODEGPE&);
};

template <typename T, typename PE>
inline
ODEGPE<T,PE>::ODEGPE()
: PESTRUCT(PE())
{
  _pODEBND_VAL = new ODEBND_VAL<T,PE>;
  _pODEBND_GSL = new ODEBND_GSL<T,PE>;
  _pFPR      = new FPRT;
  _pFPV      = new FPVT[_np];
  _pTM   = 0;
  _pTMFP = 0;
  _IVPparam = 0;
  _pP = new std::pair<double,double>[_np];
  _pTMOpt = new typename TModel<T>::Options;
}

template <typename T, typename PE>
inline
ODEGPE<T,PE>::~ODEGPE()
{
  delete _pODEBND_VAL;
  delete _pODEBND_GSL;
  delete _pFPR;
  delete[] _pFPV;
  delete _pTM;
  delete _pTMFP;
  delete[] _IVPparam;
  delete[] _pP;
  delete _pTMOpt;
}

template <typename T, typename PE>
inline void
ODEGPE<T,PE>::output_stacks
( std::ostream&os_open, std::ostream&os_inner, const unsigned int DPREC ) const
{
  return SetInv<T>::output_stacks( os_open, os_inner, DPREC );
}

template <typename T, typename PE>
inline std::pair<double,double>
ODEGPE<T,PE>::solve
( const unsigned int ns, const double*tk, const T*P,
  const std::list< Data >*Y, std::ostream&os )
{
  // Keep track of time stages, bounds, and measurement data
  _data.tk = tk;
  _data.ns = ns;
  _data.P  = P;
  _data.Y  = Y;

  // Keep track of execution times
  stats.cumul_ODEBND = stats.cumul_FPREL = stats.cumul_REDUC =
  stats.cumul_SETINV = 0.;

  // Determine PE problem structure
  _set_PE_structure();
  _init_IVP_relax();

  // Run set-inversion algorithm
  SetInv<T>::variables( _np, P );
  stats.cumul_SETINV = -time();
  const std::pair<double,double>& res = SetInv<T>::solve( os );
  stats.cumul_SETINV += time();
  return res;
}

template <typename T, typename PE>
inline void
ODEGPE<T,PE>::_set_PE_structure()
{
  // Initialize structure detection variables
  Structure *pS = new Structure[_np];
  for( unsigned int ip=0; ip<_np; ip++ )
    pS[ip].indep(ip);   
  Structure **xkS = new Structure*[_data.ns+1];

  // detect structure of initial conditions in IVP
  xkS[0] = new Structure[_nx];
  for( unsigned int ix=0; ix<_nx; ix++ ){
    xkS[0][ix] = PE().IC( ix, pS );
#ifdef MC__ODEGPE_DEBUG
    std::cout << "xkS[0][" << ix << "] <- " << xkS[0][ix]
              << std::endl;
#endif
  }

  // Detect structure of right-hand side in each stage in IVP
  for( unsigned int is=0; is<_data.ns; is++ ){
    xkS[is+1] = new Structure[_nx];
    for( unsigned int ix=0; ix<_nx; ix++ )
      xkS[is+1][ix] = PE().RHS( ix, pS, xkS[is], _data.tk[is], is ) + xkS[is][ix];
    bool iterate = true;
    while( iterate ){
      iterate = false;
      for( unsigned int ix=0; ix<_nx; ix++ ){
        Structure xiSnew = PE().RHS( ix, pS, xkS[is+1], _data.tk[is], is ) + xkS[is+1][ix];
        if( xkS[is+1][ix] != xiSnew ){
	  xkS[is+1][ix] = xiSnew; iterate = true;
	}
      }
    }
#ifdef MC__ODEGPE_DEBUG
    for( unsigned int ix=0; ix<_nx; ix++ )
      std::cout << "xkS[" << is+1 << "][" << ix << "] <- " << xkS[is+1][ix]
                << std::endl;
#endif
  }

  // Set IVP variable dependence in _IVPdepend
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

  // Detect structure of constraints (parameter dependence)
  _PEstruct[0] = 0.;
  for( unsigned int ic=0; ic<_nc; ic++ )
    _PEstruct[0] += PE().CTR( ic, pS, xkS, _data.ns ).first;
  typename std::list< Data >::const_iterator cit = _data.Y->begin();
  for( ; cit != _data.Y->end(); ++cit ){
    if( (*cit).index >= _ny || (*cit).stage > _data.ns )
      throw Exceptions( Exceptions::DATA );
    _PEstruct[0] += PE().OUT( (*cit).index, pS, xkS[(*cit).stage], (*cit).stage );
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << "PE <- " << _PEstruct[0] << std::endl;
  {int dum; std::cin >> dum;}
#endif

  // Variables to exclude from branching
  SetInv<T>::_exclude_vars.clear();
  std::map<int,bool>& PEdep = _PEstruct[0].dep();
  for( unsigned int ip=0; ip<_np; ip++ ){
    std::map<int,bool>::iterator it = PEdep.find(ip);
    if( it == PEdep.end() ) SetInv<T>::_exclude_vars.insert(ip);
  }

  // Structure of constraints (parameter & state dependence)
  for( unsigned int is=0, ixs=_np; is<=_data.ns; is++ )
    for( unsigned int ix=0; ix<_nx; ix++, ixs++)
      xkS[is][ix].indep(ixs);
  _PEstruct[1] = 0.;
  for( unsigned int ic=0; ic<_nc; ic++ )
    _PEstruct[1] += PE().CTR( ic, pS, xkS, _data.ns ).first;
  cit = _data.Y->begin();
  for( ; cit != _data.Y->end(); ++cit ){
    if( (*cit).index >= _ny || (*cit).stage > _data.ns )
      throw Exceptions( Exceptions::DATA );
    _PEstruct[1] += PE().OUT( (*cit).index, pS, xkS[(*cit).stage], (*cit).stage );
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << "PE <- " << _PEstruct[1] << std::endl;
  {int dum; std::cin >> dum;}
#endif

  // clean up
  delete[] pS;
  for( unsigned int is=0; is<=_data.ns; is++ )
    delete[] xkS[is];
  delete[] xkS;
}

template <typename T, typename PE>
inline void
ODEGPE<T,PE>::_init_IVP_relax
()
{
  delete _pTM;   _pTM   = 0;
  delete _pTMFP; _pTMFP = 0;

  // Match IVP parameters to PE problem in array _IVPparam
  delete[] _IVPparam;
  _IVPparam = new unsigned int[_IVPdepend.size()];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ )
    if( _IVPdepend.find(ip) != _IVPdepend.end() ) _IVPparam[isub++] = ip;

  switch( options.IVP_BOUNDING_TYPE ){
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

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_bound_IVP_IA
( const T*P, FPVT**FPVxk, typename Options::IVPMETH meth )
{
  // Set interval bounds for states at stage times
  stats.cumul_ODEBND -= time();
  T* Ixk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    Ixk[is] = new T[_nx];
  E Exk[_data.ns+1];
  
  // Compute state bounds at stage times
  bool failed = true;
  switch( meth ){
  case Options::VERIF:
    if( _pODEBND_VAL->bounds( _data.ns, _data.tk, P, Ixk, Exk )
        == ODEBND_VAL<T,PE>::NORMAL ) failed = false; break;
  case Options::DINEQ: default:
    if( _pODEBND_GSL->bounds( _data.ns, _data.tk, P, Ixk )
        == ODEBND_GSL<T,PE>::NORMAL ) failed = false; break;
  }
  if( failed ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] Ixk[is];
    stats.cumul_ODEBND += time();
    return SetInv<T>::FAILURE;
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << *_pODEBND_VAL;
#endif
  stats.cumul_ODEBND += time();

  // Set current independent variables (_pFPV) and state variables (FPVxk)
  //in FPRelax
  stats.cumul_FPREL -= time();
  typename SetInv<T>::STATUS flag = _bound_output( P, Ixk );
  if( flag == SetInv<T>::UNDETERMINED && options.USE_DOMAIN_REDUCTION ){   
    _relaxation_setup( P );
    std::map<int,bool>& PEdep = _PEstruct[1].dep();
    failed = false;
    for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
      for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
        if( Op<T>::diam(Ixk[is][ix]) >= PESTRUCT::INF ){
          failed = true; break;
        }
        std::map<int,bool>::iterator it = PEdep.find(ixs);
        if( it!=PEdep.end() )
          FPVxk[is][ix] = Ixk[is][ix];
      }
  }
  stats.cumul_FPREL += time();
  
  // Clean-up
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] Ixk[is];

  return flag;
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_bound_output
( const T*Ip, T* const*Ixk )
{
  bool inner = true;
  typename std::list< Data >::const_iterator cit = _data.Y->begin();
  for( ; cit != _data.Y->end(); ++cit ){
    switch( _bound_output( PE().OUT( (*cit).index, Ip, Ixk[(*cit).stage],
                           (*cit).stage ), (*cit).lower, (*cit).upper ) ){
      case SetInv<T>::OUTER:        return SetInv<T>::OUTER;
      case SetInv<T>::UNDETERMINED: inner = false; continue;
      case SetInv<T>::INNER:        default: continue;
    }
  }
  return( inner? SetInv<T>::INNER: SetInv<T>::UNDETERMINED );
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_bound_IVP_TM
( const T*P, FPVT**FPVxk, typename Options::IVPMETH meth )
{
  // Set Taylor model for parameters
  stats.cumul_ODEBND -= time();
  TVT TVp[_np];
  for( unsigned int ip=0, isub=0; ip<_np; ip++ ){
    if( _IVPdepend.find(ip) != _IVPdepend.end() ){
      TVp[ip] = TVT( _pTM, isub, P[ip] );
      isub++;
    }
    else
      TVp[ip] = P[ip];
  }

  // Set Taylor model for states at stage times
  TVT* TVxk[_data.ns+1];
  for( unsigned int is=0; is<=_data.ns; is++ )
    TVxk[is] = new TVT[_nx];
  E Exk[_data.ns+1];
  
  // Compute state bounds at stage times
  bool failed = true;
  switch( meth ){
  case Options::VERIF:
    if( _pODEBND_VAL->bounds( _data.ns, _data.tk, TVp, E(), TVxk, Exk )
        == ODEBND_VAL<T,PE>::NORMAL ) failed = false; break;
  case Options::DINEQ: default:
    if( _pODEBND_GSL->bounds( _data.ns, _data.tk, TVp, TVxk )
        == ODEBND_GSL<T,PE>::NORMAL ) failed = false; break;
  }
  if( failed ){
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];
    stats.cumul_ODEBND += time();
    return SetInv<T>::FAILURE;
  }
#ifdef MC__ODEGPE_DEBUG
  std::cout << *_pODEBND_VAL;
#endif
  stats.cumul_ODEBND += time();

  // Set current independent variables (_pFPV) and state variables (FPVxk)
  //in FPRelax
  stats.cumul_FPREL -= time();
  typename SetInv<T>::STATUS flag = _bound_output( TVp, TVxk );
  if( flag == SetInv<T>::UNDETERMINED && options.USE_DOMAIN_REDUCTION ){
    _relaxation_setup( P );
    _pFPR->update_TModel( _pTMFP, _pTM, _IVPparam );
    std::map<int,bool>& PEdep = _PEstruct[1].dep();
    for( unsigned int is=0, ixs=_np; !failed && is<=_data.ns; is++ )
      for( unsigned int ix=0; ix<_nx; ix++, ixs++ ){
        if( Op<T>::diam(TVxk[is][ix].B()) >= PESTRUCT::INF ){
          failed = true; break;
        }
        std::map<int,bool>::iterator it = PEdep.find(ixs);
        if( it!=PEdep.end() )
          FPVxk[is][ix] = FPVT( _pFPR, TVxk[is][ix], _pTMFP );
      }
  }
  stats.cumul_FPREL += time();
  
  // Clean-up
  for( unsigned int is=0; is<=_data.ns; is++ ) delete[] TVxk[is];

  return flag;
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_bound_output
( const TVT*TVp, TVT* const*TVxk )
{
  bool inner = true;
  typename std::list< Data >::const_iterator cit = _data.Y->begin();
  for( ; cit != _data.Y->end(); ++cit ){
    switch( _bound_output( PE().OUT( (*cit).index, TVp, TVxk[(*cit).stage],
                           (*cit).stage ).B(), (*cit).lower, (*cit).upper ) ){
      case SetInv<T>::OUTER:        return SetInv<T>::OUTER;
      case SetInv<T>::UNDETERMINED: inner = false; continue;
      case SetInv<T>::INNER:        default: continue;
    }
  }
  return( inner? SetInv<T>::INNER: SetInv<T>::UNDETERMINED );
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_bound_output
( const T&Iyk, const double&lo, const double&up )
{
  if( Op<T>::l(Iyk) > up || Op<T>::u(Iyk) < lo )
    return SetInv<T>::OUTER; 
  if( Op<T>::l(Iyk) < lo || Op<T>::u(Iyk) > up )
    return SetInv<T>::UNDETERMINED;
  return SetInv<T>::INNER;
}

template <typename T, typename PE>
inline void
ODEGPE<T,PE>::_relaxation_setup
( const T*P )
{
  // Set current variables in FP representation
  _pFPR->reset();
  for( unsigned int ip=0; ip<_np; ip++ )
    _pFPV[ip] = FPVT( _pFPR, ip, P[ip] );
#ifdef MC__ODEGPE_DEBUG
  std::cout << "FPR:\n" << *_pFPR;
#endif
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_relax_outputs_constraints
( FPVT*const*FPVxk )
{
  // User constraints definition
  stats.cumul_FPREL -= time();
  for( unsigned int ic=0; ic<_nc; ic++ ){
    std::pair<FPVT,PESTRUCT::t_CTR> pctr = PE().CTR( ic, _pFPV, FPVxk, _data.ns );
    switch( pctr.second ){
      case PESTRUCT::EQ: _pFPR->add_constraint( pctr.first, FPRT::EQ, 0. ); break;
      case PESTRUCT::LE: _pFPR->add_constraint( pctr.first, FPRT::LE, 0. ); break;
      default:           _pFPR->add_constraint( pctr.first, FPRT::GE, 0. ); break;
    }
  }

  // Output constraints definition
  typename std::list< Data >::const_iterator cit = _data.Y->begin();
  for( ; cit != _data.Y->end(); ++cit ){
    FPVT pyk = PE().OUT( (*cit).index, _pFPV, FPVxk[(*cit).stage], (*cit).stage );
    _pFPR->add_constraint( pyk, FPRT::LE, (*cit).upper );
    _pFPR->add_constraint( pyk, FPRT::GE, (*cit).lower );
  }
  stats.cumul_FPREL += time();

  // Node processing: Constraint propagation
  stats.cumul_REDUC -= time();
  if ( options.USE_CONSTRAINT_PROPAGATION && !_pFPR->propagate_constraints() ){
    stats.cumul_REDUC += time();
    return SetInv<T>::OUTER;
  }
  stats.cumul_REDUC += time();

  // Polyhedral relaxations
  stats.cumul_FPREL -= time();
  _pFPR->generate_cuts();
  if( options.USE_RLT_TMODEL ) _pFPR->generate_RLT_TModel( _pTMFP );
  stats.cumul_FPREL += time();

  return SetInv<T>::UNDETERMINED;
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::_domain_reduction
( FPVT*const*FPVxk, T*P, double&maxred )
{
  maxred = 0.;
  for( unsigned int ip=0; ip<_np; ip++ ){
    T P0 = P[ip];
    for( unsigned int ipb=0; ipb<2; ipb++ ){
      _pFPR->set_objective( (ipb? FPRT::MIN: FPRT::MAX), _pFPV[ip] );
#ifdef MC__ODEGPE_DEBUG
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
            T(-PESTRUCT::INF,_pFPR->get_objective()) ) ){
            P[ip] = Op<T>::l( P[ip] );
            //for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
            //throw Exceptions( Exceptions::REDUC );
            //return SetInv<T>::INFEASIBLE;
          }
          break;
        case 1: // range lower bound
          if( !Op<T>::inter(  P[ip], P[ip],
              T(_pFPR->get_objective(),PESTRUCT::INF) ) ){
            P[ip] = Op<T>::u( P[ip] );
            //for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
            //throw Exceptions( Exceptions::REDUC );
            //return SetInv<T>::INFEASIBLE;
          }
          break;
        }
        _pFPR->update_bounds( _pFPV[ip], P[ip] );
#ifdef MC__ODEGPE_DEBUG
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
        return SetInv<T>::OUTER;
      default:
#ifdef MC__ODEGPE_SHOW_FAILURE
        std::cout << "FAILURE IN DOMAIN REDUCTION LP\n"; 
#endif
        return SetInv<T>::FAILURE;
      }
    }
    maxred = std::max( 1.-Op<T>::diam(P[ip])/Op<T>::diam(P0), maxred );
  }

  return SetInv<T>::UNDETERMINED;
}

template <typename T, typename PE>
inline typename SetInv<T>::STATUS
ODEGPE<T,PE>::assess
( const unsigned int np, T*P )
{
  // Allocate variables for states in polyhedral relaxation
  stats.cumul_FPREL -= time();
  FPVT* FPVxk[_data.ns];
  for( unsigned int is=0; is<=_data.ns; is++ )
    FPVxk[is] = new FPVT[_nx];
  stats.cumul_FPREL += time();

  // Main loop for relaxation and domain reduction
  for( unsigned int ired=0; !ired || ired < options.DOMAIN_REDUCTION_MAXLOOPS; ired++ ){

    // State relaxation, inclusion test, and polyhedral relaxation
    typename SetInv<T>::STATUS status;
    try{
      switch( options.IVP_BOUNDING_TYPE ){
      // State polyhedral relaxations via Intervals
      case Options::IA:
        status = _bound_IVP_IA( P, FPVxk, options.IVP_BOUNDING_PROPAGATION ); break;
      // State polyhedral relaxations via Taylor models
      case Options::TM:
        status = _bound_IVP_TM( P, FPVxk, options.IVP_BOUNDING_PROPAGATION ); break;
      default:
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        throw Exceptions( Exceptions::UNDEF );
      }

     // Return (i) on failure; (ii) if status is established;
     // or (iii) if domain reduction is not selected
     if( status != SetInv<T>::UNDETERMINED || !options.USE_DOMAIN_REDUCTION ){
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        return status;
      }
    }
    catch(...){
      return SetInv<T>::FAILURE;
    }

    // Add user and output constraints to polyhedral relaxation
    status = _relax_outputs_constraints( FPVxk );
    if( status != SetInv<T>::UNDETERMINED ){
      for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
      return status;
    }

    // Node processing: Domain reduction
    if( options.USE_DOMAIN_REDUCTION && ired < options.DOMAIN_REDUCTION_MAXLOOPS ){
      double maxred;
      stats.cumul_REDUC -= time();
      status = _domain_reduction( FPVxk, P, maxred );
      stats.cumul_REDUC += time();
      if( status != SetInv<T>::UNDETERMINED ){
        for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
        return status;
      }

      // Reconstruct lower bounding problem when sufficient tightening was achieved
#ifdef MC__ODEGPE_SHOW_REDUCTION
      std::cerr << "max. domain reduction: " << maxred << std::endl;
#endif
      if( maxred > options.DOMAIN_REDUCTION_THRESHOLD
       && ired+1 < options.DOMAIN_REDUCTION_MAXLOOPS ) continue;
    }

    // Clean-up
    for( unsigned int is=0; is<=_data.ns; is++ ) delete[] FPVxk[is];
    return status;
  }

#ifdef MC__ODEGPE_SHOW_FAILURE
  std::cout << "FAILURE AT EXIT\n"; 
#endif
  return SetInv<T>::FAILURE;
}

template <typename T, typename PE>
inline void
ODEGPE<T,PE>::Options::display
( std::ostream&out ) const
{
  // Display ODEGPE Options
  out << std::setw(60) << "  BOUNDER TYPE FOR PARAMETRIC IVP";
  switch( IVP_BOUNDING_TYPE ){
  case IA:   out << "IA\n";   break;
  case TM:   out << "TM\n";   break;
  }
  out << std::setw(60) << "  BOUNDING STRATEGY FOR PARAMETRIC IVP";
  switch( IVP_BOUNDING_PROPAGATION ){
  case DINEQ:   out << "DINEQ\n";   break;
  case VERIF:   out << "VERIF\n";   break;
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

template <typename T, typename PE>
inline std::ostream&
operator <<
( std::ostream&out, const ODEGPE<T,PE>&pPE )
{
  out << std::endl
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ')
      << std::setw(52) << "GUARANTEED PARAMETER ESTIMATION IN MC++\n"
      << std::setfill('_') << std::setw(72) << "\n\n" << std::setfill(' ');

  // Display ODEGPE Options
  out << std::left << "SET INVERSION ALGORITHM OPTIONS:\n\n";
  pPE.SetInv<T>::options.display( out );
  pPE.ODEGPE<T,PE>::options.display( out );

  out << std::setfill('_') << std::setw(72) << " " << std::setfill(' ')
      << std::endl << std::endl;
  return out;
}

} // end namescape mc

#endif
