// Copyright (C) 2012, 2014 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

/*!
\page page_FPRELAX Polyhedral Relaxations for Factorable Programs
\author Benoit C. Chachuat
\version 0.1
\date 2011
\bug No known bugs.

\section sec_FPREL_bkg What Is Meant by Polyhedral Relaxations?

Consider the mathematical program \f$\mathcal{P}:\min\{f({\bf x}): {\bf x}^{\rm L}\leq{\bf x}\leq{\bf x}^{\rm L}, g_j({\bf x})\leq 0, j=1,\ldots,n_g\}\f$, where \f$f\f$ and \f$g_j, j=1,\ldots,n_g\f$, are factorable, potentially nonconvex, real-valued functions. A mathematical program \f$\mathcal{R}\f$ is said to be a <I>convex relaxation</I> of \f$\mathcal{P}\f$ if (i) it is convex, and (ii) its optimal solution value underestimates the optimal solution value of \f$\mathcal{P}\f$. Moreover, if all the functions (objective and constraints) in \f$\mathcal{R}\f$ are affine, then \f$\mathcal{R}\f$ is said to be a <I>polyhedral relaxation</I> of \f$\mathcal{P}\f$.

The basic procedure to generate a polyhedral relaxation of \f$\mathcal{P}\f$ follows a three-step procedure: (i) reformulate \f$\mathcal{P}\f$ by introducing intermediate variables and extra constraints so that all the constraints in the equivalent program \f$\widetilde{\mathcal{P}}\f$ are either unary and binary operations; (ii) generate convex underestimators and concave overestimators for every constraint in \f$\widetilde{\mathcal{P}}\f$; and (iii) outer-approximate the nonlinear under- or over-estimators at given, well-chosen points. Each step is detailed further below.

- <B>Step 1. Decomposition</B>
\f{align*}
\mathcal{P}:\min_{\bf x}\ & f(\bf x) & \xrightarrow{\displaystyle\text{decomp.}} && \widetilde{\mathcal{P}}:\min_{\bf v}\ & v_{\rm obj} \\
\text{s.t.}\ & {\bf g}(\bf x) \leq {\bf 0} &&& \text{s.t.}\ & {\bf A}\,{\bf v} = {\bf b} \\
& {\bf x}^{\rm L}\leq {\bf x} \leq {\bf x}^{\rm U} &&& & v_k = v_iv_j, \quad \forall (i,j,k)\in\mathcal{B}\\
& &&& & v_k = \frac{v_i}{v_j}, \quad \forall (i,j,k)\in\mathcal{F}\\
& &&& & v_k = \varphi(v_i), \quad \forall (i,k)\in\mathcal{U}\\
& &&& & {\bf v}^{\rm L}\leq {\bf v} \leq {\bf v}^{\rm U}
\f}
The advantage of this decomposition is that it can be applied to any factorable program, but also that it detects and accounts for common subexpressions, which leads to tighter relaxations. Its main drawback, on the other hand, is that it may introduce a large number of extra variables and constraints.

- <B>Step 2. Relaxation</B>
  - The bilinear terms \f$v_k = v_iv_j\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, \f$v_j^{\rm L}\leq v_j\leq v_j^{\rm U}\f$, can be replaced by their polyhedral envelopes as
\f{align*}
v_k=v_iv_j \quad \xrightarrow{\displaystyle\text{relax.}}\quad \left\{\begin{array}{l}
v_k \geq v_i^{\rm L}v_j+v_j^{\rm L}v_i-v_i^{\rm L}v_j^{\rm L}\\
v_k \geq v_i^{\rm U}v_j+v_j^{\rm U}v_i-v_i^{\rm U}v_j^{\rm U}\\
v_k \leq v_i^{\rm U}v_j+v_j^{\rm L}v_i-v_i^{\rm U}v_j^{\rm L}\\
v_k \leq v_i^{\rm L}v_j+v_j^{\rm U}v_i-v_i^{\rm L}v_j^{\rm U}
\end{array}\right.
\f}
  - The fractional terms \f$v_k = \frac{v_i}{v_j}\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, \f$v_j^{\rm L}\leq v_j\leq v_j^{\rm U}\f$, \f$0\notin[v_j^{\rm L},v_j^{\rm U}]\f$, can be rewritten as bilinear terms \f$v_i = v_jv_k\f$ and relaxed as indicated above, with the following bounds for the variables \f$v_k\f$:
\f{align*}
\min\left\{\frac{v_j^{\rm L}}{v_k^{\rm L}}, \frac{v_j^{\rm L}}{v_k^{\rm U}}, \frac{v_j^{\rm U}}{v_k^{\rm L}}, \frac{v_j^{\rm U}}{v_k^{\rm U}}\right\} =: v_k^{\rm L} \leq v_k \leq v_k^{\rm U} := \max\left\{\frac{v_j^{\rm L}}{v_k^{\rm L}}, \frac{v_j^{\rm L}}{v_k^{\rm U}}, \frac{v_j^{\rm U}}{v_k^{\rm L}}, \frac{v_j^{\rm U}}{v_k^{\rm U}}\right\}
\f}
It should be noted, however, that this approach, although being straightforward, does not generally yield the convex/concave envelopes for fractional terms (which turn out to be quite complicated nonlinear expressions -- see \ref sec_FPREL_refs).
\n
  - The univariate terms \f$v_k = \varphi(v_i)\f$, \f$v_i^{\rm L}\leq v_i\leq v_i^{\rm U}\f$, are relaxed differently depending on whether the function \f$\varphi\f$ is convex, concave, convexo-concave, etc., on \f$[v_i^{\rm L},v_i^{\rm U}]\f$.
\f{align*}
\text{convex case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l} v_k \geq \varphi(v_i)\\ v_k \leq \varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_i^{\rm L})}{v_i^{\rm U}-v_i^{\rm L}}(v_i-v_i^{\rm L}) \end{array}\right.\\
\text{concave case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l} v_k \geq \varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_i^{\rm L})}{v_i^{\rm U}-v_i^{\rm L}}(v_i-v_i^{\rm L})\\ v_k \leq \varphi(v_i) \end{array}\right.\\
\text{convexo-concave case:} & \quad v_k=\varphi(v_i) \quad \xrightarrow{\displaystyle\text{relax}} \quad \left\{\begin{array}{l}
v_k \geq \left\{\begin{array}{ll}
\varphi(v_i), & \text{if $v_i\leq v_{\rm m}^{\rm cv}$}\\
\varphi(v_i^{\rm U}) + \frac{\varphi(v_i^{\rm U})-\varphi(v_{\rm m}^{\rm cv})}{v_i^{\rm U}-v_{\rm m}^{\rm cv}}(v_i-v_i^{\rm U}), & \text{otherwise}
\end{array}\right.\\
v_k \leq \left\{\begin{array}{ll}
\varphi(v_i), & \text{if $v_i\geq v_{\rm m}^{\rm cc}$}\\
\varphi(v_i^{\rm L}) + \frac{\varphi(v_i^{\rm L})-\varphi(v_{\rm m}^{\rm cc})}{v_i^{\rm L}-v_{\rm m}^{\rm cc}}(v_i-v_i^{\rm L}), & \text{otherwise}
\end{array}\right.
\end{array}\right.\\
& \quad \text{with:}\ v_{\rm m}^{\rm cv},v_{\rm m}^{\rm cc}\in[v_i^{\rm L},v_i^{\rm U}]:\ \varphi'(v_{\rm m}^{\rm cv}) = \textstyle\frac{\varphi(v_i^{\rm U})-\varphi(v_{\rm m}^{\rm cv})}{v_i^{\rm U}-v_{\rm m}^{\rm cv}} \text{ and } \varphi'(v_{\rm m}^{\rm cc}) = \textstyle\frac{\varphi(v_i^{\rm L})-\varphi(v_{\rm m}^{\rm cc})}{v_i^{\rm L}-v_{\rm m}^{\rm cc}}
\f}
\image html exm_uni.png
\n
  .

- <B>Step 3. Polyhedral Outer-Approximation</B> Every convex, nonlinear univariate constraint generated in Step 2 is outer-approximated by constructing supporting cuts at a number of well-chosen points. Although the resulting polyhedral relaxations are inherently weaker than the nonlinear relaxations, LP solvers are currently more robust and faster than NLP solvers.\n
An iterative scheme (a.k.a. sandwich algorithm) can be applied that adds linearization points in such a way that the maximum distance \f$\delta^{\rm max}\f$ between the nonlinear constraint and its polyhedral approximation decreases as the inverse of the square of the number \f$n\f$ of linearization points; that is, \f$\delta^{\rm max}\propto \frac{1}{n^2}\f$. This algorithm proceeds as follows:
  -# Construct cuts at both interval end-points \f$v_i^{\rm L}\f$ and \f$v_i^{\rm U}\f$
  -# <B>REPEAT</B>
    - Identify an interval \f$[v_i^\ell,v_i^{\ell+1}]\f$ with maximum outer-approximation error
    - Subdivide \f$[v_i^\ell,v_i^{\ell+1}]\f$ at a suitably chosen point \f$v_i^{\rm new}\f$
    .
    <B>UNTIL</B> \f$\text{maximum outer-approximation error} < \varepsilon^{\rm tol}\f$
  .
  In particular, diverse strategies have been proposed for selecting a new linearization point \f$v_i^{\rm new}\f$. The <I>interval bisection rule</I> and the <I>maximum error rule</I> are depicted below.

\image html  OAcvx_strategy.png


\section sec_FPREL_impl How Are Polyhedral Relaxations Implemented in MC++?

The class mc::FPRelax in MC++ provides an implementation of the abovementionned three-step procedure for the construction polyhedral relaxations for general factorable programs. Both the classes mc::FPRCplex (header file <tt>fprcplex.h</tt>) and mc::FPRGurobi (header file <tt>fprgurobi.h</tt>) are derived from mc::FPRelax and solve the resulting relaxed programs using, repectively <A href="http://www-01.ibm.com/software/integration/optimization/cplex-optimizer/">CPLEX</A> and <A href="http://www.gurobi.com/">GUROBI</A> as the optimizer, in addition to generate the relaxations.

In constructing the polyhedral relaxations, interval bounds must be calculated for the intermediate factors (the results of both unary and binary operations). This is done in mc::FPRelax via standard interval arithmetic, and the template parameter of mc::FPRelax is used precisely to specify the desired interval type. This makes it convenient for the user to switch between different interval types. For convenience, MC++ comes with a default interval type named mc::Interval (header file <tt>interval.h</tt>) that can be used inside mc::FPRelax, although it should be noted that mc::Interval is not a validated implementation. For validated interval analysis, it is recommended to use third-party libraries such as <A href="http://www.ti3.tu-harburg.de/Software/PROFILEnglisch.html">PROFIL</A> (header file <tt>mcprofil.h</tt>) or <A href="http://www.math.uni-wuppertal.de/~xsc/software/filib.html">FILIB++</A> (header file <tt>mcfilib.h</tt>). The selection of alternative interval types is also possible, yet it involves instantiating the templated structure mc::Op for the corresponding types.

Also note that mc::FPRelax is <b>not</b> a verified implementation per se. Indeed, round-off errors are not accounted during the construction of the polyhedral relaxations (no matter whether or not validated interval types are used to propagate bound for the intermediate factors).

\section sec_FPREL_val How to Generate and Solve a Polyhedral Relaxation of a Factorable Program?

For illustration, suppose we want to construct a polyhedral relaxation for the following NLP and solve it in order to obtain a lower bound on the global solution:
\f{align*}
  \max_{\bf x}\ & x_1 + \exp(1-x_2) \\
  \text{s.t.} \ & x_1\,\exp(1-x_2) \leq 4 \\
  & 3 \leq x_1 \leq 6 \\
  & 0 \leq x_2 \leq 4.
\f}
The optimizer GUROBI is considered subsequently, along with the default interval type mc::Interval for simplicity:

\code
      #include "interval.h"
      typedef mc::Interval I;

      #include "fprgurobi.h"
      typedef mc::FPRGurobi<I> FPREL;
      typedef mc::FPVar<I> FPVAR;
\endcode

First, the factorable program along with the two participating variables (indexes, ranges) are defined as follows:

\code
      FPREL FP;
      FPVAR X1( &FP, 1, I(3.,6.) );
      FPVAR X2( &FP, 2, I(0.,4.) );
\endcode

Next, the objective function and inequality constraint are defined:

\code
      FP.set_objective( FPREL::MAX, X1+exp(1-X2) );
      FP.add_constraint( X1*exp(1-X2), FPREL::LE, 4. );
      std::cout << FP;
\endcode

The last line displays the following information about the factorable program:

\verbatim
    VARIABLES:
      X1 <= [  3.00000e+00 :  6.00000e+00 ]         
      X2 <= [  0.00000e+00 :  4.00000e+00 ]         

    AUXILIARY:
      Z1 <= [ -3.00000e+00 :  1.00000e+00 ]                 Z0 - X2
      Z2 <= [  4.97871e-02 :  2.71828e+00 ]                 EXP( Z1 )
      Z3 <= [  3.04979e+00 :  8.71828e+00 ]                 X1 + Z2
      Z4 <= [  1.49361e-01 :  1.63097e+01 ]                 X1 * Z2
      Z0 <= (I) 1                                           CONSTANT
      Z5 <= (D) 4                                           CONSTANT

    NO CUTS
\endverbatim

It is seen that 6 auxiliary variables are introduced during the decomposition (step 1), \f$z_0,\ldots,z_5\f$, which correspond to the various unary and binary operations in the objective and contraint functions, as well as the (integer and real) constants. Observe, in particular, that the common sub-expression \f$\exp(1-x_2)\f$ is detected. Bounds are also computed for these auxiliary variables, as indicated.

At this stage, however, the cuts in the polyhedral relaxation (steps 2 and 3) are yet to be generated. This is done as follows:

\code
      FP.generate_cuts();
      std::cout << FP;
\endcode

The following information is now displayed about the factorable program:

\verbatim
    VARIABLES:
      X1 <= [  3.00000e+00 :  6.00000e+00 ]         
      X2 <= [  0.00000e+00 :  4.00000e+00 ]         
    
    AUXILIARY:
      Z1 <= [ -3.00000e+00 :  1.00000e+00 ]                 Z0 - X2
      Z2 <= [  4.97871e-02 :  2.71828e+00 ]                 EXP( Z1 )
      Z3 <= [  3.04979e+00 :  8.71828e+00 ]                 X1 + Z2
      Z4 <= [  1.49361e-01 :  1.63097e+01 ]                 X1 * Z2
      Z0 <= (I) 1                                           CONSTANT
      Z5 <= (D) 4                                           CONSTANT
    
    CUTS:
      + 1.00000e+00Z3 - 1.00000e+00X1 - 1.00000e+00Z2 = 0.00000e+00
      + 1.00000e+00Z1 + 1.00000e+00X2 = 1.00000e+00
      + 1.00000e+00Z4 <= 4.00000e+00
      + 1.00000e+00Z4 - 4.97871e-02X1 - 6.00000e+00Z2 <= -2.98722e-01
      + 1.00000e+00Z4 - 2.71828e+00X1 - 3.00000e+00Z2 <= -8.15485e+00
      - 2.00855e+01Z2 + 1.00000e+00Z1 <= -4.00000e+00
      - 3.67879e-01Z2 + 1.00000e+00Z1 <= 0.00000e+00
      - 9.28087e-01Z2 + 1.00000e+00Z1 <= -9.25371e-01
      - 2.17368e+00Z2 + 1.00000e+00Z1 <= -1.77642e+00
      - 5.44615e-01Z2 + 1.00000e+00Z1 <= -3.92324e-01
      + 1.00000e+00Z4 - 2.71828e+00X1 - 6.00000e+00Z2 >= -1.63097e+01
      + 1.00000e+00Z4 - 4.97871e-02X1 - 3.00000e+00Z2 >= -1.49361e-01
      - 4.00000e+00Z2 + 2.66849e+00Z1 >= -8.20463e+00
\endverbatim

By default, the cuts in the polyhedral relaxation of a convex univariate terms are generated according to the <I>maximum error rule</I> (see above), and a maximum of 5 cuts are generated for each nonlinear constraint. The cut generation is also controlled by the absolute and relative tolerances on the maximum outer-approximation error, which are both set to \f$10^{-3}\f$ by default. All these default values can be altered via the options mechanism of mc::FPRelax, as explained below in the section \ref sec_FPREL_opt below.

Having generated the cuts, a solution to the resulting polyhedral relaxation can be computed as:

\code
      FP.solve();
      std::cout << "\n  RELAXATION VALUE:   " << FP.get_objective()
                << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                              << FP.get_variable(2)
                << "\n  CPU TIME [sec]:     " << FP.get_runtime()
                << std::endl;
\endcode

The following results are obtained:

\verbatim
    RELAXATION VALUE:   6.66667e+00
    RELAXATION OPTIMUM: 6.00000e+00, 3.07531e+00
    CPU TIME [sec]:     1.20902e-03
\endverbatim

Note that the solution of the polyhedral relaxation corresponds to the global optimum for this simple problem.


Prior to generating the polyhedral relaxation, <b>constraint propagation</b> can be applied as a preprocessing step to tighten the variable bounds. This is done as follows:

\code
      FP.propagate_constraints();
      FP.generate_cuts( true );
      std::cout << FP;
\endcode

where the argument `<a>true</a>' passed to the method <a href="classmc_1_1FPRelax.html#a90277409a6e357ca8ff819a3f9201289">mc::FPRelax::generate_cuts</a> indicates that the polyhedral relaxations are to be rebuilt. The following result is displayed:

\verbatim
    VARIABLES:
      X1 <= [  3.00000e+00 :  6.00000e+00 ]         
      X2 <= [  7.12318e-01 :  4.00000e+00 ]         

    AUXILIARY:
      Z1 <= [ -3.00000e+00 :  2.87682e-01 ]                 Z0 - X2
      Z2 <= [  4.97871e-02 :  1.33333e+00 ]                 EXP( Z1 )
      Z3 <= [  3.04979e+00 :  8.71828e+00 ]                 X1 + Z2
      Z4 <= [  1.49361e-01 :  4.00000e+00 ]                 X1 * Z2
      Z0 <= (I) 1                                           CONSTANT
      Z5 <= (D) 4                                           CONSTANT

    CUTS:
      + 1.00000e+00Z3 - 1.00000e+00X1 - 1.00000e+00Z2 = 0.00000e+00
      + 1.00000e+00Z1 + 1.00000e+00X2 = 1.00000e+00
      + 1.00000e+00Z4 <= 4.00000e+00
      + 1.00000e+00Z4 - 4.97871e-02X1 - 6.00000e+00Z2 <= -2.98722e-01
      + 1.00000e+00Z4 - 1.33333e+00X1 - 3.00000e+00Z2 <= -4.00000e+00
      - 2.00855e+01Z2 + 1.00000e+00Z1 <= -4.00000e+00
      - 7.50000e-01Z2 + 1.00000e+00Z1 <= -7.12318e-01
      - 1.79462e+00Z2 + 1.00000e+00Z1 <= -1.58479e+00
      - 3.84904e+00Z2 + 1.00000e+00Z1 <= -2.34782e+00
      - 1.08971e+00Z2 + 1.00000e+00Z1 <= -1.08591e+00
      + 1.00000e+00Z4 - 1.33333e+00X1 - 6.00000e+00Z2 >= -8.00000e+00
      + 1.00000e+00Z4 - 4.97871e-02X1 - 3.00000e+00Z2 >= -1.49361e-01
      - 3.28768e+00Z2 + 1.28355e+00Z1 >= -4.01432e+00
\endverbatim

Observe, in particular, that constraint propagation could tighten the bounds for variable <a>X2</a>, which results in tighter polyhedral relaxations in turn.

\section sec_FPREL_link How to Incorporate McCormick, Taylor and McCormick-Taylor Models in the Polyhedral Relaxation of a Factorable Program?

Polyhedral relaxations can also be constructed for optimization problems that contain terms for which McCormick, Taylor or McCormick-Taylor models are known. One special case where this is particularly useful is for the polyhedral relaxation of dynamic systems described by differential equations. 

\subsection ssec_FPREL_link_TM Incorporation of Taylor and McCormick-Taylor Models in a Polyhedral Relaxation

For illustration, suppose we want to construct and solve a polyhedral relaxation for the same NLP as before:
\f{align*}
  \max_{\bf x}\ & x_1 + \exp(1-x_2) \\
  \text{s.t.} \ & x_1\,\exp(1-x_2) \leq 4 \\
  & 3 \leq x_1 \leq 6 \\
  & 0 \leq x_2 \leq 4,
\f}
but with the use of a McCormick-Taylor model of order 4 for the subexpression \f$\exp(1-x_2)\f$ (see Page \ref page_TAYLOR for a description of how to compute Taylor Model for factorable functions). An implementation of this problem is the following:

\code
      #include "interval.h"
      #include "tmodel.h"
      #include "mccormick.h"
      #include "fprgurobi.h"
      typedef mc::Interval I;
      typedef mc::McCormick<I> MC;
      typedef mc::FPRGurobi<I> FPREL;
      typedef mc::FPVar<I> FPVAR;
      typedef mc::TModel<MC> TMMC;
      typedef mc::TVar<MC> TVMC;
      typedef mc::TModel<FPVAR> TMFPVAR;

      FPREL FP;
      FPVAR X1( &FP, 1, I(3.,6.) );
      FPVAR X2( &FP, 2, I(0.,4.) );

      TMMC TMod( 1, 4 );
      double X2_Ref[1] = { 2. };
      MC X2_MC( I(0.,4.), X2_Ref[0] );
      X2_MC.sub( 1, 0 );
      TVMC X2_TMMC( &TMod, 0, X2_MC );
      unsigned int TMod2FP[1] = { 2 };
      TVMC X3_TMMC = exp( 1 - X2_TMMC );
      TMFPVAR*TModFP = FP.create_TModel( &TMod, TMod2FP );
      FPVAR X3( &FP, X3_TMMC, TModFP, X2_Ref, TMod2FP );
      delete TModFP;
\endcode

The last line initializes the variable <A>X3</A> as a polyhedral enclosure of the McCormick-Taylor model <A>X3_TMMC</A>. Note that this constructor also requires to pass:
- (i) a pointer [<A>TModFP</A>] to a mc::TModel<mc::FPVar> object corresponding to the actual McCormick-taylor model mc::TModel<mc::McCormick>. This pointer can be created with the member function TModFP in mc::FPRelax, can be reused to append polyhedral relaxations of multiple McCormick-Taylor models in the same variables (e.g., for the solutions of parametric ODE systems), and must be deleted before coming out of scope to avoid memory leaks;
- (ii) a pointer [<A>X2_Ref</A>] to the point at which the McCormick relaxation of the remainder term was computed; and
- (iii) a pointer [<A>TMod2FP</A>] to an array matching the indexes of the variable that participate in the McCormick relaxation of the remainder term (subgradient components) to those of the factorable program.
.

At this point, the definition of the objective function and inequality constraint, the construction of the polyhedral relaxation, and its solution are the same as previously:

\code
      FP.set_objective( FPREL::MAX, X1+X3 );
      FP.add_constraint( X1*X3, FPREL::LE, 4. );
      FP.generate_cuts();
      FP.solve();
      std::cout << "\n  RELAXATION VALUE:   " << FP.get_objective()
                << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                              << FP.get_variable(2)
                << "\n  CPU TIME [sec]:     " << FP.get_runtime()
                << std::endl;
      std::cout << FP;
\endcode

The following results are obtained:

\verbatim
    RELAXATION VALUE:   6.66667e+00
    RELAXATION OPTIMUM: 6.00000e+00, 6.47938e-01
    CPU TIME [sec]:     1.82295e-03
\endverbatim

\subsection ssec_FPREL_link_MC Incorporation of McCormick Relaxations in a Polyhedral Relaxation

Suppose instead that a McCormick relaxation is known for the subexpression \f$\exp(1-x_2)\f$ (see Page \ref page_MCCORMICK for a description of how to compute McCormick relaxations for factorable functions). The implementation is changed as follows:

\code
      #include "interval.h"
      #include "mccormick.h"
      #include "fprgurobi.h"
      typedef mc::Interval I;
      typedef mc::McCormick<I> MC;
      typedef mc::FPRGurobi<I> FPREL;
      typedef mc::FPVar<I> FPVAR;

      FPREL FP;
      FPVAR X1( &FP, 1, I(3.,6.) );
      FPVAR X2( &FP, 2, I(0.,4.) );

      double X2_REF[1] = { 2. };
      MC X2_MC( I(0.,4.), X2_Ref[0] );
      X2_MC.sub( 1, 0 );
      unsigned int TMod2FP[1] = { 2 };
      MC X3_MC = exp(1-X2_MC);
      FPVAR X3( &FP, X3_MC, X2_Ref, TMod2FP );

      FP.set_objective( FPREL::MAX, X1+X3 );
      FP.add_constraint( X1*X3, FPREL::LE, 4. );
      FP.generate_cuts();
      FP.solve();
      std::cout << "\n  RELAXATION VALUE:   " << FP.get_objective()
                << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                              << FP.get_variable(2)
                << "\n  CPU TIME [sec]:     " << FP.get_runtime()
                << std::endl;
      std::cout << FP;
\endcode

In this case, the variable <A>X3</A> is initialized with a polyhedral enclosure of the McCormick relaxation <A>X3_MC</A>, as computed from the available subgradient information. As with McCormick-Taylor models, the corresponding constructor requires to pass: (i) a pointer to the point at which the McCormick relaxation was computed [<A>X2_REF</A>]; and (ii) a pointer to an array matching the indexes of the variable that participate in the McCormick relaxation to those of the factorable program [<A>X2_TMMC</A>].

The following results are obtained:

\verbatim
  RELAXATION VALUE:   6.66667e+00
  RELAXATION OPTIMUM: 6.00000e+00, 3.07531e+00
  CPU TIME [sec]:     1.19686e-03
\endverbatim

\section sec_FPREL_opt What are the Options Relative to the Generation of a Polyhedral Relaxation?

All the options are defined in the structure mc::FPRelax::Options. The current options are:

<TABLE border="1">
<CAPTION><EM>Options in mc::FPRelax::Options: name, type and
description</EM></CAPTION>
     <TR><TH><b>Name</b>  <TD><b>Type</b><TD><b>Default</b>
         <TD><b>Description</b>
     <TR><TH><tt>SANDWICH_ATOL</tt> <TD><tt>double</tt> <TD>1e-3
         <TD>Absolute tolerance in the outer approximation of convex/concave univariate terms using the sandwich algorithm
     <TR><TH><tt>SANDWICH_RTOL</tt> <TD><tt>double</tt> <TD>1e-3
         <TD>Relative tolerance in the outer approximation of convex/concave univariate terms using the sandwich algorithm
     <TR><TH><tt>SANDWICH_MAXCUT</tt> <TD><tt>integer</tt> <TD>5
         <TD>Maximum number of cuts in the outer approximation of convex/concave univariate terms using the sandwich algorithm
     <TR><TH><tt>SANDWICH_RULE</tt> <TD><tt>mc::FPRelax::Options::SANDWICH</tt> <TD>mc::FPRelax::Options::MAXERR
         <TD>Strategy for outer approximation of convex/concave univariate terms using the sandwich algorithm: bisection (mc::FPRelax::Options::BISECT); maximum error rule (mc::FPRelax::Options::MAXERR)
     <TR><TH><tt>NEWTON_USE</tt> <TD><tt>bool</tt> <TD>true
         <TD>Whether to use the Newton/secant method to calculate convex/concave envelopes for concavo-convex univariate functions (including odd monomial terms, sine, cosine, erf, erfc, etc.)
     <TR><TH><tt>NEWTON_TOL</tt> <TD><tt>double</tt> <TD>1e-10
         <TD>Termination tolerance in the Newton/secant method
     <TR><TH><tt>NEWTON_MAXIT</tt> <TD><tt>int</tt> <TD>100
         <TD>Maximum number of iterations in the Newton/secant method
     <TR><TH><tt>RRLT_STRATEGY</tt> <TD><tt>mc::FPRelax::Options::RRLT</tt> <TD>mc::FPRelax::Options::PRIMRRLT
         <TD>Strategy for reduced (R)RLT constraints; no cuts (mc::FPRelax::Options::NORRLT); cuts involved primary variables as multipliers only (mc::FPRelax::Options::PRIMRRLT); cuts involved primary and auxiliary variables as multipliers (mc::FPRelax::Options::ALLRRLT)
     <TR><TH><tt>RRLT_DISPLAY</tt> <TD><tt>unsigned int</tt> <TD>0
         <TD>Display level during reduced (R)RLT constraint generation; 0-no display; 1-RLT result only; >1-full display
</TABLE>

The mc::FPRelax class has a public static member called mc::FPRelax::options that can be used to set/modify the options; for example:

\code
      FP.options.SANDWICH_ATOL = 1e-3;
      FP.options.SANDWICH_RTOL = 1e-3;
      FP.options.SANDWICH_MAXCUT = 5;
      FP.options.NEWTON_USE = true;
      FP.options.NEWTON_TOL = 1e-12;
      FP.options.NEWTON_MAXIT = 100;
      FP.options.SANDWICH_RULE = mc::FPRelax::Options::MAXERR;
      FP.options.RRLT_STRATEGY = mc::FPRelax::Options::PRIMRRLT;
\endcode

\section sec_FPREL_err What Errors Can Be Encountered during the Calculation of an LP Relaxation of a Factorable Program?

Errors are managed based on the exception handling mechanism of the C++ language. Each time an error is encountered, an instance of the class mc::FPRelax::Exceptions is thrown, which contains the type of error. It is the user's responsibility to test whether an exception was thrown during an LP relaxation, and then make changes where appropriate. Should an exception be thrown and not caught by the calling program, execution will abort.

Possible errors encountered during the calculation of an LP relaxation are:

<TABLE border="1">
<CAPTION><EM>Errors during Calculation of an LP relaxation</EM></CAPTION>
     <TR><TH><b>Number</b> <TD><b>Description</b>
     <TR><TH><tt>1</tt> <TD>Error during the calculation of the convex or concave envelope of a univariate functions (Newton/secant method failed)
     <TR><TH><tt>2</tt> <TD>Error during intersection of two terms (empty intersection)
     <TR><TH><tt>3</tt> <TD>Error in constraint definition (e.g., number = number)
     <TR><TH><tt>4</tt> <TD>Error in objective definition (e.g., minimize interval)
     <TR><TH><tt>5</tt> <TD>Error during initialization with a variable of class mc::McCormick or mc::Taylor
     <TR><TH><tt>-1</tt> <TD>Error due to an operation between variables participating in different factorable programs
     <TR><TH><tt>-2</tt> <TD>Error due to calling a function/feature not yet implemented in MC++
     <TR><TH><tt>-3</tt> <TD>Internal error
</TABLE>

\section sec_FPREL_refs References

- Tawarmalani, M., and N.V. Sahinidis, <A href="http://dx.doi.org/10.1007/s10107-003-0467-6">Global optimization of mixed-integer nonlinear programs: A theoretical and computational study</A>, <I>Mathematical Programming</I>, <B>99</B>(3):563-591, 2004.
- Smith, E.M.B, and C.C. Pantelides, <A href="http://dx.doi.org/10.1016/S0098-1354(98)00286-5">A symbolic reformulation/spatial branch-and-bound algorithm for the global optimisation of nonconvex MINLPs</A>, <I>Computers & Chemical Engineering</I>, <B>23</B>(4-5):457-478, 1999.
.
*/

// TO DO:
// - [DONE] Factor out the sandwich alorithm!
// - Implement remaining simple univariate functions: sqr [DONE], sqrt [DONE], log [DONE], fabs [DONE], pow (even) [DONE], fstep [DONE], bstep [DONE], min [DONE], max [DONE]
// - Implement more complex univariate functions with Newton method: pow (odd) [DONE], sin [DONE], cos [DONE], asin [DONE], acos [DONE], tan [DONE], atan [DONE], erf [DONE], erfc [DONE], hyperbolic functions
// - Also implement functions like: bilin [DONE], xlog, arh, monod, haldane
// - [DONE] Implement a header file "lprgurobi.h" for constructing LP relaxations in Gurobi format
// - [DONE] Floating point exception sent by Gurobi w/ Profil! (but not mc::Interval...)
// - [DONE] Implement CPLEX version too

// - [BAD IDEA] Define Cuts in terms of the variables themselves, not variable indexes -> not a good idea because of Outer Approximation
// - [DONE] Factorize variables/operations in operators/functions
// - [DONE] Introduce a _type field for the variables instead of _constant and negative indices, such as { VARCONT, VARBIN, AUXCONT, AUXBIN, AUXINT, AUXREAL }
// - [DONE] Make link with McCormick, Taylor and McCormick-Taylor for variable initialization

// - [DONE] Separate factorable decomposition from cut generation
// - [DONE] Add/remove particular cuts in Gurobi/Cplex
// - [DONE] Keep track of operations <-> cuts
// - Implement trilinear terms relaxation
// - [DONE] Change name from "LP" to "FP" for "Factorable Program"
// - [DONE] Implement constraint propagation

// - [DONE] Write down documentation
// - Check that variables participating in binary operations share the same FPRelax object


#ifndef MC__FPRELAX_HPP
#define MC__FPRELAX_HPP

#include <iostream>
#include <string>
#include <sstream>
#include <stdarg.h>
#include <set>
#include <list>
#include <queue>
#include <utility>
#include <algorithm>

#include "mcfunc.hpp"
#include "mcop.hpp"
#include "mccormick.hpp"
#include "tmodel.hpp"

#define MC__FPRELAX_SCALE_SEPARABLE
#undef  MC__FPRELAX_DEBUG_PROPAGATION

namespace mc
{
template <typename T> class FPOp;
template <typename T> class lt_FPOp;
template <typename T> class FPCut;
template <typename T> class lt_FPCut;
template <typename T> class FPRLT;
template <typename T> class FPRelax;
template <typename T> class FPRRLT;

//! @brief C++ template structure for defining the numeric field of a factorable program variable
////////////////////////////////////////////////////////////////////////
//! mc::FPVarNum<T> is a C++ structure for defining  the numeric field
//! of a variable in a factorable program, which can either a bound (T),
//! a real value (double), or an integer value (int).
////////////////////////////////////////////////////////////////////////
template <typename T>
struct FPVarNum
{
  //! @brief Enumeration type for numeric variables in factorable program
  enum TYPE{
    INT=0,	//!< Integer value
    REAL,	//!< Real value
    RANGE	//!< Range
  };
  //! @brief Variable type
  TYPE t;
  //! @brief Integer/real variable value
  union{
    int n;
    double x;
  }; 
  //! @brief Variable range
  T I;

  //! @brief Constructor for an integer variable
  FPVarNum( const int i=0 ):
    t(INT), n(i), I(i)
    {}
  //! @brief Constructor for a real variable
  FPVarNum( const double d ):
    t(REAL), x(d), I(d)
    {}
  //! @brief Constructor for a variable range
  FPVarNum( const T&B ):
    t(RANGE), n(0), I(B)
    {}

  //! @brief Assignment of an integer variable
  FPVarNum<T>& operator=
    ( const int i )
    { t = INT; n = i; I = i; return *this; }
  //! @brief Assignment of a real variable
  FPVarNum<T>& operator=
    ( const double d )
    { t = REAL; x = d; I = d; return *this; }
  //! @brief Assignment of a variable range
  FPVarNum<T>& operator=
    ( const T&B )
    { t = RANGE; n = 0; I = B; return *this; }
  //! @brief Assignment of another FPVarNum variable
  FPVarNum<T>& operator=
    ( const FPVarNum<T>&num )
    { t = num.t; t==REAL? x=num.x: n; I = num.I; return *this; }
};

//! @brief C++ structure for comparing values of factorable program variables
////////////////////////////////////////////////////////////////////////
//! mc::eq_FPVarNum<T> is a C++ structure for comparing the numeric
//! field of variables in factorable programs.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct eq_FPVarNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPVarNum<T>*Num1, const FPVarNum<T>*Num2 ) const
    {
      if( Num1->t != Num2->t ) return false;
      switch( Num1->t ){
        case FPVarNum<T>::INT:   return Num1->n==Num2->n? true: false;
        case FPVarNum<T>::REAL:  return isequal( Num1->x, Num2->x );
        case FPVarNum<T>::RANGE: return Op<T>::eq( Num1->I, Num2->I );
      }
    }
};

//! @brief C++ structure for comparing values of factorable program variables
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPVarNum<T> is a C++ structure for comparing the numerical
//! field of variables in factorable programs.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_FPVarNum
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPVarNum<T>*Num1, const FPVarNum<T>*Num2 ) const
    {
      if( Num1->t < Num2->t ) return true;
      if( Num1->t > Num2->t ) return false;
      switch( Num1->t ){
        case FPVarNum<T>::INT:
          return Num1->n<Num2->n? true: false;
        case FPVarNum<T>::REAL:
          return !isequal( Num1->x, Num2->x ) && Num1->x<Num2->x? true: false;
        case FPVarNum<T>::RANGE:
          return Op<T>::lt( Num1->I, Num2->I );
      }
      return false;
    }
};

//! @brief C++ template class for defining variables in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::FPVar<T> is a C++ template class for defining variables in the
//! factored form of a factorable program.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPVar
////////////////////////////////////////////////////////////////////////
{
  // friends of this class with other classes/structures/operators
  template <typename U> friend class FPVar;
  //template <typename U> friend FPVar<U>* FPRelax<U>::_auxiliary_variable
  //  ( const U&, FPOp<U>* );
  //template <typename U> friend FPVar<U>* FPRelax<U>::_auxiliary_variable
  //  ( FPOp<U>* );
  template <typename U> friend class FPRelax;
  template <typename U> friend class FPRRLT;
  template <typename U> friend struct lt_FPVar;
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPRelax<U>& );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPOp<U>& );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPCut<U>& );
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPRLT<U>& );

  // friends of this class for operator overloading
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator+
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator+
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator-
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator-
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator*
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator*
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator*
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator*
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator*
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator/
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator/
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator/
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator/
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator/
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator^
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> max
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> max
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> max
    ( const unsigned int, const FPVar<U>* );
  template <typename U> friend FPVar<U> min
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> min
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> min
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> min
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> min
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> min
    ( const unsigned int, const FPVar<U>* );

  // friends of class FPVar for function overloading
  template <typename U> friend FPVar<U> inv
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sqr
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> exp
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> log
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sqrt
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> fabs
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> cos
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sin
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> tan
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> acos
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> asin
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> atan
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> xlog
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> arh
    ( const FPVar<U>&, const double );
  template <typename U> friend FPVar<U> erf
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> erfc
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> fstep
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> bstep
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> pow
    ( const FPVar<U>&, const int );
  template <typename U, typename V> friend FPVar<U> pow
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> pow
    ( const double, const FPVar<U>& );
  template <typename U> friend FPVar<U> ltcond
    ( const FPVar<U>&, const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> ltcond
    ( const U&, const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> gtcond
    ( const FPVar<U>&, const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> gtcond
    ( const U&, const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> bilin
    ( const FPVar<U>&, const FPVar<U>& );

public:

  // other operator overloadings
  FPVar<T>& operator=
    ( const FPVar<T>& );
  FPVar<T>& operator=
    ( const T& );
  FPVar<T>& operator=
    ( const int );
  FPVar<T>& operator=
    ( const double );
  template <typename U> FPVar<T>& operator+=
    ( const U& );
  template <typename U> FPVar<T>& operator-=
    ( const U& );
  template <typename U> FPVar<T>& operator*=
    ( const U& );
  template <typename U> FPVar<T>& operator/=
    ( const U& );

  // operator overloadings (for constraint definition)
  template <typename U> bool operator ==
    ( const U&Var )
    { _FP->add_constraint( *this, FPRelax<T>::EQ, Var ); return true; }
  template <typename U> bool operator <=
    ( const U&Var )
    { _FP->add_constraint( *this, FPRelax<T>::LE, Var ); return true; }
  template <typename U> bool operator >=
    ( const U&Var )
    { _FP->add_constraint( *this, FPRelax<T>::GE, Var ); return true; }

  /** @defgroup FP Polyhedral Relaxations for Factorable Programs
   *  @{
   */
  //! @brief Index for 'free' variables in factorable program
  static const int NOREF = -33;
  //! @brief Enumeration type for variables in factorable program
  enum TYPE{
    VARCONT=0,	//!< Native continuous variable
    VARBIN,	//!< Native binary variable
    AUXCONT,	//!< Auxiliary continuous variable
    AUXBIN,	//!< Auxiliary binary variable   
    AUXINT,	//!< Auxiliary integer constant
    AUXREAL	//!< auxiliary real constant
  };
  //! @brief Typedef for variable identifier in factorable program
  typedef std::pair< TYPE, long > pt_idVar;
  /** @} */

private:
  //! @brief Variable identifier (type and index)
  pt_idVar _id;
  //! @brief Variable numeric value (integer, real or bound)
  FPVarNum<T> _num;
  //! @brief Pointer to underlying factorable program
  FPRelax<T> *_FP;
  //! @brief Pointer to defining operation
  FPOp<T> *_Op;

public:

  /** @ingroup FP
   *  @{
   */
  //! @brief Constructor for the continuous variable <a>ix</a>, with range <a>X</a>, in the factorable program <a>FP</a>
  FPVar
    ( FPRelax<T>*FP, const unsigned long ix, const T&X )
    : _id( VARCONT, ix ), _num( X ), _FP( FP ), _Op( 0 )
    { 
      if( !FP )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );
      // Keep track of variable count
      if( FP->_nvar <= ix ) FP->_nvar = ix+1;
      // Search for existing variable with same identification in _Vars,
      // and erase it from _Vars if found along with the corresponding
      // operation in _Ops
      FPVar<T>* pVar = new FPVar<T>( *this );
      typename FPRelax<T>::it_Vars iVar = FP->_Vars.find( pVar );
      if( iVar!=FP->_Vars.end() ){
        FP->_Ops.erase( (*iVar)->_Op ); delete (*iVar)->_Op;
        delete *iVar; FP->_Vars.erase( iVar );
      }
      // Insert new variable in _Vars and corresponding operation in _Ops
      FPOp<T>* Op = new FPOp<T>( FPOp<T>::VAR, 0, 0, pVar );
      FP->_Ops.insert( Op );
      pVar->_Op = _Op = Op;
      FP->_append_var( pVar );
   }

  //! @brief Set the continuous variable <a>ix</a>, with range <a>X</a>, in the factorable program <a>FP</a>.
  FPVar<T>& set
    ( FPRelax<T>*FP, const unsigned long ix, const T&X )
    { *this = FPVar( FP, ix, X ); return *this; }

  //! @brief Constructor for the binary variable <a>ix</a>, in the factorable program <a>FP</a>
  FPVar
    ( FPRelax<T>*FP, const unsigned long ix )
    : _id( VARBIN, ix ), _num( Op<T>::zeroone() ), _FP( FP ), _Op( 0 )
    { 
      if( !FP )
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );
      // Keep track of variable count
      if( FP->_nvar <= ix ) FP->_nvar = ix+1;
      // Search for existing variable with same identification in _Vars,
      // and erase it from _Vars if found along with the corresponding
      // operation in _Ops
      FPVar<T>* pVar = new FPVar<T>( *this );
      typename FPRelax<T>::it_Vars iVar = FP->_Vars.find( pVar );
      if( iVar!=FP->_Vars.end() ){
        FP->_Ops.erase( (*iVar)->_Op ); delete (*iVar)->_Op;
        delete *iVar; FP->_Vars.erase( iVar );
      }
      // Insert new variable in _Vars and corresponding operation in _Ops
      FPOp<T>* Op = new FPOp<T>( FPOp<T>::VAR, 0, 0, pVar );
      FP->_Ops.insert( Op );
      pVar->_Op = _Op = Op;
      FP->_append_var( pVar );
   }

  //! @brief Set the binary variable <a>ix</a>, in the factorable program <a>FP</a>.
  FPVar<T>& set
    ( FPRelax<T>*FP, const unsigned long ix )
    { *this = FPVar( FP, ix ); return *this; }

  //! @brief Constructor for integer parameter
  FPVar
    ( const int i=0 )
    : _id( AUXINT, NOREF ), _num( i ), _FP( 0 ), _Op( 0 )
    {}

  //! @brief Constructor for real parameter
  FPVar
    ( const double d )
    : _id( AUXREAL, NOREF ), _num( d ), _FP( 0 ), _Op( 0 )
    {}

  //! @brief Constructor for range
  FPVar
    ( const T&X )
    : _id( AUXCONT, NOREF ), _num( X ), _FP( 0 ), _Op( 0 )
    {}

  //! @brief Constructor for McCormick convex/concave bounds
  FPVar
    ( FPRelax<T>*FP, const mc::McCormick<T>&MCX, const double*xref,
      const unsigned int*isub=0 );

  //! @brief Constructor for Taylor models
  FPVar
    ( FPRelax<T>*FP, const mc::TVar<T>&TVX, mc::TModel< FPVar<T> >*TMFP );

  //! @brief Constructor for McCormick-Taylor models
  FPVar
    ( FPRelax<T>*FP, const mc::TVar< mc::McCormick<T> >&TVMCX,
      mc::TModel< FPVar<T> >*TMFP, const double*MCXref,
      const unsigned int*isub );

  //! @brief Copy constructor
  FPVar
    ( const FPVar<T>&Var )
    : _id( Var._id ), _num( Var._num ), _FP( Var._FP ), _Op( Var._Op )
    {}
  /** @} */

private:

  //! @brief Constructor for an intermediate continuous variable, with bounds <a>bounds</a>, in the factorable program <a>FP</a>
  FPVar
    ( FPRelax<T>*FP, const T&X, FPOp<T>*Op=0 )
    : _id( AUXCONT, FP->_naux++ ), _num( X ), _FP( FP ), _Op( Op )
    {}

  //! @brief Constructor for an intermediate binary variable in the factorable program <a>FP</a>
  FPVar
    ( FPRelax<T>*FP, FPOp<T>*Oper=0 )
    : _id( AUXBIN, FP->_naux++ ), _num( Op<T>::zeroone() ), _FP( FP ), _Op( Oper )
    {}

  //! @brief Constructor for a variable with identity <a>id</a> in the factorable program <a>FP</a>
  FPVar
    ( FPRelax<T>*FP, const pt_idVar&id )
    : _id( id ), _num( 0 ), _FP( FP ), _Op( 0 )
    {}
    
public:

  /** @ingroup FP
   *  @{
   */
  //! @brief Get variable identifier
  const std::pair<TYPE,long> id() const
    { return _id; }

  //! @brief Get reference to variable identifier
  std::pair<TYPE,long>& id()
    { return _id; }

  //! @brief Get const reference to variable numeric field
  const FPVarNum<T>& num() const
    { return _num; }

  //! @brief Get variable bounds
  const T& I() const
    { return _num.I; }

  //! @brief Get/set variable bounds
  T& I()
    { return _num.I; }

  //! @brief Get const pointer to defining operation
  const FPOp<T>* Oper() const
    { return _Op; }

  //! @brief Get pointer to defining operation
  FPOp<T>*& Oper()
    { return _Op; }

  //! @brief Get const pointer to factorable program
  const FPRelax<T>* FP() const
    { return _FP; }

  //! @brief Get pointer to factorable program
  FPRelax<T>*& FP()
    { return _FP; }

  //! @brief Get variable name
  std::string name() const
    { return _name(_id); }
  /** @} */

private:

  //! @brief Get name of variable with identifier <a>id</a>
  static std::string _name
    ( const std::pair<TYPE,long> id )
    { std::ostringstream ovar;
      (id.first==VARCONT? ovar<<"X": (id.first==VARBIN? ovar<<"Y": ovar<<"Z" ));
      ovar<<id.second;
      return ovar.str(); }
};

//! @brief C++ structure for comparing variables in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPVar<T> is a C++ structure for comparing variables in a
//! factorable program.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_FPVar
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPVar<T>*Var1, const FPVar<T>*Var2 ) const
    {
      // Order variables/constants w.r.t. their types first
      if( Var1->_id.first < Var2->_id.first ) return true;
      if( Var1->_id.first > Var2->_id.first ) return false;
      // If variables, order w.r.t. their index next
      switch( Var1->_id.first ){
        case FPVar<T>::VARCONT: case FPVar<T>::VARBIN:
        case FPVar<T>::AUXCONT: case FPVar<T>::AUXBIN:
          if( Var1->_id.second < Var2->_id.second ) return true;
          if( Var1->_id.second > Var2->_id.second ) return false;
          break;
        case FPVar<T>::AUXINT: case FPVar<T>::AUXREAL:
          lt_FPVarNum<T> ltNum;
          return ltNum( &Var1->_num, &Var2->_num );
          break;
      }
      return false;
    }
};

//! @brief C++ template class for defining cuts in the relaxation of factorable programs
////////////////////////////////////////////////////////////////////////
//! mc::FPCut<T> is a C++ template class defining cuts in the relaxation
//! of factorable programs.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPCut
////////////////////////////////////////////////////////////////////////
{
  // friends of class FPCut for operator overloading
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPCut<U>& );

public:
  /** @ingroup FP
   *  @{
   */
  //! @brief Enumeration type for cuts and special sets
  enum TYPE{
    EQ=0,	//!< Equality constraint Ax=b
    LE,		//!< Inequality constraint Ax<=b
    GE,		//!< Inequality constraint Ax>=b
    SOS1,	//!< SOS1-type constraint 
    SOS2	//!< SOS2-type constraint 
  };
  /** @} */

private:
  //! @brief Pointer to defining operation
  FPOp<T> *_op;
  //! @brief Type of cut
  TYPE _type;
  //! @brief Right-hand side
  double _rhs;
  //! @brief Number of participating variables
  unsigned int _nvar;
  //! @brief Coefficients
  double* _coef;
  //! @brief Indices
  typename FPVar<T>::pt_idVar* _idvar;

public:
  /** @ingroup FP
   *  @{
   */
  //! @brief Retreive type of cut
  TYPE type() const
    { return _type; }
  //! @brief Retreive right-hand side
  double rhs() const
    { return _rhs; }
  //! @brief Retreive number of participating variables
  unsigned int nvar() const
    { return _nvar; }
  //! @brief Retreive variable coefficients
  const double* coef() const 
    { return _coef; }
  //! @brief Retreive variable indices
  const typename FPVar<T>::pt_idVar* idvar() const 
    { return _idvar; }
  //! @brief Retreive pointer to defining operation
  FPOp<T>*& op()
    { return _op; }
  //! @brief Retreive pointer to defining operation
  const FPOp<T>* op() const
    { return _op; }

  //! @brief Default constructor for cut
  FPCut
    ( FPOp<T>*op, TYPE type=EQ ):
    _op(op), _type(type), _nvar(0), _coef(0), _idvar(0)
    {}

  //! @brief Constructor for cut w/ 1 participating variable
  FPCut
    ( FPOp<T>*op, TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1 ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar  = 1;
      _coef  = new double[_nvar];
      _idvar = new typename FPVar<T>::pt_idVar[_nvar];
      _coef[0] = a1;
      _idvar[0] = id1;
    }
  //! @brief Constructor for cut w/ 2 participating variables
  FPCut
    ( FPOp<T>*op, TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1,
      const typename FPVar<T>::pt_idVar id2, const double a2 ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar  = 2;
      _coef  = new double[_nvar];
      _idvar = new typename FPVar<T>::pt_idVar[_nvar];
      _coef[0] = a1;
      _coef[1] = a2;
      _idvar[0] = id1;
      _idvar[1] = id2;
    }
  //! @brief Constructor for cut w/ 3 participating variables
  FPCut
    ( FPOp<T>*op, TYPE type, const double b,
      const typename FPVar<T>::pt_idVar id1, const double a1,
      const typename FPVar<T>::pt_idVar id2, const double a2,
      const typename FPVar<T>::pt_idVar id3, const double a3 ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar  = 3;
      _coef  = new double[_nvar];
      _idvar = new typename FPVar<T>::pt_idVar[_nvar];
      _coef[0] = a1;
      _coef[1] = a2;
      _coef[2] = a3;
      _idvar[0] = id1;
      _idvar[1] = id2;
      _idvar[2] = id3;
    }
  //! @brief Constructor for cut w/ <a>n</a> participating variables
  FPCut
    ( FPOp<T>*op, TYPE type, const double b,
      const unsigned int n, const typename FPVar<T>::pt_idVar*id,
      const double*a ):
    _op(op), _type(type), _rhs(b)
    {
      _nvar  = n;
      _coef  = new double[_nvar];
      _idvar = new typename FPVar<T>::pt_idVar[_nvar];
      for( unsigned int ivar=0; ivar<_nvar; ivar++ ){
        _coef[ivar] = a[ivar]; _idvar[ivar] = id[ivar];
      }
    }
  //! @brief Destructor
  ~FPCut()
    {
      delete[] _coef;
      delete[] _idvar;
    }
  /** @} */
};

//! @brief C++ structure for cuts comparison
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPCut is a C++ structure for cuts comparison based on their
//! types and operation.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_FPCut
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPCut<T>*Cut1, const FPCut<T>*Cut2 ) const
    {
      // Order cuts w.r.t. their types first
      if( Cut1->type() < Cut2->type() ) return true;
      // then w.r.t. their defining operation next
      return ( Cut1->op() && Cut2->op() ?
        lt_FPOp<T>()( Cut1->op(), Cut2->op() ):
        false );
    }
};

//! @brief C++ template class for RLT constraints in the relaxation of factorable programs involving Taylor models
////////////////////////////////////////////////////////////////////////
//! mc::FPCut<T> is a C++ template class defining RLT constraints in
//! the relaxation of factorable programs involving Taylor models.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPRLT
////////////////////////////////////////////////////////////////////////
{
  // friends of class FPCut for operator overloading
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPRLT<U>& );

private:
  //! @brief Structure containing information for terms in RLT constraint LHS
  struct RLTTerm
  {
    //! @brief Constructor
    RLTTerm
      ( double a, const unsigned int n, const unsigned int*e):
      coef(a), nexp(n)
      { iexp = new unsigned int[nexp];
        for(unsigned int i=0; i<nexp; i++) iexp[i]=e[i]; }
    //! @brief Constructor (copy)
    RLTTerm
      ( const RLTTerm&term ):
      coef(term.coef), nexp(term.nexp)
      { iexp = new unsigned int[nexp];
        for(unsigned int i=0; i<nexp; i++) iexp[i]=term.iexp[i]; }
    //! @brief Destructor
    ~RLTTerm()
      { delete[] iexp; }
    //! @brief Term coefficient
    double coef;
    //! @brief Term exponent array
    unsigned int *iexp;
    //! @brief Term exponent size
    unsigned int nexp;
  };

public:
  /** @ingroup FP
   *  @{
   */
  //! @brief Typedef for left-hand side entry
  typedef std::pair< const FPVar<T>*, RLTTerm > p_lhs;
  //! @brief Typedef for left-hand side map
  typedef std::map< const FPVar<T>*, RLTTerm, lt_FPVar<T> > t_lhs;
  //! @brief Typedef for left-hand side map iterator
  typedef typename t_lhs::iterator it_lhs;
  //! @brief Typedef for left-hand side map const_iterator
  typedef typename t_lhs::const_iterator cit_lhs;
  //! @brief Enumeration type for RLT constraints
  enum TYPE{
    LE=0,		//!< Less-than-or-equal RLT constraint
    GE			//!< Greater-than-or-equal RLT constraint
  };
  /** @} */

private:
  //! @brief Pointer to defining operation
  TModel< FPVar<T> > *_pTM;
  //! @brief Type of RLT constraint
  TYPE _type;
  //! @brief Order of RLT constraint
  unsigned int _order;
  //! @brief RLT constraint left-hand side
  t_lhs _lhs;  
  //! @brief RLT constraint right-hand side
  double _rhs;
  //! @brief RLT constraint flag
  std::pair<unsigned int, bool> _flag;

public:
  /** @ingroup FP
   *  @{
   */
  //! @brief Retreive pointer to underyling Taylor model
  TModel< FPVar<T> >*& TM()
    { return _pTM; }
  //! @brief Retreive const pointer to underyling Taylor model
  const TModel< FPVar<T> >* TM() const
    { return _pTM; }
  //! @brief Retreive type of RLT constraint
  TYPE type() const
    { return _type; }
  //! @brief Retreive order of RLT constraint
  unsigned int order() const
    { return _order; }
  //! @brief Retreive reference to left-hand side map
  const t_lhs& lhs() const 
    { return _lhs; }
  //! @brief Retreive right-hand side value
  double rhs() const
    { return _rhs; }
  //! @brief Retreive flag
  std::pair<unsigned int, bool> flag() const
    { return _flag; }

  //! @brief Default constructor for RLT constraint
  FPRLT
    ( TModel< FPVar<T> >*pTM ):
    _pTM(pTM), _type(GE), _order(0), _lhs(), _rhs(0.), _flag(0,true)
    {}
  //! @brief Destructor
  ~FPRLT()
    {}

  //! @brief Generate new RLT constraint by multiplying with (X-XL) or (X-XU) - <a>flag</a> indicates the current variable index in TModel and whether the lifting is w.r.t. lower or upper bound
  FPRLT<T>* product
    ( const FPVar<T>*Var, const unsigned int nexp, const unsigned int*iexp,
      std::pair<unsigned int, bool> flag ) const;
  //! @brief Generate and return pair of variables indices (pt_idVar*) and coefficient (double*) arrays for RLT constraint
  std::pair<typename FPVar<T>::pt_idVar*, double*> lhs2cut
    () const;
  //! @brief Populate pair of variables indices (pt_idVar*) and coefficient (double*) arrays for RLT constraint - returns whether array size <a>lhsmax</a> is big enough (true) or not (false)
  bool lhs2cut
    ( std::pair<typename FPVar<T>::pt_idVar*, double*>&lhs,
      const unsigned int lhsmax ) const;
  /** @} */
};

//! @brief C++ structure for storing subintervals in an outer approximation
////////////////////////////////////////////////////////////////////////
//! mc::FPOuter is a C++ structure for storing subintervals in the
//! outer approximation of a univariate convex/concave function.
////////////////////////////////////////////////////////////////////////
class FPOuter
////////////////////////////////////////////////////////////////////////
{

public:
  //! @brief Constructor
  FPOuter( const double LB, const double UB, const double MID, 
    const double GAP ):
    _xL( LB ), _xU( UB ), _xM( MID ), _gap( GAP )
    {}
  //! @brief Retreive interval lower bound
  const double xL() const
    { return _xL; }
  //! @brief Retreive interval upper bound
  const double xU() const
    { return _xU; }
  //! @brief Retreive bisection point
  const double xM() const
    { return _xM; }
  //! @brief Retreive maximum gap
  const double gap() const
    { return _gap; }

private:
  //! @brief Interval lower bound
  double _xL;
  //! @brief Interval upper bound
  double _xU;
  //! @brief Bisection point
  double _xM;
  //! @brief Maximum gap
  double _gap;
};

//! @brief C++ structure for comparing subintervals in an outer approximation
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPOuter is a C++ structure for comparing subintervals in the
//! outer approximation of a univariate convex/concave function.
////////////////////////////////////////////////////////////////////////
struct lt_FPOuter
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPOuter&Dom1, const FPOuter&Dom2 ) const
    {
      return( Dom1.gap() < Dom2.gap() );
    }
};

//! @brief C++ class for defining operations in a factorable program, and generating linear relaxations for such operations
////////////////////////////////////////////////////////////////////////
//! mc::FPOp<T> is a C++ template class for defining operations in a
//! factorable program, and generating linear relaxations for such
//! operations.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPOp
////////////////////////////////////////////////////////////////////////
{
public:

  typedef double (*p_Univ)( const double, const double*, const int* );
  typedef std::pair<double,double> (*p_dUniv)( const double, const double*,
    const int* );
  typedef void (FPOp<T>::*p_Cut)( FPRelax<T>*,
    const typename FPVar<T>::pt_idVar, const double, const double,
    const double, const typename FPVar<T>::pt_idVar, const double,
    const double );
  typedef void (FPOp<T>::*p_Cutpar)( FPRelax<T>*, const double,
    const typename FPVar<T>::pt_idVar, const double, const double,
    const typename FPVar<T>::pt_idVar, const double, const double,
    const double*, const int* );
  typedef std::priority_queue< FPOuter, std::vector<FPOuter>, lt_FPOuter > t_Outer;

  /** @ingroup FP
   *  @{
   */
  //! @brief Enumeration type for unary and binary operations
  enum TYPE{
    CNST=0, RANGE, VAR, MC, TM, TMMC,
    EQ, LE, MIN, MAX, RLT,
    PLUS, NEG, MINUS, TIMES, SCALE, DIV, INTER, BILIN,
    EXP, LOG, SQRT, SQR, IPOW, POW, COS, TAN, ASIN, ATAN,
    FABS, XLOG, ARH, ERF, FSTEP, MINF, MAXF
  };

  //! @brief Constructor
  FPOp( TYPE op, FPVar<T>*lop=0, FPVar<T>*rop=0, FPVar<T>*res=0 ):
    type( op ), pres( res ), _visited(false), _pMC(0), _pTM(0), _pTMMC(0)
    {
      // Order operands for commutative binary operators
      if( !lop
        || ( op != PLUS && op != TIMES && op != SCALE && op != INTER && op != EQ )
        || lt_FPVar<T>()( lop, rop ) )
        { plop = lop; prop = rop; }
      else
        { plop = rop; prop = lop; }
    }

  //! @brief Destructor
  ~FPOp()
    {
      delete _pMC;
      delete _pTM;
      delete _pTMMC;
    }

  //! @brief Type of operation
  TYPE type;
  //! @brief Pointer to operation result
  FPVar<T>* pres;
  //! @brief Pointer to left operand
  FPVar<T>* plop;
  //! @brief Pointer to right operand
  FPVar<T>* prop;

  //! @brief Generate polyhedral cuts for all operations in the factorable program <a>pFP</a>
  void generate_cuts
    ( FPRelax<T>*pFP );
  //! @brief Append polyhedral cuts for the current operation in the factorable program <a>pFP</a>
  void append_cuts
    ( FPRelax<T>*pFP );
  //! @brief Propagate bounds for all operations in the factorable program <a>pFP</a>
  bool propagate_bounds
    ( FPRelax<T>*pFP );
  //! @brief Tighten variable bounds for the current operation in the factorable program <a>pFP</a> via constraint propagation; on returns, <a>false</a> is an indication of an infeasible program
  bool tighten_bounds
    ( FPRelax<T>*pFP );

  //! @brief Attach Taylor model to variable to be handled during cut generation
  void attach
    ( const McCormick<T>&MC )
    { _pMC = new McCormick<T>( MC ); }
  //! @brief Attach Taylor model to variable to be handled during cut generation
  void attach
    ( const TVar<T>&TV, TModel< FPVar<T> >*TMFP )
    { _pTM = new std::pair< TVar<T>, TModel< FPVar<T> >* >( TV, TMFP ); }
  //! @brief Attach McCormick-Taylor model to variable to be handled during cut generation
  void attach
    ( const TVar< McCormick<T> >&TVMC, TModel< FPVar<T> >*TMFP )
    { _pTMMC = new std::pair< TVar< McCormick<T> >, TModel< FPVar<T> >* >( TVMC, TMFP ); }

  //! @brief Set operation status (visited or not)
  void set
    ( const bool visited=true )
    { _visited = visited; }
  //! @brief Retreive operation status (visited or not)
  const bool status() const
    { return _visited; }
  //! @brief Retreive/set operation status (visited or not)
  const bool& status()
    { return _visited; }
  /** @} */

private:
  //! @brief Whether a cut has been visited (e.g. in a tree navigation)
  bool _visited;

  //! @brief Pointer to McCormick variable
  McCormick<T> *_pMC;
  //! @brief Pointer to Taylor variable
  std::pair< TVar<T>, TModel< FPVar<T> >* > *_pTM;
  //! @brief Pointer to McCormick-Taylor variable
  std::pair< TVar< McCormick<T> >, TModel< FPVar<T> >* > *_pTMMC;

  //! @brief Variable bounds propagation for equality constraints
  bool _EQ_bounds
    ( FPRelax<T>*pFP, FPVar<T>*Var1, FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for less than equal constraints
  bool _LE_bounds
    ( FPRelax<T>*pFP, FPVar<T>*Var1, FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for  operation
  bool _PLUS_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1,
      FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for unary - operation
  bool _NEG_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1 );
  //! @brief Variable bounds propagation for - operation
  bool _MINUS_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1,
      FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for * operation
  bool _TIMES_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1,
      FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for / operation
  bool _DIV_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1,
      FPVar<T>*Var2 );
  //! @brief Variable bounds propagation for univariate exp term
  bool _EXP_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var );
  //! @brief Variable bounds propagation for univariate log term
  bool _LOG_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var );
  //! @brief Variable bounds propagation for univariate sqr term
  bool _SQR_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var );
  //! @brief Variable bounds propagation for univariate sqrt term
  bool _SQRT_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var );
  //! @brief Variable bounds propagation for univariate fabs term
  bool _FABS_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var );
  //! @brief Variable bounds propagation for univariate pow term
  bool _IPOW_bounds
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var, const FPVar<T>*Exp );

  //! @brief Append linear cuts for equality constraints
  void _EQ_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*Var1, const FPVar<T>*Var2 );
  //! @brief Append linear cuts for less than equal constraints
  void _LE_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*Var1, const FPVar<T>*Var2 );
  //! @brief Append linear cuts for + operation
  void _PLUS_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for unary - operation
  void _NEG_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1 );
  //! @brief Append linear cuts for - operation
  void _MINUS_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for * operation
  void _TIMES_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for / operation
  void _DIV_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for intersection (^) operation
  void _INTER_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for bivariate min operation
  void _MINF_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for bivariate max operation
  void _MAXF_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for bilinear product operation
  void _BILIN_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Perform variable subdivision (equal subintervals)
  FPVar<T>** _subdivide
    ( FPRelax<T>*pFP, const FPVar<T>*Var, const unsigned int Ndiv,
      double*a, typename FPVar<T>::pt_idVar*id );
  //! @brief Append semi-linear cuts at equal subintervals
  void _semilinear_cut
    ( FPRelax<T>*pFP, const unsigned int Ndiv, const double XL,
      const double XU, FPVar<T>**dVarX, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, p_Univ f, const double*rpar=0,
      const int*ipar=0 );
  //! @brief Perform variable subdivision using SOS2 sets
  FPVar<T>** _subdivide_SOS
    ( FPRelax<T>*pFP, const FPVar<T>*Var, const unsigned int Ndiv,
      double*a, typename FPVar<T>::pt_idVar*id );
  //! @brief Append semi-linear cuts for bilinear product operation using SOS2 sets
  void _semilinear_SOS
    ( FPRelax<T>*pFP, const unsigned int Ndiv, const double XL,
      const double XU, FPVar<T>**dVarX, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, p_Univ f, const double*rpar=0,
      const int*ipar=0 );
  //! @brief Computes max distance between function and outer-approximation
  std::pair< double, double > _distmax
    ( const FPRelax<T>*pFP, p_dUniv f, const double xL, const double xU,
      const double*rpar=0, const int*ipar=0 ) const;
  //! @brief Append cuts for nonlinear operation using outer-approximatio (sandwich algorithm)
  void _sandwich_cuts
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, p_Cutpar Cut, p_dUniv f, const double*rpar=0,
      const int*ipar=0 );
  //! @brief Computes solution of scalar nonlinear equation using Newton's method
  double _newton
    ( const FPRelax<T>*pFP, const double x0, const double xL,
      const double xU, p_dUniv f, const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Computes solution of scalar nonlinear equation using the secant method
  double _secant
    ( const FPRelax<T>*pFP, const double x0, const double x1, const double xL,
      const double xU, p_Univ f, const double*rusr=0, const int*iusr=0 ) const;
  //! @brief Append linear cuts for univariate exp function
  void _EXP_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts for univariate log function
  void _LOG_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of univariate exp function
  void _addcut_expcv
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar=0, const int*ipar=0 );
  //! @brief Append linear cuts corresponding to concave relaxation of univariate exp function
  void _addcut_expcc
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, const double*rpar=0, const int*ipar=0 );
  //! @brief Append linear cuts for univariate sqr function
  void _SQR_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts for univariate sqrt function
  void _SQRT_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of univariate sqr function
  void _addcut_sqrcv
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar=0, const int*ipar=0 );
  //! @brief Append linear cuts corresponding to concave relaxation of univariate sqr function
  void _addcut_sqrcc
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, const double*rpar=0, const int*ipar=0 );
  //! @brief Append linear cuts for univariate fabs function
  void _FABS_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts for univariate power function w/ integer exponents
  void _IPOW_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var,
      const FPVar<T>*Exp );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_powcv_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*iExp );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_powcc_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*iExp );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_powcv_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*iExp );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_powcc_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*iExp );
  //! @brief Append linear cuts for univariate cosine function
  void _COS_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_coscv_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_coscc_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_coscv_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_coscc_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts for univariate asin function
  void _ASIN_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of univariate asin function
  void _addcut_asin1cv
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of univariate asin function
  void _addcut_asin1cc
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to convex relaxation of univariate asin function
  void _addcut_asin2cc
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of univariate asin function
  void _addcut_asin2cv
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to convex relaxation of univariate asin function
  void _addcut_asin3cc
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of univariate asin function
  void _addcut_asin3cv
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
      const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts for univariate atan function
  void _ATAN_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts for univariate tan function
  void _TAN_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_atancv_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_atancc_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_atancv_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_atancc_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts for univariate atan function
  void _ERF_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_erfcv_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to convex relaxation of power function
  void _addcut_erfcc_lin
    ( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
      const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
      const double YL, const double YU, const double*rpar, const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_erfcv_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts corresponding to concave relaxation of power function
  void _addcut_erfcc_sec
    ( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
      const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
      const int*ipar );
  //! @brief Append linear cuts for univariate atan function
  void _FSTEP_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var );
  //! @brief Append linear cuts for bivariate max function
  void _FMAX_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
  //! @brief Append linear cuts for bivariate min function
  void _FMIN_cuts
    ( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
      const FPVar<T>*Var2 );
};

//! @brief C++ structure for comparing operations in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPOp is a C++ structure for comparing operations in a
//! factorable program based on their types and operands.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_FPOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPOp<T>*Op1, const FPOp<T>*Op2 ) const
    {
      // Sort by type of operation first
      if( Op1->type < Op2->type ) return true;
      if( Op1->type > Op2->type ) return false;
      if( Op1->type == FPOp<T>::RLT ) return false;

      // Sort by variable type next
      lt_FPVar<T> ltVar;
      if( !Op1->plop ) return ltVar( Op1->pres, Op2->pres );
      if( ltVar( Op1->plop, Op2->plop ) ) return true;
      if( ltVar( Op2->plop, Op1->plop ) ) return false;
      if( Op1->prop ) return ltVar( Op1->prop, Op2->prop );
      return false;
    }
};

//! @brief C++ structure for comparing operations in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_FPOp is a C++ structure for comparing operations in a
//! factorable program based on their types only.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct range_FPOp
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const FPOp<T>*Op1, const FPOp<T>*Op2 ) const
    { return ( Op1->type < Op2->type ); }
};

//! @brief C++ template base class for the definition and linear relaxation of factorable programs
////////////////////////////////////////////////////////////////////////
//! mc::FPRelax<T> is a C++ base class that allows definition of
//! factorable programs, and generation of linear relaxations for such
//! programs.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPRelax
////////////////////////////////////////////////////////////////////////
{
  // friends of class FPRelax
  template <typename U> friend class FPVar;
  template <typename U> friend class FPOp;
  //template <typename U> friend void FPOp<U>::append_cuts( FPRelax<U>*pFP );

  // friends of class FPVar for operator overloading
  template <typename U> friend std::ostream& operator<<
    ( std::ostream&, const FPRelax<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator+
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator+
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator+
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator-
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator-
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator-
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator*
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator*
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator*
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator*
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator*
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator/
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator/
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> operator/
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> operator/
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> operator/
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> operator^
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> max
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> max
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> max
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> max
    ( const unsigned int, const FPVar<U>* );
  template <typename U> friend FPVar<U> min
    ( const FPVar<U>&, const FPVar<U>& );
  template <typename U> friend FPVar<U> min
    ( const U&, const FPVar<U>& );
  template <typename U> friend FPVar<U> min
    ( const FPVar<U>&, const U& );
  template <typename U, typename V> friend FPVar<U> min
    ( const V&, const FPVar<U>& );
  template <typename U, typename V> friend FPVar<U> min
    ( const FPVar<U>&, const V& );
  template <typename U> friend FPVar<U> min
    ( const unsigned int, const FPVar<U>* );

  // friends of class FPVar for function overloading
  template <typename U> friend FPVar<U> inv
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sqr
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> exp
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> log
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sqrt
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> fabs
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> pow
    ( const FPVar<U>&, const int );
  template <typename U> friend FPVar<U> cos
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> sin
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> tan
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> acos
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> asin
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> atan
    ( const FPVar<U>& );
//   template <typename U> friend FPVar<U> xlog
//     ( const FPVar<U>& );
//   template <typename U> friend FPVar<U> arh
//     ( const FPVar<U>&, const double );
  template <typename U> friend FPVar<U> erf
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> erfc
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> fstep
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> bstep
    ( const FPVar<U>& );
  template <typename U> friend FPVar<U> bilin
    ( const FPVar<U>&, const FPVar<U>& );

public:
  /** @ingroup FP
   *  @{
   */
  typedef std::multimap< const FPVar<T>*, FPOp<T>*, lt_FPVar<T> > t_Vars_Ops;
  typedef typename t_Vars_Ops::iterator it_Vars_Ops;
  typedef typename t_Vars_Ops::const_iterator cit_Vars_Ops;
  typedef std::multiset< FPCut<T>*, lt_FPCut<T> > t_Cuts;
  typedef std::set< FPVar<T>*, lt_FPVar<T> > t_Vars;
  typedef std::set< FPOp<T>*,  lt_FPOp<T> > t_Ops;
  typedef typename t_Cuts::iterator it_Cuts;
  typedef typename t_Cuts::const_iterator cit_Cuts;
  typedef typename t_Vars::iterator it_Vars;
  typedef typename t_Vars::const_iterator cit_Vars;
  typedef typename t_Ops::iterator  it_Ops;
  typedef typename t_Ops::const_iterator  cit_Ops;
  typedef typename std::pair< it_Vars, bool > pt_Vars;
  typedef typename std::pair< it_Ops, bool > pt_Ops;
  typedef typename FPVar<T>::pt_idVar pt_idVar;
  /** @} */

protected:
  //! @brief Number of native variables in factorable program
  unsigned long _nvar;
  //! @brief Number of auxiliary variables in factorable program
  unsigned long _naux;

  //! @brief Set of variables in factorable program
  t_Vars _Vars;
  //! @brief Set of operations in factorable program
  t_Ops  _Ops;
  //! @brief Set of cuts in relaxed factorable program
  t_Cuts _Cuts;

public:
  /** @ingroup FP
   *  @{
   */
  //! @brief Default Constructor
  FPRelax():
    _nvar( 0 ), _naux( 0 ) 
    {}

  //! @brief Destructor
  virtual ~FPRelax()
    { reset(); }

  //! @brief FPRelax exceptions
  class Exceptions
  {
  public:
    //! @brief Enumeration type for exception handling
    enum TYPE{
      NEWTON=1, 	//!< Failure of Newton's method during the calculation of the convex or concave envelope of a univariate term
      INTER, 	//!< Error during intersection of two terms (terms may not intersect)
      CONSTRAINT,	//!< Error in constraint definition (e.g., number = number)
      OBJECTIVE,	//!< Error in objective definition (e.g., minimize interval)
      INIT, 	//!< Error during initialization with a variable of class mc::McCormick or mc::Taylor
      FP=-1, 	//!< Error due to an operation between variables participating in different factorable programs
      UNDEF=-2, 	//!< Error due to calling a function/feature not yet implemented in MC++
      INTERNAL=-3	//!< Internal error
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions( TYPE ierr ) : _ierr( ierr ){}
    //! @brief Inline function returning the error flag
    int ierr(){ return _ierr; }

  private:
    TYPE _ierr;
  };

  //! @brief FPRelax options
  struct Options
  {
    //! @brief Constructor
    Options():
      NEWTON_USE(true), NEWTON_MAXIT(100), NEWTON_TOL(1e-10),
      SANDWICH_ATOL(1e-3), SANDWICH_RTOL(1e-3), SANDWICH_MAXCUT(5),
      SANDWICH_RULE(MAXERR), FRACTIONAL_ATOL(machprec()),
      FRACTIONAL_RTOL(machprec()), BILINEAR_RULE(UNIVARIATE),
      BILINEAR_SUBDIV(1), RLT_MAXORDER(4), RLT_DISPLAY(0),
      RRLT_STRATEGY(PRIMRRLT), RRLT_DISPLAY(0), LP_FEASTOL(1e-9),
      LP_OPTIMTOL(1e-9), MILP_RELGAP(0e0), MILP_ABSGAP(0e0),
      SOLVER_PRESOLVE(true), SOLVER_DISPLAY(false), SOLVER_OUTPUT_FILE("")
      {}
    //! @brief Enumeration type for sandwich strategy
    enum SANDWICH{
      BISECT=0,	//!< Interval bisection
      MAXERR	//!< Maximum error rule
    };
    //! @brief Enumeration type for reduced (R)RLT strategy
    enum RRLT{
      //NORRLT=0,	//!< No reduced RLT constraints
      PRIMRRLT=0,	//!< RRLT with primary variables as operands only
      ALLRRLT		//!< RRLT with primary and auxiliary variables as operands
    };
    //! @brief Enumeration type for bilinear term relaxation strategy
    enum BILINEAR{
      UNIVARIATE=0,	//!< Bisection along the first variable in bilinear term
      SEPARABLE,		//!< Reformulation as separable terms before relaxation w/o SOS
      SEPARSOS		//!< Reformulation as separable terms before relaxation w/ SOS
    };
    //! @brief Flag to indicate whether Newton's method is to be used
    bool NEWTON_USE;
    //! @brief Maximal number of iterations in Newton's method
    unsigned int NEWTON_MAXIT;
    //! @brief Termination tolerance in Newton's method
    double NEWTON_TOL;
    //! @brief Absolute tolerance in outer approximation of nonlinear convex/concave terms
    double SANDWICH_ATOL;
    //! @brief Relative tolerance in outer approximation of nonlinear convex/concave terms
    double SANDWICH_RTOL;
    //! @brief Maximal number of cuts in outer approximation of nonlinear convex/concave terms
    unsigned int SANDWICH_MAXCUT;
    //! @brief Rule for outer approximation of nonlinear convex/concave terms
    SANDWICH SANDWICH_RULE;
    //! @brief Absolute tolerance in (real) fractional terms to prevent division by zero
    double FRACTIONAL_ATOL;
    //! @brief Relative tolerance in (real) fractional terms to prevent division by zero
    double FRACTIONAL_RTOL;
    //! @brief Rule for relaxation of bilinear terms
    BILINEAR BILINEAR_RULE;
    //! @brief Subdivision for semi-linear relaxation of bilinear terms
    unsigned int BILINEAR_SUBDIV;
    //! @brief Maximum order for RLT constraints
    unsigned int RLT_MAXORDER;
    //! @brief Option for RLT display
    unsigned int RLT_DISPLAY;
    //! @brief Option for RRLT strategy
    RRLT RRLT_STRATEGY;
    //! @brief Option for RRLT display
    unsigned int RRLT_DISPLAY;
    //! @brief Feasibility tolerance of LP solver
    double LP_FEASTOL;
    //! @brief Optimality tolerance of LP solver
    double LP_OPTIMTOL;
    //! @brief Relative optimality gap of MILP solver
    double MILP_RELGAP;
    //! @brief Absolute optimality gap of MILP solver
    double MILP_ABSGAP;
    //! @brief Flag to indicate whether presolve is to be used by solver
    bool SOLVER_PRESOLVE;
    //! @brief Flag to indicate whether solver output is to be displayed
    bool SOLVER_DISPLAY;
    //! @brief Name of output file for optimization model
    std::string SOLVER_OUTPUT_FILE;
  } options;
  
  //! @brief Retreive number of original variables in factorable program
  unsigned long nvar() const
    { return _nvar; }
  
  //! @brief Retreive number of auxiliary variables in factorable program
  unsigned long naux() const
    { return _naux; }
  
  //! @brief Retreive reference to set of cuts in factorable program relaxation
  const t_Cuts& Cuts() const
    { return _Cuts; }
  
  //! @brief Retreive reference to set of (all) variables in factorable program
  const t_Vars& Vars() const
    { return _Vars; }

  //! @brief Clear the factorable program (all variables and operations)
  void reset()
    { _erase_cuts(); _erase_vars(); _erase_operations(); _naux=_nvar=0; }

  //! @brief Enumeration for objective function type in factorable program
  enum OBJTYPE{
    MIN=0,	//!< Minimization
    MAX		//!< Maximization
  };

  //! @brief Set objective function in factorable program
  virtual void set_objective
    ( const typename FPRelax<T>::OBJTYPE type, const FPVar<T>&objVar );

  //! @brief Enumeration for constraint type in factorable program
  enum CTRTYPE{
    EQ=0,	//!< Equality constraint
    LE,		//!< Inequality constraint
    GE		//!< Inequality constraint
  };

  //! @brief Add constraint to factorable program
  virtual FPOp<T>* add_constraint
    ( const FPVar<T>&lhsVar, const typename FPRelax<T>::CTRTYPE type,
      const FPVar<T>&rhsVar );
  //! @brief Remove constraint from factorable program
  virtual bool remove_constraint
    ( FPOp<T>* constr )
    { return _erase_operation( constr ); }

  //! @brief Generate polyhedral cuts for the objective and all constraints in the factorable program, after removing or not the existing cuts
  void generate_cuts
    ( bool reset=false );
  //! @brief Generate polyhedral cuts for the constraint <a>constr</a> in the factorable program
  void generate_cuts
    ( FPOp<T>*constr );
  //! @brief Propagate constraints in the factorable program to reduce variable range, for the constraint <a>constr</a> or all constraints in case <a>constr=0</a>; on return, <a>false</a> is an indication of an infeasible program
  bool propagate_constraints
    ( FPOp<T>*constr=0 );
  //! @brief Generate reduction constraints via Liberti & Pantelides's RRLT approach
  void generate_reduction_constraints();
  //! @brief Generate tightening constraints via Sherali's RLT approach for variables participating in the Taylor model <a>TMFP</a>
  void generate_RLT_TModel
    ( TModel< FPVar<T> >*TMFP );

  //! @brief Build Taylor model environment mc::TModel< FPVar<U> > corresponding to environment mc::TModel<U> <a>TM</a>, with matching indices given in <a>isub</a>
  template <typename U> mc::TModel< FPVar<T> >* create_TModel
    ( const mc::TModel<U>*TM, const unsigned int*isub=0 );
  //! @brief Update Taylor model environment mc::TModel< FPVar<U> > corresponding to environment mc::TModel<U> <a>TM</a>, with matching indices given in <a>isub</a>
  template <typename U> void update_TModel
    ( mc::TModel< FPVar<T> >* TMFP, const mc::TModel<U>*TM,
      const unsigned int*isub=0 );
  /** @} */
   
protected:
  //! @brief Erase operation <a>op</a> in set <a>_Ops</a>
  bool _erase_operation
    ( FPOp<T>* op );
  //! @brief Erase all operations in set <a>_Ops</a>
  void _erase_operations()
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) delete *ito;
      _Ops.clear(); }
  //! @brief Reset all operations in set <a>_Ops</a>
  void _reset_operations()
    { it_Ops ito = _Ops.begin();
      for( ; ito != _Ops.end(); ++ito ) (*ito)->set( false ); }
  //! @brief Looks for the operation of type <a>op</a> with left and right operands <a>lop</a>, <a>rop</a> in set <a>_Ops</a> and adds it if not found
  FPOp<T>*_operation
    ( const typename FPOp<T>::TYPE op, FPVar<T>*lop, FPVar<T>*rop=0 );

  //! @brief Adds the auxiliary continuous variable with bounds <a>X</a> from operation <a>pOp</a>
  FPVar<T>* _auxiliary_variable
    ( const T&X, FPOp<T>*pOp );
  //! @brief Adds the auxiliary binary variable from operation <a>pOp</a>
  FPVar<T>* _auxiliary_variable
    ( FPOp<T>*pOp );
  //! @brief Looks for the real auxiliary constant <a>x</a> and adds it if not found
  FPVar<T>* _auxiliary_constant
    ( const double x );
  //! @brief Looks for the integer auxiliary constant <a>n</a> and adds it if not found
  FPVar<T>* _auxiliary_constant
    ( const int n );
  //! @brief Adds the interval bounds <a>I</a>
  FPVar<T>* _auxiliary_constant
    ( const T&I );

  //! @brief Erase all variables in _Vars
  void _erase_vars()
    { it_Vars itv = _Vars.begin();
      for( ; itv != _Vars.end(); ++itv ) delete *itv;
      _Vars.clear(); }
  //! @brief Appends the auxiliary variable <a>pAux</a> and define it in _Ops with type <a>tOp</a>
  void _append_aux
    ( FPVar<T>*pAux, typename FPOp<T>::TYPE tOp );
  //! @brief Appends new auxiliary variable
  virtual void _append_aux
    ( FPVar<T>*pAux );
  //! @brief Appends new original variable
  virtual void _append_var
    ( FPVar<T>*pVar );
  //! @brief Update the bounds of the variable <a>Var</a> as <a>Bnd</a> 
  virtual void _update_var
    ( FPVar<T>*pVar, const T&Bnd )
    { pVar->I() = Bnd; }
  //! @brief Search for the variable with identify <a>id</a> in <a>_Vars</a>
  FPVar<T>* _find_var
    ( const typename FPVar<T>::pt_idVar&id );

  //! @brief Erase all relaxation cuts in _Cuts
  void _erase_cuts()
    { it_Cuts itc = _Cuts.begin();
      for( ; itc != _Cuts.end(); ++itc ) delete *itc;
      _Cuts.clear(); }
  //! @brief Erase relaxation cuts corresponding to operation <a>op</a> in _Cuts
  void _erase_cuts
    ( FPOp<T>* op );
  //! @brief Appends new relaxation cut in _Cuts w/ 1 variable
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type,
      const double b, const typename FPVar<T>::pt_idVar id1,
      const double a1 );
  //! @brief Appends new relaxation cut in _Cuts w/ 2 variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type,
      const double b, const typename FPVar<T>::pt_idVar id1,
      const double a1, const typename FPVar<T>::pt_idVar id2,
      const double a2 );
  //! @brief Appends new relaxation cut in _Cuts w/ 3 variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type,
      const double b, const typename FPVar<T>::pt_idVar id1,
      const double a1, const typename FPVar<T>::pt_idVar id2,
      const double a2, const typename FPVar<T>::pt_idVar id3,
      const double a3 );
  //! @brief Appends new relaxation cut in _Cuts w/ <a>n</a> variables
  virtual FPCut<T>* _append_cut
    ( FPOp<T>*op, const typename FPCut<T>::TYPE type,
      const double b, const unsigned int n,
      const typename FPVar<T>::pt_idVar* id, const double*a );

  //! @brief Lift RLT constraints in <a>RLT_lst</a> by 1 order for variables participating in the Taylor model <a>TMFP</a>
  void _lift_RLT_TModel
    ( TModel< FPVar<T> >*TMFP, std::list< FPRLT<T>* >& RLT_lst,
      const bool disp );
  //! @brief Append cuts for RLT constraints in <a>RLT_lst</a>
  void _append_cuts_RLT_TModel
    ( std::list< FPRLT<T>* >& RLT_lst, const unsigned int nord );
};

template <typename T> class lt_RLTVar;

//! @brief C++ template class for reduced RLT constraint generation
////////////////////////////////////////////////////////////////////////
//! mc::FPRRLT<T> is a C++ class implementing Liberti & Pantelides'
//! approach to generate valid reduced RLT constraints .
////////////////////////////////////////////////////////////////////////
template <typename T>
class FPRRLT
{
public:

  typedef typename FPVar<T>::pt_idVar pt_idVar;
  typedef std::pair< FPOp<T>*, unsigned int > pt_Op;

  typedef std::set< FPVar<T>*, lt_FPVar<T> > t_Vars;
  typedef typename t_Vars::iterator it_Vars;
  typedef typename t_Vars::const_iterator cit_Vars;

  typedef std::set< FPOp<T>*,  lt_FPOp<T> > t_Ops;
  typedef typename t_Ops::iterator  it_Ops;
  typedef typename t_Ops::const_iterator  cit_Ops;

  typedef std::multimap< const FPVar<T>*, FPOp<T>*, lt_FPVar<T> > t_Vars_Ops;
  typedef typename t_Vars_Ops::iterator it_Vars_Ops;
  typedef typename t_Vars_Ops::const_iterator cit_Vars_Ops;

  typedef std::pair< FPVar<T>*, typename FPOp<T>::TYPE > t_RLTVar;
  typedef std::multimap< const t_RLTVar, FPOp<T>*, lt_RLTVar<T> > t_RLTMap;
  typedef typename t_RLTMap::iterator it_RLTMap;
  typedef typename t_RLTMap::const_iterator cit_RLTMap;
  
private:
  //! @brief Reference to set of variables in linear relaxation
  const t_Vars& _Vars;
  //! @brief Reference to set of operations in linear relaxation
  const t_Ops& _Ops;

  //! @brief Subset of linear operations
  t_Ops _linTerms;
  //! @brief Subset of variables participating in linear operations
  t_Vars _linVars;
  //! @brief Map of linear operations and participating variables
  t_Vars_Ops _linEdges;
  //! @brief Subset of variables for RLT constraint search
  t_Vars _RLTVars;
  //! @brief Subset of edges for RLT constraint search
  t_Vars_Ops _RLTEdges;

  //! @brief Size of variable subset for RLT constraint search
  unsigned int _nRLTVars;
  //! @brief Size of linear operation subset
  unsigned int _nlinTerms;

  // subsets of bilinear terms
  t_Ops _bilTerms;
  // subsets of fractional terms
  t_Ops _divTerms;
  // subsets of square terms
  t_Ops _sqrTerms;
  // subsets of square root terms
  t_Ops _sqrtTerms;

  //! @brief Map of RLT constraints and corresponding multiplier variable
  t_RLTMap _RLTMap;

  //! @brief Whether a given RLT variable has been assigned to a linear term and which term
  pt_Op* _VarAssigned;
  //! @brief Whether a given RLT variable has been visited
  bool* _VarVisited;
  //! @brief Whether a given linear term has been visited
  bool* _TermVisited;  //! @brief Display level during reduction constraint generation
  unsigned int _display;

public:
  
  //! @brief Default Constructor
  FPRRLT
    ( const t_Vars& Vars, const t_Ops& Ops ):
    _Vars( Vars ), _Ops( Ops ), _nRLTVars( 0 ), _nlinTerms( 0 ),
    _VarAssigned( 0 ), _VarVisited(0), _TermVisited(0)
    {}

  //! @brief Destructor
  ~FPRRLT()
    { delete[] _VarAssigned; delete[] _VarVisited; delete[] _TermVisited; }

  //! @brief Identify valid reduction constraints with the possible candidate multiplier variables
  t_RLTMap& map_RLT
    ( typename FPRelax<T>::Options::RRLT option, unsigned int display );

private:
  //! @brief Define subsets and map of linear operations and participating variables
  void _define_linear();
  //! @brief Define subsets of bilinear, fractional, square and square root terms
  void _define_bilinear();
  //! @brief Identify variables not yet participating in any nonlinear terms with variable <a>idMult</a>, and create corresponding bigraph edges
  void _bigraph_RLT
    ( const typename FPVar<T>::pt_idVar& idMult, const typename FPOp<T>::TYPE& idOp );
  //! @brief Identify variables not yet participating in any nonlinear terms multiplied by variable <a>idMult</a>, and create corresponding bigraph edges
  void _bigraph_RLT_mult
    ( typename FPRRLT<T>::t_Vars &linVars, typename FPRRLT<T>::t_Vars_Ops &linEdges,
      const typename FPVar<T>::pt_idVar& idMult );
  //! @brief Identify variables not yet participating in any nonlinear terms divided by variable <a>idDiv</a>, and create corresponding bigraph edges
  void _bigraph_RLT_div
    ( typename FPRRLT<T>::t_Vars &linVars, typename FPRRLT<T>::t_Vars_Ops &linEdges,
      const typename FPVar<T>::pt_idVar& idDiv );
  //! @brief Identify valid reduction constraints with the candidate reduction variable <a>VarRed</a>
  void _reduction_RLT
    ( const t_RLTVar& VarRed );
  //! @brief Determine whether an augmented path can be found that emanates from the linear term <a>pOp<\a>
  bool _augpath_RLT
    ( const pt_Op& pOp );

public:
  //! @brief Create subset of operations of given type
  t_Ops subset_op
    ( const unsigned int nOp, const typename FPOp<T>::TYPE*typeOp ) const;
  //! @brief Create subset of operations of given type
  t_Ops subset_op
    ( const typename FPOp<T>::TYPE&typeOp ) const;
  //! @brief Create subset of variables participating in operations Ops
  static t_Vars subset_var
    ( const t_Ops& Ops );
  //! @brief Create subset of variables participating in operations Ops and corresponding edges
  static t_Vars_Ops submap_var
    ( const t_Ops& Ops );
};

//! @brief C++ structure for comparing variables in a factorable program
////////////////////////////////////////////////////////////////////////
//! mc::lt_RLTVar is a C++ structure for ordering variables in an RLT
//! map.
////////////////////////////////////////////////////////////////////////
template <typename T>
struct lt_RLTVar
////////////////////////////////////////////////////////////////////////
{
  bool operator()
    ( const typename FPRRLT<T>::t_RLTVar& Var1,
      const typename FPRRLT<T>::t_RLTVar& Var2 ) const
    {
      // Order RLT variables w.r.t. their RLT operation types first
      if( Var1.second < Var2.second ) return true;
      if( Var1.second > Var2.second ) return false;
      // Order RLT variables w.r.t. to their index next
      return lt_FPVar<T>()( Var1.first, Var2.first );
    }
};

///////////////////////////////// FPOuter //////////////////////////////////////

inline std::ostream&
operator <<
( std::ostream&out, const FPOuter&Int )
{
  const int iprec = 5;
  out << std::right << std::scientific << std::setprecision(iprec)
      << "Bounds: [ " << std::setw(iprec+7) << Int.xL() << " : "
      << std::setw(iprec+7) << Int.xU() << "]    OA Point: "
      << std::setw(iprec+7) << Int.xM() << "    OA Gap: "
      << std::setw(iprec+7) << Int.gap();
  return out;
}

////////////////////////////////// FPCut ///////////////////////////////////////

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPCut<T>&Cut )
{
  const int iprec = 5;
  out << std::right << std::scientific << std::setprecision(iprec);
  
  switch( Cut._type ){
    case FPCut<T>::EQ: case FPCut<T>::LE: case FPCut<T>::GE:
      for( unsigned int k=0; k<Cut.nvar(); k++ ){
        if( isequal( Cut._coef[k], 0. ) )
          out << " + " << std::setw(iprec+6) << 0.;
        else if( Cut._coef[k] > 0. )
          out << " + " << std::setw(iprec+6) << Cut._coef[k];
        else
          out << " - " << std::setw(iprec+6) << -Cut._coef[k];
        out << FPVar<T>::_name( Cut._idvar[k] );  
      }
      break;

    case FPCut<T>::SOS1: case FPCut<T>::SOS2:
      out << " {";
      for( unsigned int k=0; k<Cut.nvar(); k++ )
        out << " " << FPVar<T>::_name( Cut._idvar[k] );
      out << " }";
  }
  
  switch( Cut._type ){
    case FPCut<T>::EQ: out << " = "; break;
    case FPCut<T>::LE: out << " <= "; break;
    case FPCut<T>::GE: out << " >= "; break;
    case FPCut<T>::SOS1: out << " SOS1"; return out;
    case FPCut<T>::SOS2: out << " SOS2"; return out;
  }
  
  out << std::setw(iprec+6) << Cut._rhs;
  return out;
}

////////////////////////////////// FPVar ///////////////////////////////////////

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPVarNum<T>&VarNum )
{
  switch( VarNum.t ){
    case FPVarNum<T>::INT:
      out << "(I) " << VarNum.n ; break;
    case FPVarNum<T>::REAL:
      out << "(D) " << VarNum.x; break;
    case FPVarNum<T>::RANGE:
      out << VarNum.I; break;
  }
  return out;
}

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPVar<T>&Var)
{
  //std::ostringstream num; num << Var._num;
  //out << Var.name() << " <= " << std::left << std::setw(40) << num.str();
  // // << "\t(" << Var._FP << ")";
  //return out;

  out << Var.name()
      << " <= " << std::left << Var._num;// << "\t(" << Var._FP << ")";
  return out;
}

template <typename T> inline 
FPVar<T>::FPVar
( FPRelax<T>*FP, const McCormick<T>&MCX, const double*pref,
  const unsigned int*isub )
: _id( AUXCONT, FP->_naux++ ), _num( MCX.I() ), _FP( FP ), _Op(0)
{
  // Insert new auxiliary in _Vars and corresponding operation in _Ops
  FPVar<T>* pAux = new FPVar<T>( *this );
  _Op = new FPOp<T>( FPOp<T>::MC, 0, 0, pAux );
  //_Op->attach( MCX );
  pAux->_Op = _Op;
  FP->_Ops.insert( _Op );
  FP->_append_aux( pAux );

  // Build up polyhedral relaxation
  if( !pref || !MCX.nsub() ) return;
  FPVar<T> LAff = MCX.cv();
  FPVar<T> UAff = MCX.cc();
  for( unsigned int ip=0; ip<MCX.nsub(); ip++ ){
    FPVar<T>*pVar = FP->_find_var( pt_idVar( VARCONT, isub?isub[ip]:ip ) );
    if( !pVar ) pVar = FP->_find_var( pt_idVar( VARBIN, isub?isub[ip]:ip ) );
    // Throw exception if variable ix+ioff is undefined
    if( !pVar ) throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );
    LAff += MCX.cvsub(ip) * ( *pVar - pref[ip] );
    UAff += MCX.ccsub(ip) * ( *pVar - pref[ip] );
  }
  FP->add_constraint( *pAux, FPRelax<T>::GE, LAff );
  FP->add_constraint( *pAux, FPRelax<T>::LE, UAff );
}

template <typename T> inline 
FPVar<T>::FPVar
( FPRelax<T>*FP, const TVar<T>&TVX, TModel< FPVar<T> >*TMFP )
: _id( AUXCONT, FP->_naux++ ), _num( TVX.B() ), _FP( FP ), _Op(0)
{
  // Insert new auxiliary in _Vars and corresponding operation in _Ops
  FPVar<T>* pAux = new FPVar<T>( *this );
  _Op = new FPOp<T>( FPOp<T>::TM, 0, 0, pAux );
  //_Op->attach( TVX, TMFP );
  pAux->_Op = _Op;
  FP->_Ops.insert( _Op );
  FP->_append_aux( pAux );

  // Build up polyhedral relaxation
  typedef mc::TVar< FPVar<T> > TVFP;
  TVFP TPFPX( TMFP, TVX );
  FP->add_constraint( *pAux, FPRelax<T>::EQ, TPFPX.B() );
}

template <typename T> inline 
FPVar<T>::FPVar
( FPRelax<T>*FP, const TVar< McCormick<T> >&TVMCX, TModel< FPVar<T> >*TMFP,
  const double*MCXref, const unsigned int*isub )
: _id( AUXCONT, FP->_naux++ ), _num( TVMCX.B().I() ), _FP( FP ), _Op(0)
{
  // Insert new auxiliary in _Vars and corresponding operation in _Ops
  FPVar<T>* pAux = new FPVar<T>( *this );
  _Op = new FPOp<T>( FPOp<T>::TMMC, 0, 0, pAux );
  //_Op->attach( TVMCX, TMFP );
  pAux->_Op = _Op;
  FP->_Ops.insert( _Op );
  FP->_append_aux( pAux );

  // Build up polyhedral relaxation by separating remainder term and
  // polynomial part
  typedef mc::McCormick<T> MC;
  typedef mc::TVar<T> TV;
  struct loc{ static FPVar<T> conv( const MC&X ){ return FPVar<T>( X.I() ); } };
  typedef mc::TVar< FPVar<T> > TVFP;
  TVFP TPFPX( TMFP, TVMCX.P(), loc::conv );
  FPVar<T> TRFPX( FP, TVMCX.R(), MCXref, isub );
  FP->add_constraint( *pAux, FPRelax<T>::EQ, TPFPX.B() + TRFPX );
}

template <typename T> inline FPVar<T>&
FPVar<T>::operator=
( const FPVar<T>&Var )
{
  if( this == &Var )
    return *this;

  _id  = Var._id;
  _num = Var._num;
  _FP  = Var._FP;
  _Op  = Var._Op;
  return *this;
}

template <typename T> inline FPVar<T>&
FPVar<T>::operator=
( const T&I )
{
  _id  = std::make_pair(AUXCONT,NOREF);
  _num = I;
  _FP  = 0;
  _Op  = 0;
  return *this;
}

template <typename T> inline FPVar<T>&
FPVar<T>::operator=
( const int n )
{
  _id  = std::make_pair(AUXINT,NOREF);
  _num = n;
  _FP  = 0;
  _Op  = 0;
  return *this;
}

template <typename T> inline FPVar<T>&
FPVar<T>::operator=
( const double x )
{
  _id  = std::make_pair(AUXREAL,NOREF);
  _num = x;
  _FP  = 0;
  _Op  = 0;
  return *this;
}

template <typename T> inline FPVar<T>
operator+
( const FPVar<T>&Var )
{
  return Var;
}

template <typename T> template <typename U> inline FPVar<T>&
FPVar<T>::operator+=
( const U&Var )
{
  FPVar<T> VarNew = *this + Var;
  *this = VarNew;
  return *this;
}

template <typename T> inline FPVar<T>
operator+
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{ 
  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.n + Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.n + Var2._num.x );
          case FPVarNum<T>::RANGE: return( (double)Var1._num.n + Var2._num.I );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.x + Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.x + Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.x + Var2._num.I );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.I + (double)Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.I + Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.I + Var2._num.I );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( Var2 + (double)Var1._num.n );
      case FPVarNum<T>::REAL:  return( Var2 + Var1._num.x );
      case FPVarNum<T>::RANGE: return( Var2 + Var1._num.I );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Var1 + (double)Var2._num.n );
      case FPVarNum<T>::REAL:  return( Var1 + Var2._num.x );
      case FPVarNum<T>::RANGE: return( Var1 + Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::PLUS, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I+Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
operator+
( const FPVar<T>&Var1, const T&Bnd2 )
{
  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( (double)Var1._num.n + Bnd2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x + Bnd2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I + Bnd2 );
    }
  }

  // Append new interval Bnd2
  FPRelax<T>* pFP = Var1._FP;
  FPVar<T>* pBnd2 = pFP->_auxiliary_constant( Bnd2 );
  return( Var1 + *pBnd2 );
}

template <typename T> inline FPVar<T>
operator+
( const T&Bnd1, const FPVar<T>&Var2 )
{
  return( Var2 + Bnd1 );
}

template <typename T,typename U> inline FPVar<T>
operator+
( const FPVar<T>&Var1, const U&Cst2 )
{
  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( Var1._num.n + Cst2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x + Cst2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I + (double)Cst2 );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var1._FP;
  FPVar<T>* pCst2 = pFP->_auxiliary_constant( Cst2 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::PLUS, Var1._Op->pres, pCst2->_Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I+(double)Cst2, pOp );
  return *pOp->pres;
}

template <typename T,typename U> inline FPVar<T>
operator+
( const U&Cst1, const FPVar<T>&Var2 )
{
  return( Var2 + Cst1 );
}

template <typename T> inline FPVar<T>
operator-
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( -Var._num.n );
      case FPVarNum<T>::REAL:  return( -Var._num.x );
      case FPVarNum<T>::RANGE: return( -Var._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::NEG, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( -Var._num.I, pOp );
  return *pOp->pres;
}

template <typename T> template <typename U> inline FPVar<T>&
FPVar<T>::operator-=
( const U&Var )
{
  FPVar<T> VarNew = *this - Var;
  *this = VarNew;
  return *this;
}

template <typename T> inline FPVar<T>
operator-
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return 0.;

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.n - Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.n - Var2._num.x );
          case FPVarNum<T>::RANGE: return( (double)Var1._num.n - Var2._num.I );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.x - Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.x - Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.x - Var2._num.I );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.I - (double)Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.I - Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.I - Var2._num.I );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( (double)Var1._num.n - Var2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x - Var2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I - Var2 );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Var1 - (double)Var2._num.n );
      case FPVarNum<T>::REAL:  return( Var1 - Var2._num.x );
      case FPVarNum<T>::RANGE: return( Var1 - Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MINUS, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I-Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
operator-
( const FPVar<T>&Var1, const T&Bnd2 )
{
  return( Var1 + (-Bnd2) );
}

template <typename T> inline FPVar<T>
operator-
( const T&Bnd1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Bnd1 - (double)Var2._num.n );
      case FPVarNum<T>::REAL:  return( Bnd1 - Var2._num.x );
      case FPVarNum<T>::RANGE: return( Bnd1 - Var2._num.I );
    }
  }

  // Append new interval Bnd1
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pBnd1 = pFP->_auxiliary_constant( Bnd1 );
  return( *pBnd1 - Var2 );
}

template <typename T,typename U> inline FPVar<T>
operator-
( const FPVar<T>&Var1, const U&Cst2 )
{
  return( Var1 + (-Cst2) );
}

template <typename T,typename U> inline FPVar<T>
operator-
( const U&Cst1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Cst1 - Var2._num.n );
      case FPVarNum<T>::REAL:  return( Cst1 - Var2._num.x );
      case FPVarNum<T>::RANGE: return( Cst1 - Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pCst1 = pFP->_auxiliary_constant( Cst1 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MINUS, pCst1->_Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Cst1-Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> template <typename U> inline FPVar<T>&
FPVar<T>::operator*=
( const U&Var )
{
  FPVar<T> VarNew = *this * Var;
  *this = VarNew;
  return *this;
}

template <typename T> inline FPVar<T>
operator*
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return sqr(Var1);

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.n * Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.n * Var2._num.x );
          case FPVarNum<T>::RANGE: return( (double)Var1._num.n * Var2._num.I );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.x * Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.x * Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.x * Var2._num.I );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.I * (double)Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.I * Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.I * Var2._num.I );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( (double)Var1._num.n * Var2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x * Var2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I * Var2 );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Var1 * (double)Var2._num.n );
      case FPVarNum<T>::REAL:  return( Var1 * Var2._num.x );
      case FPVarNum<T>::RANGE: return( Var1 * Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::TIMES, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I*Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
operator*
( const T&Bnd1, const FPVar<T>&Var2 )
{
  return( Var2 * Bnd1 );
}

template <typename T> inline FPVar<T>
operator*
( const FPVar<T>&Var1, const T&Bnd2 )
{
  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( (double)Var1._num.n * Bnd2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x * Bnd2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I * Bnd2 );
    }
  }

  // Append new interval Bnd2
  FPRelax<T>* pFP = Var1._FP;
  FPVar<T>* pBnd2 = pFP->_auxiliary_constant( Bnd2 );
  return( Var1 * *pBnd2 );
}

template <typename T,typename U> inline FPVar<T>
operator*
( const FPVar<T>&Var1, const U&Cst2 )
{
  return( Cst2 * Var1 );
}

template <typename T,typename U> inline FPVar<T>
operator*
( const U&Cst1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Cst1 * Var2._num.n );
      case FPVarNum<T>::REAL:  return( Cst1 * Var2._num.x );
      case FPVarNum<T>::RANGE: return( (double)Cst1 * Var2._num.I );
    }
  }
  //if( mc::isequal( Cst1, 0. ) ) return( 0. );

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pCst1 = pFP->_auxiliary_constant( Cst1 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::SCALE, pCst1->_Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( (double)Cst1*Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> template <typename U> inline FPVar<T>&
FPVar<T>::operator/=
( const U&Var )
{
  FPVar<T> VarNew = *this / Var;
  *this = VarNew;
  return *this;
}

template <typename T> inline FPVar<T>
operator/
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return 1.;

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.n / Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.n / Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.n / Var2._num.I );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.x / Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.x / Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.x / Var2._num.I );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Var1._num.I / Var2._num.n );
          case FPVarNum<T>::REAL:  return( Var1._num.I / Var2._num.x );
          case FPVarNum<T>::RANGE: return( Var1._num.I / Var2._num.I );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( Var1._num.n / Var2 );
      case FPVarNum<T>::REAL:  return( Var1._num.x / Var2 );
      case FPVarNum<T>::RANGE: return( Var1._num.I / Var2 );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Var1 / Var2._num.n );
      case FPVarNum<T>::REAL:  return( Var1 / Var2._num.x );
      case FPVarNum<T>::RANGE: return( Var1 / Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::DIV, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I/Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
operator/
( const FPVar<T>&Var1, const T&Int2 )
{
  return( ( 1 / Int2 ) * Var1 );
}

template <typename T> inline FPVar<T>
operator/
( const T&Bnd1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Bnd1 / Var2._num.n );
      case FPVarNum<T>::REAL:  return( Bnd1 / Var2._num.x );
      case FPVarNum<T>::RANGE: return( Bnd1 / Var2._num.I );
    }
  }

  // Append new interval Bnd1
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pBnd1 = pFP->_auxiliary_constant( Bnd1 );
  return( *pBnd1 / Var2 );
}

template <typename T, typename U> inline FPVar<T>
operator/
( const FPVar<T>&Var1, const U&Cst2 )
{
  return( ( 1 / Cst2 ) * Var1 );
}

template <typename T, typename U> inline FPVar<T>
operator/
( const U&Cst1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( Cst1 / Var2._num.n );
      case FPVarNum<T>::REAL:  return( Cst1 / Var2._num.x );
      case FPVarNum<T>::RANGE: return( Cst1 / Var2._num.I );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pCst1 = pFP->_auxiliary_constant( Cst1 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::DIV, pCst1->_Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Cst1/Var2._num.I, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
inv
( const FPVar<T>&Var )
{
  return( 1 / Var );
}

template <typename T> inline FPVar<T>
operator^
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    T Inter;
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:
            if ( Var1._num.n == Var2._num.n ) return( Var1._num.n );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::REAL:
            if( isequal( Var1._num.n, Var2._num.x ) ) return( Var1._num.n );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::RANGE:
            if( Op<T>::inter( Inter, Var1._num.n, Var2._num.I ) ) return( Var1._num.n );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:
            if( isequal( Var1._num.x, Var2._num.n ) ) return( Var2._num.n );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::REAL:
            if( isequal( Var1._num.x, Var2._num.x ) ) return( Var1._num.x );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::RANGE:
            if( Op<T>::inter( Inter, Var1._num.x, Var2._num.I ) ) return( Var1._num.x );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:
            if( Op<T>::inter( Inter, Var1._num.I, Var2._num.n ) ) return( Var2._num.n );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::REAL:
            if( Op<T>::inter( Inter, Var1._num.I, Var2._num.x ) ) return( Var2._num.x );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
          case FPVarNum<T>::RANGE:
            if( Op<T>::inter( Inter, Var1._num.I, Var2._num.I ) ) return( Inter );
            else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    T Inter;
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        if( Op<T>::inter( Inter, Var1._num.n, Var2._num.I ) ) return( Var1._num.n );
        else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
      case FPVarNum<T>::REAL:
        if( Op<T>::inter( Inter, Var1._num.x, Var2._num.I ) ) return( Var1._num.x );
        else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
      case FPVarNum<T>::RANGE:
        if( Op<T>::inter( Inter, Var1._num.I, Var2._num.I ) ) return( Inter );
        else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
    }
  }

  if( Var2._id.second == FPVar<T>::NOREF ){
    T Inter;
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:
        if( Op<T>::inter( Inter, Var1._num.I, Var2._num.n ) ) return( Var2._num.n );
        else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
      case FPVarNum<T>::REAL:
        if( Op<T>::inter( Inter, Var1._num.I, Var2._num.x ) ) return( Var2._num.x );
        else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
      case FPVarNum<T>::RANGE:
        if( Op<T>::inter( Inter, Var1._num.I, Var2._num.I ) ) return( Inter );
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::INTER, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  T Inter;
  if( !Op<T>::inter( Inter, Var1._num.I, Var2._num.I ) )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTER );
  pOp->pres = pFP->_auxiliary_variable( Inter, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
max
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( max(Var1._num.n,Var2._num.n) );
          case FPVarNum<T>::REAL:  return( max(Var1._num.n,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::max((double)Var1._num.n,Var2._num.I) );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( max(Var1._num.x,Var2._num.n) );
          case FPVarNum<T>::REAL:  return( max(Var1._num.x,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::max(Var1._num.x,Var2._num.I) );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Op<T>::max(Var1._num.I,(double)Var2._num.n) );
          case FPVarNum<T>::REAL:  return( Op<T>::max(Var1._num.I,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::max(Var1._num.I,Var2._num.I) );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( max((double)Var1._num.n,Var2) );
      case FPVarNum<T>::REAL:  return( max(Var1._num.x,Var2) );
      case FPVarNum<T>::RANGE: return( max(Var1._num.I,Var2) );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( max((double)Var2._num.n,Var1) );
      case FPVarNum<T>::REAL:  return( max(Var2._num.x,Var1) );
      case FPVarNum<T>::RANGE: return( max(Var2._num.I,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MAXF, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::max(Var1._num.I,Var2._num.I), pOp );
  return *pOp->pres;
}

template <typename T,typename U> inline FPVar<T>
max
( const FPVar<T>&Var1, const U&Cst2 )
{
  return( max( Cst2, Var1 ) );
}

template <typename T,typename U> inline FPVar<T>
max
( const U&Cst1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( max(Cst1,Var2._num.n) );
      case FPVarNum<T>::REAL:  return( max(Cst1,Var2._num.x) );
      case FPVarNum<T>::RANGE: return( Op<T>::max((double)Cst1,Var2._num.I) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pCst1 = pFP->_auxiliary_constant( Cst1 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MAXF, pCst1->_Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::max((double)Cst1,Var2._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
max
( const T&Bnd1, const FPVar<T>&Var2 )
{
  return( max( Var2, Bnd1 ) );
}

template <typename T> inline FPVar<T>
max
( const FPVar<T>&Var1, const T&Bnd2 )
{
  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( Op<T>::max((double)Var1._num.n,Bnd2) );
      case FPVarNum<T>::REAL:  return( Op<T>::max(Var1._num.x,Bnd2) );
      case FPVarNum<T>::RANGE: return( Op<T>::max(Var1._num.I,Bnd2) );
    }
  }

  // Append new interval Bnd2
  FPRelax<T>* pFP = Var1._FP;
  FPVar<T>* pBnd2 = pFP->_auxiliary_constant( Bnd2 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MAXF, Var1._Op->pres, pBnd2->_Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::max(Var1._num.I,Bnd2), pOp );
  return *pOp->pres;
  //return( max(Var1,*pBnd2) );
}

template <typename T> inline FPVar<T>
max
( const unsigned nVar, const FPVar<T>*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FPVar<T> VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = max( VarR, pVar[i] );
  return( VarR );
}

template <typename T> inline FPVar<T>
min
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  if( &Var1 == &Var2 ) return Var1;

  if( Var1._id.second == FPVar<T>::NOREF
   && Var2._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( min(Var1._num.n,Var2._num.n) );
          case FPVarNum<T>::REAL:  return( min(Var1._num.n,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::min((double)Var1._num.n,Var2._num.I) );
        }
      case FPVarNum<T>::REAL:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( min(Var1._num.x,Var2._num.n) );
          case FPVarNum<T>::REAL:  return( min(Var1._num.x,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::min(Var1._num.x,Var2._num.I) );
        }
      case FPVarNum<T>::RANGE:
        switch( Var2._num.t ){
          case FPVarNum<T>::INT:   return( Op<T>::min(Var1._num.I,(double)Var2._num.n) );
          case FPVarNum<T>::REAL:  return( Op<T>::min(Var1._num.I,Var2._num.x) );
          case FPVarNum<T>::RANGE: return( Op<T>::min(Var1._num.I,Var2._num.I) );
        }
    }
  }

  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( min((double)Var1._num.n,Var2) );
      case FPVarNum<T>::REAL:  return( min(Var1._num.x,Var2) );
      case FPVarNum<T>::RANGE: return( min(Var1._num.I,Var2) );
    }
  }
  
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( min((double)Var2._num.n,Var1) );
      case FPVarNum<T>::REAL:  return( min(Var2._num.x,Var1) );
      case FPVarNum<T>::RANGE: return( min(Var2._num.I,Var1) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  if( Var1._FP != Var2._FP )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::FP );
  FPRelax<T>* pFP = Var1._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MINF, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::min(Var1._num.I,Var2._num.I), pOp );
  return *pOp->pres;
}

template <typename T,typename U> inline FPVar<T>
min
( const FPVar<T>&Var1, const U&Cst2 )
{
  return( min( Cst2, Var1 ) );
}

template <typename T,typename U> inline FPVar<T>
min
( const U&Cst1, const FPVar<T>&Var2 )
{
  if( Var2._id.second == FPVar<T>::NOREF ){
    switch( Var2._num.t ){
      case FPVarNum<T>::INT:   return( min(Cst1,Var2._num.n) );
      case FPVarNum<T>::REAL:  return( min(Cst1,Var2._num.x) );
      case FPVarNum<T>::RANGE: return( Op<T>::min((double)Cst1,Var2._num.I) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append constant Cst2 if not defined
  FPRelax<T>* pFP = Var2._FP;
  FPVar<T>* pCst1 = pFP->_auxiliary_constant( Cst1 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MINF, pCst1->_Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::min((double)Cst1,Var2._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
min
( const T&Bnd1, const FPVar<T>&Var2 )
{
  return( min( Var2, Bnd1 ) );
}

template <typename T> inline FPVar<T>
min
( const FPVar<T>&Var1, const T&Bnd2 )
{
  if( Var1._id.second == FPVar<T>::NOREF ){
    switch( Var1._num.t ){
      case FPVarNum<T>::INT:   return( Op<T>::min((double)Var1._num.n,Bnd2) );
      case FPVarNum<T>::REAL:  return( Op<T>::min(Var1._num.x,Bnd2) );
      case FPVarNum<T>::RANGE: return( Op<T>::min(Var1._num.I,Bnd2) );
    }
  }

  // Append new interval Bnd2
  FPRelax<T>* pFP = Var1._FP;
  FPVar<T>* pBnd2 = pFP->_auxiliary_constant( Bnd2 );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::MINF, Var1._Op->pres, pBnd2->_Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::min(Var1._num.I,Bnd2), pOp );
  return *pOp->pres;
  //return( Var1 * *pBnd2 );
}

template <typename T> inline FPVar<T>
min
( const unsigned nVar, const FPVar<T>*pVar )
{
  if( nVar<=0 || !pVar ) return( 0 );
  
  FPVar<T> VarR = pVar[0];
  for( unsigned int i=1; i<nVar; i++ ) VarR = min( VarR, pVar[i] );
  return( VarR );
}

template <typename T> inline FPVar<T>
exp
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::exp( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::exp( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::exp( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::EXP, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::exp(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
log
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::log( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::log( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::log( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::LOG, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::log(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
sqr
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( mc::sqr( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( mc::sqr( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::sqr( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::SQR, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::sqr(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
sqrt
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::sqrt( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::sqrt( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::sqrt( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::SQRT, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::sqrt(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
fabs
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::fabs( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::fabs( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::fabs( Var._num.I ) );
    }
  }

  if( Op<T>::l(Var._num.I) >= 0. ) return( Var );
  if( Op<T>::u(Var._num.I) <= 0. ) return( -Var );

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::FABS, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::fabs(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
pow
( const FPVar<T>&Var, const int iExp )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::pow( Var._num.n, iExp ) );
      case FPVarNum<T>::REAL:  return( std::pow( Var._num.x, iExp ) );
      case FPVarNum<T>::RANGE: return( Op<T>::pow( Var._num.I, iExp ) );
    }
  }

  if( iExp == 0 ) return( 1. );
  if( iExp == 1 ) return Var;
  if( iExp == 2 ) return sqr(Var);
  if( iExp >= 3 ){
    FPRelax<T>* pFP = Var._FP;
    if( !(iExp%2) || pFP->options.NEWTON_USE ){ 
      // Append new intermediate variable and corresponding operation
      // (only if necessary); Also append integer exponent iExp if not defined
      FPVar<T>* piExp = pFP->_auxiliary_constant( iExp );
      FPOp<T>* pOp = pFP->_operation( FPOp<T>::IPOW, Var._Op->pres, piExp->_Op->pres );
      if( pOp->pres ) return *pOp->pres;
      pOp->pres = pFP->_auxiliary_variable( Op<T>::pow(Var._num.I,iExp), pOp );
      return *pOp->pres;
    }
    return pow( Var, iExp-1 ) * Var;
  }
  // Append new intermediate variable and corresponding operation
  // (only if necessary); Also append integer exponent iExp if not defined
  FPRelax<T>* pFP = Var._FP;
  FPVar<T>* piExp = pFP->_auxiliary_constant( iExp );
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::IPOW, Var._Op->pres, piExp->_Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::pow(Var._num.I,iExp), pOp );
  return *pOp->pres;
  //return inv( pow( Var, -iExp ) );
}

template <typename T, typename U> inline FPVar<T>
pow
( const FPVar<T>&Var1, const U&Var2 )
{
  return exp( Var2 * log( Var1 ) );
}

template <typename T> inline FPVar<T>
pow
( const double Var1, const FPVar<T>&Var2 )
{
  return exp( Var2 * std::log( Var1 ) );
}

template <typename T> inline FPVar<T>
cos
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::cos( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::cos( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::cos( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::COS, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::cos(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
sin
( const FPVar<T>&Var )
{
  return cos( Var - PI/2. );
}

template <typename T> inline FPVar<T>
tan
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::tan( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::tan( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::tan( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::TAN, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::tan(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
asin
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::asin( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::asin( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::asin( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::ASIN, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::asin(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
acos
( const FPVar<T>&Var )
{
  return( asin( -Var ) + PI/2. );
}

template <typename T> inline FPVar<T>
atan
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( std::atan( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( std::atan( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::atan( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::ATAN, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::atan(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
erf
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( ::erf( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( ::erf( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::erf( Var._num.I ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::ERF, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Op<T>::erf(Var._num.I), pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
erfc
( const FPVar<T>&Var )
{
  return ( 1. - erf( Var ) );
}

template <typename T> inline FPVar<T>
fstep
( const FPVar<T>&Var )
{
  if( Var._id.second == FPVar<T>::NOREF ){
    switch( Var._num.t ){
      case FPVarNum<T>::INT:   return( mc::fstep( Var._num.n ) );
      case FPVarNum<T>::REAL:  return( mc::fstep( Var._num.x ) );
      case FPVarNum<T>::RANGE: return( Op<T>::l( Var._num.I )>=0? T(1.):
        ( Op<T>::u( Var._num.I )<0? T(0.): Op<T>::zeroone() ) );
    }
  }

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPRelax<T>* pFP = Var._FP;
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::FSTEP, Var._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  T Bres = Op<T>::l( Var._num.I )>=0? T(1.): ( Op<T>::u( Var._num.I )<0?
    T(0.): Op<T>::zeroone() );
  pOp->pres = pFP->_auxiliary_variable( Bres, pOp );
  return *pOp->pres;
}

template <typename T> inline FPVar<T>
bstep
( const FPVar<T>&Var )
{
  return ( fstep( -Var ) );
}

template <typename T> inline FPVar<T> bilin
( const FPVar<T>&Var1, const FPVar<T>&Var2 )
{
  FPRelax<T>* pFP = Var1._FP;
  const unsigned int Ndiv = pFP->options.BILINEAR_SUBDIV;
  if( Ndiv <= 1 || Var1._id.first == FPVar<T>::AUXINT
   || Var1._id.first == FPVar<T>::AUXREAL || Var2._id.first == FPVar<T>::AUXINT
   || Var2._id.first == FPVar<T>::AUXREAL ) return( Var1 * Var2 );

  // Append new intermediate variable and corresponding operation
  // (only if necessary)
  FPOp<T>* pOp = pFP->_operation( FPOp<T>::BILIN, Var1._Op->pres, Var2._Op->pres );
  if( pOp->pres ) return *pOp->pres;
  pOp->pres = pFP->_auxiliary_variable( Var1._num.I*Var2._num.I, pOp );
  return *pOp->pres;
}

////////////////////////////////// FPOp ////////////////////////////////////////

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPOp<T>&Op)
{
  switch( Op.type ){
  case FPOp<T>::CNST:
    out << "CONSTANT";
    break;
  case FPOp<T>::RANGE:
  case FPOp<T>::VAR:
    out << "VARIABLE";
    break;
  case FPOp<T>::MC:
    out << "MCCORMICK";
    break;
  case FPOp<T>::TM:
    out << "TAYLOR";
    break;
  case FPOp<T>::TMMC:
    out << "MCCORMICK-TAYLOR";
    break;
  case FPOp<T>::MIN:
    out << "MIN " << FPVar<T>::_name( Op.plop->id() );
    break;
  case FPOp<T>::MAX:
    out << "MAX " << FPVar<T>::_name( Op.plop->id() );
    break;
  case FPOp<T>::EQ:
    out << FPVar<T>::_name( Op.plop->id() ) << " == "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::LE:
    out << FPVar<T>::_name( Op.plop->id() ) << " <= "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::PLUS:
    out << FPVar<T>::_name( Op.plop->id() ) << " + "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::NEG:
    out << "- " << FPVar<T>::_name( Op.plop->id() );
    break;
  case FPOp<T>::MINUS:
    out << FPVar<T>::_name( Op.plop->id() ) << " - "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::TIMES: case FPOp<T>::SCALE:
    out << FPVar<T>::_name( Op.plop->id() ) << " * "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::DIV:
    out << FPVar<T>::_name( Op.plop->id() ) << " / "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::INTER:
    out << FPVar<T>::_name( Op.plop->id() ) << " ^ "
        << FPVar<T>::_name( Op.prop->id() );
    break;
  case FPOp<T>::MINF:
    out << "MIN( " << FPVar<T>::_name( Op.plop->id() ) << ", "
        << FPVar<T>::_name( Op.prop->id() ) << " )";
    break;
  case FPOp<T>::MAXF:
    out << "MIN( " << FPVar<T>::_name( Op.plop->id() ) << ", "
        << FPVar<T>::_name( Op.prop->id() ) << " )";
    break;
  case FPOp<T>::BILIN:
    out << "BILIN( " << FPVar<T>::_name( Op.plop->id() ) << ", "
        << FPVar<T>::_name( Op.prop->id() ) << " )";
    break;
  case FPOp<T>::EXP:
    out << "EXP( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::LOG:
    out << "LOG( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::SQR:
    out << "SQR( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::SQRT:
    out << "SQRT( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::FABS:
    out << "FABS( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::IPOW:
    out << "POW( " << FPVar<T>::_name( Op.plop->id() ) << ", "
        << FPVar<T>::_name( Op.prop->id() ) << " )";
    break;
  case FPOp<T>::COS:
    out << "COS( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::TAN:
    out << "TAN( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::ASIN:
    out << "ASIN( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::ATAN:
    out << "ATAN( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::ERF:
    out << "ERF( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  case FPOp<T>::FSTEP:
    out << "FSTEP( " << FPVar<T>::_name( Op.plop->id() ) << " )";
    break;
  default:;
  } 
  return out;
}

template <typename T> inline bool
FPOp<T>::propagate_bounds
( FPRelax<T>*pFP )
{
  if( !tighten_bounds( pFP ) ) return false;
  if( plop && plop->Oper() && !plop->Oper()->propagate_bounds( pFP ) )
    return false;
  if( prop && prop->Oper() && !prop->Oper()->propagate_bounds( pFP ) )
    return false;
  return true;
}

template <typename T> inline bool
FPOp<T>::tighten_bounds
( FPRelax<T>*pFP )
{
  switch( type ){
    
  case FPOp<T>::MIN:
  case FPOp<T>::MAX:   return true;
  case FPOp<T>::EQ:    return _EQ_bounds( pFP, plop, prop );
  case FPOp<T>::LE:    return _LE_bounds( pFP, plop, prop );

  case FPOp<T>::PLUS:  return _PLUS_bounds( pFP, pres, plop, prop );
  case FPOp<T>::NEG:   return _NEG_bounds( pFP, pres, plop );
  case FPOp<T>::MINUS: return _MINUS_bounds( pFP, pres, plop, prop );
  case FPOp<T>::TIMES:
  case FPOp<T>::SCALE:
  case FPOp<T>::BILIN: return _TIMES_bounds( pFP, pres, plop, prop );
  case FPOp<T>::DIV:   return _DIV_bounds( pFP, pres, plop, prop );
  case FPOp<T>::MINF: 
  case FPOp<T>::MAXF: 
  case FPOp<T>::INTER: return true;

  case FPOp<T>::EXP:   return _EXP_bounds( pFP, pres, plop );
  case FPOp<T>::LOG:   return _LOG_bounds( pFP, pres, plop );
  case FPOp<T>::SQR:   return _SQR_bounds( pFP, pres, plop );
  case FPOp<T>::SQRT:  return _SQRT_bounds( pFP, pres, plop );
  case FPOp<T>::FABS:  return _FABS_bounds( pFP, pres, plop );
  case FPOp<T>::IPOW:  return _IPOW_bounds( pFP, pres, plop, prop );

  case FPOp<T>::CNST:
  case FPOp<T>::RANGE:
  case FPOp<T>::VAR:
  case FPOp<T>::MC:
  case FPOp<T>::TM:
  case FPOp<T>::TMMC:  return true;

  default:
    return true; //throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
  }
}

template <typename T> inline bool
FPOp<T>::_EQ_bounds
( FPRelax<T>*pFP, FPVar<T>*Var1, FPVar<T>*Var2 )
{
  T IV1V2;
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: EQ" << std::endl
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  if( !Op<T>::inter( IV1V2, Var1->I(), Var2->I() ) ) return false;
  pFP->_update_var( Var1, IV1V2 );
  pFP->_update_var( Var2, IV1V2 );
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_LE_bounds
( FPRelax<T>*pFP, FPVar<T>*Var1, FPVar<T>*Var2 )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: LE" << std::endl
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  if( Op<T>::u(Var2->I()) < Op<T>::l(Var1->I()) ) return false;
  if( Op<T>::u(Var1->I()) > Op<T>::u(Var2->I()) )
    pFP->_update_var( Var1, T( Op<T>::l(Var1->I()), Op<T>::u(Var2->I()) ) );
  if( Op<T>::l(Var2->I()) < Op<T>::l(Var1->I()) )
    pFP->_update_var( Var2, T( Op<T>::l(Var1->I()), Op<T>::u(Var2->I()) ) );
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_NEG_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: NEG" << std::endl
            << "  BEF.  " << Var->name() << "=" << Var->I() << std::endl;
#endif
  T VarRef;
  if( !Op<T>::inter( VarRef, Var->I(), -VarR->I() ) ) return false;
  pFP->_update_var( Var, VarRef );
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var->name() << "=" << Var->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_PLUS_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1, FPVar<T>*Var2 )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: PLUS" << std::endl
            << "  RES.  " << VarR->name() << "=" << VarR->I() << "  "
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  T VarRef;
  if( !Op<T>::inter( VarRef, Var1->I(), VarR->I()-Var2->I() ) ) return false;
  pFP->_update_var( Var1, VarRef );
  if( !Op<T>::inter( VarRef, Var2->I(), VarR->I()-Var1->I() ) ) return false;
  pFP->_update_var( Var2, VarRef );
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_MINUS_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1, FPVar<T>*Var2 )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: MINUS" << std::endl
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  T VarRef;
  if( !Op<T>::inter( VarRef, Var1->I(), VarR->I()+Var2->I() ) ) return false;
  pFP->_update_var( Var1, VarRef );
  if( !Op<T>::inter( VarRef, Var2->I(), Var1->I()-VarR->I() ) ) return false;
  pFP->_update_var( Var2, VarRef );
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_TIMES_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1, FPVar<T>*Var2 )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: TIMES" << std::endl
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  T VarRef;
  if( !Op<T>::inter( VarRef, Var2->I(), 0. ) ){
    if( !Op<T>::inter( VarRef, Var1->I(), VarR->I()/Var2->I() ) ) return false;
    pFP->_update_var( Var1, VarRef );
  }
  if( !Op<T>::inter( VarRef, Var1->I(), 0. ) ){
    if( !Op<T>::inter( VarRef, Var2->I(), VarR->I()/Var1->I() ) ) return false;
  pFP->_update_var( Var2, VarRef );
  }
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_DIV_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var1, FPVar<T>*Var2 )
{
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nCONSTRAINT PROPAGATION: DIV" << std::endl
            << "  BEF.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  T VarRef;
  if( !Op<T>::inter( VarRef, Var1->I(), VarR->I()*Var2->I() ) ) return false;
  pFP->_update_var( Var1, VarRef );
  if( !Op<T>::inter( VarRef, VarR->I(), 0. ) ){
    if( !Op<T>::inter( VarRef, Var2->I(), Var1->I()/VarR->I() ) ) return false;
  pFP->_update_var( Var2, VarRef );
  }
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "  AFT.  " << Var1->name() << "=" << Var1->I() << "  "
                          << Var2->name() << "=" << Var2->I() << std::endl;
#endif
  return true;
}

template <typename T> inline bool
FPOp<T>::_EXP_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( !Op<T>::inter( VarRef, Var->I(), Op<T>::log(VarR->I()) ) ) return false;
  pFP->_update_var( Var, VarRef );
  return true;
}

template <typename T> inline bool
FPOp<T>::_LOG_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( !Op<T>::inter( VarRef, Var->I(), Op<T>::exp(VarR->I()) ) ) return false;
  pFP->_update_var( Var, VarRef );
  return true;
}

template <typename T> inline bool
FPOp<T>::_SQR_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( Op<T>::l(Var->I()) >= 0
   && !Op<T>::inter( VarRef, Var->I(), Op<T>::sqrt(VarR->I()) ) )
      return false;
  if( Op<T>::u(Var->I()) <= 0
   && !Op<T>::inter( VarRef, Var->I(), -Op<T>::sqrt(VarR->I()) ) )
    return false;
  if( !Op<T>::inter( VarRef, Var->I(), Op<T>::sqrt(VarR->I())
   *(2.*Op<T>::zeroone()-1.) ) )
    return false;
  pFP->_update_var( Var, VarRef );
  return true;
}

template <typename T> inline bool
FPOp<T>::_SQRT_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( !Op<T>::inter( VarRef, Var->I(), Op<T>::sqr(VarR->I()) ) ) return false;
  pFP->_update_var( Var, VarRef );
  return true;
}

template <typename T> inline bool
FPOp<T>::_FABS_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( !Op<T>::inter( VarRef, Var->I(), VarR->I() ) ) return false;
  if( Op<T>::u(VarR->I()) < Op<T>::u(Var->I()) )
    pFP->_update_var( Var, T( Op<T>::l(Var->I()), Op<T>::u(VarR->I()) ) );
  return true;
}

template <typename T> inline bool
FPOp<T>::_IPOW_bounds
( FPRelax<T>*pFP, const FPVar<T>*VarR, FPVar<T>*Var,
  const FPVar<T>*Exp )
{
  const int iExp = Exp->num().n;
  if( iExp == 0 || iExp == 1 || iExp == 2 )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  if( Var->id().first == FPVar<T>::AUXINT
   || Var->id().first == FPVar<T>::AUXREAL ) return true;

  T VarRef;
  if( iExp>0 && !(iExp%2) ){
    if( Op<T>::inter( VarRef, VarR->I(), 0. ) ){
      if( Op<T>::l(Var->I()) > 0 ){
        if( !Op<T>::inter( VarRef, Var->I(),
          Op<T>::pow(VarR->I(),1./(T)iExp) ) ) return false;
        pFP->_update_var( Var, VarRef );
        return true;
      }
      else if( Op<T>::u(Var->I()) < 0 ){
        if( !Op<T>::inter( VarRef, Var->I(),
          -Op<T>::pow(VarR->I(),1./(T)iExp) ) ) return false;
        pFP->_update_var( Var, VarRef );
        return true;
      }
    }
    if( !Op<T>::inter( VarRef, Var->I(), (2.*Op<T>::zeroone()-1.)
      *std::pow(Op<T>::u(VarR->I()),1./(double)iExp) ) ) return false;
    pFP->_update_var( Var, VarRef );
    return true;
  }

  else if( iExp > 0 ){
    if( !Op<T>::inter( VarRef, Var->I(), T(mc::sign(Op<T>::l(VarR->I()))
      *std::pow(std::fabs(Op<T>::l(VarR->I())),1./(double)iExp),
      mc::sign(Op<T>::u(VarR->I()))
      *std::pow(std::fabs(Op<T>::u(VarR->I())),1./(double)iExp)) ) )
      return false;
    pFP->_update_var( Var, VarRef );
  }

  return true;
}

template <typename T> inline void
FPOp<T>::generate_cuts
( FPRelax<T>*pFP )
{
  if( _visited ) return;
  _visited = true;
  
  append_cuts( pFP );
  if( plop && plop->Oper() ) plop->Oper()->generate_cuts( pFP );
  if( prop && prop->Oper() ) prop->Oper()->generate_cuts( pFP );
  return;
}

template <typename T> inline void
FPOp<T>::append_cuts
( FPRelax<T>*pFP )
{
  switch( type ){
    
  case FPOp<T>::MIN:
  case FPOp<T>::MAX:   break;
  case FPOp<T>::EQ:    _EQ_cuts( pFP, plop, prop );          break;
  case FPOp<T>::LE:    _LE_cuts( pFP, plop, prop );          break;

  case FPOp<T>::PLUS:  _PLUS_cuts( pFP, pres, plop, prop );  break;
  case FPOp<T>::NEG:   _NEG_cuts( pFP, pres, plop );         break;
  case FPOp<T>::MINUS: _MINUS_cuts( pFP, pres, plop, prop ); break;
  case FPOp<T>::TIMES:
  case FPOp<T>::SCALE: _TIMES_cuts( pFP, pres, plop, prop ); break;
  case FPOp<T>::DIV:   _DIV_cuts( pFP, pres, plop, prop );   break;
  case FPOp<T>::INTER: _INTER_cuts( pFP, pres, plop, prop ); break;
  case FPOp<T>::MINF:  _MINF_cuts( pFP, pres, plop, prop );  break;
  case FPOp<T>::MAXF:  _MAXF_cuts( pFP, pres, plop, prop );  break;
  case FPOp<T>::BILIN: _BILIN_cuts( pFP, pres, plop, prop ); break;

  case FPOp<T>::EXP:   _EXP_cuts( pFP, pres, plop );         break;
  case FPOp<T>::LOG:   _LOG_cuts( pFP, pres, plop );         break;
  case FPOp<T>::SQR:   _SQR_cuts( pFP, pres, plop );         break;
  case FPOp<T>::SQRT:  _SQRT_cuts( pFP, pres, plop );        break;
  case FPOp<T>::FABS:  _FABS_cuts( pFP, pres, plop );        break;
  case FPOp<T>::IPOW:  _IPOW_cuts( pFP, pres, plop, prop );  break;
  case FPOp<T>::ASIN:  _ASIN_cuts( pFP, pres, plop );        break;
  case FPOp<T>::COS:   _COS_cuts( pFP, pres, plop );         break;
  case FPOp<T>::ATAN:  _ATAN_cuts( pFP, pres, plop );        break;
  case FPOp<T>::TAN:   _TAN_cuts( pFP, pres, plop );         break;
  case FPOp<T>::ERF:   _ERF_cuts( pFP, pres, plop );         break;
  case FPOp<T>::FSTEP: _FSTEP_cuts( pFP, pres, plop );       break;

  case FPOp<T>::CNST:
  case FPOp<T>::RANGE:
  case FPOp<T>::VAR:
  case FPOp<T>::MC:
  case FPOp<T>::TM:
  case FPOp<T>::TMMC:  break;

  default:
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::UNDEF );
  }
}

template <typename T> inline void
FPOp<T>::_EQ_cuts
( FPRelax<T>*pFP, const FPVar<T>*Var1, const FPVar<T>*Var2 )
{
  if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().n, Var2->id(), 1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().x, Var2->id(), 1. );
  else if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().n, Var1->id(), 1. );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().x, Var1->id(), 1. );
  else 
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), -1. );
}

template <typename T> inline void
FPOp<T>::_LE_cuts
( FPRelax<T>*pFP, const FPVar<T>*Var1, const FPVar<T>*Var2 )
{
  if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().n, Var2->id(), 1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().x, Var2->id(), 1. );
  else if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::LE, Var2->num().n, Var1->id(), 1. );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::LE, Var2->num().x, Var1->id(), 1. );
  else 
    pFP->_append_cut( this, FPCut<T>::LE, 0., Var1->id(), 1., Var2->id(), -1. );
}

template <typename T> inline void
FPOp<T>::_PLUS_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().n, VarR->id(), 1.,
      Var1->id(), -1. );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().x, VarR->id(), 1.,
      Var1->id(), -1. );
  else if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().n, VarR->id(), 1.,
      Var2->id(), -1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().x, VarR->id(), 1.,
      Var2->id(), -1. );
  else{
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(), -1.,
      Var2->id(), -1. );
  }
}

template <typename T> inline void
FPOp<T>::_NEG_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1 )
{
  if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, -Var1->num().n, VarR->id(), 1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, -Var1->num().x, VarR->id(), 1. );
  else{
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(), 1. );
  }
}

template <typename T> inline void
FPOp<T>::_MINUS_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, -Var2->num().n, VarR->id(), 1.,
      Var1->id(), -1. );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, -Var2->num().x, VarR->id(), 1.,
      Var1->id(), -1. );
  else if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().n, VarR->id(), 1.,
      Var2->id(), 1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().x, VarR->id(), 1.,
      Var2->id(), 1. );
  else{
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(), -1.,
      Var2->id(), 1. );
  }
}

template <typename T> inline void
FPOp<T>::_TIMES_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(),
      -Var2->num().n );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(),
      -Var2->num().x );
  else if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var2->id(),
      -Var1->num().n );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var2->id(),
      -Var1->num().x );
  else{
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::u(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I),
      Var2->id(), -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::l(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I),
      Var2->id(), -Op<T>::l(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::u(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I),
      Var2->id(), -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::l(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I),
      Var2->id(), -Op<T>::l(Var1->num().I) );
  }
}

template <typename T> inline void
FPOp<T>::_DIV_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., VarR->id(),
      -Var2->num().n );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., VarR->id(),
      -Var2->num().x );
  else if( Var1->id().first == FPVar<T>::AUXINT ){
    pFP->_append_cut( this, FPCut<T>::LE, Var1->num().n+Op<T>::u(VarR->num().I)
      *Op<T>::u(Var2->num().I), VarR->id(), Op<T>::u(Var2->num().I), Var2->id(),
      Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE, Var1->num().n+Op<T>::l(VarR->num().I)
      *Op<T>::l(Var2->num().I), VarR->id(), Op<T>::l(Var2->num().I), Var2->id(),
      Op<T>::l(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().n+Op<T>::u(VarR->num().I)
      *Op<T>::l(Var2->num().I), VarR->id(), Op<T>::l(Var2->num().I), Var2->id(),
      Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().n+Op<T>::l(VarR->num().I)
      *Op<T>::u(Var2->num().I), VarR->id(), Op<T>::u(Var2->num().I), Var2->id(),
      Op<T>::l(VarR->num().I) );
  }
  else if( Var1->id().first == FPVar<T>::AUXREAL ){
    pFP->_append_cut( this, FPCut<T>::LE, Var1->num().x+Op<T>::u(VarR->num().I)
      *Op<T>::u(Var2->num().I), VarR->id(), Op<T>::u(Var2->num().I), Var2->id(),
      Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE, Var1->num().x+Op<T>::l(VarR->num().I)
      *Op<T>::l(Var2->num().I), VarR->id(), Op<T>::l(Var2->num().I), Var2->id(),
      Op<T>::l(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().x+Op<T>::u(VarR->num().I)
      *Op<T>::l(Var2->num().I), VarR->id(), Op<T>::l(Var2->num().I), Var2->id(),
      Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE, Var1->num().x+Op<T>::l(VarR->num().I)
      *Op<T>::u(Var2->num().I), VarR->id(), Op<T>::u(Var2->num().I), Var2->id(),
      Op<T>::l(VarR->num().I) );
  }
  else{
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::u(VarR->num().I)*Op<T>::u(Var2->num().I), Var1->id(), 1.,
      VarR->id(), -Op<T>::u(Var2->num().I), Var2->id(), -Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::l(VarR->num().I)*Op<T>::l(Var2->num().I), Var1->id(), 1.,
      VarR->id(), -Op<T>::l(Var2->num().I), Var2->id(), -Op<T>::l(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::u(VarR->num().I)*Op<T>::l(Var2->num().I), Var1->id(), 1.,
      VarR->id(), -Op<T>::l(Var2->num().I), Var2->id(), -Op<T>::u(VarR->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::l(VarR->num().I)*Op<T>::u(Var2->num().I), Var1->id(), 1.,
      VarR->id(), -Op<T>::u(Var2->num().I), Var2->id(), -Op<T>::l(VarR->num().I) );
  }
}

template <typename T> inline void
FPOp<T>::_MINF_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  // Introduce auxiliary variables Var3, as the difference Var2-Var1
  FPVar<T>*Var3 = pFP->_auxiliary_variable( Var2->num().I-Var1->num().I, this );
  _MINUS_cuts( pFP, Var3, Var2, Var1 );
  //pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), -1., Var2->id(), 1.,    
  //  Var3->id(), -1. );

  // Introduce auxiliary variables Var4, as the absolute value of Var3
  FPVar<T>*Var4 = pFP->_auxiliary_variable( Op<T>::fabs(Var3->num().I), this );
  _FABS_cuts( pFP, Var4, Var3 );

  // Introduce auxiliary variables Var5, as the sum Var1+Var2
  FPVar<T>*Var5 = pFP->_auxiliary_variable( Var1->num().I+Var2->num().I, this );
  _PLUS_cuts( pFP, Var5, Var1, Var2 );

  // Cuts
  pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 2., Var4->id(), 1.,
    Var5->id(), -1. );
//   pFP->_append_cut( this, FPCut<T>::LE, 0., VarR->id(), 1., Var1->id(), -1. );
//   pFP->_append_cut( this, FPCut<T>::LE, 0., VarR->id(), 1., Var2->id(), -1. );
}

template <typename T> inline void
FPOp<T>::_MAXF_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  // Introduce auxiliary variables Var3, as the difference Var1-Var2
  FPVar<T>*Var3 = pFP->_auxiliary_variable( Var1->num().I-Var2->num().I, this );
  _MINUS_cuts( pFP, Var3, Var1, Var2 );
  //pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), -1.,    
  //  Var3->id(), -1. );

  // Introduce auxiliary variables Var4, as the absolute value of Var3
  FPVar<T>*Var4 = pFP->_auxiliary_variable( Op<T>::fabs(Var3->num().I), this );
  _FABS_cuts( pFP, Var4, Var3 );

  // Introduce auxiliary variables Var5, as the sum Var1+Var2
  FPVar<T>*Var5 = pFP->_auxiliary_variable( Var1->num().I+Var2->num().I, this );
  _PLUS_cuts( pFP, Var5, Var1, Var2 );

  // Cuts
  pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 2., Var4->id(), -1.,
    Var5->id(), -1. );
//   pFP->_append_cut( this, FPCut<T>::GE, 0., VarR->id(), 1., Var1->id(), -1. );
//   pFP->_append_cut( this, FPCut<T>::GE, 0., VarR->id(), 1., Var2->id(), -1. );
}

template <typename T> inline void
FPOp<T>::_INTER_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( Var1->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().n, VarR->id(), 1. );
  else if( Var1->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var1->num().x, VarR->id(), 1. );
  else 
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var1->id(), -1. );

  if( Var2->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().n, VarR->id(), 1. );
  else if( Var2->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, Var2->num().x, VarR->id(), 1. );
  else 
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1., Var2->id(), -1. );
}

template <typename T> inline void
FPOp<T>::_BILIN_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var1,
  const FPVar<T>*Var2 )
{
  if( mc::isequal( Op<T>::diam(Var1->num().I), 0. )
   || mc::isequal( Op<T>::diam(Var2->num().I), 0. ) )
    return _TIMES_cuts( pFP, VarR, Var1, Var2 );

  const unsigned int Ndiv = pFP->options.BILINEAR_SUBDIV;

  // Append more intermediate variables and corresponding semi-linear cuts
  double*a = new double[Ndiv+3];
  typename FPVar<T>::pt_idVar*id = new typename FPVar<T>::pt_idVar[Ndiv+3];

  switch( pFP->options.BILINEAR_RULE ){
   case FPRelax<T>::Options::UNIVARIATE:{

    FPVar<T>**dVar1 = _subdivide( pFP, Var1, Ndiv, a, id );
    FPVar<T>**dVarR = new FPVar<T>*[Ndiv];
    for( unsigned int idiv=0; idiv<Ndiv; idiv++ ){
      dVarR[idiv] = pFP->_auxiliary_variable( Op<T>::zeroone()
        *Op<T>::diam(Var2->num().I), this );
      a[idiv]  = Op<T>::diam(Var1->num().I)/(double)Ndiv;
      id[idiv] = dVarR[idiv]->id();
      pFP->_append_cut( this, FPCut<T>::GE, -Op<T>::u(Var2->num().I),
        dVarR[idiv]->id(), 1., dVar1[idiv]->id(), -Op<T>::diam(Var2->num().I),
        Var2->id(), -1. );
      pFP->_append_cut( this, FPCut<T>::LE, 0., dVarR[idiv]->id(), 1.,
        dVar1[idiv]->id(), -Op<T>::diam(Var2->num().I) );
      if( !idiv )
        pFP->_append_cut( this, FPCut<T>::LE, -Op<T>::l(Var2->num().I),
          dVarR[idiv]->id(), 1., Var2->id(), -1. );
      else
        pFP->_append_cut( this, FPCut<T>::LE, 0., dVarR[idiv]->id(), 1.,
          dVarR[idiv-1]->id(), -1. );
    }
    a[Ndiv]  = -1.;
    id[Ndiv] = VarR->id();
    a[Ndiv+1]  = Op<T>::l(Var2->num().I);
    id[Ndiv+1] = Var1->id();
    a[Ndiv+2]  = Op<T>::l(Var1->num().I);
    id[Ndiv+2] = Var2->id();
    pFP->_append_cut( this, FPCut<T>::EQ,
      Op<T>::l(Var1->num().I)*Op<T>::l(Var2->num().I), Ndiv+3, id, a );

    delete [] dVar1;
    delete [] dVarR;
    break;
   }

   case FPRelax<T>::Options::SEPARABLE:{

    // Introduce auxiliary variables Var3 and Var4, and subdivide them
#ifdef MC__FPRELAX_SCALE_SEPARABLE
    FPVar<T>*Var3 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I/Op<T>::diam(Var1->num().I)
      +Var2->num().I/Op<T>::diam(Var2->num().I)), this );
    FPVar<T>*Var4 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I/Op<T>::diam(Var1->num().I)
      -Var2->num().I/Op<T>::diam(Var2->num().I)), this );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(),
      1./Op<T>::diam(Var1->num().I), Var2->id(), 1./Op<T>::diam(Var2->num().I),
      Var3->id(), -2. );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(),
      1./Op<T>::diam(Var1->num().I), Var2->id(), -1./Op<T>::diam(Var2->num().I),
      Var4->id(), -2. );
#else
    FPVar<T>*Var3 = pFP->_auxiliary_variable( 0.5*(Var1->num().I+Var2->num().I), this );
    FPVar<T>*Var4 = pFP->_auxiliary_variable( 0.5*(Var1->num().I-Var2->num().I), this );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), 1.,
      Var3->id(), -2. );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), -1.,
      Var4->id(), -2. );    
#endif
    FPVar<T>**dVar3 = _subdivide( pFP, Var3, Ndiv, a, id );
    FPVar<T>**dVar4 = _subdivide( pFP, Var4, Ndiv, a, id );

    // Introduce auxiliary variables Var5 and Var6, as the squares of Var3
    // and Var4, and append corresponding cuts
    FPVar<T>*Var5 = pFP->_auxiliary_variable( Op<T>::sqr(Var3->num().I), this );
    FPVar<T>*Var6 = pFP->_auxiliary_variable( Op<T>::sqr(Var4->num().I), this );
#ifdef MC__FPRELAX_SCALE_SEPARABLE
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var5->id(), 1., Var6->id(), -1.,
      VarR->id(), -1./Op<T>::diam(Var1->num().I)/Op<T>::diam(Var2->num().I) );
#else
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var5->id(), 1., Var6->id(), -1.,
      VarR->id(), -1. );
#endif

    // Append original McCormick cuts
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::u(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I), Var2->id(),
      -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::l(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I), Var2->id(),
      -Op<T>::l(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::u(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I), Var2->id(),
      -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::l(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I), Var2->id(),
      -Op<T>::l(Var1->num().I) );

    // Append outer-approximation cuts for Var5 and Var6
    struct dloc{ static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( mc::sqr(x), 2.*x ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_sqrcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_sqrcv;
    _sandwich_cuts( pFP, Var3->id(), Op<T>::l(Var3->num().I),
      Op<T>::u(Var3->num().I), Var5->id(), Op<T>::l(Var5->num().I),
      Op<T>::u(Var5->num().I), paddcut_sqrcv, dloc::sqr );
    _sandwich_cuts( pFP, Var4->id(), Op<T>::l(Var4->num().I),
      Op<T>::u(Var4->num().I), Var6->id(), Op<T>::l(Var6->num().I),
      Op<T>::u(Var6->num().I), paddcut_sqrcv, dloc::sqr );

    // Append piecewise-affine cuts for Var5 and Var6
    struct loc{ static double sqr
      ( const double x, const double*rusr, const int*iusr )
      { return mc::sqr(x); }
    };
    _semilinear_cut( pFP, Ndiv, Op<T>::l(Var3->num().I),
      Op<T>::u(Var3->num().I), dVar3, Var5->id(), Op<T>::l(Var5->num().I),
      Op<T>::u(Var5->num().I), loc::sqr );
    _semilinear_cut( pFP, Ndiv, Op<T>::l(Var4->num().I),
      Op<T>::u(Var4->num().I), dVar4, Var6->id(), Op<T>::l(Var6->num().I),
      Op<T>::u(Var6->num().I), loc::sqr );

    delete [] dVar3;
    delete [] dVar4;
    break;
   }

   case FPRelax<T>::Options::SEPARSOS:{

    // Introduce auxiliary variables Var3 and Var4, and subdivide them
#ifdef MC__FPRELAX_SCALE_SEPARABLE
    FPVar<T>*Var3 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I/Op<T>::diam(Var1->num().I)
      +Var2->num().I/Op<T>::diam(Var2->num().I)), this );
    FPVar<T>*Var4 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I/Op<T>::diam(Var1->num().I)
      -Var2->num().I/Op<T>::diam(Var2->num().I)), this );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(),
      1./Op<T>::diam(Var1->num().I), Var2->id(), 1./Op<T>::diam(Var2->num().I),
      Var3->id(), -2. );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(),
      1./Op<T>::diam(Var1->num().I), Var2->id(), -1./Op<T>::diam(Var2->num().I),
      Var4->id(), -2. );
#else
    FPVar<T>*Var3 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I+Var2->num().I), this );
    FPVar<T>*Var4 = pFP->_auxiliary_variable(
      0.5*(Var1->num().I-Var2->num().I), this );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), 1.,
      Var3->id(), -2. );
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var1->id(), 1., Var2->id(), -1.,
      Var4->id(), -2. );    
#endif
    FPVar<T>**dVar3 = _subdivide_SOS( pFP, Var3, Ndiv, a, id );
    FPVar<T>**dVar4 = _subdivide_SOS( pFP, Var4, Ndiv, a, id );

    // Introduce auxiliary variables Var5 and Var6, as the squares of Var3
    // and Var4, and append corresponding cuts
    FPVar<T>*Var5 = pFP->_auxiliary_variable( Op<T>::sqr(Var3->num().I), this );
    FPVar<T>*Var6 = pFP->_auxiliary_variable( Op<T>::sqr(Var4->num().I), this );
#ifdef MC__FPRELAX_SCALE_SEPARABLE
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var5->id(), 1., Var6->id(), -1.,
      VarR->id(), -1./Op<T>::diam(Var1->num().I)/Op<T>::diam(Var2->num().I) );
#else
    pFP->_append_cut( this, FPCut<T>::EQ, 0., Var5->id(), 1., Var6->id(), -1.,
      VarR->id(), -1. );
#endif

    // Append original McCormick cuts
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::u(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I), Var2->id(),
      -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE,
      -Op<T>::l(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I), Var2->id(),
      -Op<T>::l(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::u(Var1->num().I)*Op<T>::l(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::l(Var2->num().I), Var2->id(),
      -Op<T>::u(Var1->num().I) );
    pFP->_append_cut( this, FPCut<T>::LE,
      -Op<T>::l(Var1->num().I)*Op<T>::u(Var2->num().I),
      VarR->id(), 1., Var1->id(), -Op<T>::u(Var2->num().I), Var2->id(),
      -Op<T>::l(Var1->num().I) );

    // Append outer-approximation cuts for Var5 and Var6
    struct dloc{ static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( mc::sqr(x), 2.*x ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_sqrcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_sqrcv;
    _sandwich_cuts( pFP, Var3->id(), Op<T>::l(Var3->num().I),
      Op<T>::u(Var3->num().I), Var5->id(), Op<T>::l(Var5->num().I),
      Op<T>::u(Var5->num().I), paddcut_sqrcv, dloc::sqr );
    _sandwich_cuts( pFP, Var4->id(), Op<T>::l(Var4->num().I),
      Op<T>::u(Var4->num().I), Var6->id(), Op<T>::l(Var6->num().I),
      Op<T>::u(Var6->num().I), paddcut_sqrcv, dloc::sqr );

    // Append piecewise-affine SOS cuts for Var5 and Var6
    struct loc{ static double sqr
      ( const double x, const double*rusr, const int*iusr )
      { return mc::sqr(x); }
    };
    _semilinear_SOS( pFP, Ndiv, Op<T>::l(Var3->num().I),
      Op<T>::u(Var3->num().I), dVar3, Var5->id(), Op<T>::l(Var5->num().I),
      Op<T>::u(Var5->num().I), loc::sqr );
    _semilinear_SOS( pFP, Ndiv, Op<T>::l(Var4->num().I),
      Op<T>::u(Var4->num().I), dVar4, Var6->id(), Op<T>::l(Var6->num().I),
      Op<T>::u(Var6->num().I), loc::sqr );

    delete [] dVar3;
    delete [] dVar4;
    break;
   }
  }

  delete [] a;
  delete [] id;
}

template <typename T> inline FPVar<T>**
FPOp<T>::_subdivide
( FPRelax<T>*pFP, const FPVar<T>*Var, const unsigned int Ndiv, double*a,
  typename FPVar<T>::pt_idVar*id )
{
  const double divVar = Op<T>::diam(Var->num().I)/(double)Ndiv;
  FPVar<T>**dVar = new FPVar<T>*[Ndiv];
  for( unsigned int idiv=0; idiv<Ndiv; idiv++ ){
    dVar[idiv] = pFP->_auxiliary_variable( Op<T>::zeroone(), this );
    a[idiv]  = divVar;
    id[idiv] = dVar[idiv]->id();
  }
  a[Ndiv]  = -1.;
  id[Ndiv] = Var->id();
  pFP->_append_cut( this, FPCut<T>::EQ, -Op<T>::l(Var->num().I), Ndiv+1,
    id, a );

  FPVar<T>**bVar = new FPVar<T>*[Ndiv-1];
  for( unsigned int idiv=0; idiv<Ndiv-1; idiv++ ){
    bVar[idiv] = pFP->_auxiliary_variable( this );
    pFP->_append_cut( this, FPCut<T>::GE, 0., dVar[idiv]->id(), 1.,
                 bVar[idiv]->id(), -1. );
    pFP->_append_cut( this, FPCut<T>::LE, 0., dVar[idiv+1]->id(), 1.,
                 bVar[idiv]->id(), -1. );
  }

  delete[] bVar;
  return dVar;
}

template <typename T> inline void
FPOp<T>::_semilinear_cut
( FPRelax<T>*pFP, const unsigned int Ndiv, const double XL, const double XU,
  FPVar<T>**dVarX, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, p_Univ f, const double*rpar, const int*ipar )
{
  double a[Ndiv+1];
  typename FPVar<T>::pt_idVar id[Ndiv+1];

  const double fL = f(XL,rpar,ipar);
  double Xi1 = XL, Fi1 = fL, Xi2, Fi2;
  for( unsigned int idiv=0; idiv<Ndiv; idiv++ ){
    Xi2 = ( idiv<Ndiv-1? Xi1+(XU-XL)/(double)Ndiv: XU );
    Fi2 = f(Xi2,rpar,ipar);
    a[idiv]  = Fi1-Fi2;
    id[idiv] = dVarX[idiv]->id();
    Xi1 = Xi2;
    Fi1 = Fi2;
  }
  a[Ndiv]  = 1.;
  id[Ndiv] = iY;
  pFP->_append_cut( this, FPCut<T>::LE, fL, Ndiv+1, id, a );
}

template <typename T> inline FPVar<T>**
FPOp<T>::_subdivide_SOS
( FPRelax<T>*pFP, const FPVar<T>*Var, const unsigned int Ndiv, double*a,
  typename FPVar<T>::pt_idVar*id )
{
  const double divVar = Op<T>::diam(Var->num().I)/(double)Ndiv;
  FPVar<T>**dVar = new FPVar<T>*[Ndiv+1];
  for( unsigned int idiv=0; idiv<Ndiv+1; idiv++ ){
    dVar[idiv] = pFP->_auxiliary_variable( Op<T>::zeroone(), this );
    a[idiv]  = 1.;
    id[idiv] = dVar[idiv]->id();
  }
  pFP->_append_cut( this, FPCut<T>::SOS2, 1., Ndiv+1, id, a );

  double Xi = Op<T>::l(Var->num().I);
  a[0]  = -Xi;
  id[0] = dVar[0]->id();
  for( unsigned int idiv=0; idiv<Ndiv; idiv++ ){
    Xi = ( idiv<Ndiv-1? Xi+divVar: Op<T>::u(Var->num().I) );
    a[idiv+1]  = -Xi;
    id[idiv+1] = dVar[idiv+1]->id();
  }
  a[Ndiv+1]  = 1.;
  id[Ndiv+1] = Var->id();
  pFP->_append_cut( this, FPCut<T>::EQ, 0., Ndiv+2, id, a );

  return dVar;
}

template <typename T> inline void
FPOp<T>::_semilinear_SOS
( FPRelax<T>*pFP, const unsigned int Ndiv, const double XL,
  const double XU, FPVar<T>**dVarX, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, p_Univ f, const double*rpar, const int*ipar )
{
  double*a = new double[Ndiv+2];
  typename FPVar<T>::pt_idVar*id = new typename FPVar<T>::pt_idVar[Ndiv+2];

  double Xi = XL, Fi = f(XL,rpar,ipar);
  a[0]  = -Fi;
  id[0] = dVarX[0]->id();
  for( unsigned int idiv=0; idiv<Ndiv; idiv++ ){
    Xi = ( idiv<Ndiv-1? Xi+(XU-XL)/(double)Ndiv: XU );
    Fi = f(Xi,rpar,ipar);
    a[idiv+1]  = -Fi;
    id[idiv+1] = dVarX[idiv+1]->id();
  }
  a[Ndiv+1]  = 1.;
  id[Ndiv+1] = iY;
  pFP->_append_cut( this, FPCut<T>::LE, 0., Ndiv+2, id, a );

  delete [] a;
  delete [] id;
}

template <typename T> inline std::pair< double, double >
FPOp<T>::_distmax
( const FPRelax<T>*pFP, p_dUniv f, const double xL, const double xU,
  const double*rpar, const int*ipar ) const
{
  std::pair<double,double> fL = f(xL,rpar,ipar), fU = f(xU,rpar,ipar);
  double Ddf = fU.second-fL.second,
         Adf = std::max(std::fabs(fL.second),std::fabs(fU.second));

  double xmid = 0.5 * ( xL + xU );
  double xmax = ( std::fabs(Ddf) - Adf*pFP->options.FRACTIONAL_RTOL
    > pFP->options.FRACTIONAL_ATOL?
    ( fU.second * xU - fL.second * xL - fU.first + fL.first ) / Ddf: xmid );
//   double xmax = ( std::fabs(Ddf) < Adf*pFP->options.FRACTIONAL_RTOL
//     + pFP->options.FRACTIONAL_ATOL? xmid:
//     ( fU.second * xU - fL.second * xL - fU.first + fL.first ) / Ddf );
  std::pair<double,double> fmax = f(xmax,rpar,ipar);
  double dmax = std::fabs( fmax.first - fL.second * ( xmax - xL ) - fL.first );

  switch( pFP->options.SANDWICH_RULE ){
  case FPRelax<T>::Options::BISECT:
    return std::make_pair( xmid, dmax );
  case FPRelax<T>::Options::MAXERR: default:
    return std::make_pair( xmax, dmax );
  }
}

template <typename T> inline void
FPOp<T>::_sandwich_cuts
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, p_Cutpar Cut, p_dUniv f, const double*rpar, const int*ipar )
{
  // Convex cuts @xL,xU
  (this->*Cut)( pFP, XL, iX, XL, XU, iY, YL, YU, rpar, ipar );
  (this->*Cut)( pFP, XU, iX, XL, XU, iY, YL, YU, rpar, ipar );
  unsigned int NBCUTS = 2;
  typename FPOp<T>::t_Outer OA;
  std::pair<double,double> worst = _distmax( pFP, f, XL, XU, rpar, ipar );
  OA.push( FPOuter( XL, XU, worst.first, worst.second ) );
  const double dtol = pFP->options.SANDWICH_ATOL
    + pFP->options.SANDWICH_RTOL*std::max(std::fabs(YL),std::fabs(YU));
#ifdef MC__FPRELAX_DEBUG_SANDWICH
  std::cerr << "DTOL: " << dtol << std::endl;
  std::cerr << OA.top() << std::endl;
#endif

  // OA @xM
  while( OA.top().gap() > dtol && NBCUTS < pFP->options.SANDWICH_MAXCUT ){
    // x - y/exp(xref) <= xref - 1, @xref=xmax
    (this->*Cut)( pFP, OA.top().xM(), iX, OA.top().xL(), OA.top().xU(),
      iY, YL, YU, rpar, ipar );
    NBCUTS++;
    worst = _distmax( pFP, f, OA.top().xL(), OA.top().xM(), rpar, ipar );
    OA.push( FPOuter( OA.top().xL(), OA.top().xM(), worst.first, worst.second ) );
#ifdef MC__FPOP_DEBUG_SANDWICH
    std::cerr << "  Pushed: "  << FPOuter( OA.top().xL(), OA.top().xM(),
      worst.first, worst.second ) << std::endl;
#endif
    worst = _distmax( pFP, f, OA.top().xM(), OA.top().xU(), rpar, ipar );
    OA.push( FPOuter( OA.top().xM(), OA.top().xU(), worst.first, worst.second ) );
#ifdef MC__FPOP_DEBUG_SANDWICH
    std::cerr << "  Pushed: " << FPOuter( OA.top().xM(), OA.top().xU(),
      worst.first, worst.second ) << std::endl;
#endif
    OA.pop();
#ifdef MC__FPOP_DEBUG_SANDWICH
    std::cerr << OA.top() << std::endl;
#endif
  }
#ifdef MC__FPOP_DEBUG_SANDWICH
  std::cerr << "NBCUTS: " << NBCUTS << std::endl;
#endif
}

template <typename T> inline double
FPOp<T>::_newton
( const FPRelax<T>*pFP, const double x0, const double xL,
  const double xU, p_dUniv f, const double*rusr, const int*iusr ) const 
{
  double xk = std::max(xL,std::min(xU,x0)), dk;
  std::pair<double,double> fk = f(xk,rusr,iusr);
  
  for( unsigned int it=0; it<pFP->options.NEWTON_MAXIT; it++ ){
    if( std::fabs(fk.first) < pFP->options.NEWTON_TOL ) return xk;
    if( fk.second == 0 )
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::NEWTON );
    dk = fk.first/fk.second;
    if( mc::isequal(xk,xL) && dk>0 ) return xk;
    if( mc::isequal(xk,xU) && dk<0 ) return xk;
    xk = std::max(xL,std::min(xU,xk-dk));
    fk = f(xk,rusr,iusr);
  }

  throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::NEWTON );
}

template <typename T> inline double
FPOp<T>::_secant
( const FPRelax<T>*pFP, const double x0, const double x1, const double xL,
  const double xU, p_Univ f, const double*rusr, const int*iusr ) const
{
  double xkm = std::max(xL,std::min(xU,x0));
  double fkm = f(xkm,rusr,iusr);
  double xk = std::max(xL,std::min(xU,x1));
  
  for( unsigned int it=0; it<pFP->options.NEWTON_MAXIT; it++ ){
    double fk = f(xk,rusr,iusr);
    if( std::fabs(fk) < pFP->options.NEWTON_TOL ) return xk;
    double Bk = (fk-fkm)/(xk-xkm);
    if( Bk == 0 )
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::NEWTON );
    if( isequal(xk,xL) && fk/Bk>0 ) return xk;
    if( isequal(xk,xU) && fk/Bk<0 ) return xk;
    xkm = xk;
    fkm = fk;
    xk = std::max(xL,std::min(xU,xk-fk/Bk));
  }

  throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::NEWTON );
}

template <typename T> inline void
FPOp<T>::_EXP_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::exp(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::exp(Var->num().x),
      VarR->id(), 1. );
  else{
    _addcut_expcc( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I) );
    struct loc{ static std::pair<double,double> exp
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::exp(x), std::exp(x) ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_expcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_expcv;
    _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), paddcut_expcv, loc::exp );
  }
}

template <typename T> inline void
FPOp<T>::_LOG_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::log(Var->num().n),
    VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::log(Var->num().x),
    VarR->id(), 1. );
  else{
    _addcut_expcc( pFP, Var->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I) );
    struct loc{ static std::pair<double,double> exp
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::exp(x), std::exp(x) ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_expcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_expcv;
    _sandwich_cuts( pFP, VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), paddcut_expcv, loc::exp );
  }
}

template <typename T> inline void
FPOp<T>::_addcut_expcv
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x - y*exp(-xref) <= xref - 1
  pFP->_append_cut( this, FPCut<T>::LE, Xref-1., iY, -std::exp(-Xref),
    iX, 1. );
}

template <typename T> inline void
FPOp<T>::_addcut_expcc
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, const double*rpar, const int*ipar )
{
  // dy*x - dx*y >= dy*xL - dx*exp(xL)
  double dX = XU - XL, dY = YU - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*std::exp(XL), iY, -dX,
    iX, dY );
}

template <typename T> inline void
FPOp<T>::_SQR_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, mc::sqr(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, mc::sqr(Var->num().x),
      VarR->id(), 1. );
  else{
    _addcut_sqrcc( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I) );
    struct dloc{ static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( mc::sqr(x), 2.*x ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_sqrcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_sqrcv;
    _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), paddcut_sqrcv, dloc::sqr );
  }
}

template <typename T> inline void
FPOp<T>::_SQRT_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::sqrt(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::sqrt(Var->num().x),
      VarR->id(), 1. );
  else{
    _addcut_sqrcc( pFP, VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I) );
    struct dloc{ static std::pair<double,double> sqr
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( mc::sqr(x), 2.*x ); }
    };
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    void (FPOp<T>::*paddcut_sqrcv)( FPRelax<T>*, const double, const pt_idVar, 
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_sqrcv;
    _sandwich_cuts( pFP, VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), paddcut_sqrcv, dloc::sqr );
  }
}

template <typename T> inline void
FPOp<T>::_addcut_sqrcv
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // 2*xref*x - y <= xref^2
  pFP->_append_cut( this, FPCut<T>::LE, mc::sqr(Xref), iY, -1., iX, 2.*Xref );
}

template <typename T> inline void
FPOp<T>::_addcut_sqrcc
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, const double*rpar, const int*ipar )
{
  // (xL+xU)*x - y >= xL*xU
  pFP->_append_cut( this, FPCut<T>::GE, XL*XU, iY, -1., iX, XL+XU );
}

template <typename T> inline void
FPOp<T>::_FABS_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::fabs(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::fabs(Var->num().x),
      VarR->id(), 1. );
  else{
    double XL = Op<T>::l(Var->num().I), dX = Op<T>::u(Var->num().I)-XL,
      YL = std::fabs(Op<T>::l(Var->num().I)),
      dY = std::fabs(Op<T>::u(Var->num().I))-YL;
    pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*YL, VarR->id(), -dX,
      Var->id(), dY );
    pFP->_append_cut( this, FPCut<T>::GE, 0., VarR->id(), 1., Var->id(), -1. );
    pFP->_append_cut( this, FPCut<T>::GE, 0., VarR->id(), 1., Var->id(), 1. );
  }
}

template <typename T> inline void
FPOp<T>::_FSTEP_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, mc::fstep(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, mc::fstep(Var->num().x),
      VarR->id(), 1. );
  else if( Op<T>::l(Var->num().I) >= 0. )
    pFP->_append_cut( this, FPCut<T>::EQ, 1., VarR->id(), 1. );    
  else if( Op<T>::u(Var->num().I) < 0. )
    pFP->_append_cut( this, FPCut<T>::EQ, 0., VarR->id(), 1. );    
  else{
    pFP->_append_cut( this, FPCut<T>::LE, 1., VarR->id(), 1. );    
    pFP->_append_cut( this, FPCut<T>::GE, 0., VarR->id(), 1. );    
    pFP->_append_cut( this, FPCut<T>::GE, 0., Var->id(), -1.,
      VarR->id(), Op<T>::u(Var->num().I) );
    pFP->_append_cut( this, FPCut<T>::GE, Op<T>::l(Var->num().I),
      Var->id(), 1., VarR->id(), Op<T>::l(Var->num().I) );
  }
}

template <typename T> inline void
FPOp<T>::_IPOW_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var,
  const FPVar<T>*Exp )
{
  const int iExp = Exp->num().n;
  if( iExp == 0 || iExp == 1 || iExp == 2 )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::pow(Var->num().n,iExp),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::pow(Var->num().x,iExp),
      VarR->id(), 1. );
  else{
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    struct loc{ static std::pair<double,double> pow
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::pow(x,*iusr), *iusr*std::pow(x,*iusr-1) ); }
    };

    if( iExp > 0 && !(iExp%2) ){
      // Append new linear cuts for power term w/ even exponent
      _addcut_powcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
      void (FPOp<T>::*paddcut_powcv_lin)( FPRelax<T>*, const double, const pt_idVar,
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_powcv_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_powcv_lin, loc::pow, 0, &iExp );
    }

    else if( iExp > 0 && pFP->options.NEWTON_USE ){
      // Append new linear cuts for power term w/ even exponent
      // -- Convex Portion
      if( Op<T>::l(Var->num().I) >= 0. ){
        _addcut_powcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
        void (FPOp<T>::*paddcut_powcv_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_powcv_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_powcv_lin, loc::pow, 0, &iExp );
      }

      // -- Concave Portion
      else if( Op<T>::u(Var->num().I) <= 0. ){
        void (FPOp<T>::*paddcut_powcc_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_powcc_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_powcc_lin, loc::pow, 0, &iExp );
        _addcut_powcv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
      }

      // -- Nonconvex/Nonconcave Portion
      else{
        struct fct{ static std::pair<double,double> powoddfunc
          ( const double x, const double*rusr, const int*iusr )
          { return std::make_pair(
              ((*iusr-1)*x-(*iusr)*(*rusr))*std::pow(x,*iusr-1) + std::pow(*rusr,*iusr),
              (*iusr)*(*iusr-1)*(x-(*rusr))*std::pow(x,*iusr-2) ); }
        };
        double xJcc = Op<T>::u(Var->num().I);
        xJcc = _newton( pFP, Op<T>::l(Var->num().I), Op<T>::l(Var->num().I), 0.,
          fct::powoddfunc, &xJcc, &iExp );
        if( mc::isequal(xJcc,Op<T>::l(Var->num().I)) ){
          _addcut_powcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
            Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
        }
        else{
          void (FPOp<T>::*paddcut_powcc_lin)( FPRelax<T>*, const double, const pt_idVar, 
            const double, const double, const pt_idVar, const double, const double,
            const double*, const int* ) = &FPOp<T>::_addcut_powcc_lin;
          _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcc,
            VarR->id(), Op<T>::l(VarR->num().I),
            Op<T>::u(VarR->num().I), paddcut_powcc_lin, loc::pow, &xJcc, &iExp );
        }
        double xJcv = Op<T>::l(Var->num().I);
        xJcv = _newton( pFP, Op<T>::u(Var->num().I), 0., Op<T>::u(Var->num().I),
          fct::powoddfunc, &xJcv, &iExp );
        if( mc::isequal(xJcv,Op<T>::u(Var->num().I)) ){
          _addcut_powcv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
        }
        else{
          void (FPOp<T>::*paddcut_powcv_lin)( FPRelax<T>*, const double, const pt_idVar, 
            const double, const double, const pt_idVar, const double, const double,
            const double*, const int* ) = &FPOp<T>::_addcut_powcv_lin;
          _sandwich_cuts( pFP, Var->id(), xJcv, Op<T>::u(Var->num().I),
            VarR->id(), Op<T>::l(VarR->num().I),
            Op<T>::u(VarR->num().I), paddcut_powcv_lin, loc::pow, &xJcv, &iExp );
        }
      }
    }

    else if( iExp < 0 ){
      // Append new linear cuts for power term w/ negative exponent
      // -- Convex Case
      if( !(iExp%2) || Op<T>::l(Var->num().I) > 0. ){
        _addcut_powcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
        void (FPOp<T>::*paddcut_powcv_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_powcv_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_powcv_lin, loc::pow, 0, &iExp );
      }
      // -- Concave Case
      else{
        void (FPOp<T>::*paddcut_powcc_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_powcc_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_powcc_lin, loc::pow, 0, &iExp );
        _addcut_powcv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, &iExp );
      }
    }

    else throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }
}

template <typename T> inline void
FPOp<T>::_addcut_powcv_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*iExp )
{
  // n*xref^(n-1)*x - y <= (n-1)*xref^n
  const double a = std::pow(Xref,(*iExp)-1);
  pFP->_append_cut( this, FPCut<T>::LE, ((*iExp)-1)*a*Xref, iY, -1., iX,
    (*iExp)*a );
}

template <typename T> inline void
FPOp<T>::_addcut_powcc_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*iExp )
{
  // n*xref^(n-1)*x - y >= (n-1)*xref^n
  const double a = std::pow(Xref,(*iExp)-1);
  pFP->_append_cut( this, FPCut<T>::GE, ((*iExp)-1)*a*Xref, iY, -1., iX,
    (*iExp)*a );
}

template <typename T> inline void
FPOp<T>::_addcut_powcv_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*iExp )
{
  // dy*x - dx*y <= dy*xL - dx*yL
  double dX = XU - XL, YL = std::pow(XL,*iExp), dY = std::pow(XU,*iExp) - YL;
  pFP->_append_cut( this, FPCut<T>::LE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_powcc_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*iExp )
{
  // dy*x - dx*y >= dy*xL - dx*yL
  double dX = XU - XL, YL = std::pow(XL,*iExp), dY = std::pow(XU,*iExp) - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_ASIN_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::asin(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::asin(Var->num().x),
      VarR->id(), 1. );
  else{
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    struct loc{ static std::pair<double,double> asin
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::asin(x), 1./std::sqrt(1-mc::sqr(x)) ); }
    };

    // Append new linear cuts for asin univariate term
    // -- Convex Portion
    if(Op<T>::l(Var->num().I) >= 0. ){
      _addcut_asin1cc( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), 0, 0 );
      void (FPOp<T>::*paddcut_asin1cv)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_asin1cv;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_asin1cv, loc::asin, 0, 0 );
    }

    // -- Concave Portion
    else if(Op<T>::u(Var->num().I) <= 0. ){
      void (FPOp<T>::*paddcut_asin2cc)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_asin2cc;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_asin2cc, loc::asin, 0, 0 );
      _addcut_asin2cv( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), 0, 0 );
    }

    // -- Nonconvex/Nonconcave Portion (convex/cocnave envelopes)
    else if( pFP->options.NEWTON_USE ){
      struct fct{ static double asinfunc
        ( const double x, const double*rusr, const int*iusr )
        { return x-(*rusr)-std::sqrt(1.-x*x)*(std::asin(x)-std::asin(*rusr)); }
      };
      double xJcc = Op<T>::u(Var->num().I);
      xJcc = _secant( pFP, 0., Op<T>::l(Var->num().I), Op<T>::l(Var->num().I),
        0., fct::asinfunc, &xJcc, 0 );
      if( mc::isequal(xJcc,Op<T>::l(Var->num().I)) ){
        _addcut_asin1cc( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_asin2cc)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_asin2cc;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcc,
          VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_asin2cc, loc::asin, &xJcc, 0 );
      }
      double xJcv = Op<T>::l(Var->num().I);
      xJcv = _secant( pFP, 0., Op<T>::u(Var->num().I), 0., Op<T>::u(Var->num().I),
        fct::asinfunc, &xJcv, 0 );
      if( mc::isequal(xJcv,Op<T>::u(Var->num().I)) ){
        _addcut_asin2cv( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_asin1cv)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_asin1cv;
        _sandwich_cuts( pFP, Var->id(), xJcv, Op<T>::u(Var->num().I),
          VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), paddcut_asin1cv, loc::asin, &xJcv, 0 );
      }
    }

    // -- Nonconvex/Nonconcave Portion (cheap convex/concave relaxations)
    else{
      _addcut_asin3cc( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), 0, 0 );
      _addcut_asin3cv( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), 0, 0 ); 
    }
  }
}

template <typename T> inline void
FPOp<T>::_addcut_asin1cv
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x - y*sqrt(1-xref^2) <= xref - asin(xref)*sqrt(1-xref^2)
  const double a = -std::sqrt(1-mc::sqr(Xref));
  pFP->_append_cut( this, FPCut<T>::LE, Xref+std::asin(Xref)*a, iY, a, iX, 1. );
}

template <typename T> inline void
FPOp<T>::_addcut_asin1cc
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, const double*rpar, const int*ipar )
{
  // dy*x - dx*y >= dy*xL - dx*asin(xL)
  double dX = XU - XL, dY = YU - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*std::asin(XL), iY, -dX,
    iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_asin2cc
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x - y*sqrt(1-xref^2) >= xref - asin(xref)*sqrt(1-xref^2)
  const double a = -std::sqrt(1-mc::sqr(Xref));
  pFP->_append_cut( this, FPCut<T>::GE, Xref+std::asin(Xref)*a, iY, a, iX, 1. );
}

template <typename T> inline void
FPOp<T>::_addcut_asin2cv
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double YL,
  const double YU, const double*rpar, const int*ipar )
{
  // dy*x - dx*y <= dy*xL - dx*asin(xL)
  double dX = XU - XL, dY = YU - YL;
  pFP->_append_cut( this, FPCut<T>::LE, dY*XL-dX*std::asin(XL), iY, -dX,
    iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_asin3cc
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x - y >= xU - asin(xU)
  pFP->_append_cut( this, FPCut<T>::GE, XU-std::asin(XU), iY, -1., iX, 1. );
}

template <typename T> inline void
FPOp<T>::_addcut_asin3cv
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x - y <= xL - asin(xL)
  pFP->_append_cut( this, FPCut<T>::LE, XL-std::asin(XL), iY, -1., iX, 1. );
}

template <typename T> inline void
FPOp<T>::_COS_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT ){
    pFP->_append_cut( this, FPCut<T>::EQ, std::sin(Var->num().n),
      VarR->id(), 1. );
    return;
  }
  if( Var->id().first == FPVar<T>::AUXREAL ){
    pFP->_append_cut( this, FPCut<T>::EQ, std::sin(Var->num().x),
      VarR->id(), 1. );
    return;
  }

  typedef typename FPVar<T>::pt_idVar pt_idVar;
  struct loc{ static std::pair<double,double> cos
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair( std::cos(x), -std::sin(x) ); }
  };
  struct fct{ static std::pair<double,double> cosfunc
    ( const double x, const double*rusr, const int*iusr )
    { return std::make_pair(
        (x-*rusr)*std::sin(x)+std::cos(x)-std::cos(*rusr),
        (x-*rusr)*std::cos(x) ); }
  };
  
  // Convex relaxation
  int kL = std::ceil( -0.5*(1.+Op<T>::l(Var->num().I)/PI) );
  double dxL = 2.*PI*kL; 
  double xL1 = Op<T>::l(Var->num().I)+dxL, xU1 = Op<T>::u(Var->num().I)+dxL;
  assert( xL1 >= -PI && xL1 <= PI );
 
  if( xU1 >= PI ){
    const int kU = std::ceil( -0.5*(1.+Op<T>::u(Var->num().I)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(Var->num().I)+dxU;
    assert( xU2 >= -PI && xU2 <= PI );

    double xJcv1 = Op<T>::l(Var->num().I);
    if( xL1 <= PI/2. ) xJcv1 = _newton( pFP, PI-dxL,
      Op<T>::l(Var->num().I), PI-dxL, fct::cosfunc, &xJcv1, 0 );
    double xJcv2 = Op<T>::u(Var->num().I);
    if( xU2 >= -PI/2. ) xJcv2 = _newton( pFP, -PI-dxU,
      -PI-dxU, Op<T>::u(Var->num().I), fct::cosfunc, &xJcv2, 0 ) ;

    void (FPOp<T>::*paddcut_coscv_lin)( FPRelax<T>*, const double, const pt_idVar,
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_coscv_lin;
    _sandwich_cuts( pFP, Var->id(), xJcv1, PI-dxL,
      VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
      paddcut_coscv_lin, loc::cos, 0, 0 );
    _sandwich_cuts( pFP, Var->id(), -PI-dxU, xJcv2, 
      VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
      paddcut_coscv_lin, loc::cos, 0, 0 );        
  }
  
  else if( xL1 >= PI/2. ){
    void (FPOp<T>::*paddcut_coscv_lin)( FPRelax<T>*, const double, const pt_idVar,
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_coscv_lin;
    _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), paddcut_coscv_lin, loc::cos, 0, 0 );
  }
  
  else if( xL1 >= -PI/2. && xU1 <= PI/2. ){
    _addcut_coscv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
  }
  
  else if( xL1 >= -PI/2. ){
    double xJcv1 = Op<T>::l(Var->num().I);
    xJcv1 = _newton( pFP, Op<T>::u(Var->num().I),
      Op<T>::l(Var->num().I), Op<T>::u(Var->num().I), fct::cosfunc, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::u(Var->num().I) ) ){
      _addcut_coscv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }
    else{
      void (FPOp<T>::*paddcut_coscv_lin)( FPRelax<T>*, const double, const pt_idVar,
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_coscv_lin;
      _sandwich_cuts( pFP, Var->id(), xJcv1, Op<T>::u(Var->num().I),
        VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
        paddcut_coscv_lin, loc::cos, 0, 0 );
    }
  }
   
  else{
    double xJcv1 = Op<T>::u(Var->num().I);
    xJcv1 = _newton( pFP, Op<T>::l(Var->num().I),
      Op<T>::l(Var->num().I), Op<T>::u(Var->num().I), fct::cosfunc, &xJcv1, 0 );
    if( mc::isequal( xJcv1, Op<T>::l(Var->num().I) ) ){
      _addcut_coscv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }
    else{
      void (FPOp<T>::*paddcut_coscv_lin)( FPRelax<T>*, const double, const pt_idVar,
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_coscv_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcv1,
        VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
        paddcut_coscv_lin, loc::cos, 0, 0 );
    }
  }
  
  // Concave relaxation
  kL = std::ceil( -0.5*(2.+Op<T>::l(Var->num().I)/PI) );
  dxL = 2.*PI*kL;
  xL1 = Op<T>::l(Var->num().I)+dxL, xU1 = Op<T>::u(Var->num().I)+dxL;
  assert( xL1 >= -2.*PI && xL1 <= 0. );

  if( xU1 >= 0. ){
    const int kU = std::ceil( -0.5*(2.+Op<T>::u(Var->num().I)/PI) );
    const double dxU = 2.*PI*kU; 
    const double xU2 = Op<T>::u(Var->num().I)+dxU;
    assert( xU2 >= -2.*PI && xU2 <= 0. );

    double xJcc1 = Op<T>::l(Var->num().I);
    if( xL1 <= -PI/2. ) xJcc1 = _newton( pFP, -dxL,
      Op<T>::l(Var->num().I), -dxL, fct::cosfunc, &xJcc1, 0 );
    double xJcc2 = Op<T>::u(Var->num().I);
    if( xU2 >= -3.*PI/2. ) xJcc2 = _newton( pFP, -2.*PI-dxU,
      -2.*PI-dxU, Op<T>::u(Var->num().I), fct::cosfunc, &xJcc2, 0 );

    void (FPOp<T>::*paddcut_coscc_lin)( FPRelax<T>*, const double, const pt_idVar,
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_coscc_lin;
    _sandwich_cuts( pFP, Var->id(), xJcc1, -dxL,
      VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
      paddcut_coscc_lin, loc::cos, 0, 0 );
    _sandwich_cuts( pFP, Var->id(), -2.*PI-dxU, xJcc2, 
      VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
      paddcut_coscc_lin, loc::cos, 0, 0 );   
  }

  else if( xL1 >= -PI/2. ){
    void (FPOp<T>::*paddcut_coscc_lin)( FPRelax<T>*, const double, const pt_idVar,
      const double, const double, const pt_idVar, const double, const double,
      const double*, const int* ) = &FPOp<T>::_addcut_coscc_lin;
    _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
      Op<T>::u(VarR->num().I), paddcut_coscc_lin, loc::cos, 0, 0 );
  }
  
  else if( xL1 >= -3.*PI/2. && xU1 <= -PI/2. ){
    _addcut_coscc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
      Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
  }
  
  else if( xL1 >= -3.*PI/2. ){
    double xJcc1 = Op<T>::l(Var->num().I);
    xJcc1 = _newton( pFP, Op<T>::u(Var->num().I),
      Op<T>::l(Var->num().I), Op<T>::u(Var->num().I), fct::cosfunc, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::u(Var->num().I) ) ){
      _addcut_coscc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }
    else{
      void (FPOp<T>::*paddcut_coscc_lin)( FPRelax<T>*, const double, const pt_idVar,
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_coscc_lin;
      _sandwich_cuts( pFP, Var->id(), xJcc1, Op<T>::u(Var->num().I),
        VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
        paddcut_coscc_lin, loc::cos, 0, 0 );
    }
  }
   
  else{
    double xJcc1 = Op<T>::u(Var->num().I);
    xJcc1 = _newton( pFP, Op<T>::l(Var->num().I),
      Op<T>::l(Var->num().I), Op<T>::u(Var->num().I), fct::cosfunc, &xJcc1, 0 );
    if( mc::isequal( xJcc1, Op<T>::l(Var->num().I) ) ){
      _addcut_coscc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }
    else{
      void (FPOp<T>::*paddcut_coscc_lin)( FPRelax<T>*, const double, const pt_idVar,
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_coscc_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcc1,
        VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
        paddcut_coscc_lin, loc::cos, 0, 0 );
    }
  }
}
  
template <typename T> inline void
FPOp<T>::_addcut_coscv_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x*sin(xref) + y >= cos(xref) + xref*sin(xref)
  pFP->_append_cut( this, FPCut<T>::GE, std::cos(Xref)+Xref*std::sin(Xref),
    iY, 1., iX, std::sin(Xref) );
}

template <typename T> inline void
FPOp<T>::_addcut_coscc_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // x*sin(xref) + y <= cos(xref) + xref*sin(xref)
  pFP->_append_cut( this, FPCut<T>::LE, std::cos(Xref)+Xref*std::sin(Xref),
    iY, 1., iX, std::sin(Xref) );
}

template <typename T> inline void
FPOp<T>::_addcut_coscv_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y <= dy*xL - dx*yL
  double dX = XU - XL, YL = std::cos(XL), dY = std::cos(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::LE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_coscc_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y >= dy*xL - dx*yL
  double dX = XU - XL, YL = std::cos(XL), dY = std::cos(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_ATAN_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::atan(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::atan(Var->num().x),
      VarR->id(), 1. );
  else{
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    struct loc{ static std::pair<double,double> atan
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::atan(x), 1./(1.+x*x) ); }
    };

    // Append new linear cuts for univariate atan term
    // -- Convex Portion
    if( Op<T>::u(Var->num().I) <= 0. ){
      _addcut_atancc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      void (FPOp<T>::*paddcut_atancv_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_atancv_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_atancv_lin, loc::atan, 0, 0 );
    }

    // -- Concave Portion
    else if( Op<T>::l(Var->num().I) >= 0. ){
      void (FPOp<T>::*paddcut_atancc_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_atancc_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_atancc_lin, loc::atan, 0, 0 );
      _addcut_atancv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }

    // -- Concavo-convex Portion
    else{
      struct fct{ static std::pair<double,double> atanfunc
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair(
            x-(*rusr)+(1+x*x)*(std::atan(*rusr)-std::atan(x)),
            2*x*(std::atan(*rusr)-std::atan(x)) ); }
      };

      double xJcc = Op<T>::l(Var->num().I);
      xJcc = _newton( pFP, Op<T>::u(Var->num().I), 0., Op<T>::u(Var->num().I),
        fct::atanfunc, &xJcc, 0 );
      if( mc::isequal( xJcc, Op<T>::u(Var->num().I) ) ){
        _addcut_atancc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_atancc_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_atancc_lin;
        _sandwich_cuts( pFP, Var->id(), xJcc, Op<T>::u(Var->num().I),
          VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
          paddcut_atancc_lin, loc::atan, &xJcc, 0 );
      }

      double xJcv = Op<T>::u(Var->num().I);
      xJcv = _newton( pFP, Op<T>::l(Var->num().I), Op<T>::l(Var->num().I), 0.,
        fct::atanfunc, &xJcv, 0 );
      if( mc::isequal(xJcv,Op<T>::l(Var->num().I)) ){
        _addcut_atancv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_atancv_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_atancv_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcv, 
          VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
          paddcut_atancv_lin, loc::atan, &xJcv, 0 );
      }
    }
  }
}

template <typename T> inline void
FPOp<T>::_TAN_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, std::tan(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, std::tan(Var->num().x),
      VarR->id(), 1. );
  else{
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    struct loc{ static std::pair<double,double> atan
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( std::atan(x), 1./(1.+x*x) ); }
    };

    // Append new linear cuts for univariate atan term
    // -- Convex Portion
    if( Op<T>::u(VarR->num().I) <= 0. ){
      _addcut_atancc_sec( pFP, VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), Var->id(), 0, 0 );
      void (FPOp<T>::*paddcut_atancv_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_atancv_lin;
      _sandwich_cuts( pFP, VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), paddcut_atancv_lin, loc::atan, 0, 0 );
    }

    // -- Concave Portion
    else if( Op<T>::l(VarR->num().I) >= 0. ){
      void (FPOp<T>::*paddcut_atancc_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_atancc_lin;
      _sandwich_cuts( pFP, VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), paddcut_atancc_lin, loc::atan, 0, 0 );
      _addcut_atancv_sec( pFP, VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), Var->id(), 0, 0 );
    }

    // -- Concavo-convex Portion
    else{
      struct fct{ static std::pair<double,double> atanfunc
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair(
            x-(*rusr)+(1+x*x)*(std::atan(*rusr)-std::atan(x)),
            2*x*(std::atan(*rusr)-std::atan(x)) ); }
      };

      double xJcc = Op<T>::l(VarR->num().I);
      xJcc = _newton( pFP, Op<T>::u(VarR->num().I), 0., Op<T>::u(VarR->num().I),
        fct::atanfunc, &xJcc, 0 );
      if( mc::isequal( xJcc, Op<T>::u(VarR->num().I) ) ){
        _addcut_atancc_sec( pFP, VarR->id(), Op<T>::l(VarR->num().I),
          Op<T>::u(VarR->num().I), Var->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_atancc_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_atancc_lin;
        _sandwich_cuts( pFP, VarR->id(), xJcc, Op<T>::u(VarR->num().I),
          Var->id(), Op<T>::l(Var->num().I), Op<T>::u(Var->num().I),
          paddcut_atancc_lin, loc::atan, &xJcc, 0 );
      }

      double xJcv = Op<T>::u(VarR->num().I);
      xJcv = _newton( pFP, Op<T>::l(VarR->num().I), Op<T>::l(VarR->num().I), 0.,
        fct::atanfunc, &xJcv, 0 );
      if( mc::isequal(xJcv,Op<T>::l(VarR->num().I)) ){
        _addcut_atancv_sec( pFP, VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), Var->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_atancv_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_atancv_lin;
        _sandwich_cuts( pFP, VarR->id(), Op<T>::l(VarR->num().I), xJcv, 
          Var->id(), Op<T>::l(Var->num().I), Op<T>::u(Var->num().I),
          paddcut_atancv_lin, loc::atan, &xJcv, 0 );
      }
    }
  }
}
  
template <typename T> inline void
FPOp<T>::_addcut_atancv_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // -x/(1+xref^2) + y >= atan(xref) - xref/(1+xref^2)
  const double dfref = 1./(1+Xref*Xref);
  pFP->_append_cut( this, FPCut<T>::GE, std::atan(Xref)-Xref*dfref,
    iY, 1., iX, -dfref );
}

template <typename T> inline void
FPOp<T>::_addcut_atancc_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // -x/(1+xref^2) + y <= atan(xref) - xref/(1+xref^2)
  const double dfref = 1./(1+Xref*Xref);
  pFP->_append_cut( this, FPCut<T>::LE, std::atan(Xref)-Xref*dfref,
    iY, 1., iX, -dfref );
}

template <typename T> inline void
FPOp<T>::_addcut_atancv_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y <= dy*xL - dx*yL
  double dX = XU - XL, YL = std::atan(XL), dY = std::atan(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::LE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_atancc_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y >= dy*xL - dx*yL
  double dX = XU - XL, YL = std::atan(XL), dY = std::atan(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_ERF_cuts
( FPRelax<T>*pFP, const FPVar<T>*VarR, const FPVar<T>*Var )
{
  if( Var->id().first == FPVar<T>::AUXINT )
    pFP->_append_cut( this, FPCut<T>::EQ, ::erf(Var->num().n),
      VarR->id(), 1. );
  else if( Var->id().first == FPVar<T>::AUXREAL )
    pFP->_append_cut( this, FPCut<T>::EQ, ::erf(Var->num().x),
      VarR->id(), 1. );
  else{
    typedef typename FPVar<T>::pt_idVar pt_idVar;
    struct loc{ static std::pair<double,double> erf
      ( const double x, const double*rusr, const int*iusr )
      { return std::make_pair( ::erf(x), 2./std::sqrt(PI)*std::exp(-x*x) ); }
    };

    // Append new linear cuts for univariate erf term
    // -- Convex Portion
    if( Op<T>::u(Var->num().I) <= 0. ){
      _addcut_erfcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      void (FPOp<T>::*paddcut_erfcv_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_erfcv_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_erfcv_lin, loc::erf, 0, 0 );
    }

    // -- Concave Portion
    else if( Op<T>::l(Var->num().I) >= 0. ){
      void (FPOp<T>::*paddcut_erfcc_lin)( FPRelax<T>*, const double, const pt_idVar, 
        const double, const double, const pt_idVar, const double, const double,
        const double*, const int* ) = &FPOp<T>::_addcut_erfcc_lin;
      _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), Op<T>::l(VarR->num().I),
        Op<T>::u(VarR->num().I), paddcut_erfcc_lin, loc::erf, 0, 0 );
      _addcut_erfcv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
    }

    // -- Concavo-convex Portion
    else{
      struct fct{ static std::pair<double,double> erffunc
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair(
            2./std::sqrt(PI)*(x-(*rusr))+(::erf(*rusr)-::erf(x))*std::exp(x*x),
            2.*x*std::exp(x*x)*(::erf(*rusr)-::erf(x)) ); }
      };

      double xJcc = Op<T>::l(Var->num().I);
      xJcc = _newton( pFP, Op<T>::u(Var->num().I), 0., Op<T>::u(Var->num().I),
        fct::erffunc, &xJcc, 0 );
      if( mc::isequal( xJcc, Op<T>::u(Var->num().I) ) ){
        _addcut_erfcc_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
          Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_erfcc_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_erfcc_lin;
        _sandwich_cuts( pFP, Var->id(), xJcc, Op<T>::u(Var->num().I),
          VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
          paddcut_erfcc_lin, loc::erf, &xJcc, 0 );
      }

      double xJcv = Op<T>::u(Var->num().I);
      xJcv = _newton( pFP, Op<T>::l(Var->num().I), Op<T>::l(Var->num().I), 0.,
        fct::erffunc, &xJcv, 0 );
      if( mc::isequal(xJcv,Op<T>::l(Var->num().I)) ){
        _addcut_erfcv_sec( pFP, Var->id(), Op<T>::l(Var->num().I),
        Op<T>::u(Var->num().I), VarR->id(), 0, 0 );
      }
      else{
        void (FPOp<T>::*paddcut_erfcv_lin)( FPRelax<T>*, const double, const pt_idVar, 
          const double, const double, const pt_idVar, const double, const double,
          const double*, const int* ) = &FPOp<T>::_addcut_erfcv_lin;
        _sandwich_cuts( pFP, Var->id(), Op<T>::l(Var->num().I), xJcv, 
          VarR->id(), Op<T>::l(VarR->num().I), Op<T>::u(VarR->num().I),
          paddcut_erfcv_lin, loc::erf, &xJcv, 0 );
      }
    }
  }
}
  
template <typename T> inline void
FPOp<T>::_addcut_erfcv_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // -x*erf'(xref) + y >= erf(xref) - xref*erf'(xref)
  const double dfref = 2./std::sqrt(PI)*std::exp(-Xref*Xref);
  pFP->_append_cut( this, FPCut<T>::GE, ::erf(Xref)-Xref*dfref,
    iY, 1., iX, -dfref );
}

template <typename T> inline void
FPOp<T>::_addcut_erfcc_lin
( FPRelax<T>*pFP, const double Xref, const typename FPVar<T>::pt_idVar iX,
  const double XL, const double XU, const typename FPVar<T>::pt_idVar iY,
  const double YL, const double YU, const double*rpar, const int*ipar )
{
  // -x*erf'(xref) + y <= erf(xref) - xref*erf'(xref)
  const double dfref = 2./std::sqrt(PI)*std::exp(-Xref*Xref);
  pFP->_append_cut( this, FPCut<T>::LE, ::erf(Xref)-Xref*dfref,
    iY, 1., iX, -dfref );
}

template <typename T> inline void
FPOp<T>::_addcut_erfcv_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y <= dy*xL - dx*yL
  double dX = XU - XL, YL = ::erf(XL), dY = ::erf(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::LE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

template <typename T> inline void
FPOp<T>::_addcut_erfcc_sec
( FPRelax<T>*pFP, const typename FPVar<T>::pt_idVar iX, const double XL,
  const double XU, const typename FPVar<T>::pt_idVar iY, const double*rpar,
  const int*ipar )
{
  // dy*x - dx*y >= dy*xL - dx*yL
  double dX = XU - XL, YL = ::erf(XL), dY = ::erf(XU) - YL;
  pFP->_append_cut( this, FPCut<T>::GE, dY*XL-dX*YL, iY, -dX, iX, dY );
}

///////////////////////////////// FPRelax //////////////////////////////////////

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPRelax<T>&FP)
{
  typename FPRelax<T>::t_Vars Vars = FP._Vars;
  typename FPRelax<T>::it_Vars itv = Vars.begin();
  out << ( FP._nvar? "\nVARIABLES:\n": "\nNO VARIABLE\n" );
  for( ; itv!=Vars.end() && (*itv)->_id.first<=FPVar<T>::VARBIN; ++itv ){
    out << "  " << **itv;
    //if( (*itv)->_Op ) out << "\t:= " << *(*itv)->_Op;
    out << std::endl;
  }
  out << ( FP._naux? "\nAUXILIARY:\n": "\nNO AUXILIARY\n" );
  for( ; itv!=Vars.end(); ++itv ){
    out << "  " << **itv;
    if( (*itv)->_Op ) out << "\t" << *(*itv)->_Op;
    out << std::endl;
  }

  typename FPRelax<T>::t_Cuts Cuts = FP._Cuts;
  typename FPRelax<T>::it_Cuts itc = Cuts.begin();
  out << ( Cuts.size()? "\nCUTS:\n": "\nNO CUTS\n" );
  for( ; itc!=Cuts.end(); ++itc )
    out << " " << **itc << std::endl;

  return out;
}

template <typename T> inline void
FPRelax<T>::set_objective
( const typename FPRelax<T>::OBJTYPE type, const FPVar<T>&objVar )
{
  // Erase any previous objective
  static FPOp<T> refOp( FPOp<T>::MIN );
  static const unsigned int nobjOp = 2;
  static typename FPOp<T>::TYPE objOp[nobjOp] = { FPOp<T>::MIN, FPOp<T>::MAX };
  for( unsigned int iOp=0; iOp<nobjOp; iOp++ ){
    refOp.type = objOp[iOp];
    std::pair< cit_Ops, cit_Ops > rangeOp = std::equal_range(
      _Ops.begin(), _Ops.end(), &refOp, range_FPOp<T>() );
    for( cit_Ops ito = rangeOp.first; ito != rangeOp.second; ++ito )
      delete *ito;
    _Ops.erase( rangeOp.first, rangeOp.second );
  }

  // Constant objective case
  if( objVar.id().second == FPVar<T>::NOREF ){
    double vobj = 0.;
    switch( objVar.num().t ){
      case FPVarNum<T>::INT:   vobj = (double)objVar.num().n; break;
      case FPVarNum<T>::REAL:  vobj = objVar.num().x; break;
      case FPVarNum<T>::RANGE:
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::OBJECTIVE );
    }
    FPVar<T>* pCst = _auxiliary_constant( vobj );
    switch( type ){
      case MIN: _operation( FPOp<T>::MIN, pCst->Oper()->pres ); break;
      case MAX: _operation( FPOp<T>::MAX, pCst->Oper()->pres ); break;
    }
    return;
  }

  // General case
  switch( type ){
    case MIN: _operation( FPOp<T>::MIN, objVar.Oper()->pres ); break;
    case MAX: _operation( FPOp<T>::MAX, objVar.Oper()->pres ); break;
  }
  return;
}

template <typename T> inline FPOp<T>*
FPRelax<T>::add_constraint
( const FPVar<T>&lhsVar, const typename FPRelax<T>::CTRTYPE type,
  const FPVar<T>&rhsVar )
{
  // Inconsistent case
  if( lhsVar.id().second == FPVar<T>::NOREF
   && rhsVar.id().second == FPVar<T>::NOREF ){
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::CONSTRAINT );
  }

  // Constant LHS case
  if( lhsVar.id().second == FPVar<T>::NOREF ){
    switch( type ){
      case LE: return add_constraint( rhsVar, GE, lhsVar );
      case EQ: return add_constraint( rhsVar, EQ, lhsVar );
      case GE: return add_constraint( rhsVar, LE, lhsVar );
    }
  }
  
  // Constant RHS case
  if( rhsVar.id().second == FPVar<T>::NOREF ){
    double vrhs = 0.;
    switch( rhsVar.num().t ){
      case FPVarNum<T>::INT:   vrhs = (double)rhsVar.num().n; break;
      case FPVarNum<T>::REAL:  vrhs = rhsVar.num().x; break;
      case FPVarNum<T>::RANGE:
        throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::CONSTRAINT );
    }
    FPVar<T>* pCst = _auxiliary_constant( vrhs );
    FPOp<T>* pOp = 0;
    switch( type ){
      case LE: pOp = _operation( FPOp<T>::LE, lhsVar.Oper()->pres,
        pCst->Oper()->pres ); break;
      case EQ: pOp = _operation( FPOp<T>::EQ, lhsVar.Oper()->pres,
        pCst->Oper()->pres ); break;
      case GE: pOp = _operation( FPOp<T>::LE, pCst->Oper()->pres,
        lhsVar.Oper()->pres ); break;
    }
    if( !pOp->pres ) pOp->pres = lhsVar.Oper()->pres;
    return pOp;
  }

  // General case
  FPOp<T>* pOp = 0;
  switch( type ){
    case LE: pOp = _operation( FPOp<T>::LE, lhsVar.Oper()->pres,
      rhsVar.Oper()->pres ); break;
    case EQ: pOp = _operation( FPOp<T>::EQ, lhsVar.Oper()->pres,
      rhsVar.Oper()->pres ); break;
    case GE: pOp = _operation( FPOp<T>::LE, rhsVar.Oper()->pres,
      lhsVar.Oper()->pres ); break;
  }
  if( !pOp->pres ) pOp->pres = lhsVar.Oper()->pres;
  return pOp;
}

template <typename T> inline FPOp<T>*
FPRelax<T>::_operation
( const typename FPOp<T>::TYPE type, FPVar<T>*lop, FPVar<T>*rop )
{
  FPOp<T>* op = new FPOp<T>( type, lop, rop );
  typename FPRelax<T>::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){ delete op; return *itop; }
  _Ops.insert( op );
  return op;
}

template <typename T> inline bool
FPRelax<T>::_erase_operation
( FPOp<T>* op )
{
  typename FPRelax<T>::it_Ops itop = _Ops.find( op );
  if( itop!=_Ops.end() ){
    _erase_cuts( op );
    delete op;
    _Ops.erase( itop );
    return true;
  }
  return false;
}
 
template <typename T> inline void
FPRelax<T>::_erase_cuts
( FPOp<T>* op )
{
  it_Cuts itc = _Cuts.begin(), itp;
  while( itc != _Cuts.end() ){
    itp = itc; ++itc;
    if( (*itp)->op() == op ){ delete *itp; _Cuts.erase( itp ); }
  } 
}
 
template <typename T> inline FPVar<T>*
FPRelax<T>::_auxiliary_variable
( const T&X, FPOp<T>*pOp )
{
  pOp->pres = new FPVar<T>( this, X, pOp );
  _append_aux( pOp->pres );
  return pOp->pres;
}
  
template <typename T> inline FPVar<T>*
FPRelax<T>::_auxiliary_variable
( FPOp<T>*pOp )
{
  pOp->pres = new FPVar<T>( this, pOp );
  _append_aux( pOp->pres );
  return pOp->pres;
}

template <typename T> inline FPVar<T>*
FPRelax<T>::_auxiliary_constant
( const double x )
{
  // Check if real constant x already defined in _Vars
  FPVar<T>* pAux = new FPVar<T>( x );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant x
  _append_aux( pAux, FPOp<T>::CNST );
  return pAux;
}

template <typename T> inline FPVar<T>*
FPRelax<T>::_auxiliary_constant
( const int n )
{
  // Check if real constant x already defined in _Vars
  FPVar<T>* pAux = new FPVar<T>( n );
  it_Vars iAux = _Vars.find( pAux );
  if( iAux!=_Vars.end() ){ delete pAux; return *iAux; }

  // Otherwise, append constant n
  _append_aux( pAux, FPOp<T>::CNST );
  return pAux;
}

template <typename T> inline FPVar<T>*
FPRelax<T>::_auxiliary_constant
( const T&I )
{
  // Append interval bounds I
  FPVar<T>* pAux = new FPVar<T>( I );
  _append_aux( pAux, FPOp<T>::RANGE );
  return pAux;
}

template <typename T> inline void
FPRelax<T>::_append_aux
( FPVar<T>*pAux, typename FPOp<T>::TYPE tOp )
{
  FPOp<T>*pOp = new FPOp<T>( tOp, 0, 0, pAux );
  _Ops.insert( pOp );
  pAux->FP() = this;
  pAux->Oper() = pOp;
  pAux->id().second = _naux++;
  _append_aux( pAux );
}

template <typename T> inline void
FPRelax<T>::_append_aux
( FPVar<T>*pAux )
{
  _Vars.insert( pAux );
}

template <typename T> inline void
FPRelax<T>::_append_var
( FPVar<T>*pVar )
{
  _Vars.insert( pVar );
}   

template <typename T> inline FPVar<T>*
FPRelax<T>::_find_var
( const typename FPVar<T>::pt_idVar&id )
{
  // Check if real constant x already defined in _Vars
  FPVar<T>* pVar = new FPVar<T>( this, id );
  it_Vars iVar = _Vars.find( pVar );
  delete pVar;
  return( iVar==_Vars.end()? 0: *iVar );
}

template <typename T> inline FPCut<T>*
FPRelax<T>::_append_cut
( FPOp<T>*op, const typename FPCut<T>::TYPE type,
  const double b, const typename FPVar<T>::pt_idVar id1,
  const double a1 )
{
  FPCut<T>* pCut = new FPCut<T>( op, type, b, id1, a1 );
  _Cuts.insert( pCut );
  return pCut;
}

template <typename T> inline FPCut<T>*
FPRelax<T>::_append_cut
( FPOp<T>*op, const typename FPCut<T>::TYPE type,
  const double b, const typename FPVar<T>::pt_idVar id1,
  const double a1, const typename FPVar<T>::pt_idVar id2,
  const double a2 )
{
  FPCut<T>* pCut = new FPCut<T>( op, type, b, id1, a1, id2, a2 );
  _Cuts.insert( pCut );
  return pCut;
}

template <typename T> inline FPCut<T>*
FPRelax<T>::_append_cut
( FPOp<T>*op, const typename FPCut<T>::TYPE type,
  const double b, const typename FPVar<T>::pt_idVar id1,
  const double a1, const typename FPVar<T>::pt_idVar id2,
  const double a2, const typename FPVar<T>::pt_idVar id3,
  const double a3 )
{
  FPCut<T>* pCut = new FPCut<T>( op, type, b, id1, a1, id2, a2, id3, a3 );
  _Cuts.insert( pCut );
  return pCut;
}

template <typename T> inline FPCut<T>*
FPRelax<T>::_append_cut
( FPOp<T>*op, const typename FPCut<T>::TYPE type,
  const double b, const unsigned int n,
  const typename FPVar<T>::pt_idVar* id, const double*a )
{
  if( !n ) throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  FPCut<T>* pCut = new FPCut<T>( op, type, b, n, id, a );
  _Cuts.insert( pCut );
  return pCut;
}

template <typename T> inline void
FPRelax<T>::generate_cuts
( const bool reset )
{
  // Reset cuts in polyhedral relaxation?
  if( reset ){
    _erase_cuts();
    _reset_operations();
  }
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nOPERATIONS BEFORE PROPAGATION:" << std::endl;
  for( it_Ops ito = _Ops.begin(); ito != _Ops.end(); ++ito )
    std::cout << **ito << "   \t" << (*ito)->status() << std::endl;
#endif

  // Generate cuts for all objective/constraints
  static FPOp<T> refOp( FPOp<T>::EQ );
  static const unsigned int nfctOp = 4;
  static typename FPOp<T>::TYPE fctOp[nfctOp] = { FPOp<T>::EQ, FPOp<T>::LE,
    FPOp<T>::MIN, FPOp<T>::MAX };
  for( unsigned int iOp=0; iOp<nfctOp; iOp++ ){
    refOp.type = fctOp[iOp];
    std::pair< cit_Ops, cit_Ops > rangeOp = std::equal_range(
      _Ops.begin(), _Ops.end(), &refOp, range_FPOp<T>() );
    cit_Ops ito = rangeOp.first;
    for( ; ito != rangeOp.second; ++ito ) (*ito)->generate_cuts( this );
  }
 
#ifdef MC__FPRELAX_DEBUG_PROPAGATION
  std::cout << "\nOPERATIONS AFTER PROPAGATION:" << std::endl;
  for( it_Ops ito = _Ops.begin(); ito != _Ops.end(); ++ito )
    std::cout << **ito << "   \t" << (*ito)->status() << std::endl;
#endif

  return;
}

template <typename T> inline void
FPRelax<T>::generate_cuts
( FPOp<T>*constr )
{
  // Generate cuts for constr only
  if( !constr ) return;
  cit_Ops itop = _Ops.find( constr );
  if( itop == _Ops.end() ) return;
  return (*itop)->generate_cuts( this );
}

template <typename T> inline bool
FPRelax<T>::propagate_constraints
( FPOp<T>*constr )
{
  // Propagate constraint for constr only
  if( constr ){
    cit_Ops itop = _Ops.find( constr );
    if( itop == _Ops.end() ) return true;
    return (*itop)->propagate_bounds( this );
  }

  // Propagate all constraints
  static FPOp<T> refOp( FPOp<T>::EQ );
  static const unsigned int nfctOp = 4;
  static typename FPOp<T>::TYPE fctOp[nfctOp] = { FPOp<T>::EQ, FPOp<T>::LE,
    FPOp<T>::MIN, FPOp<T>::MAX };
  bool keepgoing = true;
  for( unsigned int iOp=0; iOp<nfctOp && keepgoing; iOp++ ){
    refOp.type = fctOp[iOp];
    std::pair< cit_Ops, cit_Ops > rangeOp = std::equal_range(
      _Ops.begin(), _Ops.end(), &refOp, range_FPOp<T>() );
    cit_Ops ito = rangeOp.first;
    for( ; ito != rangeOp.second && keepgoing; ++ito )
      keepgoing = (*ito)->propagate_bounds( this );
  }
  return keepgoing;
}

template <typename T> inline void
FPRelax<T>::generate_reduction_constraints()
{
  FPRRLT<T> RRLT( _Vars, _Ops );
  typename FPRRLT<T>::t_RLTMap& RLTCuts
    = RRLT.map_RLT( options.RRLT_STRATEGY, options.RRLT_DISPLAY );
  typename FPRRLT<T>::cit_RLTMap it = RLTCuts.begin();

  for( ; it != RLTCuts.end(); ++it ){
    switch( (*it).first.second ){
    // Case reduction constraint via multiplier variable
    case FPOp<T>::TIMES:
      switch( (*it).second->type ){
        case FPOp<T>::EQ:
          *(*it).first.first * *(*it).second->plop
            == *(*it).first.first * *(*it).second->prop;
          break;
        case FPOp<T>::LE:
          *(*it).first.first * *(*it).second->plop - Op<T>::l((*it).first.first->I())
            * ( *(*it).second->plop - *(*it).second->prop )
            <= *(*it).first.first * *(*it).second->prop;
          *(*it).first.first * *(*it).second->plop - Op<T>::u((*it).first.first->I())
            * ( *(*it).second->plop - *(*it).second->prop )
            >= *(*it).first.first * *(*it).second->prop;
          break;
        case FPOp<T>::PLUS:
          *(*it).first.first * *(*it).second->pres
            == *(*it).first.first * *(*it).second->plop
             + *(*it).first.first * *(*it).second->prop;
          break;
        case FPOp<T>::NEG:
          *(*it).first.first * *(*it).second->pres
            + *(*it).first.first * *(*it).second->plop == 0;
          break;
        case FPOp<T>::MINUS:
          *(*it).first.first * *(*it).second->pres
            == *(*it).first.first * *(*it).second->plop
             - *(*it).first.first * *(*it).second->prop;
          break;
        case FPOp<T>::SCALE:
          *(*it).first.first * *(*it).second->pres
            == ( *(*it).first.first * *(*it).second->plop )
             * *(*it).second->prop;
          break;
        default:
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      }
      break;
      
    // Case reduction constraint via divider variable
    case FPOp<T>::DIV:
      if( Op<T>::l((*it).first.first->I()) * Op<T>::u((*it).first.first->I()) <= 0. )
        break;
      switch( (*it).second->type ){
        case FPOp<T>::EQ:
          *(*it).second->plop / *(*it).first.first
            == *(*it).second->prop / *(*it).first.first;
          break;
        case FPOp<T>::LE:
          //*(*it).first.first * *(*it).second->plop - Op<T>::l((*it).first.first->I())
          //  * ( *(*it).second->plop - *(*it).second->prop )
          //  <= *(*it).first.first * *(*it).second->prop;
          //*(*it).first.first * *(*it).second->plop - Op<T>::u((*it).first.first->I())
          //  * ( *(*it).second->plop - *(*it).second->prop )
          //  >= *(*it).first.first * *(*it).second->prop;
          break;
        case FPOp<T>::PLUS:
          *(*it).second->pres / *(*it).first.first
            == *(*it).second->plop / *(*it).first.first
             + *(*it).second->prop / *(*it).first.first;
          break;
        case FPOp<T>::NEG:
           *(*it).second->pres / *(*it).first.first
            + *(*it).second->plop / *(*it).first.first == 0;
          break;
        case FPOp<T>::MINUS:
          *(*it).second->pres / *(*it).first.first
            == *(*it).second->plop / *(*it).first.first
             - *(*it).second->prop / *(*it).first.first;
          break;
        case FPOp<T>::SCALE:
          *(*it).second->pres / *(*it).first.first
            == ( *(*it).second->plop / *(*it).first.first )
             * *(*it).second->prop;
          break;
        default:
          throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
      }
      break;

    // Oups - Reduction constraint undefined...
    default:
      throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
    }
  }
}

template <typename T> template <typename U> inline mc::TModel< FPVar<T> >*
FPRelax<T>::create_TModel
( const mc::TModel<U>*TM, const unsigned int*isub )
{
  // Build Taylor model environment mc::TModel< FPVar<T> >
  if( !TM ) return 0;
  typedef mc::TModel< FPVar<T> > TMFPT;
  typedef mc::TVar< FPVar<T> > TVFPT;
  TMFPT *TMFP = new TMFPT( TM->nvar(), TM->nord() );
  TMFP->options = TM->options;
  for( unsigned int ix=0; ix<TM->nvar(); ix++ ){
    FPVar<T>*pVar = _find_var( pt_idVar( FPVar<T>::VARCONT, isub?isub[ix]:ix ) );
    if( !pVar ) pVar = _find_var( pt_idVar( FPVar<T>::VARBIN, isub?isub[ix]:ix ) );
    // Throw exception if variable ix+ioff is undefined
    if( !pVar ) throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );
    TVFPT TVFP( TMFP, ix, *pVar );
  }
  return TMFP;
}

template <typename T> template <typename U> inline void
FPRelax<T>::update_TModel
( mc::TModel< FPVar<T> >*TMFP, const mc::TModel<U>*TM,
  const unsigned int*isub )
{
  if( !TM || !TMFP || TMFP->nvar()!=TM->nvar() )
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );

  // Update Taylor model environment mc::TModel< FPVar<T> >
  typedef mc::TModel< FPVar<T> > TMFPT;
  typedef mc::TVar< FPVar<T> > TVFPT;
  TMFP->options = TM->options;
  TMFP->reset();
  for( unsigned int ix=0; ix<TM->nvar(); ix++ ){
    FPVar<T>*pVar = _find_var( pt_idVar( FPVar<T>::VARCONT, isub?isub[ix]:ix ) );
    if( !pVar ) pVar = _find_var( pt_idVar( FPVar<T>::VARBIN, isub?isub[ix]:ix ) );
    // Throw exception if variable ix+ioff is undefined
    if( !pVar ) throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INIT );
    TVFPT TVFP( TMFP, ix, *pVar );
  }
}

template <typename T> inline void
FPRelax<T>::generate_RLT_TModel
( mc::TModel< FPVar<T> >*TMFP )
{
  if( !TMFP || options.RLT_MAXORDER < 2 ) return;
  unsigned int MAXORD = ( options.RLT_MAXORDER < TMFP->nord()?
                          options.RLT_MAXORDER: TMFP->nord() );

  // Generate RLT constraints up to order MAXORD
  std::list< FPRLT<T>* > RLT_lst;
  RLT_lst.push_back( new FPRLT<T>( TMFP ) );
  for( unsigned int iord=1; iord<=MAXORD; iord++ ){
    const bool disp = ( options.RLT_DISPLAY >= 2
                        || ( options.RLT_DISPLAY == 1 && iord == MAXORD )?
                        true: false );
    if( disp ) std::cout << "\nRLT Constraints: Order" << iord << std::endl;
    _lift_RLT_TModel( TMFP, RLT_lst, disp );
  }

  // Append RLT cuts of order MAXORD
  _append_cuts_RLT_TModel( RLT_lst, MAXORD );
}
 
template <typename T> inline void
FPRelax<T>::_append_cuts_RLT_TModel
( std::list< FPRLT<T>* >& RLT_lst, const unsigned int nord )
{
  // Create new RLT operation - needed to order cuts
  FPOp<T>* op = new FPOp<T>( FPOp<T>::RLT );
  _Ops.insert( op );

  // Allocate variable indices/coefficient arrays (no more than 2*nord terms)
  const unsigned int lhsmax = std::pow(2,nord);
  std::pair<typename FPVar<T>::pt_idVar*, double*> lhs
    = make_pair( new typename FPVar<T>::pt_idVar[lhsmax], new double[lhsmax] );

  // Append all cuts in RLT list
  typedef typename std::list< FPRLT<T>* >::iterator it_RLT;
  it_RLT it = RLT_lst.begin();
  for( ; it != RLT_lst.end(); ++it ){
    // Populate variable indices/coefficient arrays
    assert( (*it)->lhs2cut( lhs, lhsmax ) );
    // Append cut 
    typename FPCut<T>::TYPE type = ( (*it)->type()==FPRLT<T>::LE?
                                     FPCut<T>::LE: FPCut<T>::GE );
    _append_cut( op, type, (*it)->rhs(), (*it)->lhs().size(),
                 lhs.first, lhs.second );
    // Clean up RLT entry
    delete (*it); (*it) = 0;
  }
  // Clean up index/coefficient arrays
  delete[] lhs.first; delete[] lhs.second;
}

template <typename T> inline void
FPRelax<T>::_lift_RLT_TModel
( mc::TModel< FPVar<T> >*TMFP, std::list< FPRLT<T>* >& RLT_lst,
  const bool disp )
{
  // Lift current RLT constraints by one order
  const unsigned int nlst = RLT_lst.size();
  for( unsigned int ilst=0; ilst<nlst; ilst++ ){
    const std::pair<unsigned int,bool> flag = (*RLT_lst.begin())->flag();

    // Multiply every RLT constraint by every variable, avoiding duplicates
    for( unsigned int ivar=flag.first; ivar< TMFP->nvar(); ivar++ ){
      const unsigned int ipos = TMFP->nvar() - ivar;
      unsigned int ipb = ( ivar==flag.first? 1-flag.second: 0 );
      for( ; ipb<=1; ipb++ ){
        const std::pair<unsigned int,bool> newflag( ivar, 1-ipb );
        //std::cout << pTMFP_bnd[ipos] << std::endl;
        FPRLT<T>* RLTprod = (*RLT_lst.begin())->product(
          TMFP->bndmon()+ipos, TMFP->nvar(),
          TMFP->expmon()+ipos*TMFP->nvar(), newflag
        );
        RLT_lst.push_back( RLTprod );  // insert new RLT constraint for product with X-XL or X-XU
        if( disp ) std::cout << *RLTprod << std::endl;
      }
    }

    // Remove previous RLT constraint - Only keep new lifted RLT constraints
    delete (*RLT_lst.begin());
    RLT_lst.pop_front();
  }
  if( disp ){ int dum; std::cin >> dum; }
}

////////////////////////////////// FPRLT //////////////////////////////////////

template <typename T> inline std::pair<typename FPVar<T>::pt_idVar*, double*>
FPRLT<T>::lhs2cut
() const
{
  // Allocate arrays of variable indices and coefficients
  double*coef = new double[_lhs.size()];
  typename FPVar<T>::pt_idVar*id = new typename FPVar<T>::pt_idVar[_lhs.size()];

  // Populate arrays of variable indices and coefficients
  cit_lhs cit = _lhs.begin();
  for( unsigned int ilhs=0; cit != _lhs.end(); ++cit, ilhs++ ){
    id[ilhs] = (*cit).first->id();
    coef[ilhs] = (*cit).second.coef;
  }
  return std::make_pair(id,coef);  
}

template <typename T> inline bool
FPRLT<T>::lhs2cut
( std::pair<typename FPVar<T>::pt_idVar*, double*>&lhs,
  const unsigned int lhsmax ) const
{
  //std::cout << lhsmax << "  " << _lhs.size() << std::endl;
  if( _lhs.size() > lhsmax ) return false;
  
  // Populate arrays of variable indices and coefficients
  cit_lhs cit = _lhs.begin();
  for( unsigned int ilhs=0; cit != _lhs.end(); ++cit, ilhs++ ){
    lhs.first[ilhs] = (*cit).first->id();
    lhs.second[ilhs] = (*cit).second.coef;
  }
  return true;
}

template <typename T> inline FPRLT<T>*
FPRLT<T>::product
( const FPVar<T>*Var, const unsigned int nexp, const unsigned int*iexp,
  std::pair<unsigned int, bool> flag ) const
{
  // Generate new RLT constraint via product with (X-XL) or (X-XU)
  FPRLT<T>* RLT = new FPRLT<T>( _pTM );

  // Initialize RLT constraint
  if( _lhs.empty() ){
    RLT->_lhs.insert( p_lhs( Var, RLTTerm( 1., nexp, iexp ) ) );
    RLT->_rhs = ( flag.second? Op<T>::l( Var->I() ): Op<T>::u( Var->I() ) );
    RLT->_type = ( flag.second? GE: LE );
    RLT->_flag = flag;
    return RLT;
  }

  // Account for LHS and RHS product with XL/XU
  cit_lhs cit = _lhs.begin();
  const double c = ( flag.second? -Op<T>::l( Var->I() ): Op<T>::u( Var->I() ) );
  for( ; !isequal( c, 0. ) && cit != _lhs.end(); ++cit )
    RLT->_lhs.insert(
      p_lhs( (*cit).first, RLTTerm( c * (*cit).second.coef, nexp, (*cit).second.iexp ) )
    );
  RLT->_rhs = c * _rhs;

  // Account for RHS product with X
  it_lhs it = RLT->_lhs.find( Var );
  if( !isequal( _rhs, 0. ) && it != RLT->_lhs.end() ){
    (*it).second.coef += ( flag.second? -_rhs: _rhs );
    if( isequal( (*it).second.coef, 0. ) ) RLT->_lhs.erase( it );
  }
  else if( !isequal( _rhs, 0. ) )
    RLT->_lhs.insert( p_lhs( Var, RLTTerm( ( flag.second? -_rhs: _rhs ), nexp, iexp ) ) );

  // Account for LHS product with X
  cit = _lhs.begin();
  unsigned int*iexpprod = new unsigned int[nexp];
  assert( _pTM );
  for( ; cit != _lhs.end(); ++cit ){

    // Location of new product term
    for( unsigned int i=0; i<nexp; i++ )
      iexpprod[i] = (*cit).second.iexp[i] + iexp[i];
    const FPVar<T>* Varprod = _pTM->bndmon() + _pTM->loc_expmon(iexpprod);
    it = RLT->_lhs.find( Varprod );

    // Account for new product term in LHS
    if( !isequal( (*cit).second.coef, 0. ) && it != RLT->_lhs.end() ){
      (*it).second.coef += ( flag.second? (*cit).second.coef: -(*cit).second.coef );
      if( isequal( (*it).second.coef, 0. ) ) RLT->_lhs.erase( it );
    }
    else if( !isequal( (*cit).second.coef, 0. ) )
      RLT->_lhs.insert(
        p_lhs( Varprod, RLTTerm( ( flag.second? (*cit).second.coef: -(*cit).second.coef ),
                                 nexp, iexpprod ) )
      );
  }
  delete[] iexpprod;

  // Set new RLT constraint sense and flag, then return
  RLT->_type = _type;
  RLT->_flag = flag;
  return RLT;
}

template <typename T> inline std::ostream&
operator <<
( std::ostream&out, const FPRLT<T>&RLT )
{
  const int iprec = 5;
  out << std::right << std::scientific << std::setprecision(iprec);

  // RLT constraint is uninitialized
  if( RLT._lhs.empty() ){
    out << " Uninitialized RLT constraint";
    return out;
  }

  // Display RLT constraint LHS
  typename FPRLT<T>::cit_lhs it = RLT._lhs.begin();
  for( ; it != RLT._lhs.end(); ++it ){
    if( isequal( (*it).second.coef, 0. ) )
      out << " + " << std::setw(iprec+6) << 0.;
    else if( (*it).second.coef > 0. )
      out << " + " << std::setw(iprec+6) <<  (*it).second.coef;
    else
      out << " - " << std::setw(iprec+6) << -(*it).second.coef;
    out << FPVar<T>::_name( (*it).first->id() );  
  }
  
  // Display RLT constraint sense
  switch( RLT._type ){
    case FPRLT<T>::LE: out << " <= "; break;
    case FPRLT<T>::GE: out << " >= "; break;
  }
  
  // Display RLT constraint LHS
  out << std::setw(iprec+6) << RLT._rhs;
  return out;
}

////////////////////////////////// FPRRLT //////////////////////////////////////

template <typename T> inline void
FPRRLT<T>::_define_linear()
{
  static const unsigned int nlinOp = 6;
  static const typename FPOp<T>::TYPE linOp[nlinOp] = { FPOp<T>::EQ,
    FPOp<T>::LE, FPOp<T>::PLUS, FPOp<T>::NEG, FPOp<T>::MINUS, FPOp<T>::SCALE };

  // Create subset of linear operations
  _linTerms  = subset_op( nlinOp, linOp );
  if( _display >= 2 ){
    std::cout << "\nLINEAR TERMS ACCOUNTED FOR IN RLT:\n";
    for( it_Ops it = _linTerms.begin(); it != _linTerms.end(); ++it )
      std::cout << *(*it) << std::endl;
  }
  _nlinTerms = _linTerms.size();
  delete[] _TermVisited;
  _TermVisited  = new bool[ _nlinTerms ];

  // Create subset of variables participating in linear operations
  _linVars = subset_var( _linTerms );
  if( _display >= 2 ){
    std::cout << "\nVARIABLES PARTICIPATING IN RLT TERMS:\n";
    for( it_Vars it = _linVars.begin(); it != _linVars.end(); ++it )
      std::cout << *(*it) << std::endl;
  }

  // Create map of linear operations and participating variables
  _linEdges = submap_var( _linTerms );
  if( _display >= 2 ){
    std::cout << "\nRLT TERMS <-> VARIABLES EDGES:\n";
    for( it_Vars_Ops it = _linEdges.begin(); it != _linEdges.end(); ++it )
      std::cout << *(*it).first << "  <->  " << *(*it).second << std::endl;
  }

  // Update size of linear operation arrays
  if( _nlinTerms == _linTerms.size() ) return;
  _nlinTerms = _linTerms.size();
  delete[] _TermVisited;
  _TermVisited  = new bool[ _nlinTerms ];
}

template <typename T> inline void
FPRRLT<T>::_define_bilinear()
{
  // subsets of bilinear terms
  static const unsigned int nbilOp = 2;
  static typename FPOp<T>::TYPE bilOp[nbilOp] = { FPOp<T>::TIMES,
    FPOp<T>::BILIN };
  _bilTerms  = subset_op( nbilOp, bilOp );

  // subsets of fractional terms
  _divTerms  = subset_op( FPOp<T>::DIV );

  // subsets of square terms
  _sqrTerms  = subset_op( FPOp<T>::SQR );

  // subsets of square root terms
  _sqrtTerms = subset_op( FPOp<T>::SQRT );
}

template <typename T> inline void
FPRRLT<T>::_bigraph_RLT_mult
( typename FPRRLT<T>::t_Vars &linVars, typename FPRRLT<T>::t_Vars_Ops &linEdges,
  const typename FPVar<T>::pt_idVar& idMult )
{
  // Erase those variables participating with <a>idMult</a> in bilinear terms
  for( cit_Ops it = _bilTerms.begin(); it != _bilTerms.end(); ++it ){
    if( (*it)->plop->id() == idMult ){
      linVars.erase( (*it)->prop ); linEdges.erase( (*it)->prop );
    }
    else if( (*it)->prop->id() == idMult ){
      linVars.erase( (*it)->plop ); linEdges.erase( (*it)->plop );
    }
  }

  // Erase those variables participating with <a>idMult</a> in fractional terms
  for( cit_Ops it = _divTerms.begin(); it != _divTerms.end(); ++it ){
    if( (*it)->pres->id() == idMult ){
      linVars.erase( (*it)->prop ); linEdges.erase( (*it)->prop );
    }
    else if( (*it)->prop->id() == idMult ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idMult</a> in square terms
  for( cit_Ops it = _sqrTerms.begin(); it != _sqrTerms.end(); ++it ){
    if( (*it)->plop->id() == idMult ){
      linVars.erase( (*it)->plop ); linEdges.erase( (*it)->plop );
    }
  }

  // Erase those variables participating with <a>idMult</a> in square root terms
  for( cit_Ops it = _sqrtTerms.begin(); it != _sqrtTerms.end(); ++it ){
    if( (*it)->pres->id() == idMult ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  if( _display >= 2 ){
    std::cout << "\nVARIABLES PARTICIPATING IN RLT TERMS,\n"
              << "YET NOT PARTICIPATING IN ANY BILINEAR TERM WITH "
              << FPVar<T>::_name( idMult ) << ":\n";
    for( cit_Vars itv = linVars.begin(); itv != linVars.end(); ++itv )
      std::cout << *(*itv) << std::endl;
    std::cout << "\nRLT GRAPH FOR : " << FPVar<T>::_name( idMult )
              << ": OPERATION <-> VARIABLE\n";
    for( cit_Vars_Ops itvo = linEdges.begin(); itvo != linEdges.end(); ++itvo )
      std::cout << *(*itvo).second << "  <->  " << *(*itvo).first << std::endl;
  }
}

template <typename T> inline void
FPRRLT<T>::_bigraph_RLT_div
( typename FPRRLT<T>::t_Vars &linVars, typename FPRRLT<T>::t_Vars_Ops &linEdges,
  const typename FPVar<T>::pt_idVar& idDiv )
{
  // Erase <a>idDiv</a> since division by itself yield 1
  for( cit_Vars it = linVars.begin(); it != linVars.end(); ++it )
    if( (*it)->id() == idDiv ){ linVars.erase( *it ); break; }

  // Erase those variables participating with <a>idDiv</a> in bilinear terms
  for( cit_Ops it = _bilTerms.begin(); it != _bilTerms.end(); ++it ){
    if( (*it)->plop->id() == idDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
    else if( (*it)->prop->id() == idDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in fractional terms
  for( cit_Ops it = _divTerms.begin(); it != _divTerms.end(); ++it ){
    if( (*it)->pres->id() == idDiv ){
      linVars.erase( (*it)->plop ); linEdges.erase( (*it)->plop );
    }
    else if( (*it)->prop->id() == idDiv ){
      linVars.erase( (*it)->plop ); linEdges.erase( (*it)->plop );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in square terms
  for( cit_Ops it = _sqrTerms.begin(); it != _sqrTerms.end(); ++it ){
    if( (*it)->plop->id() == idDiv ){
      linVars.erase( (*it)->pres ); linEdges.erase( (*it)->pres );
    }
  }

  // Erase those variables participating with <a>idDiv</a> in square root terms
  for( cit_Ops it = _sqrtTerms.begin(); it != _sqrtTerms.end(); ++it ){
    if( (*it)->pres->id() == idDiv ){
      linVars.erase( (*it)->plop ); linEdges.erase( (*it)->plop );
    }
  }

  if( _display >= 2 ){
    std::cout << "\nVARIABLES PARTICIPATING IN RLT TERMS,\n"
              << "YET NOT PARTICIPATING IN ANY FRACTIONAL TERM WITH "
              << FPVar<T>::_name( idDiv ) << ":\n";
    for( cit_Vars itv = linVars.begin(); itv != linVars.end(); ++itv )
      std::cout << *(*itv) << std::endl;
    std::cout << "\nRLT GRAPH FOR : " << FPVar<T>::_name( idDiv )
              << ": OPERATION <-> VARIABLE\n";
    for( cit_Vars_Ops itvo = linEdges.begin(); itvo != linEdges.end(); ++itvo )
      std::cout << *(*itvo).second << "  <->  " << *(*itvo).first << std::endl;
  }
}

template <typename T> inline void
FPRRLT<T>::_bigraph_RLT
( const typename FPVar<T>::pt_idVar& idVar, const typename FPOp<T>::TYPE& idOp )
{
  // Initialize variable nodes and edges in bigraph
  _RLTVars = _linVars;
  _RLTEdges = _linEdges;

  // Erase those variables participating in nonlinear terms with variable <a>idVar</a>
  switch( idOp ){
  case FPOp<T>::TIMES:
    _bigraph_RLT_mult( _RLTVars, _RLTEdges, idVar ); break;
  case FPOp<T>::DIV:
    _bigraph_RLT_div( _RLTVars, _RLTEdges, idVar );  break;
  default:
    throw typename FPRelax<T>::Exceptions( FPRelax<T>::Exceptions::INTERNAL );
  }

  // Update size of RLT variable arrays
  if(  _nRLTVars == _RLTVars.size() ) return;
  _nRLTVars = _RLTVars.size();
  delete[] _VarAssigned;
  delete[] _VarVisited;
  _VarAssigned = new pt_Op[ _nRLTVars ];
  _VarVisited  = new bool[ _nRLTVars ];
}

template <typename T> inline typename FPRRLT<T>::t_RLTMap&
FPRRLT<T>::map_RLT
( typename FPRelax<T>::Options::RRLT option, unsigned int display )
{
  // No RLT constraint requested
  _RLTMap.clear();
  //if( option == FPRelax<T>::Options::NORRLT ) return _RLTMap;
  _display = display;

  // Define subsets and map of linear operations and participating variables
  _define_linear();

  // Define subsets of bilinear, fractional, square and square root terms
  _define_bilinear();

  // MAIN LOOP: Valid RRLT cut for each candidate multiplier/divider variables
  for( cit_Vars itVar = _Vars.begin(); itVar != _Vars.end(); ++itVar ){
    const pt_idVar& idVar = (*itVar)->id();

    // Multiplying/dividing constraints by a constant is worthless
    if( idVar.first == FPVar<T>::AUXINT
      || idVar.first == FPVar<T>::AUXREAL ) continue;

    // Multiplier/divider variables limited to primary variables (no auxiliary)
    if( option == FPRelax<T>::Options::PRIMRRLT
    && ( idVar.first == FPVar<T>::AUXCONT
      || idVar.first == FPVar<T>::AUXBIN ) ) continue;

    // Erase those variables participating in nonlinear terms multiplied by
    // variable <a>idVar</a>
    _bigraph_RLT( idVar, FPOp<T>::TIMES );

    // Identify valid reduction constraints with the candidate multiplier
    // variable
    _reduction_RLT( std::make_pair( (*itVar), FPOp<T>::TIMES ) );

    // Identify variables not yet participating in any bilinear terms with
    // the candidate multiplier variable
    _bigraph_RLT( idVar, FPOp<T>::DIV );

    // Identify valid reduction constraints with the candidate divider
    // variable
    _reduction_RLT( std::make_pair( (*itVar), FPOp<T>::DIV ) );
  }

  if( _display >= 1 ){
    if( _RLTMap.begin() == _RLTMap.end() )
      std::cout << "\n NO VALID REDUCTION CONSTRAINTS FOUND\n";
    else{
      std::cout << "\n VALID REDUCTION CONSTRAINTS: "
                << " REDUCTION VARIABLE <-TYPE-> OPERATION\n";
      for( cit_RLTMap it = _RLTMap.begin(); it != _RLTMap.end(); ++it ){
        std::cout << "  " << *(*it).first.first << "  <";
        switch( (*it).first.second ){
        case FPOp<T>::TIMES:
          std::cout << "MULT"; break;
        case FPOp<T>::DIV:
          std::cout << "DIV";  break;
        default:
          std::cout << "???";  break;
        }
        std::cout << ">  " << *(*it).second << std::endl;
      }
    }
  }

  return _RLTMap;
}

template <typename T> inline void
FPRRLT<T>::_reduction_RLT
( const t_RLTVar& VarRed )
{
  // Matching initially empty
  FPOp<T>* NA = 0;
  for( unsigned int iv=0; iv<_nRLTVars; iv++ )
    _VarAssigned[iv] = std::make_pair( NA, 0 );

  // Try to construct an augmenting emanating from each constraint
  cit_Ops ito = _linTerms.begin();
  for( unsigned int io=0; ito != _linTerms.end(); ++ito, io++ ){
    for( unsigned int j=0; j<_nlinTerms; j++ ) _TermVisited[j] = false;
    for( unsigned int i=0; i<_nRLTVars; i++ )  _VarVisited[i] = false;
    // No augmenting emanating found from current constraint
    if( !_augpath_RLT( std::make_pair(*ito,io) ) ){
      cit_Ops jto = _linTerms.begin();
      // Append new valid reduction RLT constraints
      for( unsigned int jo=0; jto != _linTerms.end(); ++jto, jo++ )
        if( _TermVisited[jo] )
          _RLTMap.insert( std::make_pair( VarRed, *jto ) );
    }
  }
}

template <typename T> inline bool
FPRRLT<T>::_augpath_RLT
( const pt_Op& pOp )
{
  _TermVisited[pOp.second] = true;

  // Try and find immediate assignment
  cit_Vars itv = _RLTVars.begin();
  for( unsigned int iv=0; itv != _RLTVars.end(); ++itv, iv++ ){
    std::pair< cit_Vars_Ops, cit_Vars_Ops > rangeOp;
    if( _VarAssigned[iv].first ) continue;
    rangeOp = _RLTEdges.equal_range( *itv );
    for( cit_Vars_Ops itvo = rangeOp.first; itvo != rangeOp.second; ++itvo ){
      if( (*itvo).second != pOp.first ) continue;
      _VarAssigned[iv] = pOp;
      return true;
    }
  }

  // Try and find an augmenting path starting from another variable
  itv = _RLTVars.begin();
  for( unsigned int iv=0; itv != _RLTVars.end(); ++itv, iv++ ){
    std::pair< cit_Vars_Ops, cit_Vars_Ops > rangeOp;
    if( _VarVisited[iv] ) continue;
    rangeOp = _RLTEdges.equal_range( *itv );
    for( cit_Vars_Ops itvo = rangeOp.first; itvo != rangeOp.second; ++itvo ){
      if( (*itvo).second != pOp.first ) continue;
      _VarVisited[iv] = true;
      if( !_augpath_RLT( _VarAssigned[iv] ) ) continue;
      _VarAssigned[iv] = pOp;
      return true;
    }
  }

  // Failed to find an augmenting path
  return false;
}

template <typename T> inline typename FPRRLT<T>::t_Ops
FPRRLT<T>::subset_op
( const unsigned int nOp, const typename FPOp<T>::TYPE*typeOp ) const
{
  // Create subset of operations of given type
  t_Ops Ops;
  static FPOp<T> refOp( FPOp<T>::CNST );
  for( unsigned int iOp=0; iOp<nOp; iOp++ ){
    refOp.type = typeOp[iOp];
    std::pair< cit_Ops, cit_Ops > rangeOp = std::equal_range(
      _Ops.begin(), _Ops.end(), &refOp, range_FPOp<T>() );
    Ops.insert( rangeOp.first, rangeOp.second );
  }
  return Ops;
}

template <typename T> inline typename FPRRLT<T>::t_Ops
FPRRLT<T>::subset_op
( const typename FPOp<T>::TYPE&typeOp ) const
{
  // Create subset of operations of given type
  t_Ops Ops;
  static FPOp<T> refOp( FPOp<T>::CNST );
  refOp.type = typeOp;
  std::pair< cit_Ops, cit_Ops > rangeOp = std::equal_range(
    _Ops.begin(), _Ops.end(), &refOp, range_FPOp<T>() );
  Ops.insert( rangeOp.first, rangeOp.second );
  return Ops;
}

template <typename T> inline typename FPRRLT<T>::t_Vars
FPRRLT<T>::subset_var
( const typename FPRRLT<T>::t_Ops& Ops )
{
  // Create subset of variables participating in a set of operations Ops
  t_Vars Vars;
  for( cit_Ops it = Ops.begin(); it != Ops.end(); ++it ){
    if( (*it)->plop ){
      const pt_idVar& idlop = (*it)->plop->id();
      if( idlop.first != FPVar<T>::AUXINT && idlop.first != FPVar<T>::AUXREAL )
        Vars.insert( (*it)->plop );
    }
    if( (*it)->prop ){
      const pt_idVar& idrop = (*it)->prop->id();
      if( idrop.first != FPVar<T>::AUXINT && idrop.first != FPVar<T>::AUXREAL )
        Vars.insert( (*it)->prop );
    }
    if( (*it)->pres && (*it)->type != FPOp<T>::EQ && (*it)->type != FPOp<T>::LE ){
      const pt_idVar& idres = (*it)->pres->id();
      if( idres.first != FPVar<T>::AUXINT && idres.first != FPVar<T>::AUXREAL )
        Vars.insert( (*it)->pres );
    }
  }
  return Vars;
}

template <typename T> inline typename FPRRLT<T>::t_Vars_Ops
FPRRLT<T>::submap_var
( const typename FPRRLT<T>::t_Ops& Ops )
{
  // Create subset of variables participating in operations Ops
  t_Vars_Ops Vars_Ops;
  for( cit_Ops it = Ops.begin(); it != Ops.end(); ++it ){
    if( (*it)->plop ){
      const pt_idVar& idlop = (*it)->plop->id();
      if( idlop.first != FPVar<T>::AUXINT && idlop.first != FPVar<T>::AUXREAL )
        Vars_Ops.insert( std::make_pair( (*it)->plop, *it ) );
    }
    if( (*it)->prop ){
      const pt_idVar& idrop = (*it)->prop->id();
      if( idrop.first != FPVar<T>::AUXINT && idrop.first != FPVar<T>::AUXREAL )
        Vars_Ops.insert( std::make_pair( (*it)->prop, *it ) );
    }
    if( (*it)->pres && (*it)->type != FPOp<T>::EQ && (*it)->type != FPOp<T>::LE ){
      const pt_idVar& idres = (*it)->pres->id();
      if( idres.first != FPVar<T>::AUXINT && idres.first != FPVar<T>::AUXREAL )
        Vars_Ops.insert( std::make_pair( (*it)->pres, *it ) );
    }
  }
  return Vars_Ops;
}

} // namespace mc

#include "mcop.hpp"

namespace mc
{

//! @brief Specialization of the mc::Op templated structure to allow usage of type mc::FPVar as template type in other MC++ classes, such as mc::Taylor
template <> template<typename T> struct Op< mc::FPVar<T> >
{
  typedef mc::FPVar<T> FPV;
  typedef mc::FPRelax<T> FPR;
  static FPV point( const double c ) { return FPV(c); }
  static FPV zeroone() { return FPV( mc::Op<T>::zeroone() ); }
  static void I(FPV& x, const FPV&y) { x = y; }
  static double l(const FPV& x) { return mc::Op<T>::l(x.num().I); }
  static double u(const FPV& x) { return mc::Op<T>::u(x.num().I); }
  static double abs (const FPV& x) { return mc::Op<T>::abs(x.num().I);  }
  static double mid (const FPV& x) { return mc::Op<T>::mid(x.num().I);  }
  static double diam(const FPV& x) { return mc::Op<T>::diam(x.num().I); }
  static FPV inv (const FPV& x) { return mc::inv(x);  }
  static FPV sqr (const FPV& x) { return mc::sqr(x);  }
  static FPV sqrt(const FPV& x) { return mc::sqrt(x); }
  static FPV log (const FPV& x) { return mc::log(x);  }
  static FPV xlog(const FPV& x) { return mc::xlog(x); }
  static FPV fabs(const FPV& x) { return mc::fabs(x); }
  static FPV exp (const FPV& x) { return mc::exp(x);  }
  static FPV sin (const FPV& x) { return mc::sin(x);  }
  static FPV cos (const FPV& x) { return mc::cos(x);  }
  static FPV tan (const FPV& x) { return mc::tan(x);  }
  static FPV asin(const FPV& x) { return mc::asin(x); }
  static FPV acos(const FPV& x) { return mc::acos(x); }
  static FPV atan(const FPV& x) { return mc::atan(x); }
  static FPV erf (const FPV& x) { return mc::erf(x);  }
  static FPV erfc(const FPV& x) { return mc::erfc(x); }
  static FPV min (const FPV& x, const FPV& y) { return mc::min(x,y);  }
  static FPV max (const FPV& x, const FPV& y) { return mc::max(x,y);  }
  static FPV arh (const FPV& x, const double k) { return mc::arh(x,k); }
  template <typename X, typename Y> static FPV pow(const X& x, const Y& y)
    { return mc::pow(x,y); }
  static FPV hull(const FPV& x, const FPV& y)
    { return mc::Op<T>::hull(x.num().I,y.num().I); }
  static bool inter(FPV& xIy, const FPV& x, const FPV& y)
    { try{ xIy = x^y; return true; }
      catch(typename FPR::Exceptions &eObj){ return false; } }
  static bool eq(const FPV& x, const FPV& y)
    { return mc::Op<T>::eq(x.num().I,y.num().I); }
  static bool ne(const FPV& x, const FPV& y)
    { return mc::Op<T>::ne(x.num().I,y.num().I); }
  static bool lt(const FPV& x, const FPV& y)
    { return mc::Op<T>::lt(x.num().I,y.num().I); }
  static bool le(const FPV& x, const FPV& y)
    { return mc::Op<T>::le(x.num().I,y.num().I); }
  static bool gt(const FPV& x, const FPV& y)
    { return mc::Op<T>::gt(x.num().I,y.num().I); }
  static bool ge(const FPV& x, const FPV& y)
    { return mc::Op<T>::ge(x.num().I,y.num().I); }
};

} // namespace mc

#endif
