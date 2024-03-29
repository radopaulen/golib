//#define USE_PROFIL
//#define USE_FILIB
#define USE_CPLEX

#include <fstream>
#include <iomanip>

#ifdef USE_PROFIL
  #include "mcprofil.h"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.h"
    //typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
    typedef filib::interval<double> I;
  #else
    #include "interval.h"
    typedef mc::Interval I;
  #endif
#endif

#ifdef USE_CPLEX
  #include "fprcplex.h"
  typedef mc::FPRCplex<I> FPREL;
#else
  #include "fprgurobi.h"
  typedef mc::FPRGurobi<I> FPREL;
#endif
typedef mc::FPVar<I> FPVAR;
typedef mc::FPOp<I>* FPCTR;

#include "tmodel.h"
#include "mccormick.h"
typedef mc::McCormick<I> MC;
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;
typedef mc::TModel<MC> TMMC;
typedef mc::TVar<MC> TVMC;
typedef mc::TModel<FPVAR> TMFPVAR;

using namespace std;
//using namespace mc;

// AVAILABLE 1-D TEST FUNCTIONS:
//  - func_odd		range [0.1:0.9] 
//  - func_erf		range [-1:1] 	
//  - func_cubic	range [-e:e] 	
//  - func_sin		range [PI/6:PI/3]
//  - func_min		range [-1:1]	
//  - func_ltcond	range [-1:1]	
//  - func_gtcond	range [-1:1]	
//  - func_exp		range [-1:1]	
//  - func_asin		range [-0.9:0.9]
#define TEST_FUNCTION func_exp	// <-- select test function here
const int nX = 500;		// <-- select discretization here
const double Xl=0., Xu=1.;	// <-- bounds on variable X here
const double X0=0.5;			// <-- reference point
//const double Xl=mc::PI/6., Xu=mc::PI/3.;	// <-- bounds on variable X here
//const double X0=mc::PI/4.;			// <-- reference point
#define SAVE_RESULTS
#define LP_OUTPUT 1

///////////////////////////////////////////////////////////
template <class T>
T func_exp
( const T&x )
///////////////////////////////////////////////////////////
{
  using mc::sqr;
  //T f = x*2/x;
  //T f = fabs(pow(x,x));
  //T f = pow(x,3);
  //T f = exp(x);
  //T f = exp(log(x));
  //T f = x*2.*exp(-fabs(pow(x*2.,3))+3);
  T f = x*exp(-sqr(x));
  //T f = x*exp(log(x+1));
  //T f = sqrt(fabs(x));
  return f;
}

///////////////////////////////////////////////////////////
template <class T>
T func_invtrig
( const T&x )
///////////////////////////////////////////////////////////
{
  //T f = acos(x);
  T f = atan(x);
  return f;
}

///////////////////////////////////////////////////////////
template <class T>
T func_trig
( const T&x )
///////////////////////////////////////////////////////////
{
  T f = cos(pow(x,2))*sin(pow(x,-3));
  //T f = sin(pow(x,-3));
  //T f = tan(x);
  //T f = erfc(x);
  return f;
}

///////////////////////////////////////////////////////////
template <class T>
T func_nonsmooth
( const T&x )
///////////////////////////////////////////////////////////
{
  using mc::sqr;
  using mc::fstep;
  using mc::bstep;
  //T f = min(pow(x-.5,3),-x+1.);
  //T f = min(pow(x-.5,2),0.5);//pow(x+0.5,2));
  //T m[3] = { pow(x-.5,2), pow(x,2), pow(x+0.5,2) };
  //T f = min(3,m);
  //T f = min(x,-x);
  //T f = fabs(x-(-x));
  //T f = bstep( sqr(x)-0.25 ) * pow(x-0.5,3)
  //    + fstep( sqr(x)-0.25 ) * pow(-x,3);
  T f = fstep( sqr(x)-0.25 ) * (sqr(x)-0.8)
      + bstep( sqr(x)-0.25 ) * (0.1-sqr(x));
  return f;
}

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

#ifdef SAVE_RESULTS
  ofstream res( "FPREL-1D.out", ios_base::out );
  res << scientific << setprecision(5) << right;
#endif

  try{ 
    // Construct factorable program
    FPREL FP;
    FP.options.NEWTON_USE      = true;
    FP.options.NEWTON_MAXIT    = 100;
    FP.options.NEWTON_TOL      = 1e-12;
    FP.options.SANDWICH_RTOL   = 1e-5;
    FP.options.SANDWICH_ATOL   = 1e-5;
    FP.options.SANDWICH_MAXCUT = 10;
    FP.options.SANDWICH_RULE   = FPREL::Options::MAXERR;
 
    FPVAR Xrel( &FP, 0, I(Xl,Xu) );
    double Xrel_Ref[1] = { X0 };
    unsigned int MC2FP[1] = { 0 };

    //TM TMod( 1, 4 );
    //TMod.options.BOUNDER_TYPE = TM::Options::HYBRID;
    //TMod.options.PROPAGATE_BNDT = false;
    //TMod.options.SCALE_VARIABLES = false;
    //TV Xrel_TM( &TMod, 0, I(Xl,Xu) );
    //TV Zrel_TM = TEST_FUNCTION( Xrel_TM );
    //TMFPVAR*TModFP = FP.TModelFP( &TMod, MC2FP );
    //FPVAR Zrel( &FP, Zrel_TM, TModFP );

    TMMC TMod( 1, 4 );
    TMod.options.BOUNDER_TYPE = TM::Options::HYBRID;
    TMod.options.PROPAGATE_BNDT = false;
    TMod.options.SCALE_VARIABLES = false;
    TVMC Xrel_TMMC( &TMod, 0, MC(I(Xl,Xu),X0).sub(1,0) );
    TVMC Zrel_TMMC = TEST_FUNCTION( Xrel_TMMC );
    cout << Zrel_TMMC;
    TMFPVAR*TModFP = FP.TModelFP( &TMod, MC2FP );
    FPVAR Zrel( &FP, Zrel_TMMC, TModFP, Xrel_Ref, MC2FP );

    //FPVAR Zrel = TEST_FUNCTION( Xrel );

    FP.set_objective( FPREL::MIN, Zrel );
    cout << FP;
    FP.generate_cuts();
    cout << FP;

    // Solve LP repeatedly at grid point
    FPCTR CX1 = 0;
    for( int iX=0; iX<nX; iX++ ){ 

      // Set variable at grid point
      double Xval = Xl+iX*(Xu-Xl)/(nX-1.);
      CX1 = FP.add_constraint( Xrel, FPREL::EQ, Xval );
      FP.generate_cuts( CX1 );

      // Original function value
      double Zval = TEST_FUNCTION( Xval );

      // Polyhedral lower bound
      FP.set_objective( FPREL::MIN, Zrel );
      FP.solve(); // 0, LP_OUTPUT , "relax.lp" );
      double Zcv = FP.get_objective();

      // Polyhedral upper bound
      FP.set_objective( FPREL::MAX, Zrel );
      FP.solve();// 0, LP_OUTPUT , "relax.lp" );
      double Zcc = FP.get_objective();

#ifdef SAVE_RESULTS
      res << setw(14) << Xval << setw(14) << Zval
          << setw(14) << mc::Op<I>::l(Zrel.I())
          << setw(14) << mc::Op<I>::u(Zrel.I())
          << setw(14) << Zcv << setw(14) << Zcc
          << endl;
#endif

      // Remove extra constraint
      FP.remove_constraint( CX1 );
    }   
  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension." << endl
         << "Execution aborted..." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( FPREL::Exceptions& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in factorable program." << endl
         << "Execution aborted..." << endl;
    return eObj.ierr();
  }
#ifdef USE_CPLEX
  catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
   }
#else
  catch(GRBException& ex) {
    cerr << "Error code = " << ex.getErrorCode() << endl;
    cerr << ex.getMessage() << endl;
  }
#endif
  catch(...) {
    cout << "Unidentified execution error" << endl;
  }

  return 0;
}
