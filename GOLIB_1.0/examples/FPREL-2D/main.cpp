//#define USE_PROFIL
//#define USE_FILIB
//#define USE_CPLEX

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

using namespace std;
using namespace mc;

#define TEST_FUNCTION func_bilin	// <-- select test function here
const int nX = 40, nY = 40;		// <-- select discretization here
const double Xl=-1, Xu=1;		// <-- bounds on variable X here
const double Yl=0.5, Yu=2;		// <-- bounds on variable Y here
#define SAVE_RESULTS
#define LP_OUTPUT 1

///////////////////////////////////////////////////////////
template< class T >
T func_bilin
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  //T f = bilin(x,y);
  T f = x/y;
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_sqrt_fabs
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = sqrt(fabs(x-y));
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_step
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = fstep(x)*sqr(y) + bstep(x)*pow(y,3);
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_exp
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = x*y*(x*(exp(x)-exp(-x))-y*(exp(y)-exp(-y)));
  //T f = pow(x*exp(fabs(x)/y),3);
  //T f = x*exp(x+pow(y,2));
  //T f = x*exp(x+pow(y,2)) - pow(y,2);
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_trig
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y);
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_pow_inv
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = -1./(pow(x-4.,2)+pow(y-4.,2)+0.1)
        -1./(pow(x-1.,2)+pow(y-1.,2)+0.2)
        -1./(pow(x-8.,2)+pow(y-8.,2)+0.2);
//   T f = +1./(pow(x-1.,3)+pow(y-1.,3)+0.1)
//         -1./(pow(x-2.,2)+pow(y-3.,4)+0.2)
//         +1./(pow(x-3.,3)+pow(y-2.,1)+0.2);
  return f;
}

///////////////////////////////////////////////////////////
template< class T >
T func_max
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f[3] = { -pow(x+y-2.,2), -pow(x+y,2), -pow(x+y+2.,2) };
  return max( 3, f );
}

///////////////////////////////////////////////////////////
template< class T >
T func_norm
( const T&x, const T&y )
///////////////////////////////////////////////////////////
{
  T f = sqrt(pow(x,2)+pow(y,2));
  return f;
}

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

#ifdef SAVE_RESULTS
  ofstream res( "FPREL-2D.out", ios_base::out );
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
    FP.options.SANDWICH_MAXCUT = 40;
    FP.options.SANDWICH_RULE   = FPREL::Options::MAXERR;
    FP.options.DIVISION_RTOL   = 1e-5;
    FP.options.DIVISION_ATOL   = 1e-5;
    FP.options.BILINEAR_RULE   = FPREL::Options::SEPARABLE;
    FP.options.BILINEAR_SUBDIV = 4;

    FPVAR Xrel( &FP, 0, I(Xl,Xu) );
    FPVAR Yrel( &FP, 1, I(Yl,Yu) );
    FPVAR Zrel = TEST_FUNCTION( Xrel, Yrel );

    FP.set_objective( FPREL::MIN, Zrel );
    cout << FP;
    FP.generate_cuts();
    cout << FP;

    // Solve LP repeatedly at grid point
    FPCTR CX1 = 0, CX2 = 0;
    for( int iX=0; iX<nX; iX++ ){ 

     // Set variable X1 at grid point
     double Xval = Xl+iX*(Xu-Xl)/(nX-1.);
     CX1 = FP.add_constraint( Xrel, FPREL::EQ, Xval );
     FP.generate_cuts( CX1 );

     for( int iY=0; iY<nY; iY++ ){ 

      // Set variable X2 at grid point
      double Yval = Yl+iY*(Yu-Yl)/(nY-1.);
      CX2 = FP.add_constraint( Yrel, FPREL::EQ, Yval );
      FP.generate_cuts( CX2 );

      // Original function value
      double Zval = TEST_FUNCTION( Xval, Yval );

      // Convex bound
      FP.set_objective( FPREL::MIN, Zrel );
      FP.solve();// 0, LP_OUTPUT , "relax.lp" );
      double Zcv = FP.get_objective();

      // Concave bound
      FP.set_objective( FPREL::MAX, Zrel );
      FP.solve();// 0, LP_OUTPUT , "relax.lp" );
      double Zcc = FP.get_objective();

#ifdef SAVE_RESULTS
      res << setw(14) << Xval << setw(14) << Yval << setw(14) << Zval
          << setw(14) << Op<I>::l(Zrel.I()) << setw(14) << Op<I>::u(Zrel.I())
          << setw(14) << Zcv << setw(14) << Zcc
          << endl;
#endif

      // Remove extra constraint for X2
      FP.remove_constraint( CX2 );
     }

#ifdef SAVE_RESULTS
     res << endl;
#endif      
     // Remove extra constraint for X1
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
    cout << "Error during optimization" << endl;
  }

  return 0;
}
