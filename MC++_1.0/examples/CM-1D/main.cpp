#define TEST_ABS	// <-- select test function here
const int NTE = 11;	// <-- select Chebyshev expansion order here
const int NX = 500;	// <-- select X discretization here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef USE_PROFIL	// <-- specify to use PROFIL for interval arithmetic
#undef USE_FILIB	// <-- specify to use FILIB++ for interval arithmetic

////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

using namespace std;
using namespace mc;

////////////////////////////////////////////////////////////////////////

#if defined( TEST_POW )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  //return 1. - 5.*x - pow(x,2)/2. + pow(x,3)/3.;
  return x*x*x + x*x;
  //return 1. - 5.*x - x*x/2. + x*x*x/3. + x*x*x*x*x/9.;
}

#elif defined( TEST_ABS )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return sqrt(pow(x,2));
}

#elif defined( TEST_EXP )
const double XL   = -1.;	// <-- X range lower bound
const double XU   =  2.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return x*exp(-pow(x,2));
  //return x*exp(-x*x);
}

#elif defined( TEST_INV )
const double XL   = -3.;	// <-- X range lower bound
const double XU   = -1.;	// <-- X range upper bound
template <class T>
T myfunc
( const T&x )
{
  return 1./x;
}

#elif defined( TEST_TRIG )
const double XL   =  0.;	// <-- range lower bound
const double XU   =  PI/3.;	// <-- range upper bound
template <class T>
T myfunc
( const T&x )
{
  return tan(cos(x*atan(x)));
}

#elif defined( TEST_TRIG2 )
const double XL   = PI/6.;	// <-- X range lower bound
const double XU   = PI/3.;	// <-- X range upper bound
const double Xref = PI/4.;	// <-- X ref point for McCormick
template <class T>
T myfunc
( const T&x )
{
  return sin(pow(x,-3))*cos(sqrt(x));
}
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{

#ifdef SAVE_RESULTS
  ofstream res( "CM-1D.out", ios_base::out );
  res << std::scientific << std::setprecision(5) << std::right;
#endif

  try{ 

    // Define Chebyshev model environment
    CM mod( 1, NTE );

    // <-- set options here -->
    mod.options.BOUNDER_TYPE = CM::Options::LSB;

    // Define variable X, and evaluate Chebyshev model
    CV CVX( &mod, 0, I(XL,XU) );
    CV CVF = myfunc( CVX );
    std::cout << "\nChebyshev model of f(x):" << CVF;

    // Repeated calculations at grid points (for display)
    for( int iX=0; iX<NX; iX++ ){ 

      double DX = XL+iX*(XU-XL)/(NX-1.);
      double DF = myfunc( DX );
      CVX.set( &mod, 0, I(XL,XU) );
      CVF = myfunc( CVX );
      double PF = CVF.P( &DX );
      I IF = PF + CVF.R();

#ifdef SAVE_RESULTS
      res << std::setw(14) << DX << std::setw(14) << std::setw(14) << DF
          << std::setw(14) << PF << std::setw(14) << IF.l()  << setw(14) << IF.u()
          << std::setw(14) << CVF.B().l()  << std::setw(14) << CVF.B().u()
          << std::endl;
#endif
    }

  }
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( CM::Exceptions &eObj ){
    cerr << "Error " << eObj.ierr()
         << " in Chebyshev model computation:" << endl
	 << eObj.what() << endl
         << "Aborts." << endl;
    return eObj.ierr();
  }

#ifdef SAVE_RESULTS
  res.close();
#endif
  return 0;
}
