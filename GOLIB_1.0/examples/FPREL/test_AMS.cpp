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
    typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
    //typedef filib::interval<double> I;
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

#include "tmodel.h"
#include "mccormick.h"
typedef mc::McCormick<I> MC;
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;
typedef mc::TModel<MC> TMMC;
typedef mc::TVar<MC> TVMC;
typedef mc::TModel<FPVAR> TMFPVAR;


 double Ftmval(double X1,double X2, double X3, double X1ref,double X2ref, double X3ref){
       return 2.90000e+01+4.50000e+00*(X3-X3ref)+9.75000e+00*(X2-X2ref)+1.05000e+01*(X1-X1ref)
             +4.50000e+00*(X2-X2ref)*(X3-X3ref)+1.00000e+00*mc::sqr(X2-X2ref)+1.00000e+00*(X1-X1ref)*(X3-X3ref)
             +1.50000e+00*(X1-X1ref)*(X2-X2ref)+1.00000e+00*mc::sqr(X1-X1ref)+1.00000e+00*(X1-X1ref)*(X2-X2ref)*(X3-X3ref);   
   };
 double Frval(double X1,double X2, double X3){
    return X1*(X2+X2*X3)+mc::sqr(X1)+mc::sqr(X2)+X2;
   };

int main()
{  
 try{ 

  FPREL FP;
  FP.options.SANDWICH_ATOL = 1e-3;
  FP.options.SANDWICH_RTOL = 1e-3;
  FP.options.SANDWICH_MAXCUT = 0;
  FP.options.NEWTON_USE = true;
  FP.options.NEWTON_TOL = 1e-12;
  FP.options.NEWTON_MAXIT = 100;
  FP.options.SANDWICH_RULE = FPREL::Options::MAXERR;
  FP.options.RRLT_STRATEGY = FPREL::Options::PRIMRRLT;

//  I  X1_I( -3.,3. );
//  I  X2_I( -4.,4. );
//  I  X3_I( -2.,2. );
  I  X1_I( 3.,6. );
  I  X2_I( -2.,4. );
  I  X3_I( -1.,2. );
 
  double X1ref = mid(X1_I);
  double X2ref = mid(X2_I);
  double X3ref = mid(X3_I);

  FPVAR X1( &FP, 1, X1_I );
  FPVAR X2( &FP, 2, X2_I );
  FPVAR X3( &FP, 3, X3_I ); 

// Solve polyhedral relaxation w/ Taylor models

  unsigned int TM2FP[3] = { 1,2,3 };

  TM TMod( 3, 3 );
  TMod.options.BOUNDER_TYPE = TM::Options::NAIVE;
  TMod.options.PROPAGATE_BNDT = false;
  TMod.options.SCALE_VARIABLES = false;
  TMod.options.CENTER_REMAINDER = true;

  TV X1_TM( &TMod, 0, X1_I );
  TV X2_TM( &TMod, 1, X2_I );
  TV X3_TM( &TMod, 2, X3_I );

/////////////////////////////////////////////////////////////////////
//Case a
/////////////////////////////////////////////////////////////////////
   //Construct TM of f = X1*X2*X3 + X1*X2 + sqr(X1) + sqr(X2) + X2 
   TV X4_TM = (X1_TM)*(X2_TM)*(X3_TM)+(X1_TM)*(X2_TM)+sqr(X1_TM)+sqr(X2_TM)+X2_TM; 




  std::cout << "X4_TM\n" << X4_TM <<std::endl;
  
  TMFPVAR*TModFP = FP.TModelFP( &TMod, TM2FP );
  FPVAR FTM( &FP, X4_TM, TModFP );
  delete TModFP;

  std::cout <<"FTM " << FP;
  int dum;std::cin >> dum;

  //Solve TM relaxation using built-in TM2FP feature
  FP.set_objective( FPREL::MIN, FTM );

  FP.generate_cuts( true );
  std::cout <<"here\n";
  int status = FP.solve();
  std::cout << "\n  Case a:  ";
  std::cout << "\n  RELAXATION STATUS:  " << status
            << "\n  RELAXATION VALUE:   " << FP.get_objective()
            << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                          << FP.get_variable(2) << ", "
                                          << FP.get_variable(3) << ", "
            << "\n  CPU TIME [sec]:     " << FP.get_runtime()
	    << std::endl;
  std::cin >> dum;
/////////////////////////////////////////////////////////////////////
//Case b
/////////////////////////////////////////////////////////////////////
  FP.reset();
  X1.set( &FP, 1, X1_I );
  X2.set( &FP, 2, X2_I );
  X3.set( &FP, 3, X3_I ); 
  //Extract TM coeffs and hard-code the TM as a regular expression; rather than by built-in TM2FP feature
/*  a0   =  2.90000e+01        0  0  0    B0   = [  1.00000e+00 :  1.00000e+00 ]
    a1   =  4.50000e+00        0  0  1    B1   = [ -1.50000e+00 :  1.50000e+00 ]
    a2   =  9.75000e+00        0  1  0    B2   = [ -3.00000e+00 :  3.00000e+00 ]
    a3   =  1.05000e+01        1  0  0    B3   = [ -1.50000e+00 :  1.50000e+00 ]
    a5   =  4.50000e+00        0  1  1    B5   = [ -4.50000e+00 :  4.50000e+00 ]
    a6   =  1.00000e+00        0  2  0    B6   = [  0.00000e+00 :  9.00000e+00 ]
    a7   =  1.00000e+00        1  0  1    B7   = [ -2.25000e+00 :  2.25000e+00 ]
    a8   =  1.50000e+00        1  1  0    B8   = [ -4.50000e+00 :  4.50000e+00 ]
    a9   =  1.00000e+00        2  0  0    B9   = [  0.00000e+00 :  2.25000e+00 ]
    a15  =  1.00000e+00        1  1  1    B15  = [ -6.75000e+00 :  6.75000e+00 ]*/
  FPVAR Ftm = 2.90000e+01+4.50000e+00*(X3-X3ref)+9.75000e+00*(X2-X2ref)+1.05000e+01*(X1-X1ref)
             +4.50000e+00*(X2-X2ref)*(X3-X3ref)+1.00000e+00*sqr(X2-X2ref)+1.00000e+00*(X1-X1ref)*(X3-X3ref)
             +1.50000e+00*(X1-X1ref)*(X2-X2ref)+1.00000e+00*sqr(X1-X1ref)+1.00000e+00*(X1-X1ref)*(X2-X2ref)*(X3-X3ref);   

  std::cout <<"Ftm " << FP;
  std::cin >> dum;
  //Solve hard-coded TM expression
  FP.set_objective( FPREL::MIN, Ftm );
  FP.generate_cuts( true );

  status = FP.solve();
  std::cout << "\n  Case b:  ";
  std::cout << "\n  RELAXATION STATUS:  " << status
            << "\n  RELAXATION VALUE:   " << FP.get_objective()
            << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                          << FP.get_variable(2) << ", "
                                          << FP.get_variable(3) << ", "
            << "\n  CPU TIME [sec]:     " << FP.get_runtime()
	    << std::endl;
  std::cin >> dum;
/////////////////////////////////////////////////////////////////////
//Case c
/////////////////////////////////////////////////////////////////////
  FP.reset();
  X1.set( &FP, 1, X1_I );
  X2.set( &FP, 2, X2_I );
  X3.set( &FP, 3, X3_I ); 
  //Reformulate TM as follows
  FPVAR Fr = X1*(X2+X2*X3)+sqr(X1)+sqr(X2)+X2; 
  std::cout <<"Fr " << FP;


  std::cin >> dum;

  //Solve yet another reformulation of same TM expression
  std::cout << "\n  Case c:  ";
  FP.set_objective( FPREL::MIN, Fr );
  FP.generate_cuts( true );

  status = FP.solve();
  std::cout << "\n  RELAXATION STATUS:  " << status
            << "\n  RELAXATION VALUE:   " << FP.get_objective()
            << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
                                          << FP.get_variable(2) << ", "
                                          << FP.get_variable(3) << ", "
            << "\n  CPU TIME [sec]:     " << FP.get_runtime()
	    << std::endl;
//Are the three expressions equivalent?
  std::cout<<"\nAre the three expressions equivalent?\n";
  double xval[3] = {6.,1,1.};
  std::cout<<"FTM(6.,1,1.)  "<<X4_TM.P(xval) <<std::endl;
  std::cout<<"Ftm(6.,1,1.)  "<<Ftmval(6.,1,1.,X1ref,X2ref,X3ref) <<std::endl;
  std::cout<<"Fr(6.,1,1.)  "<< Frval(6.,1,1.)  <<std::endl;

 }
#ifndef USE_PROFIL
#ifndef USE_FILIB
 catch( I::Exceptions& eObj ){
   std::cerr << "Error " << eObj.ierr()
        << " in natural interval extension." << std::endl
        << "Execution aborted..." << std::endl;
   return eObj.ierr();
 }
#endif
#endif
 catch( TM::Exceptions& eObj ){
   std::cerr << "Error " << eObj.ierr()
        << " in Taylor model." << std::endl
        << "Execution aborted..." << std::endl;
   return eObj.ierr();
 }
 catch( TMMC::Exceptions& eObj ){
   std::cerr << "Error " << eObj.ierr()
        << " in McCormick-Taylor model." << std::endl
        << "Execution aborted..." << std::endl;
   return eObj.ierr();
 }
 catch( FPREL::Exceptions &eObj ){
   std::cerr << "Error " << eObj.ierr()
        << " in polyhedral relaxation." << std::endl
        << "Execution aborted..." << std::endl;
   return eObj.ierr();
 }
#ifdef USE_CPLEX
 catch (IloException& ex) {
   std::cerr << "Error: " << ex << std::endl;
  }
#else
 catch(GRBException& ex) {
   std::cerr << "Error code = " << ex.getErrorCode() << std::endl;
   std::cerr << ex.getMessage() << std::endl;
 }
#endif
 catch(...) {
   std::cout << "Error during optimization" << std::endl;
 }

 return 0;
}
