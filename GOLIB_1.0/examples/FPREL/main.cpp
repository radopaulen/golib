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

int dum; 
int status;
 
int main()
{  
 try{ 

  FPREL FP;
  FP.options.SANDWICH_ATOL = 1e-3;
  FP.options.SANDWICH_RTOL = 1e-3;
  FP.options.SANDWICH_MAXCUT = 5;
  FP.options.NEWTON_USE = true;
  FP.options.NEWTON_TOL = 1e-12;
  FP.options.NEWTON_MAXIT = 100;
  FP.options.SANDWICH_RULE = FPREL::Options::MAXERR;
  FP.options.RRLT_STRATEGY = FPREL::Options::PRIMRRLT;
  FP.options.LP_FEASTOL  = 1e-9;
  FP.options.LP_OPTIMTOL = 1e-9;
  FP.options.SOLVER_DISPLAY = false;


//   // Example TM
//   FP.reset();
//   FP.options.RRLT_STRATEGY = FPREL::Options::ALLRRLT;
//   FP.options.RRLT_DISPLAY  = 1;
// 
//   FPVAR X1( &FP, 1, I(0,1) );
//   X1.breakpts().insert(0.5);
//   std::set<double>::iterator it = X1.breakpts().begin();
//   for( ; it != X1.breakpts().end(); ++it )
//     std::cout << *it << " ";
//   std::cout << std::endl;
//   { int dum; std::cin >> dum; }
//   FPVAR X2( &FP, 2, I(1,2) );
// 
//   TM TMod( 2, 3 );
//   TMod.options.BOUNDER_TYPE = TM::Options::EIGEN;
//   TV X1_TM( &TMod, 0, X1.I() );
//   TV X2_TM( &TMod, 1, X2.I() );
//   unsigned int TM2FP[2] = { 1, 2 };
//   TV X3_TM = X1_TM*X2_TM*(X1_TM*(exp(X1_TM)-exp(-X1_TM))
//              -X2_TM*(exp(X2_TM)-exp(-X2_TM)));
//   TMFPVAR*TModFP = FP.TModelFP( &TMod, TM2FP );
//   FPVAR X4( &FP, X3_TM, TModFP );
//   FP.set_objective( FPREL::MIN, X1 );
//   //FP.add_constraint( X1, FPREL::LE, 1 );
//   delete TModFP;
// 
//   FP.generate_reduction_constraints();
//   FP.generate_cuts( true );
//   std::cout << FP;
//   std::cin >> dum;  


  // Example Splitter
  FP.reset();
  FP.options.BILINEAR_SUBDIV  = 1;
  FP.options.BILINEAR_RULE    = FPREL::Options::SEPARSOS;
  FP.options.RRLT_STRATEGY    = FPREL::Options::ALLRRLT;
  FP.options.RRLT_DISPLAY     = 2;
  FP.options.BREAKPOINTS_USE  = true;
  FP.options.BREAKPOINTS_ATOL = 1e-2;

  FPVAR F1( &FP, 1, I(0,100) );
  FPVAR F2( &FP, 2, I(0,100) );
  FPVAR F ( &FP, 3, I(0,200) );
  FPVAR x1( &FP, 4, I(0,1) );
  x1.add_breakpt(0.5);
  x1.add_breakpt(0.489);
  x1.add_breakpt(0.99);
  std::set<double>::iterator it = x1.breakpts().begin();
  for( ; it != x1.breakpts().end(); ++it )
    std::cout << *it << " ";
  std::cout << std::endl;
  { int dum; std::cin >> dum; }
  FPVAR x2( &FP, 5, I(0,1) );

  FP.set_objective( FPREL::MIN, F );
  //FP.add_constraint( F1+F2, FPREL::EQ, F );
  FP.add_constraint( x1*F,  FPREL::EQ, F1 );
  FP.add_constraint( x2*F,  FPREL::EQ, F2 );
  FP.add_constraint( x1+x2, FPREL::EQ, 1 );
  FP.add_constraint( F1, FPREL::LE, 50 );
  FP.add_constraint( F2, FPREL::LE, 50 );
  std::cout << FP;

  //FP.generate_reduction_constraints();
  FP.generate_cuts( true );
  std::cout << FP;
// 
//   status = FP.solve();
//   std::cout << "\n  RELAXATION STATUS:  " << status
//             << "\n  RELAXATION VALUE:   " << FP.get_objective()
//             << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
//                                           << FP.get_variable(2) << ", "
//                                           << FP.get_variable(3)
//             << "\n  CPU TIME [sec]:     " << FP.get_runtime()
// 	    << std::endl;
//   std::cin >> dum;  
// 
// 
//   // Solve polyhedral relaxation w/o constraint propagation
//   FP.reset();
//   FPVAR X1( &FP, 1, I(3.,6.) );
//   FPVAR X2( &FP, 2, I(0.,4.) );
//   FPVAR Y( &FP, 3 );
// 
//   FP.set_objective( FPREL::MAX, X1+exp(1-X2)-Y );
//   FP.add_constraint( X1*exp(1-X2)-Y, FPREL::LE, 4. );
//   std::cout << FP;
// 
//   FP.generate_cuts();
//   std::cout << FP;
// 
//   status = FP.solve();
//   std::cout << "\n  RELAXATION STATUS:  " << status
//             << "\n  RELAXATION VALUE:   " << FP.get_objective()
//             << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
//                                           << FP.get_variable(2) << ", "
//                                           << FP.get_variable(3)
//             << "\n  CPU TIME [sec]:     " << FP.get_runtime()
// 	    << std::endl;
// 
//   
//   std::cin >> dum;  
// 
//   // Solve polyhedral relaxation w/ constraint propagation
//   FP.propagate_constraints();
//   std::cout << FP;
// 
//   FP.generate_cuts( true );
//   std::cout << FP;
// 
//   status = FP.solve();
//   std::cout << "\n  RELAXATION STATUS:  " << status
//             << "\n  RELAXATION VALUE:   " << FP.get_objective()
//             << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
//                                           << FP.get_variable(2) << ", "
//                                           << FP.get_variable(3)
//             << "\n  CPU TIME [sec]:     " << FP.get_runtime()
// 	    << std::endl;
// 
//   std::cin >> dum;  
// 
//   // Solve polyhedral relaxation w/ McCormick and Taylor models
//   FP.reset();
//   X1.set( &FP, 1, I(3.,6.) );
//   X2.set( &FP, 2, I(0.,4.) );
// 
//   double X2_Ref[1] = { 2. };
//   unsigned int MC2FP[1] = { 2 };
// 
//   I  X2_I( 0.,4. );
//   MC X2_MC( X2_I, X2_Ref[0] );
//   X2_MC.sub( 1, 0 );
//   MC X3_MC = exp(1-X2_MC);
//   FPVAR X3( &FP, X3_MC, X2_Ref, MC2FP );
// 
//   TM TMod( 1, 4 );
//   TV X2_TM( &TMod, 0, X2_I );
//   TV X4_TM = exp( 1 - X2_TM );
//   TMFPVAR*TModFP = FP.TModelFP( &TMod, MC2FP );
//   FPVAR X4( &FP, X4_TM, TModFP );
//   //delete TModFP;
// 
//   TMMC MCTMod( 1, 4 );
//   TVMC X2_TMMC( &MCTMod, 0, X2_MC );
//   TVMC X5_TMMC = exp( 1 - X2_TMMC );
//   //TMFPVAR*TModFP = FP.TModelFP( &MCTMod, MC2FP );
//   FPVAR X5( &FP, X5_TMMC, TModFP, X2_Ref, MC2FP );
//   delete TModFP;
// 
//   FP.set_objective( FPREL::MAX, X1+X3 );
//   FP.add_constraint( X1*X3, FPREL::LE, 4. );
//   FP.add_constraint( X1*X4, FPREL::LE, 4. );
//   FP.add_constraint( X1*X5, FPREL::LE, 4. );
//   std::cout << FP;
// 
// //  FP.propagate_constraints();
// //  std::cout << FP;
// 
//   FP.generate_cuts( true );
//   std::cout << FP;
// 
//   status = FP.solve();
//   std::cout << "\n  RELAXATION STATUS:  " << status
//             << "\n  RELAXATION VALUE:   " << FP.get_objective()
//             << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1) << ", "
//                                           << FP.get_variable(2)
//             << "\n  CPU TIME [sec]:     " << FP.get_runtime()
// 	    << std::endl;
// 
//   std::cin >> dum;  
// 
//   // Example 7 [Ryoo & Sahinidis, Comp Chem Eng 1995]
//   FP.reset();
//   FP.options.BILINEAR_SUBDIV = 1;
//   FP.options.BILINEAR_RULE   = FPREL::Options::SEPARSOS;
//   FP.options.RRLT_STRATEGY   = FPREL::Options::PRIMRRLT;
//   FP.options.RRLT_DISPLAY    = 1;
// 
//   const unsigned int NP = 10;
//   FPVAR Y1( &FP, 1, I(0,300) );
//   FPVAR Y2( &FP, 2, I(0,300) );
//   FPVAR Y3( &FP, 3, I(0,100) );
//   FPVAR Y4( &FP, 4, I(0,200) );
//   FPVAR Y5( &FP, 5, I(0,100) );
//   FPVAR Y6( &FP, 6, I(0,300) );
//   FPVAR Y7( &FP, 7, I(0,100) );
//   FPVAR Y8( &FP, 8, I(0,200) );
//   FPVAR Y9( &FP, 9, I(0,200) );
//   FPVAR Y10( &FP, 10, I(1,3) );
// 
//   FP.set_objective( FPREL::MIN, -9*Y5-15*Y9+6*Y1+16*Y2+10*Y6 );
//   //FP.add_constraint( Y1+Y2, FPREL::EQ, Y3+Y4 );
//   FP.add_constraint( Y3+Y7, FPREL::EQ, Y5 );
//   FP.add_constraint( Y4+Y8, FPREL::EQ, Y9 );
//   FP.add_constraint( Y7+Y8, FPREL::EQ, Y6 );
//   FP.add_constraint( Y10*Y3+2*Y7, FPREL::LE, 2.5*Y5 );
//   FP.add_constraint( Y10*Y4+2*Y8, FPREL::LE, 1.5*Y9 );
//   FP.add_constraint( 3*Y1+Y2, FPREL::LE, Y10*(Y3+Y4) );
//   //FP.add_constraint( bilin(Y10,Y3)+2*Y7, FPREL::LE, 2.5*Y5 );
//   //FP.add_constraint( bilin(Y10,Y4)+2*Y8, FPREL::LE, 1.5*Y9 );
//   //FP.add_constraint( 3*Y1+Y2, FPREL::LE, bilin(Y10,Y3)+bilin(Y10,Y4) );
// 
//   FP.generate_reduction_constraints();
//   FP.generate_cuts( true );
//   std::cout << FP;
//   //return 0;
// 
//   status = FP.solve();
//   std::cout << "\n  RELAXATION STATUS:  " << status
//             << "\n  RELAXATION VALUE:   " << FP.get_objective()
//             << "\n  RELAXATION OPTIMUM: " << FP.get_variable(1);
//   for( unsigned int ip=1; ip<NP; ip++ )
//     std::cout << "\n                      " << FP.get_variable(ip+1);
//   std::cout << "\n  CPU TIME [sec]:     " << FP.get_runtime()
// 	    << std::endl;
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
