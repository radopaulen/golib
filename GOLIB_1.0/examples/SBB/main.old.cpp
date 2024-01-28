//#define USE_PROFIL
//#define USE_FILIB
//#define USE_CPLEX

#include <fstream>
#include <iomanip>
#include <cassert>

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
  #include "lprcplex.h"
  typedef mc::LPRCplex<I> LPRG;
#else
  #include "lprgurobi.h"
  typedef mc::LPRGurobi<I> LPRG;
#endif
typedef mc::LPVar<I> LPV;

#include "sbb.h"
#include "nlpipopt.h"
//#include "coin/IpIpoptApplication.hpp"

using namespace std;
using namespace mc;

#define LP_OUTPUT 1
#define EXAMPLE2

////////////////////////////////////////////////////////////////////////
template <typename T>
class TestProblem: public mc::NLP<T>
{
public:

// Example 1 [Ryoo & Sahinidis, Comp Chem Eng 1995]
#if defined( EXAMPLE1 )
  int NDIV;
  TestProblem(): mc::NLP<T>( 2, 1 )
    {
      NDIV = 2;
      // Default bounds
      this->Bounds( 0 ) = std::make_pair( 0., 6. );
      this->Bounds( 1 ) = std::make_pair( 0., 4. );
      // Default initial guess
      this->Guess( 0 ) = 1.;
      this->Guess( 1 ) = 2.;
    }
  std::pair<typename TestProblem<T>::OBJTYPE, T> Objective
    ( const T*p )
    {
      //return std::make_pair( TestProblem<T>::MIN, -(p[0]*p[1]) );
      return std::make_pair( TestProblem<T>::MIN, bilin(p[0],p[1],NDIV) );
    }
  std::pair<typename TestProblem<T>::CONSTRTYPE, T> Constraint
    ( const unsigned int ic, const T*p )
    {
      assert( ic < nc() );
      switch( ic ){
        //case 0: return std::make_pair( TestProblem<T>::LE, p[0]*p[1]-4 );
        case 0: return std::make_pair( TestProblem<T>::LE, bilin(p[0],p[1],NDIV)-4 );
      }
    }

// Example 7 [Ryoo & Sahinidis, Comp Chem Eng 1995]
#elif defined( EXAMPLE2 )
  int NDIV;
  TestProblem(): mc::NLP<T>( 10, 7 )
    {
      NDIV = 2;
      // Default Bounds
      this->Bounds( 0 ) = std::make_pair( 0., 300. );
      this->Bounds( 1 ) = std::make_pair( 0., 300. );
      this->Bounds( 2 ) = std::make_pair( 0., 100. );
      this->Bounds( 3 ) = std::make_pair( 0., 200. );
      this->Bounds( 4 ) = std::make_pair( 0., 100. );
      this->Bounds( 5 ) = std::make_pair( 0., 300. );
      this->Bounds( 6 ) = std::make_pair( 0., 100. );
      this->Bounds( 7 ) = std::make_pair( 0., 200. );
      this->Bounds( 8 ) = std::make_pair( 0., 200. );
      this->Bounds( 9 ) = std::make_pair( 1., 3. );
      // Default guess
      this->Guess( 0 ) = 1.;
      this->Guess( 1 ) = 1.;
      this->Guess( 2 ) = 1.;
      this->Guess( 3 ) = 1.;
      this->Guess( 4 ) = 1.;
      this->Guess( 5 ) = 1.;
      this->Guess( 6 ) = 1.;
      this->Guess( 7 ) = 1.;
      this->Guess( 8 ) = 1.;
      this->Guess( 9 ) = 1.;
    }
  std::pair<typename TestProblem<T>::OBJTYPE, T> Objective
    ( const T*p )
    {
      return std::make_pair( TestProblem<T>::MIN, -9*p[4]-15*p[8]+6*p[0]+16*p[1]+10*p[5] );
    }
  std::pair<typename TestProblem<T>::CONSTRTYPE, T> Constraint
    ( const unsigned int ic, const T*p )
    {
      assert( ic < nc() );
      switch( ic ){
        case 0: return std::make_pair( TestProblem<T>::EQ, p[0]+p[1]-p[2]-p[3] );
        case 1: return std::make_pair( TestProblem<T>::EQ, p[2]+p[6]-p[4] );
        case 2: return std::make_pair( TestProblem<T>::EQ, p[3]+p[7]-p[8] );
        case 3: return std::make_pair( TestProblem<T>::EQ, p[6]+p[7]-p[5] );
        case 4: return std::make_pair( TestProblem<T>::LE, bilin(p[9],p[2],NDIV)+2*p[6]-2.5*p[4] );
        case 5: return std::make_pair( TestProblem<T>::LE, bilin(p[9],p[3],NDIV)+2*p[7]-1.5*p[8] );
        case 6: return std::make_pair( TestProblem<T>::LE, 3*p[0]+p[1]-bilin(p[9],p[2],NDIV)-bilin(p[9],p[3],NDIV) );
      }
    }
#endif
};

////////////////////////////////////////////////////////////////////////
SBB<I>::STATUS UPB
( const unsigned int np, I*P, double*p, double&f )
////////////////////////////////////////////////////////////////////////
{
  // Construct NLP representation of TestProblem
  TestProblem<double>*NLO = new TestProblem<double>;
  for( unsigned int ip=0; ip<np; ip++ ){
    NLO->Bounds(ip) = std::make_pair( Op<I>::l(P[ip]), Op<I>::u(P[ip]) );
    NLO->Guess(ip)  = p[ip];
  }

  // Construct NLP structure and derivatives
  TestProblem<Structure>*NLOS = new TestProblem<Structure>;
  typedef fadbad::F< double > Fdouble;
  TestProblem<Fdouble>*NLOF   = new TestProblem<Fdouble>;
  //typedef fadbad::B< double > Bdouble;
  //TestProblem<Bdouble>*NLOB   = new TestProblem<Bdouble>;
  typedef fadbad::B< fadbad::F< double > > BFdouble;
  TestProblem<BFdouble>*NLOBF = new TestProblem<BFdouble>;
  

  // Solve NLP locally
  Ipopt::SmartPtr<NLPIPOPT> pNLPIPOPT = new NLPIPOPT( NLO, NLOS, NLOF, NLOBF );
  Ipopt::ApplicationReturnStatus status = pNLPIPOPT->solve( 3, 100, 1e-9 );
  //status = pNLPIPOPT->solve( 3, 100, 1e-9 );
  if( status == Ipopt::Solve_Succeeded ){
    f = pNLPIPOPT->get_objective();
    for( unsigned int ip=0; ip<np; ip++ )
      p[ip] = pNLPIPOPT->get_variable(ip);

    // Display NLP optimization results
    cout << "  Status:  " << status
         << endl
         << "NLP OPTIMIZATION SOLUTION: " << endl
         << "  F* = " << f << endl;
    for( unsigned int ip=0; ip<np; ip++ )
      cout << "  p*(" << ip << ") = " << p[ip]
//            << "  m*(" << ip << ") = " << m[ip]
           << endl;
  }

  delete NLO;
  delete NLOS;
  delete NLOF;
  //delete NLOB;
  delete NLOBF;

  return( status == Ipopt::Solve_Succeeded? SBB<I>::NORMAL:  SBB<I>::FAILURE );
}

////////////////////////////////////////////////////////////////////////
SBB<I>::STATUS UFEAS
( const unsigned int np, I*P, double*p )
////////////////////////////////////////////////////////////////////////
{
  // Construct NLP representation of TestProblem
  TestProblem<double>*NLO = new TestProblem<double>;
  for( unsigned int ip=0; ip<np; ip++ ){
    NLO->Bounds(ip) = std::make_pair( Op<I>::l(P[ip]), Op<I>::u(P[ip]) );
    NLO->Guess(ip)  = p[ip];
  }

  // Test feasibility
  const double atol = 1e-9;
  for( unsigned int ic=0; ic<NLO->nc(); ic++ ){
    std::pair< NLP<double>::CONSTRTYPE,double> g =  NLO->Constraint( ic, p );
    switch( g.first ){
    case NLP<double>::EQ:
      cout << "  G(" << ic << ") = " << g.second << " ?= 0" << endl;
      if( fabs(g.second) <= atol ) continue;
      return SBB<I>::INFEASIBLE;
    case NLP<double>::LE:
      cout << "  G(" << ic << ") = " << g.second << " ?<= 0" << endl;
      if( g.second <= atol ) continue;
      return SBB<I>::INFEASIBLE;
    case NLP<double>::GE:
      cout << "  G(" << ic << ") = " << g.second << " ?>= 0" << endl;
      if( g.second >= -atol ) continue;
      return SBB<I>::INFEASIBLE;
    }
  }

  delete NLO; 
  return SBB<I>::NORMAL;
}

////////////////////////////////////////////////////////////////////////
SBB<I>::STATUS LPB
( const unsigned int np, I*P, double*p, double&f )
////////////////////////////////////////////////////////////////////////
{  
  try{ 

    // Construct LP relaxation of TestProblem
    TestProblem<LPV> NLO;
    assert( np == NLO.np );
    LPRG LPG;

    LPV* PLP = new LPV[np];
    for( unsigned int ip=0; ip<NLO.np(); ip++ )
      PLP[ip] = LPV( &LPG, ip, P[ip] );

    std::pair<TestProblem<LPV>::OBJTYPE,LPV> F = NLO.Objective( PLP );
    switch( F.first ){
      case TestProblem<LPV>::MIN:
        LPG.set_objective( LPRG::MIN, F.second ); break;
      case TestProblem<LPV>::MAX:
        LPG.set_objective( LPRG::MAX, F.second ); break;
    }
    
    for( unsigned int ic=0; ic<NLO.nc(); ic++ ){
      std::pair<TestProblem<LPV>::CONSTRTYPE,LPV> F = NLO.Constraint( ic, PLP );
      switch( F.first ){
        case TestProblem<LPV>::EQ:
          LPG.add_constraint( F.second, LPRG::EQ, 0 ); break;
        case TestProblem<LPV>::LE:
          LPG.add_constraint( F.second, LPRG::LE, 0 ); break;
        case TestProblem<LPV>::GE:
          LPG.add_constraint( F.second, LPRG::GE, 0 ); break;
      }
    }
    
    delete [] PLP;
    
    // Solve LP relaxation
    cout << LPG;
    int status = LPG.solve( 0, LP_OUTPUT , "relax.lp" );

    // Gather LP relaxation results
    f = LPG.get_objective();
    for( unsigned int ip=0; ip<NLO.np(); ip++ ){
      p[ip] = LPG.get_variable(ip);
//       m[ip] = LPG.get_reduced_cost(ip);
    }

    // Display LP relaxation results
    cout << "  Status:  " << status
         << "  Runtime: " << LPG.get_runtime() << " CPUsec" << endl
         << endl
         << "LP RELAXATION SOLUTION: " << endl
         << "  F* = " << f << endl;
    for( unsigned int ip=0; ip<NLO.np(); ip++ )
      cout << "  p*(" << ip << ") = " << p[ip]
//            << "  m*(" << ip << ") = " << m[ip]
//            << "  Status: "  << LPG.get_basis_status(ip)
           << endl;
  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( mc::I_Excp& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension." << endl
         << "Execution aborted..." << endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( mc::LPExcp& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in LP relaxation." << endl
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

////////////////////////////////////////////////////////////////////////
SBB<I>::STATUS SUBPB
( const SBB<I>::TASK itask, const unsigned int np, I*P, double*p, double&f,
  const double UB )
////////////////////////////////////////////////////////////////////////
{
  switch( itask ){
  case SBB<I>::LBP:
    return LPB( np, P, p, f );
  
  case SBB<I>::REDUC:
    return SBB<I>::NORMAL;
  
  case SBB<I>::UBP:
    return UPB( np, P, p, f );
  
  case SBB<I>::FEAS:
    return UFEAS( np, P, p );
  }
}


////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  LPRG::options.NEWTON_USE      = true;
  LPRG::options.NEWTON_MAXIT    = 100;
  LPRG::options.NEWTON_TOL      = 1e-12;
  LPRG::options.SANDWICH_RTOL   = 1e-3;
  LPRG::options.SANDWICH_ATOL   = 1e-3;
  LPRG::options.SANDWICH_MAXCUT = 10;
  LPRG::options.SANDWICH_RULE   = mc::LPOpt::MAXERR;
  LPRG::options.DIVISION_RTOL   = 1e-5;
  LPRG::options.DIVISION_ATOL   = 1e-5;
  
  TestProblem<LPV> NLO;
  double f;
  double*p = new double[NLO.np()];
  I*P = new I[NLO.np()];
  for( unsigned int ip=0; ip<NLO.np(); ip++ ){
    std::pair<double,double> Pip = NLO.Bounds(ip);
    P[ip] = I(Pip.first, Pip.second);
    p[ip] = NLO.Guess(ip);
  }

  //UPB( NLO.np(), P, p, f );
  //LPB( NLO.np(), P, p, f );

  SBB<I> GO( SUBPB );
  GO.set_variables( NLO.np(), P );
  GO.solve( p );
  
  delete [] p;
  delete [] P;
  
//   SBB Problem;
//   const unsigned int NP = 2;
//   I P[NP] = { I(0,6), I(0,4) };
//   Problem.set_variables( NP, P );
//   Problem.options.BRANCHING_OPTION = mc::SBBOpt::BISECTION;
//   std::pair<double,double*> opt = Problem.solve();

  return 0;
}
