//#define USE_PROFIL
//#define USE_FILIB
#define USE_CPLEX

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

using namespace std;
using namespace mc;

#define NOOUTPUT
#define EXAMPLE2

////////////////////////////////////////////////////////////////////////
template <typename U>
class TestProblem: public mc::NLPPB, public mc::NLPOBJ<U>, public mc::NLPCTR<U>
{
public:

// Example 1 [Ryoo & Sahinidis, Comp Chem Eng 1995]
#if defined( EXAMPLE1 )
  int NDIV;
  TestProblem(): mc::NLPPB( 2, 1 )
    {
      NDIV = 2;
      // Default bounds
      this->Bounds( 0 ) = std::make_pair( 0., 6. );
      this->Bounds( 1 ) = std::make_pair( 0., 4. );
      // Default initial guess
      this->Guess( 0 ) = 3.;
      this->Guess( 1 ) = 2.;
    }
  std::pair<typename mc::NLPOBJ<U>::TYPE,U> operator()
    ( const U*p )
    {
      //return std::make_pair( mc::NLPOBJ<U>::MIN, -(p[0]*p[1]) );
      return std::make_pair( mc::NLPOBJ<U>::MIN, -(p[0]+p[1]) );
    }
  std::pair<typename mc::NLPCTR<U>::TYPE,U> operator()
    ( const unsigned int ic, const U*p )
    {
      assert( ic < nc() );
      switch( ic ){
        case 0: return std::make_pair( mc::NLPCTR<U>::LE, p[0]*p[1]-4 );
        //case 0: return std::make_pair( mc::NLPCTR<U>::LE, bilin(p[0],p[1],NDIV)-4 );
      }
    }

// Example 7 [Ryoo & Sahinidis, Comp Chem Eng 1995]
#elif defined( EXAMPLE2 )
  int NDIV;
  TestProblem(): mc::NLPPB( 10, 7 )
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
  std::pair<typename mc::NLPOBJ<U>::TYPE,U> operator()
    ( const U*p )
    {
      return std::make_pair( mc::NLPOBJ<U>::MIN, -9*p[4]-15*p[8]+6*p[0]+16*p[1]+10*p[5] );
      //return std::make_pair( mc::NLPOBJ<U>::MAX, 9*p[4]+15*p[8]-6*p[0]-16*p[1]-10*p[5] );
    }
  std::pair<typename mc::NLPCTR<U>::TYPE,U> operator()
    ( const unsigned int ic, const U*p )
    {
      assert( ic < nc() );
      switch( ic ){
        case 0: return std::make_pair( mc::NLPCTR<U>::EQ, p[0]+p[1]-p[2]-p[3] );
        case 1: return std::make_pair( mc::NLPCTR<U>::EQ, p[2]+p[6]-p[4] );
        case 2: return std::make_pair( mc::NLPCTR<U>::EQ, p[3]+p[7]-p[8] );
        case 3: return std::make_pair( mc::NLPCTR<U>::EQ, p[6]+p[7]-p[5] );
        case 4: return std::make_pair( mc::NLPCTR<U>::LE, bilin(p[9],p[2],NDIV)+2*p[6]-2.5*p[4] );
        case 5: return std::make_pair( mc::NLPCTR<U>::LE, bilin(p[9],p[3],NDIV)+2*p[7]-1.5*p[8] );
        case 6: return std::make_pair( mc::NLPCTR<U>::LE, 3*p[0]+p[1]-bilin(p[2]+p[3],p[9],NDIV) );
        //case 6: return std::make_pair( mc::NLPCTR<U>::LE, 3*p[0]+p[1]-bilin(p[9],p[2],NDIV)-bilin(p[9],p[3],NDIV) );
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
  TestProblem<double> NLO;
  //typedef std::pair<NLP::OBJTYPE,double> (NLO::*OBJ)   ( const double*p );
  //NLPIPOPT::OBJ pf = &TestProblem<double>::Objective;
  for( unsigned int ip=0; ip<np; ip++ ){
    NLO.Bounds(ip) = std::make_pair( Op<I>::l(P[ip]), Op<I>::u(P[ip]) );
    NLO.Guess(ip)  = p[ip];
  }

  // Construct NLP structure and derivatives
  TestProblem<Structure> NLOS;
  typedef fadbad::F< double > Fdouble;
  TestProblem<Fdouble> NLOF;
  //typedef fadbad::B< double > Bdouble;
  //TestProblem<Bdouble> NLOB;
  typedef fadbad::B< fadbad::F< double > > BFdouble;
  TestProblem<BFdouble> NLOBF;

  // Solve NLP locally
  Ipopt::SmartPtr<NLPIPOPT> pNLPIPOPT = new NLPIPOPT( NLO.np(), NLO.nc(), NLO.Bounds() );
  pNLPIPOPT->objective( NLO, NLOS, NLOF, NLOBF );
  pNLPIPOPT->constraint( NLO, NLOS, NLOF, NLOBF );
#ifndef NOOUTPUT
  Ipopt::ApplicationReturnStatus status = pNLPIPOPT->solve( NLO.Guess(), 5, 100, 1e-9 );
#else
  Ipopt::ApplicationReturnStatus status = pNLPIPOPT->solve( NLO.Guess(), -1, 100, 1e-9 );
#endif
  if( status == Ipopt::Solve_Succeeded ){
    f = pNLPIPOPT->get_objective();
    for( unsigned int ip=0; ip<np; ip++ )
      p[ip] = pNLPIPOPT->get_variable(ip);

    // Display NLP optimization results
#ifndef NOOUTPUT
    cout << "  Status:  " << status
         << endl
         << "NLP OPTIMIZATION SOLUTION: " << endl
         << "  F* = " << f << endl;
    for( unsigned int ip=0; ip<np; ip++ )
      cout << "  p*(" << ip << ") = " << p[ip]
//            << "  m*(" << ip << ") = " << m[ip]
           << endl;
  int tmp;
  cin >> tmp;
#endif
  }

  return( status == Ipopt::Solve_Succeeded? SBB<I>::NORMAL:  SBB<I>::FAILURE );
}

////////////////////////////////////////////////////////////////////////
SBB<I>::STATUS UFEAS
( const unsigned int np, I*P, double*p, double&f )
////////////////////////////////////////////////////////////////////////
{
  // Construct NLP representation of TestProblem
  TestProblem<double> NLO;
  for( unsigned int ip=0; ip<np; ip++ ){
    NLO.Bounds(ip) = std::make_pair( Op<I>::l(P[ip]), Op<I>::u(P[ip]) );
    NLO.Guess(ip)  = p[ip];
  }

  // Test feasibility
  const double atol = 1e-8;
  for( unsigned int ic=0; ic<NLO.nc(); ic++ ){
    std::pair<NLPCTR<double>::TYPE,double> g = NLO( ic, p );
    switch( g.first ){
    case NLPCTR<double>::EQ:
#ifndef NOOUTPUT
      cout << "  G(" << ic << ") = " << g.second << " ?= 0" << endl;
#endif
      if( fabs(g.second) <= atol ) continue;
      return SBB<I>::INFEASIBLE;
    case NLPCTR<double>::LE:
#ifndef NOOUTPUT
      cout << "  G(" << ic << ") = " << g.second << " ?<= 0" << endl;
#endif
      if( g.second <= atol ) continue;
      return SBB<I>::INFEASIBLE;
    case NLPCTR<double>::GE:
#ifndef NOOUTPUT
      cout << "  G(" << ic << ") = " << g.second << " ?>= 0" << endl;
#endif
      if( g.second >= -atol ) continue;
      return SBB<I>::INFEASIBLE;
    }
  }

  f = NLO( p ).second;
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

    std::pair<NLPOBJ<LPV>::TYPE,LPV> F = NLO( PLP );
    switch( F.first ){
      case NLPOBJ<LPV>::MIN:
        LPG.set_objective( LPRG::MIN, F.second ); break;
      case NLPOBJ<LPV>::MAX:
        LPG.set_objective( LPRG::MAX, F.second ); break;
    }
    
    for( unsigned int ic=0; ic<NLO.nc(); ic++ ){
      std::pair<NLPCTR<LPV>::TYPE,LPV> F = NLO( ic, PLP );
      switch( F.first ){
        case NLPCTR<LPV>::EQ:
          LPG.add_constraint( F.second, LPRG::EQ, 0 ); break;
        case NLPCTR<LPV>::LE:
          LPG.add_constraint( F.second, LPRG::LE, 0 ); break;
        case NLPCTR<LPV>::GE:
          LPG.add_constraint( F.second, LPRG::GE, 0 ); break;
      }
    }
    
    delete [] PLP;
    
    // Solve LP relaxation
#ifndef NOOUTPUT
    cout << LPG;
    int status = LPG.solve( 0, 1, "relax.lp" );
    cout << "  Status:  " << status << endl;
#else
    int status = LPG.solve( 0, 0, "relax.lp" );
#endif

#ifdef USE_CPLEX
    if( status == IloAlgorithm::Optimal ){
#else
    if( status == GRB_OPTIMAL ){
#endif
      // Gather LP relaxation results
      f = LPG.get_objective();
      for( unsigned int ip=0; ip<NLO.np(); ip++ ){
        p[ip] = LPG.get_variable(ip);
        // m[ip] = LPG.get_reduced_cost(ip);
      }

      // Display LP relaxation results
#ifndef NOOUTPUT
      cout << "  Runtime: " << LPG.get_runtime() << " CPUsec" << endl
           << endl
           << "LP RELAXATION SOLUTION: " << endl
           << "  F* = " << f << endl;
      for( unsigned int ip=0; ip<NLO.np(); ip++ )
        cout << "  p*(" << ip << ") = " << p[ip]
             << endl;
#endif
      return SBB<I>::NORMAL;
    }

#ifdef USE_CPLEX
    else if( status == IloAlgorithm::Infeasible
          || status == IloAlgorithm::InfeasibleOrUnbounded ){
#else
    else if( status == GRB_INFEASIBLE 
          || status == GRB_INF_OR_UNBD ){
#endif
      return SBB<I>::INFEASIBLE;
    }

    else{
      return SBB<I>::FAILURE;
    }
  }
  
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( mc::I_Excp& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in natural interval extension." << endl
         << "Execution aborted..." << endl;
    return SBB<I>::FAILURE;
  }
#endif
#endif
  catch( mc::LPExcp& eObj ){
    cerr << "Error " << eObj.ierr()
         << " in LP relaxation." << endl
         << "Execution aborted..." << endl;
    return SBB<I>::FAILURE;
  }
#ifdef USE_CPLEX
  catch (IloException& ex) {
    cerr << "Error: " << ex << endl;
    return SBB<I>::FAILURE;
   }
#else
  catch(GRBException& ex) {
    cerr << "Error code = " << ex.getErrorCode() << endl;
    cerr << ex.getMessage() << endl;
    return SBB<I>::FAILURE;
  }
#endif
  catch(...) {
    cout << "Error during optimization" << endl;
    return SBB<I>::FAILURE;
  }

  return SBB<I>::NORMAL;
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
    return UFEAS( np, P, p, f );
  }
}

////////////////////////////////////////////////////////////////////////
std::set<unsigned int> EXCLUDE()
////////////////////////////////////////////////////////////////////////
{
  // Construct NLP structure
  TestProblem<Structure> NLOS;
  Structure *pS = new Structure[NLOS.np()];
  for( unsigned int ip=0; ip<NLOS.np(); ip++ ) pS[ip].indep(ip);   
  Structure L = NLOS( pS ).second;
  for( unsigned int ic=0; ic<NLOS.nc(); ic++ ) L += NLOS( ic, pS ).second;
  delete[] pS;
#ifndef NOOUTPUT
  cout << "  NLP dependence: " << L << endl;
#endif

  // Build exclusion set
  std::set<unsigned int> exclude_set;
  for( unsigned int ip=0; ip<NLOS.np(); ip++ ){
    if( L.dep().find(ip) == L.dep().end() 
     || (*L.dep().find(ip)).second ){
#ifndef NOOUTPUT
      cout << "  Exclude variable " << ip << endl;
#endif
      exclude_set.insert( ip );
    }
  }
  return exclude_set;
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
  //LPRG::options.BILINEAR_RULE   = mc::LPOpt::UNIVARIATE;
  LPRG::options.BILINEAR_RULE   = mc::LPOpt::SEPARSOS;
  
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
  //cout << ( UFEAS( NLO.np(), P, p )==SBB<I>::NORMAL? "FEASIBLE": "INFEASIBLE" )
  //     << endl;
  //LPB( NLO.np(), P, p, f );

  SBB<I> GO( SUBPB );
  GO.variables( NLO.np(), P );
  //GO.options.PROBLEM_TYPE = SBBOpt::MAXIMIZE;
  GO.options.BRANCHING_STRATEGY = SBBOpt::OMEGA;
  GO.options.BRANCHING_VARIABLE_CRITERION = SBBOpt::RANGEREL;
  GO.options.BRANCHING_VARIABLE_EXCLUSION = EXCLUDE();
  GO.options.DISPLAY = true;
  GO.options.MAX_CPU_TIME = 1e3;
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

template class mc::SBB<mc::Interval>;
