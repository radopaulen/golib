#include <fstream>
#include <iomanip>
#include <stdexcept>

#include "ellipsoid_cpplapack.hpp"

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  //mc::Ellipsoid::options.PSDCHK = false;
  //mc::Ellipsoid::options.REGTOL = 1e-8;

  CPPL::dcovector c1(2);
  c1(0) = c1(1) = 1.;
  CPPL::dsymatrix Q1(2);
  //Q1(0,0) = .5; Q1(0,1) = 0.1; Q1(1,1) = 1.;
  Q1(0,0) = 5.; Q1(0,1) = 4.; Q1(1,1) = 5.;
  mc::Ellipsoid E1( Q1, c1 ); 
  std::cout << "-- E1:" << E1 << std::endl;
  std::cout << "rank: " << E1.rankQ() << std::endl;

  // Testing eignevalues
  //std::pair< const CPPL::dcovector&, const CPPL::dgematrix& > eigE1 = E1.eigQ();
  //std::cout << "-- E1 eigenvalues:\n" << eigE1.first << std::endl;
  //std::cout << "-- E1 eigenvectors:\n" << eigE1.second << std::endl;

  // Testing square-root
  const CPPL::dsymatrix sqrtE1 = E1.sqrtQ();
  std::cout << "-- E1 square-root:\n" << sqrtE1 << std::endl;

  // Testing linear transformation
  const unsigned int mA1 = 2;
  CPPL::dgematrix A1(mA1,2);
  A1(0,0) = 1.; A1(0,1) = -1.;
  A1(1,0) = 1.; A1(1,1) = 1.;
  std::cout << "-- A1:\n" << A1 << std::endl;
  mc::Ellipsoid A1E1 = mtimes( E1, A1 ); 
  std::cout << "-- A1*E1:" << A1E1 << std::endl;
  mc::Ellipsoid A1invA1E1 = mtimes( A1E1, i(A1) ); 
  std::cout << "-- A1inv*A1*E1:" << A1invA1E1 << std::endl;

  CPPL::dcovector c2(2);
  c2(0) = c2(1) = 0.;
  CPPL::dsymatrix Q2(2);
  Q2(0,0) = 1.; Q2(0,1) = 0.2; Q2(1,1) = .5;
  mc::Ellipsoid E2( Q2, c2 ); 
  std::cout << "-- E2:" << E2 << std::endl;

  // Testing Minkowski sum (outer-approximation)
  std::vector<CPPL::dcovector> D;
  CPPL::dcovector d1(2); d1(0) = 1.; d1(1) = 0.; D.push_back( d1 );
  CPPL::dcovector d2(2); d2(0) = 0.; d2(1) = 1.; D.push_back( d2 );
  CPPL::dcovector d3(2); d3(0) = 1.; d3(1) = 1.; D.push_back( d3 );
  std::vector<mc::Ellipsoid> E;
  E.push_back( E1 );
  E.push_back( E2 );
  std::vector<mc::Ellipsoid> EA = minksum_ea( E, D ); 
  typename std::vector<mc::Ellipsoid>::iterator EAi = EA.begin();
  for( unsigned int i=1; EAi!=EA.end(); ++EAi, i++ )
     std::cout << "-- E1 + E2 (dir. #" << i << "):" << *EAi << std::endl;
 
  //mc::Ellipsoid::options.RANKTOL = 1e-5;
  const unsigned int mA2 = 3;
  CPPL::dgematrix A2(mA2,2);
  A2(0,0) = 1.; A2(0,1) = 1.;
  A2(1,0) = 0.; A2(1,1) = 0.;
  A2(2,0) = 1.; A2(2,1) = 1.;
  std::cout << "-- A2:\n" << A2 << std::endl;
  CPPL::dcovector b2(mA2);
  b2(0) = 1.; b2(1) = -1.; b2(2) = 0.;
  std::cout << "-- b2:\n" << b2 << std::endl;
  mc::Ellipsoid A2E1 = mtimes( E1, A2, b2 ); 
  std::cout << "-- A2*E1+b2:" << A2E1 << std::endl;
  std::cout << "sing val: " << t(A2E1.svdQ().first) << std::endl;
  std::cout << "rank: " << A2E1.rankQ() << std::endl;
  A2E1.regQ();
  std::cout << "-- A2*E1 reg:" << A2E1 << std::endl;
  std::cout << "sing val: " << t(A2E1.svdQ().first) << std::endl;
  std::cout << "rank: " << A2E1.rankQ();
  mc::Ellipsoid A2E1inv = inv( A2E1 ); 
  std::cout << "-- A2*E1 inv:" << A2E1inv << std::endl;
  std::cout << "-- A2*E1 inv inv:" << inv(A2E1inv) << std::endl;

  // Testing ellipsoid extension
  double c1_3 = 1.;
  CPPL::drovector Q1_3(3);
  Q1_3(0) = Q1_3(1) = 0.; Q1_3(2) = 1.;
  E1.extend( Q1_3, c1_3 );
  std::cout << "-- E1 ext:" << E1 << std::endl;

  // Testing intersection with hyperplanes / halfspaces
  CPPL::dcovector c3(3);
  c3(0) = c3(1) = 0.5; c3(2) = 0.;
  CPPL::dsymatrix Q3(3); Q3.zero();
  Q3(0,0) = Q3(1,1) = Q3(2,2) = 1.;
  mc::Ellipsoid E3( Q3, c3 );
  std::cout << "-- E3:" << E3 << std::endl;
  std::cout << "rank: " << E3.rankQ() << std::endl;

  CPPL::dcovector hp0(3); hp0(0) = 1.; hp0(1) = hp0(2) = 0.;
  std::pair<CPPL::dcovector,double> HP0( hp0, .2 );//std::sqrt(3.) );
  mc::Ellipsoid E3HP0 = hpintersection( E3, HP0 );
  std::cout << "-- E3 inter HP0:" << E3HP0;
  std::cout << "sing val: " << t(E3HP0.svdQ().first);
  std::cout << "rank: " << E3HP0.rankQ() << std::endl << std::endl;

  CPPL::dcovector hp1(3); hp1(0) = 0.; hp1(1) = 1.; hp1(2) = 0.;
  std::pair<CPPL::dcovector,double> HP1( hp1, 0. );
  mc::Ellipsoid E3HP0HP1 = hpintersection( E3HP0, HP1 );
  std::cout << "-- E3 inter HP0 inter HP1:" << E3HP0HP1;
  std::cout << "sing val: " << t(E3HP0HP1.svdQ().first);
  std::cout << "rank: " << E3HP0HP1.rankQ() << std::endl << std::endl;

  mc::Ellipsoid E3HP0P = intersection_ea( E3, HP0 );
  std::cout << "-- E3 inter HP0+ :" << E3HP0P;
  std::cout << "sing val: " << t(E3HP0P.svdQ().first);
  std::cout << "rank: " << E3HP0P.rankQ() << std::endl << std::endl;

  std::pair<CPPL::dcovector,double> HP0M( -hp0, -.2 );//std::sqrt(3.) );
  mc::Ellipsoid E3HP0PHP0M = intersection_ea( E3HP0P, HP0M );
  std::cout << "-- E3 inter HP0+ inter HP0- :" << E3HP0PHP0M;
  std::cout << "sing val: " << t(E3HP0PHP0M.svdQ().first);
  std::cout << "rank: " << E3HP0PHP0M.rankQ() << std::endl << std::endl;

  mc::Ellipsoid E3HP0M = intersection_ea( E3, HP0M );
  std::cout << "-- E3 inter HP0- :" << E3HP0M;
  std::cout << "sing val: " << t(E3HP0M.svdQ().first);
  std::cout << "rank: " << E3HP0M.rankQ() << std::endl << std::endl;

  mc::Ellipsoid E3HP0MHP0P = intersection_ea( E3HP0M, HP0 );
  std::cout << "-- E3 inter HP0- inter HP0+ :" << E3HP0MHP0P;
  std::cout << "sing val: " << t(E3HP0MHP0P.svdQ().first);
  std::cout << "rank: " << E3HP0MHP0P.rankQ() << std::endl << std::endl;

/*
  CPPL::dcovector hp2(3); hp2(2) = 1.; hp2(0) = hp2(1) = 0.;
  std::pair<CPPL::dcovector,double> HP2( hp2, 0.5 );
  mc::Ellipsoid E3HP1HP2 = hpintersection( E3HP1, HP2 );
  std::cout << "-- E3 inter HP1 inter HP2:" << E3HP1HP2 << std::endl;
  std::cout << "sing val: " << t(E3HP1HP2.svdQ().first) << std::endl;
  std::cout << "rank: " << E3HP1HP2.rankQ() << std::endl;
*/
//  CPPL::dcovector hp1(3); hp1(0) = 1.; hp1(1) = hp1(2) = 0.;
//  std::pair<CPPL::dcovector,double> HP1( hp1, 1. );
//  mc::Ellipsoid E1HP1 = hpintersection( E1, HP1 );
//  std::cout << "-- E1 inter HP1:" << E1HP1 << std::endl;
//  CPPL::dcovector hp2(3); hp2(0) = 0.; hp2(1) = 1.; hp2(2) = -1.;
//  std::pair<CPPL::dcovector,double> HP2( hp2, 0. );
//  mc::Ellipsoid E1HP2 = hpintersection( E1, HP2 );
//  std::cout << "-- E1 inter HP2:" << E1HP2 << std::endl;
//  mc::Ellipsoid E1HP1HP2 = hpintersection( E1HP1, HP2 );
//  std::cout << "-- E1 inter HP1 inter HP2:" << E1HP1HP2 << std::endl;
//  mc::Ellipsoid E1HP2HP1 = hpintersection( E1HP2, HP1 );
//  std::cout << "-- E1 inter HP2 inter HP1:" << E1HP2HP1 << std::endl;

  //mc::Ellipsoid Einter = intersection_ea( E[0], E[1] );
  //std::cout << Einter << std::endl;

  return 0;
}
