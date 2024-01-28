//#define USE_PROFIL
//#define USE_FILIB

#include <fstream>
#include <iomanip>
#include "dosbb.h"

#ifdef USE_PROFIL
  #include "mcprofil.h"
  typedef INTERVAL I;
  namespace fadbad
  {
    //! @brief C++ structure for specialization of the fadbad::Op structure to allow usage of the PROFIL interval arithmetic type INTERVAL inside the classes fadbad::T, fadbad::B and fadbad::T in FADBAD++
    template <> struct Op<INTERVAL>{
      typedef double Base;
      typedef INTERVAL I;
      static Base myInteger( const int i ) { return Base(i); }
      static Base myZero() { return myInteger(0); }
      static Base myOne() { return myInteger(1);}
      static Base myTwo() { return myInteger(2); }
      static double myPI() { return ::Constant::Pi; }
      static I myPos( const I& x ) { return  x; }
      static I myNeg( const I& x ) { return -x; }
      template <typename U> static I& myCadd( I& x, const U& y ) { return x+=y; }
      template <typename U> static I& myCsub( I& x, const U& y ) { return x-=y; }
      template <typename U> static I& myCmul( I& x, const U& y ) { return x*=y; }
      template <typename U> static I& myCdiv( I& x, const U& y ) { return x/=y; }
      static I myInv( const I& x ) { return myOne()/x; }
      static I mySqr( const I& x ) { return Sqr(x); }
      template <typename X, typename Y> static I myPow( const X& x, const Y& y ) { return ::Power( x, y ); }
      static I mySqrt( const I& x ) { return ::Sqrt( x ); }
      static I myLog( const I& x ) { return ::Log( x ); }
      static I myExp( const I& x ) { return ::Exp( x ); }
      static I mySin( const I& x ) { return ::Sin( x ); }
      static I myCos( const I& x ) { return ::Cos( x ); }
      static I myTan( const I& x ) { return ::Tan( x ); }
      static I myAsin( const I& x ) { return ::ArcSin( x ); }
      static I myAcos( const I& x ) { return ::ArcCos( x ); }
      static I myAtan( const I& x ) { return ::ArcTan( x ); }
      static bool myEq( const I& x, const I& y ) { return x==y; }
      static bool myNe( const I& x, const I& y ) { return x!=y; }
      static bool myLt( const I& x, const I& y ) { return x<y; }
      static bool myLe( const I& x, const I& y ) { return x<=y; }
      static bool myGt( const I& x, const I& y ) { return y<x; }
      static bool myGe( const I& x, const I& y ) { return y<=x; }
    };
  } // namespace fadbad
#else
  #ifdef USE_FILIB
    #include "mcfilib.h"
    typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
    //typedef filib::interval<double> I;
  namespace fadbad
  {
    //! @brief C++ structure for specialization of the fadbad::Op structure to allow usage of the PROFIL interval arithmetic type INTERVAL inside the classes fadbad::T, fadbad::B and fadbad::T in FADBAD++
    template <> struct Op<filib::interval<double,filib::native_switched,filib::i_mode_extended> >{
    //template <> struct Op<filib::interval<double> >{
      typedef double Base;
      typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
      //typedef filib::interval<double> I;
      static Base myInteger( const int i ) { return Base(i); }
      static Base myZero() { return myInteger(0); }
      static Base myOne() { return myInteger(1);}
      static Base myTwo() { return myInteger(2); }
      static double myPI() { return mc::PI; }
      static I myPos( const I& x ) { return  x; }
      static I myNeg( const I& x ) { return -x; }
      template <typename U> static I& myCadd( I& x, const U& y ) { return x+=y; }
      template <typename U> static I& myCsub( I& x, const U& y ) { return x-=y; }
      template <typename U> static I& myCmul( I& x, const U& y ) { return x*=y; }
      template <typename U> static I& myCdiv( I& x, const U& y ) { return x/=y; }
      static I myInv( const I& x ) { return myOne()/x; }
      static I mySqr( const I& x ) { return filib::sqr(x); }
      template <typename X> static I myPow( const X& x, const int n ) { return filib::power( x, n ); }
      template <typename X, typename Y> static I myPow( const X& x, const Y& y ) { return filib::pow( x, y ); }
      static I mySqrt( const I& x ) { return filib::sqrt( x ); }
      static I myLog( const I& x ) { return filib::log( x ); }
      static I myExp( const I& x ) { return filib::exp( x ); }
      static I mySin( const I& x ) { return filib::sin( x ); }
      static I myCos( const I& x ) { return filib::cos( x ); }
      static I myTan( const I& x ) { return filib::tan( x ); }
      static I myAsin( const I& x ) { return filib::asin( x ); }
      static I myAcos( const I& x ) { return filib::acos( x ); }
      static I myAtan( const I& x ) { return filib::atan( x ); }
      static bool myEq( const I& x, const I& y ) { return x.seq(y); }
      static bool myNe( const I& x, const I& y ) { return x.sne(y); }
      static bool myLt( const I& x, const I& y ) { return x.slt(y); }
      static bool myLe( const I& x, const I& y ) { return x.sle(y); }
      static bool myGt( const I& x, const I& y ) { return x.sgt(y); }
      static bool myGe( const I& x, const I& y ) { return x.sge(y); }
    };
  } // namespace fadbad
  #else
    #include "interval.h"
    typedef mc::Interval I;
  #endif
#endif
typedef mc::McCormick<I> MC;
typedef mc::TModel<I> TM;
typedef mc::TModel<MC> TMMC;

#define OC8			// <- select the example here

#if defined( OC1 )
////////////////////////////////////////////////////////////////////////
// MINIMUM TIME CONTROL PROBLEM WITH TERMINAL EQUALITY CONSTRAINTS
// 
//    MIN  tf
//   tf, u, p
//     s.t.   x1(tf) = x1f
//            x2(tf) = x2f                 _
//             dx1dt = x2                   |  t in (0,tf]
//             dx2dt = (u-x1-2*x2)*t       _|
//             x1(0) = x10, x2(0) = p
//            umin <= u  <= umax
//            pmin <= p  <= pmax
//           tfmin <= tf <= tfmax
// 
//     where:  x1f   =  1.0
//             x2f   =  0.0
//             x10   =  0.0
//             umin  = -0.5
//             umax  =  1.5
//             pmin  = -1.0
//             pmax  =  1.0
//             tfmin =  0.1
//             tfmax =  2.0
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 10;
const unsigned int NP = 2+NS;
const unsigned int NX = 3;
const unsigned int NC = 2;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : mc::DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return p[0]*x[1];
        case 1: return p[0]*(p[2+is]*x[0]-2*x[1])*x[2];
        case 2: return p[0];
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 0.;
        case 1: return p[1];
        case 2: return 0.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      return std::make_pair( p[0], MIN );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        case 0: return std::make_pair( xk[ns][0]-1., EQ );
        case 1: return std::make_pair( xk[ns][1]-0., EQ );
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( OC2 )
////////////////////////////////////////////////////////////////////////
// TEST PROBLEM #1 IN FLOUDAS ET AL
// 
//    MIN  -z(1)^2
//     u
//     s.t.  dzdt = -z^2+u,  t in (0,1]
//           z(0) = 9
//           umin <= u  <= umax
// 
//     where:  umin  = -5
//             umax  =  5
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 8;
const unsigned int NP = NS;
const unsigned int NX = 1;
const unsigned int NC = 0;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : mc::DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      using mc::sqr;
      switch( ix ){
        case 0: return -sqr(x[0])+p[is];
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 9.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      using mc::sqr;
      return std::make_pair( sqr(xk[ns][0]), MAX );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( OC3 )
////////////////////////////////////////////////////////////////////////
// TEST PROBLEM #2 IN FLOUDAS ET AL (SINGULAR CONTROL PROBLEM)
// 
//    MIN  z4(1)
//     u                                                     _
//     s.t.  dz1dt = z2,                                      |
//           dz2dt = -z3*u+16*t-8,                            |
//           dz3dt = u,                                       |  t in (0,1]
//           dz4dt = z1^2+z2^2+5e-4*(z2+16*t-8-0.1*z3*u^2)^2,_|
//           z(0) = [ 0, -1, -sqrt(5), 0 ]
//           -4 <= u  <= 10
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 4;
const unsigned int NP = NS;
const unsigned int NX = 5;
const unsigned int NC = 0;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      using mc::sqr;
      switch( ix ){
        case 0: return x[1];
        case 1: return -x[2]*p[is]+16.*x[4]-8.;
        case 2: return p[is];
        case 3: return sqr(x[0])+sqr(x[1])
                       +5e-4*sqr(x[1]+16.*x[4]-8.-0.1*x[2]*sqr(p[is]));
        case 4: return 1.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 0.;
        case 1: return -1.;
        case 2: return -sqrt(5.);
        case 3: return 0.;
        case 4: return 0.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      return std::make_pair( -xk[ns][3], MAX );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( OC4 )
////////////////////////////////////////////////////////////////////////
// TEST PROBLEM #4 IN FLOUDAS ET AL (OIL SHALE PYROLYSIS)
// 
//    MAX  z2(1)
//    p,T                                                    _
//     s.t.  dz1dt = -p*(k1*z1+k3*z1*z2+k4*z1*z2+k5*z1*z2),   |
//           dz2dt = p*(k1*z1-k2*z2+k3*z1*z2),                |  t in (0,1]
//              ki = ai*exp(-bi/(R*T)), i=1,...,5            _|
//           z(0) = [ 1, 0 ],
//           698.15 <= T <= 748.15,
//                7 <= p <= 11
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 2;
const unsigned int NP = NS+1;
const unsigned int NX = 2;
const unsigned int NC = 0;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : mc::DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      double A[5] = { 8.86, 24.25, 23.67, 18.75, 20.70 };
      double B[5] = { 10215.4, 18820.5, 17008.9, 14190.8, 15599.8 };
      using mc::exp;
      using std::exp;
      switch( ix ){
        case 0: return -p[0]*(exp(A[0]-B[0]*p[1+is])+(exp(A[2]-B[2]*p[1+is])
                       +exp(A[3]-B[3]*p[1+is])+exp(A[4]-B[4]*p[1+is]))*x[1])*x[0];
        case 1: return p[0]*(exp(A[0]-B[0]*p[1+is])*x[0]-exp(A[1]-B[1]*p[1+is])*x[1]
                       +exp(A[2]-B[2]*p[1+is])*x[0]*x[1]);
        default: throw std::runtime_error("invalid size");
      }
    }

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

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      return std::make_pair( xk[ns][1], MAX );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( OC5 )
////////////////////////////////////////////////////////////////////////
// TEST PROBLEM (VAN DE VUSSE REACTIONS)
// 
//    MAX  z2(10)
//    p,T                                              _
//     s.t.  dCAdt = F/V*(CAin-CA) - k1*CA - k3*CA^2,   |
//           dCBdt = -F/V*CB + k1*CA - k2*CB           _|  t in (0,10]
//           CA(0) = 1, CB(0) = 0, CAin = 10, V = 10,
//           k1 = 0.5, k2 = 1, k3 = 0.1,
//           0 <= F <= 15
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 5;
const unsigned int NP = NS;
const unsigned int NX = 2;
const unsigned int NC = 0;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : mc::DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      using mc::sqr;
      double k1=5e-1, k2=1e0, k3=1e-1, V=1e1, CAin=1e1;
      switch( ix ){
        case 0: return p[is]/V*(CAin-x[0]) - k1*x[0] - k3*sqr(x[0]);
        case 1: return -p[is]/V*x[1] + k1*x[0] - k2*x[1];
        default: throw std::runtime_error("invalid size");
      }
    }

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

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      return std::make_pair( xk[ns][1], MAX );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( OC8 )
////////////////////////////////////////////////////////////////////////
// TEST PROBLEM: FLOW CONTROL IN SEMI-BATCH REACTOR
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 10;
const unsigned int NP = NS;
const unsigned int NX = 3;
const unsigned int NC = 2;

const double k1 = 0.053, k2 = 0.128;
const double CA0 =0.72, CB0 = 0.05, CBin = 5., V0 = 1., pl = 0, pu = 1e-3;

class DYNOPT : public mc::DOSTRUCT
{
public:
  DYNOPT()
    : mc::DOSTRUCT( NP, NX, NC )
    {}

  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      using mc::sqr;
      switch( ix ){
        case 0: return -k1*x[0]*x[1] - p[is]/x[2]*x[0];
        case 1: return -k1*x[0]*x[1] - 2*k2*sqr(x[1]) + p[is]/x[2]*(CBin - x[1]);
	case 2: return p[is];
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return CA0;
        case 1: return CB0;
        case 2: return V0;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  std::pair<T,t_OBJ> OBJ
    ( const T*p, T* const*xk, const unsigned int ns )
    {
      return std::make_pair( CA0*V0/xk[ns][2] - xk[ns][0], MAX );
    }

  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*xk, const unsigned int ns )
    {
      switch( ic ){
        case 0: return std::make_pair( xk[ns][1] - 0.02, LE );
        case 1: return std::make_pair( 0.5*(xk[ns][0]+CBin-xk[ns][1])
                                     - 0.5*(CA0+CBin-CB0)*V0/xk[ns][2] - 0.005, LE );
        default: throw std::runtime_error("invalid size");
      }
    }
};

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
#if defined( OC1 )
  double p0[NP], tk[NS+1];
  I P[NP];
  p0[0] = 1.;   P[0] = I( 0.1, 2. );
  p0[1] = 0.5;  P[1] = I( -1., 1. );
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[2+is]  = I( -0.5, 1.5 );
    p0[2+is] = 0.;
    tk[1+is] = tk[is]+1./(double)NS;
  }
#elif defined( OC2 )
  double p0[NP], tk[NS+1];
  I P[NP];
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[is]  = I( -5., 5. );
    p0[is] = 0.;
    tk[1+is] = tk[is]+1./(double)NS;
  }
#elif defined( OC3 )
  double p0[NP], tk[NS+1];
  I P[NP];
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[is]  = I( -4., 10. );
    p0[is] = 0.;
    tk[1+is] = tk[is]+1./(double)NS;
  }
#elif defined( OC4 )
  double p0[NP], tk[NS+1];
  I P[NP];
  p0[0] = 9.;   P[0] = I( 7., 11. );
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[1+is]  = 1./I( 698.15, 748.15 );
    p0[1+is] = 1./723.15;
    tk[1+is] = tk[is]+1./(double)NS;
  }
#elif defined( OC5 )
  double p0[NP], tk[NS+1], TF=1e1;
  I P[NP];
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[is]  = I( 0., 15. );
    p0[is] = 6.;
    tk[1+is] = tk[is]+TF/(double)NS;
  }
//     P[0]  = I( 0., 6.750910936 );
//     P[1]  = I( 0., 3.75 );
//     P[2]  = I( 0., 2.616014892 );
//     P[3]  = I( 0., 3.75 );
//     P[4]  = I( 0., 5.810975484 );
#elif defined( OC8 )
  double p0[NP], tk[NS+1], TF=5e1;
  I P[NP];
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    //P[is]  = I( 1e-4, 1.e-3 );
    P[is]  = I( pl, pu );
    //P[is]  = I( 1.01e-4, 1.02e-4 );
    p0[is] = 1.015e-4;
    tk[1+is] = tk[is]+TF/(double)NS;
  }
#endif
  
  mc::DOSBB<I,DYNOPT> GDO;

  GDO.options.TAYLOR_MODEL_ORDER          = 2;
  GDO.options.IVP_BOUNDING_TYPE           = mc::DOSBB<I,DYNOPT>::Options::IA;
  GDO.options.IVP_BOUNDING_PROPAGATION    = mc::DOSBB<I,DYNOPT>::Options::DINEQ;//VERIF;
  GDO.options.USE_CONSTRAINT_PROPAGATION  = false;	//try false/true
  GDO.options.USE_DOMAIN_REDUCTION        = true;	//try false/true
  GDO.options.DOMAIN_REDUCTION_MAXLOOPS   = 10;	
  GDO.options.USE_RLT_TMODEL              = true;	//try false/true
	
  GDO.options_TModel().BOUNDER_TYPE       = mc::TModel<I>::Options::LSB;
  GDO.options_TModel().SCALE_VARIABLES    = false;

  GDO.options_SBB().BRANCHING_STRATEGY    = mc::SBB<I>::Options::MIDPOINT;
//   GDO.options_SBB().BRANCHING_STRATEGY = mc::SBB<I>::Options::OMEGA;
//   GDO.options_SBB().BRANCHING_VARIABLE_CRITERION = mc::SBB<I>::Options::RELIAB;
//   GDO.options_SBB().BRANCHING_RELIABILITY_THRESHOLD = 1;
  GDO.options_SBB().DISPLAY               = 2;
  GDO.options_SBB().MAX_CPU_TIME          = 1e6;
  GDO.options_SBB().MAX_NODES             = 100000;
  GDO.options_SBB().ABSOLUTE_TOLERANCE    = 1e-3;
  GDO.options_SBB().RELATIVE_TOLERANCE    = 1e-3;

  GDO.options_ODEBND().TSORDER            = 10;
  GDO.options_ODEBND().HMIN               = 1.e-8;

  GDO.options_ODEBND_GSL().INTMETH        = mc::ODEBND_GSL<I,DYNOPT>::Options::RKF45;//RK8PD;
  GDO.options_ODEBND_GSL().WRAPMIT        = mc::ODEBND_GSL<I,DYNOPT>::Options::ELLIPS;//DINEQ;//LINTRANS;
  GDO.options_ODEBND_GSL().ORDMIT         = 2;
  GDO.options_ODEBND_GSL().RTOL           = GDO.options_ODEBND_GSL().ATOL = 1e-7;
  GDO.options_ODEBND_GSL().HMIN           = 1.e-7;
  GDO.options_ODEBND_GSL().QTOL           = 1.e-12;
  GDO.options_ODEBND_GSL().DISPLAY        = 0;

  GDO.options_DOSEQ().DISPLAY  = 0;
  GDO.options_DOSEQ().MAXITER  = 50;
  GDO.options_DOSEQ().CVTOL    = 1e-9;
  GDO.options_DOSEQ().GRADIENT = mc::DOSEQ<DYNOPT>::Options::FORWARD;
//   GDO.options_DOSEQ().SPARSE   = true;
  GDO.options_DOSEQ().HESSIAN  = mc::DOSEQ<DYNOPT>::Options::LBFGS;

//  GDO.options_FPRelax().SANDWICH_MAXCUT = 5;
  GDO.options_FPRelax().RLT_MAXORDER      = 3;
  GDO.options_FPRelax().RLT_DISPLAY       = 0;
// 
//   GDO.options_CVODES().RELTOL  = GDO.options_CVODES().ABSTOL = 1e-8;
//   GDO.options_CVODES().FSAERR  = true;

//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.DISPLAY = 0;
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.ASACHKPT = 4000;
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.FSAERR = true;
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.HESSFORMAT = mc::CVODES<DYNOPT>::Options::HESSLOWER;
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.FDATOL =
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.FDRTOL = 1e-4;
//   GDO.options_DOSEQ().CVODES<DYNOPT>::options.FDCEN  = true;

  try{ 
    std::cout << GDO;
    std::pair<double,const double*> optim;
  
//     GDO.options.TAYLOR_MODEL_ORDER = 4;
//     GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TM;
//     GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::EIGEN;
//     std::cout << "ODE BOUNDER: TM   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//       << "  TM BOUNDER: EIGEN   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//       << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
//     optim = GDO.solve( NS, tk, P, p0 );
// 
//     GDO.options.TAYLOR_MODEL_ORDER = 5;
//     GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TM;
//     GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::EIGEN;
//     std::cout << "ODE BOUNDER: TM   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//       << "  TM BOUNDER: EIGEN   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//       << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
//     optim = GDO.solve( NS, tk, P, p0 );
// 
//     GDO.options.TAYLOR_MODEL_ORDER = 3;
//     GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TM;
//     GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::EIGEN;
//     std::cout << "ODE BOUNDER: TM   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//       << "  TM BOUNDER: EIGEN   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//       << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
//     optim = GDO.solve( NS, tk, P, p0 );
//   
//     GDO.options.TAYLOR_MODEL_ORDER = 4;
//     GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TMMC;
//     GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::EIGEN;
//     std::cout << "ODE BOUNDER: TMMC   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//       << "  TM BOUNDER: EIGEN   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//       << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
//     optim = GDO.solve( NS, tk, P, p0 );
// 
//     GDO.options.TAYLOR_MODEL_ORDER = 5;
//     GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TMMC;
//     GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::EIGEN;
//     std::cout << "ODE BOUNDER: TMMC   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//       << "  TM BOUNDER: EIGEN   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//       << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
//     optim = GDO.solve( NS, tk, P, p0 );

//    GDO.options.TAYLOR_MODEL_ORDER = 4;
//    GDO.options.IVP_BOUNDING_TYPE = mc::DOSBB<I,DYNOPT>::Options::TM;
//    GDO.options_TModel().BOUNDER_TYPE = mc::TModel<I>::Options::HYBRID;
//    std::cout << "ODE BOUNDER: TMMC   ORDER: " << GDO.options.TAYLOR_MODEL_ORDER
//      << "  TM BOUNDER: HYBRID   RTOL=" << GDO.options_SBB().RELATIVE_TOLERANCE
//      << "  ATOL=" << GDO.options_SBB().ABSOLUTE_TOLERANCE << std::endl;
    optim = GDO.solve( NS, tk, P, p0 );
    std::cout << std::fixed << std::setprecision(1)
      << "\n  Total Time: " << GDO.stats.cumul_SBB
      << "\n  DOSEQ: " << GDO.stats.cumul_DOLOC/GDO.stats.cumul_SBB*1e2 << "%"
      << "    ODEBND: " << GDO.stats.cumul_ODEBND/GDO.stats.cumul_SBB*1e2 << "%"
      << "    LPREL:  " << GDO.stats.cumul_FPREL/GDO.stats.cumul_SBB*1e2 << "%"
      << "    REDUC:  " << GDO.stats.cumul_REDUC/GDO.stats.cumul_SBB*1e2 << "%"
      << "    DOREL:  " << GDO.stats.cumul_DOREL/GDO.stats.cumul_SBB*1e2 << "%"
      << std::endl;

//   if( status == 0 ){
//     std::cout << "DO (GLOBAL) SOLUTION: " << std::endl;
//     std::cout << "  f* = " << pGDO->solution().f << std::endl;
//     for( unsigned int ip=0; ip<NP; ip++ )
//       std::cout << "  p*(" << ip << ") = " << pGDO->solution().p[ip]
//                 << std::endl;
//   }
  }
#ifndef USE_PROFIL
#ifndef USE_FILIB
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension." << std::endl
              << "Execution aborted..." << std::endl;
    return eObj.ierr();
  }
#endif
#endif
  catch( TM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in Taylor model computation." << std::endl
              << "Execution aborted..." << std::endl;
    return eObj.ierr();
  }
  catch( TMMC::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in McCormick-Taylor model computation." << std::endl
              << "Execution aborted..." << std::endl;
    return eObj.ierr();
  }
  catch( mc::DOSBB<I,DYNOPT>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr() << " in DOSBB:" << std::endl
              << eObj.what() << std::endl
              << "Execution aborted..." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
