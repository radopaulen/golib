#undef USE_PROFIL   // <- change to #define for selecting PROFIL
#undef USE_FILIB    // <- change to #define for selecting FILIB++

#include <fstream>
#include <ostream>
#include <iomanip>
#include "odebnd_gsl.hpp"
#include "odebnd_val.hpp"

#ifdef USE_PROFIL
  #include "mcprofil.hpp"
  typedef INTERVAL I;
#else
  #ifdef USE_FILIB
    #include "mcfilib.hpp"
    typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
  #else
    #include "interval.hpp"
    typedef mc::Interval I;
  #endif
#endif

typedef mc::Ellipsoid E;
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;

////////////////////////////////////////////////////////////////////////

#define AM2b				// <- select the example here
const unsigned int NTM   = 3;	// <- Order of Taylor model
const unsigned int NSAMP = 40;	// <- Number of sampling points
#define SAVE_RESULTS			// <- Saving the results to file

////////////////////////////////////////////////////////////////////////

#if defined( LV )
std::string example = "LV";
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = p*x1*(1-x2)          |  t in (0,tf]
//             dx2dt = p*x2*(x1-1)         _|
//             x1(0) = 1.2, x2(0) = 1.1
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 300;
const unsigned int NP = 1;
const unsigned int NX = 2;
const unsigned int NI = 1;
double T0 = 0., TF = 15.;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX, NI )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return p[0]*x[0]*(1.-x[1]);
        case 1: return p[0]*x[1]*(x[0]-1.);
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      using mc::sqr;
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.2+sqr(p[0]-3.);
        case 1: return 1.1+(p[0]-3.);
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename TX, typename TP, typename TT>
  TX INV
    ( const unsigned int ii, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ii < ni() );
      switch( ii ){
        case 0: return x[0]-IC(0,p)-log(x[0])+log(IC(0,p))
                      +x[1]-IC(1,p)-log(x[1])+log(IC(1,p));
        default: throw std::runtime_error("invalid index");
      }
    }
};

#elif defined( OSCIL )
std::string example = "OSCIL";
////////////////////////////////////////////////////////////////////////
//                                   _
//             dx1dt = x2 + p         |  t in (0,tf]
//             dx2dt = -x1 - p       _|
//             x1(0) = 0., x2(0) = 0.
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 500;
const unsigned int NP = 1;
const unsigned int NX = 2;
double T0 = 0., TF = 10.;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return  x[1] + p[0];
        case 1: return -x[0] - p[0];
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.;//p[0];
        case 1: return 0.;
        default: throw std::runtime_error("invalid index");
      }
    }
};

#elif defined( COMPART )
std::string example = "COMPART";
////////////////////////////////////////////////////////////////////////
//  COMPARTMENTAL MODEL FROM THE PAPER BY WALTER AND KIEFFER (2011)
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 100;
const unsigned int NP = 3;
const unsigned int NX = 2;
double T0 = 0., TF = (double)NS;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
	case 0: return -(p[0]+p[2])*x[0]+p[1]*x[1];
	case 1: return p[0]*x[0]-p[1]*x[1];
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
};

#elif defined( LTV )
std::string example = "LTV";
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = x2                   |  t in (0,tf]
//             dx2dt = (u-x1-2*x2)*t       _|
//             x1(0) = x10, x2(0) = x20
//            umin <= u  <= umax
//           tfmin <= tf <= tfmax
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 3;
const unsigned int NP = 2+NS;
const unsigned int NX = 3;
double T0 = 0., TF = 1.;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return p[0]*x[1];
        case 1: return p[0]*(p[2+is]*x[0]-2.*x[1])*x[2];
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
};

#elif defined( RXN1 )
std::string example = "RXN1";
////////////////////////////////////////////////////////////////////////
//                                         _
//             r1 = k1 * x1 * x2            |
//             r2 = 2e1* x1 * x3            |
//             dx1dt = - r1 - r2            |  t in (0,0.1]
//             dx2dt = - r1                 |
//             dx3dt =  r1 - r2             |
//             dx4dt =  r2                 _|
//             x1(0) = 1, x2(0) = x20
//             50 <= k1  <= 500
//           0.95 <= x20 <= 1.05
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 500;
const unsigned int NP = 2;
const unsigned int NX = 4;
const unsigned int NI = 2;
double T0 = 0., TF = 5e-1;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX, NI )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX r1 = p[0] * x[0] * x[1];
      TX r2 = 2e1  * x[0] * x[2];
      double S1[4] = { -1., -1.,  1.,  0. };
      double S2[4] = { -1,   0., -1.,  1. };
      return r1 * S1[ix] + r2 * S2[ix];
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.;
        case 1: return p[1];
        case 2: case 3: return 0.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename TX, typename TP, typename TT>
  TX INV
    ( const unsigned int ii, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ii < ni() );
      switch( ii ){
        case 0: return          x[1] + x[2] + x[3] - IC(1,p);
        case 1: return   x[0] - x[1]        + x[3] - IC(0,p) + IC(1,p);
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename TX, typename TP, typename TT>
  TX BND
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX xi;
      switch( ix ){
        case 0: mc::Op<TX>::inter( xi, x[0], TX(0.,1.) ); break;
        case 1: mc::Op<TX>::inter( xi, x[1], TX(0.,mc::Op<TX>::u(p[1])) ); break;
        case 2: mc::Op<TX>::inter( xi, x[2], TX(0.,1.) ); break;
        case 3: mc::Op<TX>::inter( xi, x[3], TX(0.,1.) ); break;
        default: throw std::runtime_error("invalid size");
      }
      return xi;
    }
};

#elif defined( RXN1B )
std::string example = "RXN1B";
////////////////////////////////////////////////////////////////////////
//                                         _
//             r1 = k1 * x1 * x2            |
//             r2 = 2e1* x1 * x3            |
//             dx1dt = - r1 - r2            |  t in (0,0.1]
//             dx2dt = - r1                 |
//             dx3dt =  r1 - r2             |
//             dx4dt =  r2                 _|
//             x1(0) = 1, x2(0) = x20
//             50 <= k1  <= 500
//           0.95 <= x20 <= 1.05
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 500;
const unsigned int NP = 2;
const unsigned int NX = 2;
double T0 = 0., TF = 5e1;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX r1 = p[0] * x[0] * x[1];
      TX r2 = 2e1  * x[0] * ( -IC(0,p) + 2.*IC(1,p) + x[0] - 2.*x[1] );
      double S1[4] = { -1., -1. };
      double S2[4] = { -1,   0. };
      return r1 * S1[ix] + r2 * S2[ix];
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.;
        case 1: return p[1];
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( RXN1C )
std::string example = "RXN1C";
////////////////////////////////////////////////////////////////////////
//                                            _
//             r1 = k1 * x1 * x2 - 1e-1 * x3   |
//             r2 = 2e1* x1 * x3 - 1e-1 * x4   |
//             dx1dt = - r1 - r2               |  t in (0,0.1]
//             dx2dt = - r1                    |
//             dx3dt =  r1 - r2                |
//             dx4dt =  r2                    _|
//             x1(0) = 1, x2(0) = x20
//             50 <= k1  <= 500
//           0.95 <= x20 <= 1.05
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 1000;
const unsigned int NP = 2;
const unsigned int NX = 4;
const unsigned int NI = 2;
double T0 = 0., TF = 5e0;

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX, NI )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX r1 = p[0] * x[0] * x[1] - 1e0 * x[2];
      TX r2 = 2e1  * x[0] * x[2] - 1e0 * x[3];
      double S1[4] = { -1., -1.,  1.,  0. };
      double S2[4] = { -1,   0., -1.,  1. };
      return r1 * S1[ix] + r2 * S2[ix];
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.;
        case 1: return p[1];
        case 2: case 3: return 0.;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename TX, typename TP, typename TT>
  TX INV
    ( const unsigned int ii, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ii < ni() );
      switch( ii ){
        case 0: return          x[1] + x[2] + x[3] - IC(1,p);
        case 1: return   x[0] - x[1]        + x[3] - IC(0,p) + IC(1,p);
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename TX, typename TP, typename TT>
  TX BND
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX xi;
      switch( ix ){
        case 0: mc::Op<TX>::inter( xi, x[0], TX(0.,1.) ); break;
        case 1: mc::Op<TX>::inter( xi, x[1], TX(0.,mc::Op<TX>::u(p[1])) ); break;
        case 2: mc::Op<TX>::inter( xi, x[2], TX(0.,1.) ); break;
        case 3: mc::Op<TX>::inter( xi, x[3], TX(0.,1.) ); break;
        default: throw std::runtime_error("invalid size");
      }
      return xi;
    }
};

#elif defined( RXN2 )
std::string example = "RXN2";
////////////////////////////////////////////////////////////////////////
//                                                   _
//             r1 = p1 * k10 * exp(-E1/RT) * x1 * x2 |
//             r2 = p2 * k20 * exp(-E2/RT) * x3      |
//             dx1dt = - r1 - F/V*x1                 |  t in (0,5]
//             dx2dt = - r1 + F/V*(cBin-x2)          |
//             dx3dt = r1 - r2 - F/V*x3              |
//             dx4dt = F                            _|
//             x1(0) = cA0, x2(0) = x3(0) = 0, x4(0) = V0
//             0.8 <= p1, p2 <= 1.2
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 100;
const unsigned int NP = 2;
const unsigned int NX = 4;
double T0 = 0., TF = 3e0;

const double k10 = 4e0;		// k10         [L/mol.h]
const double E1 = 6e3;		// E1          [J/mol]
const double k20 = 8e2;		// k20         [L/h]
const double E2 = 2e4;		// E2          [J/mol]
const double R = 8.314e0;	// R           [J/mol.K]
const double cA0 = 10.;		// cA0         [mol/L]
const double cB0 = 0.;		// cB0         [mol/L]
const double cC0 = 0.;		// cC0         [mol/L]
const double V0 = 1.;		// V0          [L]
const double cBin = 20.;	// cBin        [mol/L]
const double T = 273.15+35.;	// T	       [K]
const double F = 0.3;		// F	       [L/h]

class IVP : public mc::ODESTRUCT
{
public:
  IVP()
    : mc::ODESTRUCT( NP, NX )
    {}

  template <typename TX, typename TP, typename TT>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is )
    {
      assert( ix < nx() );
      TX r1 = p[0] * k10*exp(-E1/R/T) * x[0] * x[1];
      TX r2 = p[1] * k20*exp(-E2/R/T) * x[2];
      switch( ix ){
        case 0: return -r1-F/x[3]*x[0];
        case 1: return -r1+F/x[3]*(cBin-x[1]);
        case 2: return r1-r2-F/x[3]*x[2];
        case 3: return F;
        default: throw std::runtime_error("invalid size");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return cA0;
        case 1: return cB0;
        case 2: return cC0;
	case 3: return V0;
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( REPR )
std::string example = "Repressilator";
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = p*x1*(1-x2)          |  t in (0,tf]
//             dx2dt = p*x2*(x1-1)         _|
//             x1(0) = 1.2, x2(0) = 1.1
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 5000;  // Stages?
const unsigned int NP = 3;    // Number of parameters
const unsigned int NX = 3;    // Number of states
double T0 = 0., TF = 5e2;     // initial and final time

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
         (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
         {
         assert(ix < nx() );
         switch(ix)
             {
             //using mc::Op<TX>::pow;
             using mc::pow;
             using std::pow;
             case 0: return p[0]/(1.+pow(x[2],3))-p[1]*x[0]+p[2];
             case 1: return p[0]/(1.+pow(x[0],3))-p[1]*x[1]+p[2];
             case 2: return p[0]/(1.+pow(x[1],3))-p[1]*x[2]+p[2];
             default: throw std::runtime_error("invalid size");
             }
         }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
            switch(ix)
            {
            case 0:  return 5.5;
            case 1:  return 3.5;
            case 2:  return 4.5;
            default: throw std::runtime_error("invalid size");
            }
        }
};

#elif defined( CUBIC )
std::string example = "Cubic Oscillator";
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = p*x1*(1-x2)          |  t in (0,tf]
//             dx2dt = p*x2*(x1-1)         _|
//             x1(0) = 1.2, x2(0) = 1.1
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 1000;    // Stages?
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 2;    // Number of states
double T0 = 0., TF = 200.;     // initial and final time

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
         (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
         {
         assert(ix < nx() );
         const double q = 1e-1;
         switch(ix)
             {
             using mc::sqr;
             case 0: return q*x[0]*(1.-x[0]*x[0]) + x[1]*(1.-q*x[0]*x[1]);
             case 1: return q*x[1]*(1.-x[1]*x[1]) - x[0]*(1.+q*x[0]*x[1]);
             default: throw std::runtime_error("invalid size");
             }
         }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
            switch(ix)
            {
            case 0:  return p[0];
            case 1:  return p[1];
            default: throw std::runtime_error("invalid size");
            }
        }
};

#elif defined( CUBIC2 )
std::string example = "Cubic Oscillator with Dissipative Term";
////////////////////////////////////////////////////////////////////////
//                                         _
//             dx1dt = p*x1*(1-x2)          |  t in (0,tf]
//             dx2dt = p*x2*(x1-1)         _|
//             x1(0) = 1.2, x2(0) = 1.1
//             pmin <= p <= pmax
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 2000;    // Stages?
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 2;    // Number of states
double T0 = 0., TF = 200.;     // initial and final time

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
         (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
         {
         assert(ix < nx() );
         const double q = 1e-1;
         switch(ix)
             {
             using mc::sqr;
             case 0: return q*x[0]*(1.-x[0]*x[0]) + x[1]*(1.-q*x[0]*x[1]);
             case 1: return q*x[1]*(1.-x[1]*x[1]) - x[0]*(1.+q*x[0]*x[1]) - 0.2*x[1];
             default: throw std::runtime_error("invalid size");
             }
         }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
            switch(ix)
            {
            case 0:  return p[0];
            case 1:  return p[1];
            default: throw std::runtime_error("invalid size");
            }
        }
};

#elif defined( AD )
std::string example = "Anaerobic Digester";
////////////////////////////////////////////////////////////////////////
//  ANAEROBIC DIGESTION MODEL FROM THE PAPER BY BERNARD ET AL (2001)
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 1000;  // Stages?
const unsigned int NP = 1;    // Number of parameters
const unsigned int NX = 6;    // Number of states
const unsigned int NI = 2;    // Number of invariants
double T0 = 0., TF = 100.;     // initial and final time

namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 7.1e0;	// g(COD)/L
const double mumax2  = 0.74e0;	// /day
const double KS2     = 9.28e0;	// mmol/L
const double KI2     = 256e0;	// mmol/L
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 42.14e0;	// g(COD)/g(VSS)
const double k2      = 116.5e0;	// mmol/g(VSS)
const double k3      = 268.0e0;	// mmol/g(VSS)
const double k4      =  50.6e0;	// mmol/g(VSS)
const double k5      = 343.6e0;	// mmol/g(VSS)
const double k6      = 453.0e0;	// mmol/g(VSS)
const double kVFA    = 64e-3;	// g(COD)/mmol
const double Ka      = 1.5e-2;	// mmol/L
const double Kb      = 6.5e-4;	// mmol/L
const double KH      = 16e0;	// mmol/L/atm
const double S1in    = 5e0;	// g(COD)/L
const double S2in    = 80e0;	// mmol/L
const double Cin     = 0e0;	// mmol/L
const double Zin     = 50e0;	// mmol/L
const double PT      = 1e0;	// atm
const double X10     = .5e0;	// g(VSS)/L
const double X20     = 1e0;	// g(VSS)/L
const double S10     = 1e0;	// g(COD)/L
const double S20     = 5e0;	// mmol/L
const double C0      = 40e0;	// mmol/L
const double Z0      = 50e0;	// mmol/L
}
using namespace AD_param;

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
        (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
        {
          assert(ix < nx() );
          TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
          TX phi  = x[5] + x[3] - x[4] + KH * PT + k6 / kLa * mu2 * x[1];
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * ( x[5] + x[3] - x[4] ) ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( x[5] + x[3] - x[4] - KH * PCO2 );
          switch( ix ){
            case 0: return ( mu1 - alpha * p[0] ) * x[0];
            case 1: return ( mu2 - alpha * p[0] ) * x[1];
            case 2: return p[0] * ( S1in - x[2] ) - k1 * mu1 * x[0];
            case 3: return p[0] * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * mu2 * x[1];
            case 4: return p[0] * ( Zin - x[4] );
            case 5: return p[0] * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
            default: throw std::runtime_error("invalid size");
          }
        }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
          switch(ix){
            case 0: return X10;
            case 1: return X20;
            case 2: return S10;
            case 3: return S20;
            case 4: return Z0;
            case 5: return C0;
            default: throw std::runtime_error("invalid size");
          }
        }
};

#elif defined( AD2 )
std::string example = "Anaerobic Digester";
////////////////////////////////////////////////////////////////////////
//  ANAEROBIC DIGESTION MODEL FROM THE PAPER BY BERNARD ET AL (2001)
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 4000;  // Stages?
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 6;    // Number of states
double T0 = 0., TF = 40.;     // initial and final time

namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 7.1e0;	// g(COD)/L
const double mumax2  = 0.74e0;	// /day
const double KS2     = 9.28e0;	// mmol/L
const double KI2     = 256e0;	// mmol/L
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 42.14e0;	// g(COD)/g(VSS)
const double k2      = 116.5e0;	// mmol/g(VSS)
const double k3      = 268.0e0;	// mmol/g(VSS)
const double k4      =  50.6e0;	// mmol/g(VSS)
const double k5      = 343.6e0;	// mmol/g(VSS)
const double k6      = 453.0e0;	// mmol/g(VSS)
const double kVFA    = 64e-3;	// g(COD)/mmol
const double Ka      = 1.5e-2;	// mmol/L
const double Kb      = 6.5e-4;	// mmol/L
const double KH      = 16e0;	// mmol/L/atm
const double S1in    = 5e0;	// g(COD)/L
const double S2in    = 80e0;	// mmol/L
const double Cin     = 0e0;	// mmol/L
const double Zin     = 50e0;	// mmol/L
const double PT      = 1e0;	// atm
const double X10     = .5e0;	// g(VSS)/L
const double X20     = 1e0;	// g(VSS)/L
const double S10     = 1e0;	// g(COD)/L
const double S20     = 5e0;	// mmol/L
const double C0      = 40e0;	// mmol/L
const double Z0      = 50e0;	// mmol/L
const double D       = 0.4e0;	// /day
}
using namespace AD_param;

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
        (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
        {
          assert(ix < nx() );
          TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
          TX phi  = x[5] + x[3] - x[4] + KH * PT + k6 / kLa * mu2 * x[1];
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * ( x[5] + x[3] - x[4] ) ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( x[5] + x[3] - x[4] - KH * PCO2 );
          switch( ix ){
            case 0: return ( mu1 - alpha * D ) * x[0];
            case 1: return ( mu2 - alpha * D ) * x[1];
            case 2: return D * ( S1in - x[2] ) - k1 * mu1 * x[0];
            case 3: return D * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * mu2 * x[1];
            case 4: return D * ( Zin - x[4] );
            case 5: return D * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
            default: throw std::runtime_error("invalid size");
          }
        }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
          switch(ix){
            case 0: return X10 * p[0];
            case 1: return X20 * p[1];
            case 2: return S10;
            case 3: return S20;
            case 4: return Z0;
            case 5: return C0;
            default: throw std::runtime_error("invalid size");
          }
        }
};

#elif defined( AD3 )
std::string example = "Anaerobic Digester";
////////////////////////////////////////////////////////////////////////
//  ANAEROBIC DIGESTION MODEL FROM THE PAPER BY BERNARD ET AL (2001)
////////////////////////////////////////////////////////////////////////
const unsigned int NS = 1;  // Stages?
const unsigned int NP = 3;    // Number of parameters
const unsigned int NX = 6;    // Number of states
double T0 = 0., TF = 100.;     // initial and final time

namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 7.1e0;	// g(COD)/L
const double mumax2  = 0.74e0;	// /day
const double KS2     = 9.28e0;	// mmol/L
const double KI2     = 256e0;	// mmol/L
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 42.14e0;	// g(COD)/g(VSS)
const double k2      = 116.5e0;	// mmol/g(VSS)
const double k3      = 268.0e0;	// mmol/g(VSS)
const double k4      =  50.6e0;	// mmol/g(VSS)
const double k5      = 343.6e0;	// mmol/g(VSS)
const double k6      = 453.0e0;	// mmol/g(VSS)
const double KH      = 16e0;	// mmol/L/atm
const double PT      = 1e0;	// atm
const double S1in    = 5e0;	// g(COD)/L
const double S2in    = 80e0;	// mmol/L
const double Cin     = 0e0;	// mmol/L
const double Zin     = 50e0;	// mmol/L
const double X10     = .5e0;	// g(VSS)/L
const double X20     = 1e0;	// g(VSS)/L
const double S10     = 1e0;	// g(COD)/L
const double S20     = 5e0;	// mmol/L
const double C0      = 40e0;	// mmol/L
const double Z0      = 50e0;	// mmol/L
const double D       = 0.4e0;	// /day
}
using namespace AD_param;

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
        (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
        {
          assert(ix < nx() );
          TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
          TX phi  = x[5] + x[3] - x[4] + KH * PT + k6 / kLa * mu2 * x[1];
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * ( x[5] + x[3] - x[4] ) ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( x[5] + x[3] - x[4] - KH * PCO2 );
          switch( ix ){
            case 0: return ( mu1 - alpha * D ) * x[0];
            case 1: return ( mu2 - alpha * D ) * x[1];
            case 2: return D * ( S1in - x[2] ) - k1 * mu1 * x[0];
            case 3: return D * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * mu2 * x[1];
            case 4: return D * ( Zin - x[4] );
            case 5: return D * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
            default: throw std::runtime_error("invalid size");
          }
        }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
          switch(ix){
            case 0: return X10 * p[0];
            case 1: return X20 * p[1];
            case 2: return S10;
            case 3: return S20;
            case 4: return Z0;
            case 5: return C0 * p[2];
            default: throw std::runtime_error("invalid size");
          }
        }
};

#elif defined( AM2 )
std::string example = "Anaerobic Digester Model AM2";
////////////////////////////////////////////////////////////////////////
//  MODIFIED ANAEROBIC DIGESTION MODEL AM2 DEVELOPPED IN THE SCOPE
//  OF THE EUROPEAN PROJECT TELEMAC (DELIVERABLE 3.2)
////////////////////////////////////////////////////////////////////////
namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 1e0;	// g(COD)/L       [7.1]
const double mumax2  = 1.169e0;	// /day           [0.74e0]
const double KS2     = 50e0;	// mmol/L         [9.28]
const double KI2     = 1e10;	// mmol/L         [256.]
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 16.342e0;// g(COD)/g(VSS)  [27.4]
const double k2      = 47.3e0;	// mmol/g(VSS)
const double k3      = 86.496e0;// mmol/g(VSS)    [83.1]
const double k4      = 57.5e0;	// mmol/g(VSS)
const double k5      = 125.9e0;	// mmol/g(VSS)
const double k6      = 154.2e0;	// mmol/g(VSS)
const double KH      = 16e0;	// mmol/atm
const double PT      = 1.03e0;	// atm
const double X10     = 1.3441;	// g(VSS)/L
const double X20     = 2.9595;	// g(VSS)/L
const double S10     = 0.85;	// g(COD)/L   [0.12548]
const double S20     = 6.5;	// mmol/L     [6.4618]
const double C0      = 95.;	// mmol/L     [99.368]
const double Z0      = 90.;	// mmol/L     [97.845]
const double V       = 720;     // L          [717.3]
const double kVFA    = 64e-3;	// g(COD)/mmol
const double Kb      = 6.5e-4;	// mmol/L
}
using namespace AD_param;

const unsigned int NT = 4;    // Number of time stages
double T[NT+1]   = { 0., 0.266, 1.266, 2.362, 4. }; // Time stages [day]
double Qin[NT]   = { 8., 47., 47., 9. };            // Influent flowrate [L/h]
double CODin[NT] = { 18.5, 18.5, 37., 18.5 };       // Influent COD [g(COD)/L]
double VFAin[NT] = { 103., 103., 206., 103. };      // Influent VFA [mmol/L]
double TICin[NT] = { 9., 9., 18., 9. };             // Influent TIC [mmol/L]
double ALKin[NT] = { 98., 98., 195., 98. };         // Influent alk [mmol/L]
//double ALKin[NT] = { 90., 90., 180., 90. };         // Influent alk [mmol/L]
const unsigned int NS = NT*20;  // Total number of stages
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 6;    // Number of states

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
        (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
        {
          assert(ix < nx() );
          // Influent values
          unsigned int it = is*NT/NS;
	  double D = Qin[it]*24./V;
	  double S1in = CODin[it]-VFAin[it]*kVFA;
	  double S2in = VFAin[it];
	  double Zin = ALKin[it];
	  double Cin = TICin[it];
	  //std::cout << t << "  " << it << "  " << D << "  " << S1in << std::endl;
	  // Transfers and reactions
          TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
	  TX CO2  = x[5] + x[3] - x[4];
          TX qCH4 = k6 * mu2 * x[1];
          TX phi  = CO2 + KH * PT + qCH4 / kLa;
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * CO2 ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( CO2 - KH * PCO2 );
          switch( ix ){
            case 0: return ( mu1 - alpha * D ) * x[0];
            case 1: return ( mu2 - alpha * D ) * x[1];
            case 2: return D * ( S1in - x[2] ) - k1 * p[0] * mu1 * x[0];
            case 3: return D * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * p[1] * mu2 * x[1];
            case 4: return D * ( Zin - x[4] );
            case 5: return D * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
            default: throw std::runtime_error("invalid size");
          }
        }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
          switch(ix){
            case 0: return X10;
            case 1: return X20;
            case 2: return S10;
            case 3: return S20;
            case 4: return Z0;
            case 5: return C0;
            default: throw std::runtime_error("invalid size");
          }
        }
};

#elif defined( AM2b )
std::string example = "Anaerobic Digester Model AM2";
////////////////////////////////////////////////////////////////////////
//  MODIFIED ANAEROBIC DIGESTION MODEL AM2 DEVELOPPED IN THE SCOPE
//  OF THE EUROPEAN PROJECT TELEMAC (DELIVERABLE 3.2)
////////////////////////////////////////////////////////////////////////
namespace AD_param{
const double mumax1  = 1.2e0;	// /day
const double KS1     = 7.1e0;	// g(COD)/L
const double mumax2  = 0.74e0;	// /day
const double KS2     = 9.28e0;	// mmol/L
const double KI2     = 256e0;	// mmol/L
const double alpha   = 0.5;	// -
const double kLa     = 19.8e0;	// /day
const double k1      = 42.14e0;	// g(COD)/g(VSS)
const double k2      = 116.5e0;	// mmol/g(VSS)
const double k3      = 268.0e0;	// mmol/g(VSS)
const double k4      =  50.6e0;	// mmol/g(VSS)
const double k5      = 343.6e0;	// mmol/g(VSS)
const double k6      = 453.0e0;	// mmol/g(VSS)
const double kVFA    = 64e-3;	// g(COD)/mmol
const double Ka      = 1.5e-2;	// mmol/L
const double Kb      = 6.5e-4;	// mmol/L
const double KH      = 16e0;	// mmol/L/atm
const double PT      = 1e0;	// atm
const double X10     = .5e0;	// g(VSS)/L
const double X20     = 1e0;	// g(VSS)/L
const double S10     = 1e0;	// g(COD)/L
const double S20     = 5e0;	// mmol/L
const double C0      = 40e0;	// mmol/L
const double Z0      = 50e0;	// mmol/L
}
using namespace AD_param;

const unsigned int NT = 4;    // Number of time stages
double T[NT+1]   = { 0., 1., 2., 3., 4. };          // Time stages [day]
double Din[NT]   = { 0.25, 1., 1., 0.25 };          // Dilution rate [/day]
double CODin[NT] = { 7.5, 7.5, 15., 7.5 };             // Influent COD [g(COD)/L]
double VFAin[NT] = { 80., 80., 160., 80. };         // Influent VFA [mmol/L]
double TICin[NT] = { 5., 5., 10., 5. };             // Influent TIC [mmol/L]
double ALKin[NT] = { 50., 50., 100., 50. };         // Influent alk [mmol/L]
const unsigned int NS = NT*6;  // Total number of stages
const unsigned int NP = 2;    // Number of parameters
const unsigned int NX = 6;    // Number of states

class IVP : public mc::ODESTRUCT
{
    public:
        IVP(): mc::ODESTRUCT( NP, NX )
        {}

        template <typename TX, typename TP, typename TT>
        TX RHS
        (const unsigned int ix, const TP*p, const TX*x, const TT&t, const unsigned int is)
        {
          assert(ix < nx() );
          // Influent values
          unsigned int it = is*NT/NS;
	  double D = Din[it];
	  double S1in = CODin[it]-VFAin[it]*kVFA;
	  double S2in = VFAin[it];
	  double Zin = ALKin[it];
	  double Cin = TICin[it];
	  //std::cout << t << "  " << it << "  " << D << "  " << S1in << std::endl;
	  // Transfers and reactions
          TX mu1  = mumax1 * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
	  TX CO2  = x[5] + x[3] - x[4];
          TX qCH4 = k6 * mu2 * x[1];
          TX phi  = CO2 + KH * PT + qCH4 / kLa;
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * CO2 ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( CO2 - KH * PCO2 ) * p[1];
          switch( ix ){
            case 0: return ( mu1 - alpha * D * p[0] ) * x[0];
            case 1: return ( mu2 - alpha * D * p[0] ) * x[1];
            case 2: return D * ( S1in - x[2] ) - k1 * mu1 * x[0];
            case 3: return D * ( S2in - x[3] ) + k2 * mu1 * x[0] - k3 * mu2 * x[1];
            case 4: return D * ( Zin - x[4] );
            case 5: return D * ( Cin - x[5] ) + k4 * mu1 * x[0] + k5 * mu2 * x[1] - qCO2;
            default: throw std::runtime_error("invalid size");
          }
        }

        template <typename T>
        T IC
        ( const unsigned int ix, const T*p )
        {
          switch(ix){
            case 0: return X10;
            case 1: return X20;
            case 2: return S10;
            case 3: return S20;
            case 4: return Z0;
            case 5: return C0;
            default: throw std::runtime_error("invalid size");
          }
        }
};
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  // Time stages
  double tk[NS+1];
#if defined( AM2 ) || defined( AM2b )
  for( unsigned int it=0, k=0; it<NT; it++ ){
    tk[k] = T[it]; k++;
    for( unsigned int is=1; is<NS/NT; is++, k++ )
      tk[k] = tk[k-1]+(T[it+1]-T[it])*(double)NT/(double)NS;
  }
  tk[NS] = T[NT];
  for( unsigned int is=0; is<=NS; is++ )
    std::cout << tk[is] << std::endl;
#else
  tk[0] = T0;
  for( unsigned int is=0; is<NS; is++ )
    tk[1+is] = tk[is]+(TF-T0)/(double)NS;
  
#endif

  I Ip[NP];
#if defined( LV )
  Ip[0] = I( 2.95, 3.05 );
#elif defined( OSCIL )
  Ip[0] = I( -0.1, 0.1 );
#elif defined( LTV )
  Ip[0] = I( 0.1, 2. );
  Ip[1] = I( -1., 1. );
  for( unsigned int is=0; is<NS; is++ )
    Ip[2+is] = I( -0.5, 1.5 );
#elif defined( COMPART )
  for( unsigned int ip=0; ip<NP; ip++ )
    Ip[ip] = I( 4e-1, .5e0 );
#elif defined( RXN1 )
  //Ip[0] = I( 5e1, 5e2 );
  Ip[0] = I( 5e1, 6e1 );
  //Ip[1] = I( 1.0, 1.0 );
  Ip[1] = I( 0.95, 1.05 );
#elif defined( RXN1B )
  //Ip[0] = I( 5e1, 5e2 );
  Ip[0] = I( 5e1, 15e1 );
  //Ip[1] = I( 1.0, 1.0 );
  Ip[1] = I( 0.95, 1.05 );
#elif defined( RXN1C )
  //Ip[0] = I( 5e1, 5e2 );
  Ip[0] = I( 5e1, 6e1 );
  //Ip[1] = I( 1.0, 1.0 );
  Ip[1] = I( 0.95, 1.05 );
#elif defined( RXN2 )
  Ip[0] = Ip[1] = I( 0.8, 1.2 );
#elif defined( REPR )
  Ip[0] = I(215.45, 215.55);
  Ip[1] = I(0.999,1.001);
  Ip[2] = I(1.499,1.501);
  //Ip[0] = I(215., 216.);
  //Ip[1] = I(0.995,1.005);
  //Ip[2] = I(1.495,1.505);
#elif defined( CUBIC )
  //Ip[0] = I(1.99,2.01);
  //Ip[1] = I(-.01,.01);
  Ip[0] = I(1.5,2.5);
  Ip[1] = I(-.1,.1);
#elif defined( CUBIC2 )
  //Ip[0] = I(1.99,2.01);
  //Ip[1] = I(-.01,.01);
  Ip[0] = I(1.5,2.5);
  Ip[1] = I(-.1,.1);
#elif defined( AD )
  //Ip[0] = I(0.3,0.4);
  Ip[0] = I(0.4,0.401);
#elif defined( AD2 )
  Ip[0] = I(0.98,1.02);
  Ip[1] = I(0.98,1.02);
#elif defined( AD3 )
  Ip[0] = I(0.98,1.02);
  Ip[1] = I(0.98,1.02);
  Ip[2] = I(0.98,1.02);
  //Ip[0] = I(0.94,1.06);
  //Ip[1] = I(0.94,1.06);
  //Ip[2] = I(0.94,1.06);
  //Ip[0] = I(0.96,1.04);
  //Ip[1] = I(0.96,1.04);
  //Ip[2] = I(0.96,1.04);
#elif defined( AM2 ) || defined( AM2b )
  Ip[0] = I(0.999,1.001);
  Ip[1] = I(0.999,1.001);
#endif

  // Parameters / state variables
  I *Ixk[NS+1];
  E Exk[NS+1];
  double *Hxk[NS+1];
  for( unsigned int is=0; is<=NS; is++ ){
    Ixk[is] = new I[NX];
    Hxk[is] = new double[NX];
  }

  TM TMenv( NP, NTM );
  TMenv.options.SCALE_VARIABLES = false;
  TMenv.options.CENTER_REMAINDER = true;
  TV TMp[NP];
  for( unsigned int ip=0; ip<NP; ip++ ){
    TMp[ip].set( &TMenv, ip, Ip[ip] );
  }
  TV *TMxk[NS+1];
  E RExk[NS+1];
  double *HBxk[NS+1], *HRxk[NS+1];
  for( unsigned int is=0; is<=NS; is++ ){
    TMxk[is]  = new TV[NX];
    HBxk[is] = new double[NX];
    HRxk[is] = new double[NX];
  }

  //////////////////////////////////////////////////////////////////////
  // USING: ODESLV_GSL<I,IVP>

  mc::ODESLV_GSL<I,IVP>* pSLVGSL = new mc::ODESLV_GSL<I,IVP>();
  pSLVGSL->ODESLV_GSL<I,IVP>::options.DISPLAY = 1;
  pSLVGSL->ODESLV_GSL<I,IVP>::options.RTOL    = 1e-12;
  pSLVGSL->ODESLV_GSL<I,IVP>::options.ATOL    = 1e-12;
  pSLVGSL->ODESLV_GSL<I,IVP>::options.HMIN    = 1e-10;
  pSLVGSL->ODESLV_GSL<I,IVP>::options.H0      = 1e-3;
  pSLVGSL->ODESLV_GSL<I,IVP>::options.INTMETH = mc::ODEBND_GSL<I,IVP>::Options::RK8PD; //RKF45;
#if defined( SAVE_RESULTS )
  pSLVGSL->ODESLV_GSL<I,IVP>::options.RESRECORD = true;
#endif

  // Test: Sampled bounds
  std::cout << "\nCOMPUTE APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  pSLVGSL->bounds( NS, tk, Ip, Ixk, NSAMP );   
#if defined( SAVE_RESULTS )
  std::ofstream bndrec( "exact.out", std::ios_base::out );
  pSLVGSL->record( bndrec );
#endif

  delete pSLVGSL;

  //////////////////////////////////////////////////////////////////////
  // TESTING: ODEBND_VAL<I,IVP>
/*
  mc::ODEBND_VAL<I,IVP>* pBNDVAL = new mc::ODEBND_VAL<I,IVP>();
  pBNDVAL->ODEBND_VAL<I,IVP>::options.DISPLAY = 1;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.WRAPMIT = mc::ODEBND_VAL<I,IVP>::Options::ELLIPS;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.SCALING = false;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.TOL     = 1e-6;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.ATOL    = 1e-6;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.QTOL    = 1e-10;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.TSORDER = 5;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.HREDUC  = 0.8;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.HMIN    = 1e-6;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.HMAX    = 1e5;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.TSTOP   = false;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.ORDMIT  = 1;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.USEINV  = false;
  pBNDVAL->ODEBND_VAL<I,IVP>::options.USENBND = false;
#if defined( SAVE_RESULTS )
  pBNDVAL->ODEBND_VAL<I,IVP>::options.RESRECORD = true;
#endif

  // Test: Validated integration
  std::cout << "\nCOMPUTE INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  pBNDVAL->bounds( NS, tk, Ip, Ixk, Exk );  
  std::cout << " FINAL TIME\t" << pBNDVAL->final_time() << std::endl;
#if defined( SAVE_RESULTS )
  std::ostringstream namef4;
  switch( pBNDVAL->ODEBND_VAL<I,IVP>::options.WRAPMIT ){
    case mc::ODEBND_VAL<I,IVP>::Options::NONE:   namef4 << "valint_INT"; break;
    case mc::ODEBND_VAL<I,IVP>::Options::ELLIPS: namef4 << "valint_ELL"; break;
  }
  if( pBNDVAL->ODEBND_VAL<I,IVP>::options.USEINV ) namef4 << "_inv.out"; 
  namef4 << ".out";
  std::ofstream bndrec4( namef4.str().c_str(), std::ios_base::out );
  pBNDVAL->record( bndrec4 );
#endif

  //std::cout << std::endl;
  //pBNDVAL->hausdorff( NS, tk, Ip, Hxk, NSAMP, pSLVGSL->options );

  // Test: Validated integration w/ Taylor models
  std::cout << "\nCOMPUTE TAYLOR MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  //for( unsigned int it=0; it<100; it++ )
    pBNDVAL->bounds( NS, tk, TMp, E(), TMxk, RExk );
  std::cout << " FINAL TIME\t" << pBNDVAL->final_time() << std::endl;
#if defined( SAVE_RESULTS )
  std::ostringstream namef5;
  switch( pBNDVAL->ODEBND_VAL<I,IVP>::options.WRAPMIT ){
    case mc::ODEBND_VAL<I,IVP>::Options::NONE:   namef5 << "valint_TM" << NTM << "+INT"; break;
    case mc::ODEBND_VAL<I,IVP>::Options::ELLIPS: namef5 << "valint_TM" << NTM << "+ELL"; break;
  }
  if( pBNDVAL->ODEBND_VAL<I,IVP>::options.USEINV ) namef5 << "_inv.out"; 
  namef5 << ".out";
  std::ofstream bndrec5( namef5.str().c_str(), std::ios_base::out );
  pBNDVAL->record( bndrec5 );
#endif

  //std::cout << std::endl;
  //pBNDGSL->hausdorff( NS, tk, TMp, HBxk, HRxk, NSAMP );

  delete pBNDVAL;
*/
  //////////////////////////////////////////////////////////////////////
  // TESTING: ODEBND_GSL<I,IVP>
/*
  mc::ODEBND_GSL<I,IVP>* pBNDGSL = new mc::ODEBND_GSL<I,IVP>();
  pBNDGSL->ODEBND_GSL<I,IVP>::options.DISPLAY = 1;
  pBNDGSL->ODEBND_GSL<I,IVP>::options.RTOL    = 1e-6;
  pBNDGSL->ODEBND_GSL<I,IVP>::options.ATOL    = 1e-8;
  pBNDGSL->ODEBND_GSL<I,IVP>::options.H0      = 1e-9;
  //pBNDGSL->ODEBND_GSL<I,IVP>::options.QTOL    = 1e-10;
  pBNDGSL->ODEBND_GSL<I,IVP>::options.INTMETH = mc::ODEBND_GSL<I,IVP>::Options::RKF45;
  pBNDGSL->ODEBND_GSL<I,IVP>::options.WRAPMIT = mc::ODEBND_GSL<I,IVP>::Options::ELLIPS;//DINEQ
  pBNDGSL->ODEBND_GSL<I,IVP>::options.ORDMIT = 1;
#if defined( SAVE_RESULTS )
  pBNDGSL->ODEBND_GSL<I,IVP>::options.RESRECORD = true;
#endif

  // Test: Differential Inequalities
  std::cout << "\nCOMPUTE INTERVAL ENCLOSURE OF REACHABLE SET:\n\n";
  for( unsigned int it=0; it<1; it++ )
  pBNDGSL->bounds( NS, tk, Ip, Ixk );  
  std::cout << " FINAL TIME\t" << pBNDGSL->final_time() << std::endl;
#if defined( SAVE_RESULTS )
  std::ostringstream namef2;
  switch( pBNDGSL->ODEBND_GSL<I,IVP>::options.WRAPMIT ){
    case mc::ODEBND_GSL<I,IVP>::Options::NONE:   namef2 << "gslint_INT"; break;
    case mc::ODEBND_GSL<I,IVP>::Options::DINEQ:  namef2 << "gslint_DI";  break;
    case mc::ODEBND_GSL<I,IVP>::Options::ELLIPS: namef2 << "gslint_ELL" << pBNDGSL->ODEBND_GSL<I,IVP>::options.ORDMIT; break;
    default: namef2 << "gslint"; break;
  }
  namef2 << ".out";
  std::ofstream bndrec2( namef2.str().c_str(), std::ios_base::out );
  pBNDGSL->record( bndrec2 );
#endif  

  //std::cout << std::endl;
  //pBNDGSL->hausdorff( NS, tk, Ip, Hxk, NSAMP );

  // Test: Differential Inequalities with Taylor models
  std::cout << "\nCOMPUTE TAYLOR MODEL ENCLOSURE OF REACHABLE SET:\n\n";
  //for( unsigned int it=0; it<5; it++ )
  pBNDGSL->bounds( NS, tk, TMp, TMxk );
  std::cout <<  " FINAL TIME\t" << pBNDGSL->final_time() << std::endl;
#if defined( SAVE_RESULTS )
  std::ostringstream namef3;
  switch( pBNDGSL->ODEBND_GSL<I,IVP>::options.WRAPMIT ){
    case mc::ODEBND_GSL<I,IVP>::Options::NONE:   namef3 << "gslint_TM" << NTM << "+INT"; break;
    case mc::ODEBND_GSL<I,IVP>::Options::DINEQ:  namef3 << "gslint_TM" << NTM << "+DI";  break;
    case mc::ODEBND_GSL<I,IVP>::Options::ELLIPS: namef3 << "gslint_TM" << NTM << "+ELL" << pBNDGSL->ODEBND_GSL<I,IVP>::options.ORDMIT; break;
    default: namef3 << "gslint_TM" << NTM; break;
  }
  namef3 << ".out";
  std::ofstream bndrec3( namef3.str().c_str(), std::ios_base::out );
  pBNDGSL->record( bndrec3 );
#endif

  //std::cout << std::endl;
  //pBNDGSL->hausdorff( NS, tk, TMp, HBxk, HRxk, NSAMP );

  delete pBNDGSL;
*/
  for( unsigned int is=0; is<=NS; is++ ){
    delete[] Ixk[is];
    delete[] Hxk[is];
    delete[] HBxk[is];
    delete[] HRxk[is];
    delete[] TMxk[is];
  }

  return 0;
}
