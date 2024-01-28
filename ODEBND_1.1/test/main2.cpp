#undef USE_PROFIL   // <- change to #define for selecting PROFIL
#undef USE_FILIB    // <- change to #define for selecting FILIB++

#include <fstream>
#include <iomanip>
#include "odeslv_gsl.hpp"

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

////////////////////////////////////////////////////////////////////////

#define REPR			// <- select the example here
const unsigned int NSAMP = 5;	// <- Number of sampling points
#define SAVE_RESULTS		// <- Saving the results to file

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

const unsigned int NS = 1;
const unsigned int NP = 1;
const unsigned int NX = 2;
double T0 = 0., TF = 2.;

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
        case 0: return p[0]*x[0]*(1.-x[1]);
        case 1: return p[0]*x[1]*(x[0]-1.);
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.2;
        case 1: return 1.1;
        default: throw std::runtime_error("invalid index");
      }
    }
};

#elif defined( OSCIL)
std::string example = "OSCIL";
////////////////////////////////////////////////////////////////////////
//                                    _
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
        case 0: return p[0];
        case 1: return 0.;
        default: throw std::runtime_error("invalid index");
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

const unsigned int NS = 100000;
const unsigned int NP = 2;
const unsigned int NX = 4;
double T0 = 0., TF = 5e3;

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
const unsigned int NS = 100;  // Stages?
const unsigned int NP = 3;    // Number of parameters
const unsigned int NX = 3;    // Number of states
double T0 = 0., TF = 300.;     // initial and final time

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
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  // Time stages
  double tk[NS+1];
  tk[0] = T0;
  for( unsigned int is=0; is<NS; is++ )
    tk[1+is] = tk[is]+(TF-T0)/(double)NS;

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
#elif defined( RXN1 )
  Ip[0] = I( 5e1, 5.01e1 );
  //Ip[1] = I( 1.0, 1.0 );
  //Ip[0] = I( 5e1, 5e2 );
  Ip[1] = I( 0.995, 1.005 );
#elif defined( RXN2 )
  Ip[0] = Ip[1] = I( 0.8, 1.2 );
#elif defined( REPR )
  Ip[0] = I(215.5, 216.5);
  Ip[1] = I(0.995,1.005);
  Ip[2] = I(1.495,1.505);
#endif

  //////////////////////////////////////////////////////////////////////
  // TESTING: ODESLV_GSL<I,IVP>

  mc::ODESLV_GSL<I,IVP>* pIVP = new mc::ODESLV_GSL<I,IVP>();
  pIVP->ODESLV_GSL<I,IVP>::options.DISPLAY = 1;
  pIVP->ODESLV_GSL<I,IVP>::options.RTOL    = 1e-8;
  pIVP->ODESLV_GSL<I,IVP>::options.ATOL    = 1e-8;
  pIVP->ODESLV_GSL<I,IVP>::options.H0      = 1e-3;
  pIVP->ODESLV_GSL<I,IVP>::options.INTMETH = mc::ODESLV_GSL<I,IVP>::Options::RKF45;
#if defined( SAVE_RESULTS )
  pIVP->ODESLV_GSL<I,IVP>::options.RESRECORD = true;
#endif

  I *Ixk[NS+1];
  for( unsigned int is=0; is<=NS; is++ )
    Ixk[is]  = new I[NX];

  // Test: Sampled bounds
  std::cout << "\nCOMPUTE APPROXIMATE ENCLOSURE OF REACHABLE SET:\n\n";
  pIVP->bounds( NS, tk, Ip, Ixk, NSAMP );   
#if defined( SAVE_RESULTS )
  std::ofstream bndrec( "exact.out", std::ios_base::out );
  pIVP->record( bndrec );
#endif  

  for( unsigned int is=0; is<=NS; is++ )
    delete[] Ixk[is];

  delete pIVP;

  return 0;
}
