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

const unsigned int NTM   = 3;	// <- Order of Taylor model
const unsigned int NSAMP = 40;	// <- Number of sampling points
#define SAVE_RESULTS			// <- Saving the results to file

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
          TX mu1  = mumax1 * p[0] * x[2] / ( x[2] + KS1 );
          TX mu2  = mumax2 * p[1] * x[3] / ( x[3] + KS2 + x[3]*x[3] / KI2 );
	  TX CO2  = x[5] + x[3] - x[4];
          TX qCH4 = k6 * mu2 * x[1];
          TX phi  = CO2 + KH * PT + qCH4 / kLa;
          TX PCO2 = ( phi - sqrt( phi*phi - 4e0 * KH * PT * CO2 ) ) / ( 2e0 * KH );
          TX qCO2 = kLa * ( CO2 - KH * PCO2 );
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

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  // Time stages
  double tk[NS+1];
  for( unsigned int it=0, k=0; it<NT; it++ ){
    tk[k] = T[it]; k++;
    for( unsigned int is=1; is<NS/NT; is++, k++ )
      tk[k] = tk[k-1]+(T[it+1]-T[it])*(double)NT/(double)NS;
  }
  tk[NS] = T[NT];
  for( unsigned int is=0; is<=NS; is++ )
    std::cout << tk[is] << std::endl;

  I Ip[NP];
  Ip[0] = I(0.98,1.02);
  Ip[1] = I(0.98,1.02);

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

  for( unsigned int is=0; is<=NS; is++ ){
    delete[] Ixk[is];
    delete[] Hxk[is];
    delete[] HBxk[is];
    delete[] HRxk[is];
    delete[] TMxk[is];
  }

  return 0;
}
