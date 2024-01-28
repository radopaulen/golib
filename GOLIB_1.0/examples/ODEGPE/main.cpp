//#define USE_PROFIL
//#define USE_FILIB

#include <fstream>
#include <iomanip>
#include "odegpe.hpp"
#include "main.hpp"

#define PE1			// <- select the example here

#if defined( PE1 )
////////////////////////////////////////////////////////////////////////
// PROBLEM FROM KIEFFER & WALTER (2011)
////////////////////////////////////////////////////////////////////////

// Number of time stages
const unsigned int NS = 15;
// Number of model parameters
const unsigned int NP = 3;
// Number of state variables
const unsigned int NX = 2;
// Number of output variables
const unsigned int NY = 1;
// Number of constraints
const unsigned int NC = 0;

class PE : public virtual mc::PESTRUCT
{
public:
  PE()
    : mc::PESTRUCT( NP, NX, NY, NC )
    {}

  // ODEs right-hand side
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

  // Initial conditions
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

  // Initial conditions
  template <typename TX, typename TP>
  TX OUT
    ( const unsigned int iy, const TP*p, const TX*x, const unsigned int is )
    {
      assert( iy < ny() );
      switch( iy ){
	case 0: return x[1];
	default: throw std::runtime_error("invalid size");
      }
    }

  // Constraints
  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*x, const unsigned int ns )
    {
      throw std::runtime_error("invalid size");
    }
};

#elif defined( PE1b )
////////////////////////////////////////////////////////////////////////
// PROBLEM FROM KIEFFER & WALTER (2011)
////////////////////////////////////////////////////////////////////////

// Number of time stages
const unsigned int NS = 15;
// Number of model parameters
const unsigned int NP = 4;
// Number of state variables
const unsigned int NX = 2;
// Number of output variables
const unsigned int NY = 1;
// Number of constraints
const unsigned int NC = 0;

class PE : public virtual mc::PESTRUCT
{
public:
  PE()
    : mc::PESTRUCT( NP, NX, NY, NC )
    {}

  // ODEs right-hand side
  template <typename TX, typename TP>
  TX RHS
    ( const unsigned int ix, const TP*p, const TX*x, const unsigned int is )
    {
      assert( ix < nx() );
      switch( ix ){
	case 0: return -(p[0]+p[2])*x[0]+p[1]*x[1];
	case 1: return p[0]*x[0]-p[1]*x[1];
	default: throw std::runtime_error("invalid size");
      }
    }

  // Initial conditions
  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
	case 0: return 1.;
	case 1: return p[3];
	default: throw std::runtime_error("invalid size");
      }
    }

  // Initial conditions
  template <typename TX, typename TP>
  TX OUT
    ( const unsigned int iy, const TP*p, const TX*x, const unsigned int is )
    {
      assert( iy < ny() );
      switch( iy ){
	case 0: return x[1];
	default: throw std::runtime_error("invalid size");
      }
    }

  // Constraints
  template <typename T>
  std::pair<T,t_CTR> CTR
    ( const unsigned int ic, const T*p, T* const*x, const unsigned int ns )
    {
      throw std::runtime_error("invalid size");
    }
};

#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
#if defined( PE1 )
  I P[NP];
  for( unsigned int ip=0; ip<NP; ip++ )
    P[ip] = I( 1e-2, 1e0 );

  double tk[NS+1];
  for( unsigned int is=0; is<NS+1; is++ )
    tk[is] = is;

  std::list< mc::ODEGPE<I,PE>::Data > Y;
  double yk[NS] = { .36, .47, .49, .48, .46, .44, .42, .4, .38, .36, .35, .33, .31, .3, .28 };
  double e = 5e-3;
  for( unsigned int is=0; is<NS; is++ )
    Y.push_back( mc::ODEGPE<I,PE>::Data( is+1, 0, yk[is]-e, yk[is]+e ) );

#elif defined( PE1b )
  I P[NP];
  for( unsigned int ip=0; ip<3; ip++ )
    P[ip] = I( 1e-2, 1e0 );
  P[3] = I( 0e0, 1e-1 );

  double tk[NS+1];
  for( unsigned int is=0; is<NS+1; is++ )
    tk[is] = is;

  std::list< mc::ODEGPE<I,PE>::Data > Y;
  double yk[NS] = { .36, .47, .49, .48, .46, .44, .42, .4, .38, .36, .35, .33, .31, .3, .28 };
  double e = 5e-3;
  for( unsigned int is=0; is<NS; is++ )
    Y.push_back( mc::ODEGPE<I,PE>::Data( is+1, 0, yk[is]-e, yk[is]+e ) );
#endif

  mc::ODEGPE<I,PE> GPE;

  GPE.options.TAYLOR_MODEL_ORDER          = 4;
  GPE.options.IVP_BOUNDING_TYPE           = mc::ODEGPE<I,PE>::Options::TM;//IA
  GPE.options.IVP_BOUNDING_PROPAGATION    = mc::ODEGPE<I,PE>::Options::DINEQ;//VERIF;
  GPE.options.USE_CONSTRAINT_PROPAGATION  = false;	//try false/true
  GPE.options.USE_DOMAIN_REDUCTION        = true;	//try false/true
  GPE.options.DOMAIN_REDUCTION_MAXLOOPS   = 10;	
  GPE.options.USE_RLT_TMODEL              = false;	//try false/true
  GPE.options.REUSE_TMODEL_THRESHOLD      = 1e-7;//0.e0
	
  GPE.options_TModel().BOUNDER_TYPE       = mc::TModel<I>::Options::LSB;
  GPE.options_TModel().SCALE_VARIABLES    = true;

  GPE.options_SetInv().DISPLAY            = 1;
  GPE.options_SetInv().MAX_CPU_TIME       = 1e4;
  GPE.options_SetInv().MAX_NODES          = 1000000;
  GPE.options_SetInv().ABSOLUTE_TOLERANCE = 8e-6; //1e-5;
  GPE.options_SetInv().RELATIVE_TOLERANCE = 8e-6; //1e-5;

  GPE.options_ODEBND_VAL().TSORDER        = 10;
  GPE.options_ODEBND_VAL().HMIN           = 1.e-8;
  GPE.options_ODEBND_VAL().WRAPMIT        = mc::ODEBND_VAL<I,PE>::Options::ELLIPS;//NONE;
  GPE.options_ODEBND_VAL().ORDMIT         = 1;

  GPE.options_ODEBND_GSL().INTMETH        = mc::ODEBND_GSL<I,PE>::Options::RKF45;//RK8PD;
  GPE.options_ODEBND_GSL().WRAPMIT        = mc::ODEBND_GSL<I,PE>::Options::ELLIPS;//DINEQ;
  GPE.options_ODEBND_GSL().ORDMIT         = 0;
  GPE.options_ODEBND_GSL().RTOL           =
  GPE.options_ODEBND_GSL().ATOL           = 1.e-7;
  GPE.options_ODEBND_GSL().HMIN           = 1.e-7;
  GPE.options_ODEBND_GSL().QTOL           = 1.e-12;
  GPE.options_ODEBND_GSL().DISPLAY        = 0;

  // GPE.options_FPRelax().SANDWICH_MAXCUT = 5;
  GPE.options_FPRelax().RLT_MAXORDER      = 4;
  GPE.options_FPRelax().RLT_DISPLAY       = 0;

  try{
    // Send results to file?
    std::ostream&os_iter = std::cout;
    //std::ofstream os_iter(  "ODEGPE_TM4_VERIF_DOMRED_TMSTORE_1e-6", std::ios_base::out );
    //std::ofstream os_open(  "ODEGPE_TM2_VERIF_NODOMRED_OPEN_1e-6", std::ios_base::out );
    //std::ofstream os_inner( "ODEGPE_TM2_VERIF_NODOMRED_INNER_1e-6", std::ios_base::out );

    // Solve GPE
    os_iter << GPE;
    std::pair<double,double> gpe_res = GPE.solve( NS, tk, P, &Y, os_iter );
    std::cout << std::fixed << std::setprecision(1)
      << "\n  Total Time: " << GPE.stats.cumul_SETINV
      << "\n  ODEBND: " << GPE.stats.cumul_ODEBND/GPE.stats.cumul_SETINV*1e2 << "%"
      << "    FPREL:  " << GPE.stats.cumul_FPREL/GPE.stats.cumul_SETINV*1e2 << "%"
      << "    REDUC:  " << GPE.stats.cumul_REDUC/GPE.stats.cumul_SETINV*1e2 << "%"
      << std::endl;
    //GPE.output_stacks( os_open, os_inner );

    //os_open.close();
    //os_inner.close();
    //os_iter.close();
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
  catch( mc::ODEGPE<I,PE>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr() << " in DOSBB:" << std::endl
              << eObj.what() << std::endl
              << "Execution aborted..." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
