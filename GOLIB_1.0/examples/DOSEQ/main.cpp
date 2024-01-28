#include <fstream>
#include <iomanip>
//#include "mccvodes.hpp"
#include "doseq.hpp"

////////////////////////////////////////////////////////////////////////
// MINIMUM TIME CONTROL PROBLEM WITH TERMINAL EQUALITY CONSTRAINTS
// 
//    MIN  tf
//   tf, u, p
//     s.t.   x1(tf) = x1f
//            x2(tf) = x2f                 _
//             dx1dt = x2                   |  t in (0,tf]
//             dx2dt = (u-x1-2*x2)*t       _|
//             x1(0) = x10, x2(0) = x20
//            umin <= u  <= umax
//           tfmin <= tf <= tfmax
// 
//     where:  x1f   = 1.0
//             x2f   = p
//             x10   = 0.0
//             x20   = 0.0
//             umin  = -0.5
//             umax  =  1.5
//             tfmin = 0.1
//             tfmax = 10.0
////////////////////////////////////////////////////////////////////////

const unsigned int NS = 3;
const unsigned int NP = 2+NS;
const unsigned int NX = 3;
const unsigned int NC = 2;

class DYNOPT : public virtual mc::DOSTRUCT
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

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  double p0[NP], tk[NS+1];
  std::pair<double,double> P[NP];
  p0[0] = 6.;   P[0] = std::make_pair( 0.1, 2. );
  p0[1] = 0.5;  P[1] = std::make_pair( -1., 1. );
  tk[0] = 0.;
  for( unsigned int is=0; is<NS; is++ ){
    P[2+is]  = std::make_pair( -0.5, 1.5 );
    p0[2+is] = -1.+2*(double)is/(double)NS;
    tk[1+is] = tk[is]+1./(double)NS;
  }
  
/*
  //////////////////////////////////////////////////////////////////////
  // TESTING: CVODES<DYNOPT>

  mc::CVODES<DYNOPT>* pIVP = new mc::CVODES<DYNOPT>();
  pIVP->CVODES<DYNOPT>::options.DISPLAY = 1;
  pIVP->CVODES<DYNOPT>::options.ASACHKPT = 1000;
  pIVP->CVODES<DYNOPT>::options.FSAERR = true;
  pIVP->CVODES<DYNOPT>::options.HESSFORMAT = mc::CVODES<DYNOPT>::OPTIONS::HESSLOWER;
  pIVP->CVODES<DYNOPT>::options.FDATOL = pIVP->CVODES<DYNOPT>::options.FDRTOL = 1e-4;
  pIVP->CVODES<DYNOPT>::options.FDCEN  = true;

  double *xk[NS+1], *xkp[NS+1], *lk[NS+1], *lkp[NS+1], q[NP*(1+NC)], qp[NP*NP];
  for( unsigned int is=0; is<=NS; is++ ){
    tk[is]  = (is==0? 0.: tk[is-1]+1./(double)NS);
    xk[is]  = new double[NX];
    xkp[is] = new double[NX*NP];
    lk[is]  = new double[NX*(1+NC)];
    lkp[is] = new double[2*NX*NP];
  }

  double f, fp[NP], g[NC], *gp[NC], mu0[NC], Lpp[NP*NP];
  for( unsigned int ic=0; ic<NC; ic++ ){
    gp[ic] = new double[NP];
    mu0[ic] = 1.;
  }

  if( pIVP->states( NS, tk, p0, xk ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "IVP SOLUTION: " << std::endl;
    for( unsigned int is=0; is<=NS; is++ ){
      std::cout << "  Time: t[" << is << "] = " << tk[is]
                << std::endl;
      for( unsigned int ix=0; ix<NX; ix++ ){
        std::cout << "         xk(" << ix << ") = " << xk[is][ix] << std::endl;
      }
    }
  }

  if( pIVP->states_ASA( NS, tk, p0, xk, lk, q ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "IVP SOLUTION: " << std::endl;
    for( unsigned int is=0; is<=NS; is++ ){
      std::cout << "  Time: t[" << is << "] = " << tk[is]
                << std::endl;
      for( unsigned int ix=0; ix<NX; ix++ ){
        std::cout << "         xk(" << ix << ") = " << xk[is][ix];
        for( unsigned int ic=0; ic<1+NC; ic++ )
          std::cout << "  lk(" << ic << "," << ix << ") = " << lk[is][ix+ic*NX];
        std::cout << std::endl;
      }
    }
    std::cout << "  QUADRATURES:" << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ ){
      std::cout << "       ";
      for( unsigned int ic=0; ic<1+NC; ic++ )
        std::cout << "  q(" << ic << "," << ip << ") = " << q[ip+ic*NP];
      std::cout << std::endl;
    }
  }

  if( pIVP->states_FSA( NS, tk, p0, xk, xkp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "IVP SOLUTION: " << std::endl;
    for( unsigned int is=0; is<=NS; is++ ){
      std::cout << "  Time: t[" << is << "] = " << tk[is]
                << std::endl;
      for( unsigned int ix=0; ix<NX; ix++ ){
        std::cout << "         xk(" << ix << ") = " << xk[is][ix];
        for( unsigned int ip=0; ip<NP; ip++ )
          std::cout << "  xkp(" << ip << "," << ix << ") = " << xkp[is][ix+ip*NX];
        std::cout << std::endl;
      }
    }
  }
  
  if( pIVP->states_DSOASA( NS, tk, p0, mu0, xk, xkp, lkp, qp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "IVP SOLUTION: " << std::endl;
    for( unsigned int is=0; is<=NS; is++ ){
      std::cout << "  Time: t[" << is << "] = " << tk[is]
                << std::endl;
      for( unsigned int ix=0; ix<NX; ix++ ){
        std::cout << "         xk(" << ix << ") = " << xk[is][ix];
        for( unsigned int ip=0; ip<NP; ip++ )
          std::cout << "  xkp(" << ip << "," << ix << ") = " << xkp[is][ix+ip*NX];
        for( unsigned int jp=0; jp<NP; jp++ )
          std::cout << "  lkp(" << jp << "," << ix << ") = " << lkp[is][ix+jp*2*NX]
	            << "  lkp(" << jp << "," << NX+ix << ") = " << lkp[is][NX+ix+jp*2*NX];
        std::cout << std::endl;
      }
    }
    std::cout << "  QUADRATURES:" << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ ){
      std::cout << "       ";
      for( unsigned int jp=0; jp<NP; jp++ )
        std::cout << "  qp(" << jp << "," << ip << ") = " << qp[ip+jp*NP];
      std::cout << std::endl;
    }
  }

  if( pIVP->states_FDSA( NS, tk, p0, xk, xkp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "IVP SOLUTION: " << std::endl;
    for( unsigned int is=0; is<=NS; is++ ){
      std::cout << "  Time: t[" << is << "] = " << tk[is]
                << std::endl;
      for( unsigned int ix=0; ix<NX; ix++ ){
        std::cout << "         xk(" << ix << ") = " << xk[is][ix];
        for( unsigned int ip=0; ip<NP; ip++ )
          std::cout << "  xkp(" << ip << "," << ix << ") = " << xkp[is][ix+ip*NX];
        std::cout << std::endl;
      }
    }
  }

  if( pIVP->functions( NS, tk, p0, f, g ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "\nFUNCTIONS: " << std::endl;
    std::cout << "  f = " << f << std::endl;
    for( unsigned int ic=0; ic<NC; ic++ ){
      std::cout << "  g(" << ic << ") = " << g[ic] << std::endl;
    }
  }

  if( pIVP->functions_FSA( NS, tk, p0, f, fp, g, gp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "\nFUNCTIONS (FSA): " << std::endl;
    std::cout << "  f = " << f;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  fp(" << ip << ") = " << fp[ip];
    std::cout << std::endl;
    for( unsigned int ic=0; ic<NC; ic++ ){
      std::cout << "  g(" << ic << ") = " << g[ic];
      for( unsigned int ip=0; ip<NP; ip++ )
        std::cout << "  gp(" << ic << "," << ip << ") = " << gp[ic][ip];
      std::cout << std::endl;
    }
  }

  if( pIVP->functions_ASA( NS, tk, p0, f, fp, g, gp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "\nFUNCTIONS (ASA): " << std::endl;
    std::cout << "  f = " << f;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  fp(" << ip << ") = " << fp[ip];
    std::cout << std::endl;
    for( unsigned int ic=0; ic<NC; ic++ ){
      std::cout << "  g(" << ic << ") = " << g[ic];
      for( unsigned int ip=0; ip<NP; ip++ )
        std::cout << "  gp(" << ic << "," << ip << ") = " << gp[ic][ip];
      std::cout << std::endl;
    }
  }

  if( pIVP->functions_DSOASA( NS, tk, p0, mu0, f, fp, g, gp, Lpp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << std::right << std::scientific << std::setprecision(5);
    std::cout << "\nFUNCTIONS (DSOASA): " << std::endl;
    std::cout << "  f = " << f;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  fp(" << ip << ") = " << std::setw(12) << fp[ip];
    std::cout << std::endl;
    for( unsigned int ic=0; ic<NC; ic++ ){
      std::cout << "  g(" << ic << ") = " << std::setw(12) << g[ic];
      for( unsigned int ip=0; ip<NP; ip++ )
        std::cout << "  gp(" << ic << "," << ip << ") = " << std::setw(12) << gp[ic][ip];
      std::cout << std::endl;
    }
    switch( pIVP->CVODES<DYNOPT>::options.HESSFORMAT ){
    case mc::CVODES<DYNOPT>::OPTIONS::HESSFULL:
      for( unsigned int ip=0, ie=0; ip<NP; ip++ ){
        for( unsigned int jp=0; jp<NP; jp++, ie++ )
          std::cout << "  Lpp(" << ip << "," << jp << ") = " << std::setw(12) << Lpp[ie];
        std::cout << std::endl;
      }
      break;
    case mc::CVODES<DYNOPT>::OPTIONS::HESSLOWER:
      for( unsigned int ip=0, ie=0; ip<NP; ip++ ){
        for( unsigned int jp=0; jp<=ip; jp++, ie++ )
          std::cout << "  Lpp(" << ip << "," << jp << ") = " << std::setw(12) << Lpp[jp*NP-(jp*(jp-1))/2+ip-jp];
        std::cout << std::endl;
      }
      break;
    case mc::CVODES<DYNOPT>::OPTIONS::HESSUPPER:
      for( unsigned int ip=0; ip<NP; ip++ ){
        for( unsigned int jp=ip; jp<NP; jp++ )
          std::cout << "  Lpp(" << ip << "," << jp << ") = " << std::setw(12) << Lpp[jp*(jp+1)/2+ip];
        std::cout << std::endl;
      }
      break;
    }
  }

  if( pIVP->functions_FDSA( NS, tk, p0, f, fp, g, gp ) == mc::CVODES<DYNOPT>::NORMAL ){
    std::cout << "\nFUNCTIONS (FD): " << std::endl;
    std::cout << "  f = " << f;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  fp(" << ip << ") = " << fp[ip];
    std::cout << std::endl;
    for( unsigned int ic=0; ic<NC; ic++ ){
      std::cout << "  g(" << ic << ") = " << g[ic];
      for( unsigned int ip=0; ip<NP; ip++ )
        std::cout << "  gp(" << ic << "," << ip << ") = " << gp[ic][ip];
      std::cout << std::endl;
    }
  }
  
  for( unsigned int is=0; is<=NS; is++ ){
    delete[] xk[is];
    delete[] xkp[is];
    delete[] lk[is];
    delete[] lkp[is];
  }
  for( unsigned int ic=0; ic<NC; ic++ ){
    delete[] gp[ic];
  }
  delete pDO;
*/
  //////////////////////////////////////////////////////////////////////
  // TESTING: DOSEQ<DYNOPT>
  
  Ipopt::SmartPtr< mc::DOSEQ<DYNOPT> > pDO = new mc::DOSEQ<DYNOPT>();

  pDO->options.DISPLAY  = 5;
  pDO->options.MAXITER  = 1000;
  pDO->options.CVTOL    = 1e-9;
  pDO->options.GRADIENT = mc::DOSEQ<DYNOPT>::Options::FORWARD;
  pDO->options.SPARSE   = true;
  pDO->options.HESSIAN  = mc::DOSEQ<DYNOPT>::Options::LBFGS;

  pDO->CVODES<DYNOPT>::options.DISPLAY = 0;
  pDO->CVODES<DYNOPT>::options.ASACHKPT = 4000;
  pDO->CVODES<DYNOPT>::options.FSAERR = true;
  pDO->CVODES<DYNOPT>::options.HESSFORMAT = mc::CVODES<DYNOPT>::Options::HESSLOWER;
  pDO->CVODES<DYNOPT>::options.FDATOL =
  pDO->CVODES<DYNOPT>::options.FDRTOL = 1e-4;
  pDO->CVODES<DYNOPT>::options.FDCEN  = true;

  Ipopt::ApplicationReturnStatus status = pDO->solve( NS, tk, P, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "DO (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << pDO->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pDO->solution().p[ip]
                << std::endl;
  }

  return 0;
}
