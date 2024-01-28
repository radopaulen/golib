#define USE_PROFIL
//#define USE_FILIB

#include <fstream>
#include <iomanip>
#include "nlpsbb.h"

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

#define NLP1			// <- select the example here

#if defined( NLP1 )
////////////////////////////////////////////////////////////////////////
// EXAMPLE 1 IN RYOO & SAHINIDIS, C&CE, 1995
////////////////////////////////////////////////////////////////////////
const unsigned int NP = 2;
const unsigned int NC = 1;

class NLP: public virtual mc::NLPSTRUCT
{
public:
  NLP(): mc::NLPSTRUCT( NP, NC )
    {}

  template <typename U>
  std::pair<U,t_OBJ> OBJ
    ( const U*p )
    {
      return std::make_pair( p[0]+p[1], MAX );
    }

  template <typename U>
  std::pair<U,t_CTR> CTR
    ( const unsigned int ic, const U*p )
    {
      switch( ic ){
        case 0: return std::make_pair( p[0]*p[1]-4, LE );
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( NLP2 )
////////////////////////////////////////////////////////////////////////
// EXAMPLE 7 IN RYOO & SAHINIDIS, C&CE, 1995
////////////////////////////////////////////////////////////////////////
const unsigned int NP = 10;
const unsigned int NC = 7;

class NLP: public virtual mc::NLPSTRUCT
{
public:
  NLP(): mc::NLPSTRUCT( NP, NC )
    {}

  template <typename U>
  std::pair<U,t_OBJ> OBJ
    ( const U*p )
    {
      return std::make_pair( -9*p[4]-15*p[8]+6*p[0]+16*p[1]+10*p[5], MIN );
    }

  template <typename U>
  std::pair<U,t_CTR> CTR
    ( const unsigned int ic, const U*p )
    {
      switch( ic ){
        case 0: return std::make_pair( p[0]+p[1]-p[2]-p[3], EQ );
        case 1: return std::make_pair( p[2]+p[6]-p[4], EQ );
        case 2: return std::make_pair( p[3]+p[7]-p[8], EQ );
        case 3: return std::make_pair( p[6]+p[7]-p[5], EQ );
        case 4: return std::make_pair( p[9]*p[2]+2*p[6]-2.5*p[4], LE );
        case 5: return std::make_pair( p[9]*p[3]+2*p[7]-1.5*p[8], LE );
        //case 6: return std::make_pair( 3*p[0]+p[1]-p[9]*p[2]-p[9]*p[3], LE );
        case 6: return std::make_pair( 3*p[0]+p[1]-p[9]*(p[2]+p[3]), LE );
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( NLP3 )
////////////////////////////////////////////////////////////////////////
// EXAMPLE 1 IN BEN-TAL ET AL., MATH PROG, 1994
////////////////////////////////////////////////////////////////////////
const unsigned int NP = 7;
const unsigned int NC = 6;

class NLP: public virtual mc::NLPSTRUCT
{
public:
  NLP(): mc::NLPSTRUCT( NP, NC )
    {}

  template <typename U>
  std::pair<U,t_OBJ> OBJ
    ( const U*p )
    {
      const U *q = p, *y = q+3, *z = y+2;  
      U f = -z[0]+5*z[1]
           + ( 9-6*q[0]-16*q[1]-15*q[2])*y[0]
           + (15-6*q[0]-16*q[1]-15*q[2])*y[1];
      return std::make_pair( f, MAX );
    }

  template <typename U>
  std::pair<U,t_CTR> CTR
    ( const unsigned int ic, const U*p )
    {
      const U *q = p, *y = q+3, *z = y+2;  
      switch( ic ){
        case 0:
          return std::make_pair( q[2]*y[0]+q[2]*y[1]-50, LE );
        case 1:
          return std::make_pair( y[0]+z[0]-100, LE );
          break;
        case 2:
          return std::make_pair( y[1]+z[1]-200, LE );
          break;
        case 3:
          return std::make_pair( (3*q[0]+q[1]+q[2]-2.5)*y[0]-0.5*z[0], LE );
          break;
        case 4:
          return std::make_pair( (3*q[0]+q[1]+q[2]-1.5)*y[1]+0.5*z[1], LE );
          break;
        case 5:
          return std::make_pair( q[0]+q[1]+q[2]-1, EQ );
          break;
        default: throw std::runtime_error("invalid size");
      }
    }
};

#elif defined( NLP4 )
////////////////////////////////////////////////////////////////////////
// EXAMPLE 2 IN BEN-TAL ET AL., MATH PROG, 1994
////////////////////////////////////////////////////////////////////////
const unsigned int NP = 32;
const unsigned int NC = 19;

class NLP: public virtual mc::NLPSTRUCT
{
public:
  NLP(): mc::NLPSTRUCT( NP, NC )
    {}

  template <typename U>
  std::pair<U,t_OBJ> OBJ
    ( const U*p )
    {
      const U *q = p, *y = q+12, *z = y+15;  
      U f = 8*z[0]+5*z[1]+9*z[2]+6*z[3]+4*z[4]
           + (18-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[0]
           + (18-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[1]
           + (18-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[2]
           + (15-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[3]
           + (15-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[4]
           + (15-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[5]
           + (19-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[6]
           + (19-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[7]
           + (19-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[8]
           + (16-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[9]
           + (16-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[10]
           + (16-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[11]
           + (14-6*q[0]-16*q[1]-15*q[2] -12*q[3]) *y[12]
           + (14-6*q[4]-16*q[5]-15*q[6] -12*q[7]) *y[13]
           + (14-6*q[8]-16*q[9]-15*q[10]-12*q[11])*y[14];
      return std::make_pair( f, MAX );
    }

  template <typename U>
  std::pair<U,t_CTR> CTR
    ( const unsigned int ic, const U*p )
    {
      const U *q = p, *y = q+12, *z = y+15;  
      switch( ic ){
        case 0:
          return std::make_pair(
            q[2]*y[0] +q[2]*y[3] +q[2]*y[6] +q[2]*y[9]  +q[2]*y[12]
           +q[6]*y[1] +q[6]*y[4] +q[6]*y[7] +q[6]*y[10] +q[6]*y[13]
           +q[10]*y[2]+q[10]*y[5]+q[10]*y[8]+q[10]*y[11]+q[10]*y[14]-50, LE );

        case 1:
          return std::make_pair( y[0]+y[1]+y[2]+z[0]-100, LE );
          break;
        case 2:
          return std::make_pair( y[3]+y[4]+y[5]+z[1]-200, LE );
          break;
        case 3:
          return std::make_pair( y[6]+y[7]+y[8]+z[2]-100, LE );
          break;
        case 4:
          return std::make_pair( y[9]+y[10]+y[11]+z[3]-100, LE );
          break;
        case 5:
          return std::make_pair( y[12]+y[13]+y[14]+z[4]-100, LE );
          break;

        case 6:
          return std::make_pair(
            (3*q[0]+q[1]+q[2]+1.5*q[3]-2.5)*y[0]
           +(3*q[4]+q[5]+q[6]+1.5*q[7]-2.5)*y[1]
           +(3*q[8]+q[9]+q[10]+1.5*q[11]-2.5)*y[2] - 0.5*z[0], LE );
          break;
        case 7:
          return std::make_pair(
            (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[0]
           +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[1]
           +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[2] + 0.5*z[0], LE );
          break;
        case 8:
          return std::make_pair(
            (3*q[0]+q[1]+q[2]+1.5*q[3]-1.5)*y[3]
           +(3*q[4]+q[5]+q[6]+1.5*q[7]-1.5)*y[4]
           +(3*q[8]+q[9]+q[10]+1.5*q[11]-1.5)*y[5] + 0.5*z[1], LE );
          break;
        case 9:
          return std::make_pair(
            (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.5)*y[3]
           +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.5)*y[4]
           +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.5)*y[5], LE );
          break;
        case 10:
          return std::make_pair(
            (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[6]
           +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[7]
           +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[8], LE );
          break;
        case 11:
          return std::make_pair(
            (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2.6)*y[6]
           +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2.6)*y[7]
           +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2.6)*y[8] - 0.1*z[2], LE );
          break;
        case 12:
          return std::make_pair(
            (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[9]
           +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[10]
           +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[11], LE );
          break;
        case 13:
          return std::make_pair(
            (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[9]
           +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[10]
           +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[11] + 0.5*z[3], LE );
          break;
        case 14:
          return std::make_pair(
            (3*q[0]+q[1]+q[2]+1.5*q[3]-2)*y[12]
           +(3*q[4]+q[5]+q[6]+1.5*q[7]-2)*y[13]
           +(3*q[8]+q[9]+q[10]+1.5*q[11]-2)*y[14], LE );
          break;
        case 15:
          return std::make_pair(
            (q[0]+3*q[1]+2.5*q[2]+2.5*q[3]-2)*y[12]
           +(q[4]+3*q[5]+2.5*q[6]+2.5*q[7]-2)*y[13]
           +(q[8]+3*q[9]+2.5*q[10]+2.5*q[11]-2)*y[14] + 0.5*z[4], LE );
          break;

        case 16:
          return std::make_pair( q[0]+q[1]+q[2]+q[3]-1, EQ );
          break;
        case 17:
          return std::make_pair( q[4]+q[5]+q[6]+q[7]-1, EQ );
          break;
        case 18:
          return std::make_pair( q[8]+q[9]+q[10]+q[11]-1, EQ );
          break;
        default: throw std::runtime_error("invalid size");
      }
    }
};
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
#if defined( NLP1 )
  I P[NP] = { I(0,6), I(0,4) };
  double p0[NP] = { 1, 2 };
#elif defined( NLP2 )
  I P[NP] = { I(0,300), I(0,300), I(0,100), I(0,200), I(0,100), I(0,300),
              I(0,100), I(0,200), I(0,200), I(1,3) };
  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };
#elif defined( NLP3 )
  I P[NP] = { I(0,1), I(0,1), I(0,1), I(0,100), I(0,200), I(0,100), I(0,200) };
  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0 };
#elif defined( NLP4 )
  I P[NP] = { I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,1), I(0,1), I(0,1), I(0,1), I(0,1), I(0,1),
              I(0,100), I(0,100), I(0,100), I(0,200), I(0,200), I(0,200),
              I(0,100), I(0,100), I(0,100), I(0,100), I(0,100), I(0,100),
              I(0,100), I(0,100), I(0,100),
              I(0,100), I(0,200), I(0,100), I(0,100), I(0,100) };
//  double p0[NP] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
//                    0, 0, 0, 0, 0 };
  double p0[NP] = { 1, 0, 0, 0, 0, 0, 0, 1, 0.2757, 0, 0, 0.7243,
                    51.5527, 0, 17.9504, 0, 200, 0, 0, 0, 0,
                    7.9609, 0, 92.0391, 15.7866, 20.5622, 63.6512,
                    30.4969, 0, 100, 0, 0 };
#endif

  std::pair<double,const double*> optim;
  mc::NLPSBB<I,NLP> GNLO;

/*   GNLO.options.TAYLOR_MODEL_ORDER         = 4;
  GNLO.options.TAYLOR_MODEL_CUTS          = mc::NLPSBB<I,NLP>::Options::NONE;
  GNLO.options.USE_CONSTRAINT_PROPAGATION = false;
  GNLO.options.USE_REDUCTION_CONSTRAINTS  = true;
  GNLO.options.USE_DOMAIN_REDUCTION       = false;

  GNLO.options_SBB().BRANCHING_VARIABLE_CRITERION = mc::SBB<I>::Options::RGREL;
//   GNLO.options_SBB().BRANCHING_RELIABILITY_THRESHOLD = 1;
//   GNLO.options_SBB().DISPLAY = 2;
//   GNLO.options_SBB().MAX_CPU_TIME = 1e6;
  GNLO.options_SBB().MAX_NODES = 0;
//   GNLO.options_SBB().ABSOLUTE_TOLERANCE = 1e-4;
//   GNLO.options_SBB().RELATIVE_TOLERANCE = 1e-4;

  GNLO.options_IPOPT().DISPLAY  = 0;
  GNLO.options_IPOPT().MAXITER  = 100;
  GNLO.options_IPOPT().CVTOL    = 1e-7;
  GNLO.options_IPOPT().GRADIENT = mc::IPOPT<NLP>::Options::FORWARD;
  GNLO.options_IPOPT().HESSIAN  = mc::IPOPT<NLP>::Options::LBFGS;

  GNLO.options_FPRelax().BILINEAR_SUBDIV = 1;
  GNLO.options_FPRelax().BILINEAR_RULE   = mc::FPRelax<I>::Options::SEPARSOS;
  GNLO.options_FPRelax().RRLT_STRATEGY   = mc::FPRelax<I>::Options::PRIMRRLT;
  GNLO.options_FPRelax().RRLT_DISPLAY    = 0;
 */
  std::cout << GNLO;
  optim = GNLO.solve( P, p0 );
  std::cout << std::fixed << std::setprecision(1)
    << "\n  Total Time: " << GNLO.stats.cumul_SBB
    << "\n  NLPLOC: " << GNLO.stats.cumul_NLPLOC/GNLO.stats.cumul_SBB*1e2 << "%"
    << "    FPREL:  " << GNLO.stats.cumul_FPREL /GNLO.stats.cumul_SBB*1e2 << "%"
    << "    REDUC:  " << GNLO.stats.cumul_REDUC /GNLO.stats.cumul_SBB*1e2 << "%"
    << "    NLPREL: " << GNLO.stats.cumul_NLPREL/GNLO.stats.cumul_SBB*1e2 << "%"
    << std::endl;

  return 0;
}
