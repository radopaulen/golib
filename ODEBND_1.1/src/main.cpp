#include <fstream>
#include <iomanip>
#include <stdexcept>
#include "odetaylor.hpp"

#include "interval.hpp"
typedef mc::Interval I;
typedef mc::TModel<I> TM;
typedef mc::TVar<I> TV;

////////////////////////////////////////////////////////////////////////

#define REPR			// <- select the example here
const unsigned int NTM   = 5;	// <- Order of Taylor model
const unsigned int NSAMP = 50;	// <- Number of sampling points
#define SAVE_RESULTS		// <- Saving the results to file

////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
//
//             dxdt = cos(t)    t in (0,tf]
//             x(0) = p
////////////////////////////////////////////////////////////////////////

const unsigned int NE = 10;
const unsigned int NS = 1;
const unsigned int NP = 0;
const unsigned int NX = 1;
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
        case 0: return cos(x[0])*t; //*p[0]; //cos(t);
        default: throw std::runtime_error("invalid index");
      }
    }

  template <typename T>
  T IC
    ( const unsigned int ix, const T*p )
    {
      assert( ix < nx() );
      switch( ix ){
        case 0: return 1.; //p[0];
        default: throw std::runtime_error("invalid index");
      }
    }
};

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  

  mc::ODETAYLOR<double,IVP> pIVP;
  double t0 = 0.;
  double p = 1.;
  double d = 2.;
  double x0 = IVP().IC( 0, &p );

  double**f = new double*[NE+1];
  for( unsigned int k=0; k<=NE; k++ ) f[k] = new double[NX];

  pIVP.T_expand( 0, t0, &x0, &p, NE, f );

  double**dfdx = new double*[NE+1];
  for( unsigned int k=0; k<=NE; k++ ) dfdx[k] = new double[NX*NX];
  double**dfdp = new double*[NE+1];
  for( unsigned int k=0; k<=NE; k++ ) dfdp[k] = new double[NX*NP];

  pIVP.TF_expand( 0, t0, &x0, &p, NE, f, dfdx ); //, dfdp );

  double**d2fdx2dd = new double*[NE+1];
  for( unsigned int k=0; k<=NE; k++ ) d2fdx2dd[k] = new double[NX];

  pIVP.TFFd_expand( 0, t0, &x0, &p, NE, f, &d, dfdx );

  for( unsigned int k=0; k<=NE; k++ ){
    delete[] f[k];
    delete[] dfdx[k];
    delete[] dfdp[k];
    delete[] d2fdx2dd[k];
  }
  delete[] f;
  delete[] dfdx;
  delete[] dfdp;
  delete[] d2fdx2dd;
  
  return 0;
}
