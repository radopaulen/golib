#include <fstream>
#include <iomanip>
#include "mcipopt.hpp"

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

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  std::pair<double,double> P[NP];
  P[0] = std::make_pair( 0., 6. );
  P[1] = std::make_pair( 0., 4. );
  double p0[NP];
  p0[0] = 1.;
  p0[1] = 2.;

  Ipopt::SmartPtr< mc::IPOPT<NLP> > pNLP = new mc::IPOPT<NLP>();
  pNLP->options.DISPLAY = 5;
  pNLP->options.MAXITER = 100;
  pNLP->options.CVTOL = 1e-9;
  pNLP->options.GRADIENT = mc::IPOPT<NLP>::Options::BACKWARD;
  pNLP->options.HESSIAN  = mc::IPOPT<NLP>::Options::EXACT;

  Ipopt::ApplicationReturnStatus status = pNLP->solve( P, p0 );

  if( status == Ipopt::Solve_Succeeded ){
    std::cout << "NLP (LOCAL) SOLUTION: " << std::endl;
    std::cout << "  f* = " << pNLP->solution().f << std::endl;
    for( unsigned int ip=0; ip<NP; ip++ )
      std::cout << "  p*(" << ip << ") = " << pNLP->solution().p[ip]
                << std::endl;
  }

  return 0;
}
