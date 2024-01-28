#include "ffunc.hpp"
#include <fstream>
#include <iomanip>

int main()
{
  const int NX = 4, NF = 2;
  mc::FFTree FF;
  mc::FFVar X[NX];
  for( int i=0; i<NX; i++ ) X[i].set( &FF, i );

  mc::FFVar F[NF] = { X[2]*X[3]-X[0],
                      X[0]*pow(exp(X[2]*X[3]*2)+3.,3)+X[1] };
  std::cout << FF;
  std::cout << "Variable dependence of F[0]: " << F[0].dep() << std::endl;
  std::cout << "Variable dependence of F[1]: " << F[1].dep() << std::endl;

  std::list<mc::FFOp*> f1_op = FF.subtree( 1, F );     FF.output( f1_op );
  std::list<mc::FFOp*> f2_op = FF.subtree( 1, F+1 );   FF.output( f2_op );
  std::list<mc::FFOp*> f1_f2_op = FF.subtree( NF, F ); FF.output( f1_f2_op );

  std::ofstream of1f2( "f1_f2.dot", std::ios_base::out );
  FF.dot_script( NF, F, of1f2 );
  of1f2.close();

  std::ofstream of1( "f1.dot", std::ios_base::out );
  FF.dot_script( 1, F, of1 );
  of1.close();

  return 0;
}
