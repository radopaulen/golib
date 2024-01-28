
#include "EA.hpp"
#include <fstream>
#include <iomanip>


int main()
{
  const int NX = 2, NF = 2;
  mc::FFTree FF;
  mc::FFVar X[NX];
  for( int i=0; i<NX; i++ ) X[i].set( &FF, i );

  mc::FFVar F[NF] = { 3*X[0]-2*X[1] + X[0]*X[1] ,
                      5*(X[0]*X[1]) - pow(X[0],4) };

  mc::Ellipsoidal_Image( &FF, NF, F );

  return 0;
}
