#include <fstream>
#include <iostream>
#include <iomanip>

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  std::string str_open;
  std::ifstream if_open;
  while( if_open.good() ){
    std::cout << "Name of result file (open nodes)?\n";
    std::cin >> str_open;
    std::cout << str_open.c_str() << std::endl;
    if_open.open( str_open.c_str() );
  }
  //if_open.open( "ODEGPE_TM2_VERIF_DOMRED_OPEN_1e-6" );

  // Read line-by-line
  std::string line;
  while( getline( if_open, line ) )
    std::cout << line.c_str() << std::endl;

  return 0;
}
