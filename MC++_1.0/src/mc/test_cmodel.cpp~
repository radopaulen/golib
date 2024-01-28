#include <iostream>
#include "interval.hpp"
#include "cmodel.hpp"
#include "tmodel.hpp"

using namespace mc;

typedef CModel<Interval> CM;
typedef CVar<Interval> CV;
typedef TModel<Interval> TM;
typedef TVar<Interval> TV;
typedef Interval I;

int main()
{
  unsigned int order(3);

  CM Cmodel(2,order);

  CV X(&Cmodel, 0, I(-2,2));
  CV Y(&Cmodel, 1, I(-1,1));

  std::cout <<"X = " <<  X;
  std::cout <<"Y = " <<  Y;
  std::cout <<"X+Y = " << X+Y;
  std::cout <<"X-Y = " << X-Y;
  std::cout <<"X*Y = " << X*Y;
  std::cout <<"Y*Y = " << Y*Y;
  std::cout <<"X*X*X = " << X*X*X;

  Cmodel.options.BOUNDER_TYPE = CM::Options::EIGEN;
  std::cout <<"Y*Y+X*Y-X*X+X = " << Y*Y+X*Y-X*X+X;

  TM Tmodel(2,order);
  Tmodel.options.SCALE_VARIABLES = true;

  TV TMX(&Tmodel, 0, I(-2,2));
  TV TMY(&Tmodel, 1, I(-1,1));

  std::cout <<"X = " <<  TMX;
  std::cout <<"Y = " <<  TMY;
  std::cout << "X*X = " << TMX*TMX;
  
Tmodel.options.BOUNDER_TYPE = TM::Options::EIGEN;
  std::cout <<"Y*Y+X*Y-X*X+X = " << TMY*TMY+TMX*TMY-TMX*TMX+TMX;

}
