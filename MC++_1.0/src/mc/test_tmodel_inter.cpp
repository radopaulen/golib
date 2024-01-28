#include <iostream>
#include "interval.hpp"
#include "tmodel.hpp"

using namespace mc;

typedef TModel<Interval> TM;
typedef TVar<Interval> TV;
typedef Interval I;

int main()
{
  unsigned int order(3);

  TM model(1,order);
  model.options.CENTER_REMAINDER = true;
  model.options.REF_MIDPOINT = true;

  TV p, x, y;
  p.set(&model);
  p.set(I(-.5,1)); 

  I x1(-1,.5);

  x = x1;

  std::cout <<"x = " <<  x;
  std::cout <<"p = " <<  p;
  std::cout <<"Intersect? " << inter(x,p,x) << std::endl;
  std::cout <<"Result = " <<  x;
}
