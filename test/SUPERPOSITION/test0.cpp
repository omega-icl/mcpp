#include <iostream>
#include "interval.hpp"
#include "asmodel.hpp"

int main(int argc, char* argv[])
{
  using namespace mc;
  ASModel<Interval> mod( 2, 8 );

  ASVar<Interval> X( &mod, 0, Interval(-2,1) );
  std::cout << "ASM of X:\n" << X << std::endl;

  ASVar<Interval> Y( &mod, 1, Interval(-1,2) );
  std::cout << "ASM of Y:\n" << Y << std::endl;

  ASVar<Interval> Z = relu( X+Y ) + relu( X );
  std::cout << "ASM of ReLU(X+Y):\n" << Z << std::endl;

  return 0;
}



