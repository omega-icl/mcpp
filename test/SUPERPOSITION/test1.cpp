#include <iostream>
#include "interval.hpp"
#include "ismodel.hpp"

int main(int argc, char* argv[])
{
  using namespace mc;
  ISModel<Interval> mod( 2, 8 );

  ISVar<Interval> X( &mod, 0, Interval(-1.,1) );
  std::cout << "ISM of X:\n" << X << std::endl;
  //std::cout << X.B() << std::endl;

  ISVar<Interval> Y( &mod, 1, Interval(3.,3.1) );
  std::cout << "ISM of Y:\n" << Y << std::endl;

  ISVar<Interval> Z( -X );
  std::cout << "ISM of Z(-X):\n" << Z << std::endl;

  Z *= Y;
  std::cout << "ISM of Z*=Y:\n" << Z << std::endl;

  Z += Y;
  std::cout << "ISM of Z+=Y:\n" << Z << std::endl;

  Z = 1 - X + Y;
  std::cout << "ISM of Z=1-X+Y:\n" << Z << std::endl;

  Z = 0.;
  std::cout << "ISM of Z=0.:\n" << Z << std::endl;

  Z = sqr(X);
  std::cout << "ISM of Z=sqr(X):\n" << Z << std::endl;
  Z.display();

  Z = sqr(X+Y);
  std::cout << "ISM of Z=sqr(X+Y):\n" << Z << std::endl;

  Z *= Y;
  std::cout << "ISM of Z*=Y:\n" << Z << std::endl;

  return 0;
}



