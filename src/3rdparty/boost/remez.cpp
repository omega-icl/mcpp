#include <iostream>
#include <iomanip>
#include "remez.hpp"

int main()
{
  double a = 0.6;
  boost::math::tools::remez_minimax<double>
  minimax
  //( [=](const double& x) { return std::log(x+2.); }, // double precision required
  ( [=](const double& x) { return std::exp(x); }, // double precision required
  //( [=](const double& x) { return std::atan(x); }, // double precision required
  //( [=](const double& x) { return std::pow(x,a); }, // double precision required
    4, 4, // 3rd degree polynomial
    0, 1, // -1 to 1 range
    //0.1, 2, // -1 to 1 range
    false, false, 0, 64 ); // other params: bool pin = true, bool rel_err = false, int sk = 0, int bits = 0

  for( unsigned iter=0; iter<10; ++iter ){
    minimax.iterate();
    const boost::math::tools::polynomial<double> num = minimax.numerator();
    const boost::math::tools::polynomial<double> den = minimax.denominator();
    const double e = minimax.max_error();
    const double de = minimax.max_change();
    std::cout << iter << ": [ " << std::right << std::scientific << std::setprecision(15);
    for( unsigned k=0; k<num.size(); ++k )
      std::cout << std::setw(23) << num[k];
    std::cout << " ] / [ ";
    for( unsigned k=0; k<den.size(); ++k )
      std::cout << std::setw(23) << den[k];
    std::cout << " ] +/- " << std::setw(23) << e << std::endl;
    if( de < 1e-5 ) break;
  }
  
  return 0;
}
