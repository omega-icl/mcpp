#include <iostream>
#include <iomanip>
#include "remez.hpp"

int main()
{
  double a = 0.6;
  boost::math::tools::remez_minimax<double>
  approximate_atan
  //( [=](const double& x) { return std::atan(x); }, // double precision required
  ( [=](const double& x) { return std::pow(x,a); }, // double precision required
    3, 0, // 3rd degree polynomial
  //  -1, 1, // -1 to 1 range
    0.1, 2, // -1 to 1 range
    false, false, 0, 64 ); // other params: bool pin = true, bool rel_err = false, int sk = 0, int bits = 0

  for( unsigned iter=0; iter<50; ++iter ){
    approximate_atan.iterate();
    const boost::math::tools::polynomial<double> a = approximate_atan.numerator();
    const double e = approximate_atan.max_error();
    const double de = approximate_atan.max_change();
    std::cout << iter << ": [ " << std::right << std::scientific << std::setprecision(15);
    for( unsigned k=0; k<a.size(); ++k )
      std::cout << std::setw(23) << a[k];
    std::cout << " ] +/- " << std::setw(23) << e << std::endl;
    if( de < 1e-7 ) break;
  }
  
  return 0;
}
