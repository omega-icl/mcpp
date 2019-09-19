#include <fstream>
#include "interval.hpp"
#include "ismodel.hpp"
#include "ffunc.hpp"

#undef MC__ISMODEL_TRACE

typedef mc::Interval I;
typedef mc::ISModel<I> ISM;
typedef mc::ISVar<I> ISV;

int main()
{
  ISM mod( 2, 8 );
  I IX(1.,2.);
  I IY(0.,1.);
  ISV ISX( &mod, 0, IX );
  ISV ISY( &mod, 1, IY );

  auto IF  = IX*exp(IX+sqr(IY))-sqr(IY); 
  std::cout << "Interval inclusion of f:\n" << IF << std::endl;
  
  auto ISF = ISX*exp(ISX+sqr(ISY))-sqr(ISY);
  std::cout << "Interval superposition model of f:\n" << ISF << std::endl;
  auto&& mat = ISF.C();

  std::ofstream ofile( "test2.dat", std::ofstream::out );
  ISF.display( ofile );
  ofile.close();

  ///////////////////////////

  mc::FFGraph DAG;
  mc::FFVar X( &DAG), Y(&DAG);
  mc::FFVar F = X*exp(X+sqr(Y))-sqr(Y);

  // Evaluate in interval arithmetic
  std::vector<I> IWK;
  DAG.eval( IWK, 1, &F, &IF, 1, &X, &IX, 1, &Y, &IY );
  std::cout << "Interval inclusion of f:\n" << IF << std::endl;

  // Evaluate in interval superposition arithmetic
  std::vector<ISV> ISWK;
  DAG.eval( ISWK, 1, &F, &ISF, 1, &X, &ISX, 1, &Y, &ISY );
  std::cout << "Interval superposition model of f:\n" << ISF << std::endl;
  //std::cout << "Interval inclusion of f:\n" << IF << std::endl;

  return 0;
}
