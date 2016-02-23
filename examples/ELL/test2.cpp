#include <fstream>
#include "ffunc.hpp"
#include "ellimage.hpp"
#include "interval.hpp"
typedef mc::Interval    I     ;

int main()
{
  /* Define enviroment to store factorable function FF */
  mc::FFGraph FF;
  /* Define independent variables X participating in FF */
  const unsigned int NX = 2; // No. of independent variables
  mc::FFVar X[NX];
  /* Set independent variables X in factorable function FF */
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &FF );
  /* Define dependent variables F participating in FF */ 
  const unsigned int NF = 2; // No. of dependent variables
  mc::FFVar F[NF] = { mc::sqrt(X[0]+X[1]) + X[0]*X[1], 
                     mc::pow(X[0]-X[1],2) + 3*X[1] };
  /* Output FF decomposition represented as code list */
  std::cout << FF << std::endl;
  /* Output FF decomposition represented as DAG */
  std::ofstream output_F( "F.dot", std::ios_base::out );
  FF.dot_script( NF, F, output_F ); // Generates DOT graph
  output_F.close();
  
  CPPL::dcovector q0(NF); CPPL::dsymatrix Q0(NF);
  q0(0) = 3.;     Q0(0,0) = 2.;
  q0(1) = 4.;     Q0(1,0) = 1.; Q0(1,1) = 2.; 
  
  //CPPL::dssmatrix depmap = FF.depmap(NF,F,NX,X);
  // Construct ellipsoidal domain
  mc::EllImg<I> Img;
  mc::Ellipsoid::options.PSDCHK = false;
  Img.set( Q0, q0 );//, depmap );
  Img.options.CHEBUSE  = false;

  // Set independent variables in the ellipsoidal model
  mc::EllVar<I> EX[NX];
  for( long i=0; i<NX; ++i ) EX[i].set( Img, i ); 
  mc::EllVar<I> EF[NF];
  FF.eval( NF, F, EF, NX, X, EX );
  Img.output();

  return 0;
}





