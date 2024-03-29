/*****************************************************************************/
/*                                 noname                                    */
/*****************************************************************************/

//=============================================================================
#include <iostream>
#include <cstdlib>
#include <ctime>
#include "cpplapack.h"
using namespace std;

//=============================================================================
/*! main */
int main(int argc, char** argv)
{
  srand(unsigned(time(NULL)));
  //int M(5), N(4), KL(2), KU(1);
  int M(10), N(10), KL(1), KU(1);
  
  CPPL::dgbmatrix A(M,N,KL,KU);
  for(int i=0; i<A.m; i++){ for(int j=0; j<A.n; j++){
    //if(!( i-j>A.kl || j-i>A.ku )){ A(i,j) =double( rand() /(RAND_MAX/10) ); }
    if(!( i-j>A.kl || j-i>A.ku )){ A(i,j) =10.*double(i)+double(j); }
  }}
  
  cout << "A =\n" << A << endl;
  for(int i=0; i<A.m; i++){ for(int j=0; j<A.n; j++){
    if(!( i-j>A.kl || j-i>A.ku )){
      cout << "A(" << i << "," << j << ") =" << A(i,j) << endl;
    }
  }}
  
  
  
  const CPPL::dgbmatrix B(A);
  cout << "B =\n" << B << endl;
  
  //B*=10.; //compile error
  //B(0,0)=0.; //compile error
  
  //cout << "A+B=\n" << A+B << endl;
  
  //// write/read ////
  B.write( "tmp.txt" );
  CPPL::dgbmatrix C;
  C.read( "tmp.txt" );
  cout << "C-B =\n" << C-B << "<-Should be zero." << endl;

  return 0;
}

/*****************************************************************************/
