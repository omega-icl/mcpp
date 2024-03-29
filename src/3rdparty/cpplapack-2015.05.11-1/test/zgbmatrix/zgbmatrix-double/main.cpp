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
  int M(4), N(3), KL(1), KU(2);
  
  CPPL::zgbmatrix A(M,N,KL,KU);
  for(int i=0; i<A.m; i++){ for(int j=0; j<A.n; j++){
    if(!( i-j>A.kl || j-i>A.ku )){
      A(i,j) =complex<double>(rand()/(RAND_MAX/10), rand()/(RAND_MAX/10));
    }
  }}
  
  cout << "A =\n" << A << endl;
  cout << "A*10. =\n" << A*10. << endl;
  cout << "A/10. =\n" << A/10. << endl;
  
  cout << "#### A*=10.; ####" << endl;
  A*=10.;
  cout << "A =\n" << A << endl;
  cout << "#### A/=10.; ####" << endl;
  A/=10.;
  cout << "A =\n" << A << endl;
  
  return 0;
}

/*****************************************************************************/
