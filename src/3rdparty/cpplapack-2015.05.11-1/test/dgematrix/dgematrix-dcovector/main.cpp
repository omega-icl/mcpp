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
  int M(4), N(3);
  
  CPPL::dgematrix A(M,N);
  CPPL::dcovector x(N);
  for(int i=0; i<A.m; i++){ for(int j=0; j<A.n; j++){
	A(i,j) =double( rand() /(RAND_MAX/10) );
  }}
  for(int i=0; i<x.l; i++){
	x(i) =double( rand() /(RAND_MAX/10) );
  }
  
  cout << "A =\n" << A << endl;
  cout << "x =\n" << x << endl;
  cout << "A*x =\n" << A*x << endl;
  
  return 0;
}

/*****************************************************************************/
