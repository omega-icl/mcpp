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
  int N(3), M(4);
  
  CPPL::drovector x(M), y;
  CPPL::dgematrix A(M,N);
  for(int i=0; i<x.l; i++){
    x(i) =double( rand() /(RAND_MAX/10) );
  }
  for(int i=0; i<A.m; i++){ for(int j=0; j<A.n; j++){
    A(i,j) =double( rand() /(RAND_MAX/10) );
  }}
  
  cout << "x =\n" << x << endl;
  cout << "A =\n" << A << endl;
  
  cout << "#### y=x*A; ####" << endl;
  y=x*A;
  cout << "y =\n" << y << endl;
  
  return 0;
}

/*****************************************************************************/
