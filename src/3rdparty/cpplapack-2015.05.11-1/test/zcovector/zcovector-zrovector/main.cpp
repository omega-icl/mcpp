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
  int M(4);
  
  CPPL::zcovector x(M);
  CPPL::zrovector y(M);
  
  for(int i=0; i<x.l; i++){
    x(i) =complex<double>(rand()/(RAND_MAX/10), rand()/(RAND_MAX/10));
  }
  for(int i=0; i<y.l; i++){
    y(i) =complex<double>(rand()/(RAND_MAX/10), rand()/(RAND_MAX/10));
  }
  
  cout << "x =\n" << x << endl;
  cout << "y =\n" << y << endl;
  cout << "x*y =\n" << x*y << endl;
  
  return 0;
}

/*****************************************************************************/
