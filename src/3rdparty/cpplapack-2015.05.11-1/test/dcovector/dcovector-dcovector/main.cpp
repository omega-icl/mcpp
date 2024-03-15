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
  
  CPPL::dcovector x(M), y(M), z;
  for(int i=0; i<x.l; i++){
	x(i) =double( rand() /(RAND_MAX/10) );
  }
  for(int i=0; i<y.l; i++){
	y(i) =double( rand() /(RAND_MAX/10) );
  }
  
  cout << "x =\n" << x << endl;
  cout << "y =\n" << y << endl;
  
  cout << "x+x =\n" << x+x << endl;
  cout << "x-x =\n" << x-x << endl;
  cout << "x%y =\n" << x%y << endl;
  
  cout << "#### z=x; ####" << endl;
  z=x;
  cout << "z =\n" << z << endl;
  cout << "#### z=x+x-x; ####" << endl;
  z=x+x-x;
  cout << "z =\n" << z << endl;
  cout << "#### z+=x; ####" << endl;
  z+=x;
  cout << "z =\n" << z << endl;
  cout << "#### z-=x; ####" << endl;
  z-=x;
  cout << "z =\n" << z << endl;
  
  cout << "hadamerd(x,x) =\n" << hadamerd(x,x) << endl;
  cout << "hadamerd(x,x) =\n" << hadamerd(x,hadamerd(x,x)) << endl;
  return 0;
}

/*****************************************************************************/
