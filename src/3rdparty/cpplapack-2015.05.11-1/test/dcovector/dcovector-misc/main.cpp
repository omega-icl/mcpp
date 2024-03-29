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
  int L(4);
  
  CPPL::dcovector x(L), y;
  for(int i=0; i<x.l; i++){
	x(i) =double( rand() /(RAND_MAX/10) );
  }
  
  cout << "x =\n" << x << endl;
  cout << "#### y.copy(x) ####" << endl;
  y.copy(x);
  cout << "y =\n" << y << endl;
  
  cout << "#### y.clear() ####" << endl;
  y.clear();
  cout << "y =\n" << y << endl;
  
  cout << "#### y.resize(2) & y.zero() ####" << endl;
  y.resize(2);
  y.zero();
  cout << "y =\n" << y << endl;
  
  return 0;
}

/*****************************************************************************/
