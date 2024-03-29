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
  
  CPPL::dcovector x(L);
  for(int i=0; i<x.l; i++){
	x(i) =double( rand() /(RAND_MAX/10) );
  }
  
  cout << "x =\n" << x << endl;
  cout << "x*10. =\n" << x*10. << endl;
  cout << "x/10. =\n" << x/10. << endl;
  
  cout << "#### x*=10.; ####" << endl;
  x*=10.;
  cout << "x =\n" << x << endl;
  cout << "#### x/=10.; ####" << endl;
  x/=10.;
  cout << "x =\n" << x << endl;
  
  return 0;
}

/*****************************************************************************/
