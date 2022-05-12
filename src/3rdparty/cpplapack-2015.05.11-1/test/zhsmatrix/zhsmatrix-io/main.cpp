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
  int N(5), CAP(4);
  
  CPPL::zhsmatrix A(N,CAP);
  A.put(0,0, complex<double>(1.,0.) );
  A.put(3,2, complex<double>(3.,4.) );
  A.put(1,2, complex<double>(5.,6.) );
  A.put(4,1, complex<double>(7.,8.) );
  cout << "A =\n" << A << endl;
  cout << "A(0,0) = " << A(0,0) << endl;
  cout << "A(3,2) = " << A(3,2) << endl;
  cout << "A(1,2) = " << A(1,2) << endl;
  cout << "A(4,1) = " << A(4,1) << endl;
  
  //A.put(1,2, 4.5);
  //A.add(1,2, 0.1);
  //A.sub(1,2, 0.1);
  //A.mult(1,2, 10.);
  //A.div(1,2, 10.);
  //A.del(1,2);
  A.del(1);
  cout << "A =\n" << A << endl;
  
  //// write/read ////
  const CPPL::zhsmatrix B(A);
  cout << "B =\n" << B << endl;
  B.write( "tmp.txt" );
  
  CPPL::zhsmatrix C;
  C.read( "tmp.txt" );
  cout << "C =\n" << C << endl;
  cout << "C-B =\n" << C-B << "<-Should be zero." << endl;
  
  return 0;
}

/*****************************************************************************/
