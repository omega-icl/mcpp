#include <cpplapack.h>
#include <ctime>

//=============================================================================
inline void ended()
{
  system("@pause");
}

//=============================================================================
int main(int /*argc*/, char** /*argv*/)
{
  //////// to keep the command prompt open ////////
  atexit(ended);
  
  //////// declare A and y ////////
  const size_t size(3);
  CPPL::dgematrix A(size,size); //3x3 matrix
  CPPL::dcovector y(size); //3 column vector
  
  //////// set components of A and y at random ////////
  srand(unsigned(time(NULL)));
  for(size_t i=0; i<size; i++){
    for(size_t j=0; j<size; j++){
      A(i,j) =double(rand())/double(RAND_MAX);
    }
    y(i) =double(rand())/double(RAND_MAX);
  }
  std::cout << "A=\n" << A << std::endl; //print A
  std::cout << "y=\n" << y << std::endl; //print y
  
  //////// solve A*x=y ////////
  A.dgesv(y); //call dgesv of LAPACK (y becomes x)
  std::cout << "x=\n" << y << std::endl; //print x
  
  return 0;
}
