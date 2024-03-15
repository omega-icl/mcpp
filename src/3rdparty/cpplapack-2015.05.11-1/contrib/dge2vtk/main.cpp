#include <cpplapack.h>

//=============================================================================
int main(int argc, char** argv)
{
  //////// argc check ////////
  if(argc!=2){
    std::cerr << "[ERROR] invalid usage" << std::endl;
    std::cerr << "USAGE: " << argv[0] << " xxxx.dgematrix" << std::endl;
    exit(1);
  }
  
  //////// open ////////
  CPPL::dgematrix mat(argv[1]);
  std::cerr << "mat size = " << mat.m << "x" << mat.n << std::endl;
  const long mn(mat.m*mat.n);
  
  
  /////// write header ///////
  std::cout << "# vtk DataFile Version 2.0" << std::endl;
  std::cout << "made by dge2vtk" << std::endl;
  std::cout << "ASCII" << std::endl;
  std::cout << "DATASET UNSTRUCTURED_GRID" << std::endl;
  std::cout << std::endl;
  
  /////// write POINTS ///////
  std::cout << "POINTS " << mn << " double" << std::endl;
  for(long i=0; i<mat.m; i++){
    for(long j=0; j<mat.n; j++){
      std::cout << j << " " << i << " " << "0" << std::endl;
    }
  }
  std::cout << std::endl;
  
  /////// write CELLS ///////
  std::cout << "CELLS " << mn << " " << 2*mn << std::endl;
  for(long i=0; i<mn; i++){
    std::cout << "1 " << i << std::endl;
  }
  std::cout << std::endl;
  
  /////// write CELL_TYPES ///////
  std::cout.unsetf(std::ios::showpos);
  std::cout << "CELL_TYPES " << mn << std::endl;
  std::cout.setf(std::ios::showpos);
  for(long i=0; i<mn; i++){
    std::cout << "1" << std::endl; // VTK_VERTEX = 1
  }
  std::cout << std::endl;
  
  //////// write CELL_DATA ////////
  std::cout.unsetf(std::ios::showpos);
  std::cout << "CELL_DATA " << mn << std::endl;
  std::cout.setf(std::ios::showpos);
  std::cout << "SCALARS value double 1" << std::endl;// SCALAR = 1
  std::cout << "LOOKUP_TABLE default" << std::endl;
  for(long i=0; i<mat.m; i++){
    for(long j=0; j<mat.n; j++){
      std::cout << mat(i,j) << std::endl;
    }
  }
  std::cout << std::endl;
  
  return 0;
}
