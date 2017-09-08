//=============================================================================
/*! dsysv_check */
void dsysv_check_vector()
{
  cout << "############ check dsysv vector ############" << endl;
  
  srand(unsigned(time(NULL)));
  int N(3);
  
  //// make dsymatrix A  ////
  CPPL::dsymatrix A(N);
  for(int i=0; i<A.n; i++){
    for(int j=0; j<A.n; j++){
      A(i,j) =double( rand() /(RAND_MAX/10) );
    }
  }
  
  //// make dcovector b ////
  CPPL::dcovector b(N);
  for(int i=0; i<b.l; i++){
    b(i) =double( rand() /(RAND_MAX/10) );
  }
  
  //// make A_original and b_original ////
  CPPL::dsymatrix A_original(A);
  CPPL::dcovector b_original(b);
  cout << "A_original=\n" << A_original << endl;
  cout << "b_original=\n" << b_original << endl;
  
  //// solve Ax=b ////
  A.dsysv(b);
  
  //// print A, b and A_original*y ////
  cout << "modified A=\n" << A << endl;
  cout << "x=\n" << b << endl;
  cout << "b_original -A_original*x = (should be zero)\n" << b_original -A_original*b << endl;
}

//=============================================================================
void dsysv_check_matrix()
{
  cout << "############ check dsysv matrix ############" << endl;
  
  srand(unsigned(time(NULL)));
  int N(3);
  

  //// make dsymatrix A  ////
  CPPL::dsymatrix A(N);
  for(int i=0; i<A.n; i++){
    for(int j=0; j<A.n; j++){
      A(i,j) =double( rand() /(RAND_MAX/10) );
    }
  }
   
  //// make dsymatrix B  ////
  CPPL::dgematrix B(N,N);
  for(int i=0; i<B.m; i++){
    for(int j=0; j<B.n; j++){
      B(i,j) =double( rand() /(RAND_MAX/10) );
    }
  }
   
  //// make A_original and B_original ////
  CPPL::dsymatrix A_original(A);
  CPPL::dgematrix B_original(B);
  cout << "A_original=\n" << A_original << endl;
  cout << "B_original=\n" << B_original << endl;
  
  //// solve A*X=B ////
  A.dsysv(B);
  
  //// print A, B and A_original*B ////
  cout << "modified A=\n" << A << endl;
  cout << "X=\n" << B << endl;
  cout << "B_original -A_original*X = (should be zero)\n" << B_original -A_original*B << endl;
}
