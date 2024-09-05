#undef  USE_CHEB       // <-- whether to perform the decomposition in Chebyshev basis
#define USE_DAG        // <-- whether to define a DAG of the expressions
////////////////////////////////////////////////////////////////////////

#include "spoly.hpp"

#if defined( USE_DAG )
 #include "ffunc.hpp"
 typedef mc::SMon<mc::FFVar const*,mc::lt_FFVar>  t_SMon;
 typedef mc::SPoly<mc::FFVar const*,mc::lt_FFVar> t_SPoly;
#else
 typedef mc::SMon<>  t_SMon;
 typedef mc::SPoly<> t_SPoly;
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  try{

    // Construction of multivariate polynomials
#if defined( USE_CHEB )
    t_SPoly::options.BASIS = t_SPoly::Options::CHEB;
#else
    t_SPoly::options.BASIS = t_SPoly::Options::MONOM;
#endif
    t_SPoly::options.DISPLINE = true;

    size_t const NX = 3;
    t_SPoly X[NX], P;
#if defined( USE_DAG )
    mc::FFGraph DAG;
    mc::FFVar DAGX[NX];
    for( size_t i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
#else
    for( size_t i=0; i<NX; i++ ) X[i].var( i );
#endif
    //P = pow( X[0], 4 ) + pow( X[0], 6 );
    P = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
    //P = pow( X[0], 3 );//2 * sqr( X[0] ) - 1;

    std::cout << "\nSparse multivariate polynomial:\n";
    std::cout << "P = " << P << std::endl;

    // Evaluate polynomial
#if defined( USE_DAG )
    std::map<mc::FFVar const*,double,mc::lt_FFVar> Xval{ { &DAGX[0], 0.5 }, { &DAGX[1], -0.5 }, { &DAGX[2], 0.75 } };
#else
    std::map<unsigned,double> Xval( { { 0, 0.5 }, { 1, -0.5 }, { 2, 0.75 } } );
#endif
    std::cout << "\nSparse multivariate polynomial value:\n";
    std::cout << "P = " << P.eval( Xval ) << std::endl;
    
    // Differentiate polynomial
#if defined( USE_DAG )
    t_SPoly dPdX0 = P.diff( &DAGX[0] );
#else
    t_SPoly dPdX0 = P.diff( 0 );
#endif
    std::cout << "\nSparse multivariate polynomial derivative:\n";
    std::cout << "dPdX0 = " << dPdX0 << std::endl;

    // Evaluate polynomial derivative
    std::cout << "\nSparse multivariate polynomial value:\n";
    std::cout << "dPdX0 = " << dPdX0.eval( Xval ) << std::endl;

#if defined( USE_DAG )
    t_SPoly dPdX1 = P.diff( &DAGX[1] );
#else
    t_SPoly dPdX1 = P.diff( 1 );
#endif
    std::cout << "\nSparse multivariate polynomial derivative:\n";
    std::cout << "dPdX1 = " << dPdX1 << std::endl;

    // Evaluate polynomial derivative
    std::cout << "\nSparse multivariate polynomial value:\n";
    std::cout << "dPdX1 = " << dPdX1.eval( Xval ) << std::endl;

    // Factor polynomial with respect to X0

#if defined( USE_DAG )
    auto&& spmap = P.factor( &DAGX[0] );
#else
    auto&& spmap = P.factor( 0 );
#endif
    std::cout << "\nSparse multivariate polynomial derivative:\n";
    std::cout << "P = ";
    bool first=true;
    for( auto const& [ord,spol] : spmap ){
#if defined( USE_DAG )
      std::cout << (first?"( ":" + ( ") << spol << (ord?" ) "+t_SMon(&DAGX[0],ord).display(t_SPoly::options.BASIS):" )");
#else
      std::cout << (first?"( ":" + ( ") << spol << (ord?" ) "+t_SMon(0,ord).display(t_SPoly::options.BASIS):" )");
#endif
      first = false;
    }
    std::cout << std::endl;

  }

  catch( t_SPoly::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in polynomial expression:" << std::endl
	      << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
