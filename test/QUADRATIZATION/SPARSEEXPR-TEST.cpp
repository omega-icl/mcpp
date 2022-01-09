#include <iomanip>

#include "sparseexpr.hpp"

///////////////////////////////////////////////////////////////////////////////
int main()
///////////////////////////////////////////////////////////////////////////////
{
  try{
    mc::FFGraph DAG;
    const unsigned NX = 4, NF = 1;
    mc::FFVar X[NX];
    for( unsigned i(0); i<NX; i++ ) X[i].set( &DAG );
    mc::FFVar F[NF];
    //F[0] = pow( X[0] + X[1], 3 );
    //F[0] = pow( X[0] + 1 / sqr( X[1] ), 3 );
    //F[0] = 4*sqr(X[0]) - 2.1*pow(X[0],4) + pow(X[0],6)/3 + X[0]*X[1] - 4*sqr(X[1]) + 4*pow(X[1],4);
    F[0] = 250*exp(X[2])*X[0] + 250*pow(X[3],0.6)*X[1];
    //F[1] = exp( 2 * sqr( X[1] ) - 1 );
    std::cout << DAG;

    mc::SparseExpr::SPolyExpr::options.BASIS = mc::SparseExpr::SPolyExpr::Options::MONOM;
    mc::SparseEnv SPE( &DAG );
    SPE.options.LIFTDIV = true;//false;//
    SPE.options.LIFTIPOW = true;//false;//

    SPE.process( NF, F );

    std::cout << std::endl << SPE.Var().size() << " participating variables: ";
    for( auto&& var : SPE.Var() ) std::cout << var << " ";
    std::cout << std::endl;
    std::cout << std::endl << SPE.Aux().size() << " auxiliary variables: ";
    for( auto&& aux : SPE.Aux() ) std::cout << *aux.first << "->" << *aux.second << " ";
    std::cout << std::endl;
    std::cout << std::endl << SPE.Poly().size() << " polynomial constraints: " << std::endl;
    for( auto&& expr : SPE.Poly() ) DAG.output( DAG.subgraph( 1, &expr ) );
    //std::cout << std::endl;
    std::cout << std::endl << SPE.Trans().size() << " transcendental constraints: " << std::endl;
    for( auto&& expr : SPE.Trans() ) DAG.output( DAG.subgraph( 1, &expr ) );

//    std::vector<mc::SPolyExpr> SPVar;
//    for( auto&& var : SPE.Var() ) SPVar.push_back( var );
//    std::vector<mc::SPolyExpr> SPPoly( SPE.Poly().size() );
//    DAG.eval( SPE.Poly().size(), SPE.Poly().data(), SPPoly.data(), SPE.Var().size(), SPE.Var().data(), SPVar.data() );
//    unsigned int count = 0;
//    for( auto&& expr : SPPoly )
//      std::cout << std::endl << "Polynomial constraint #" << ++count << ":" << expr << std::endl;
//    std::cout << std::endl;

//    mc::QuadEnv QF( &DAG );
//    mc::QuadEnv::options.REDUC = mc::QuadEnv::Options::ALL;
//    QF.process( SPPoly.size(), SPPoly.data() );
//    std::cout << std::endl << "Sparse quadratic form: " << QF << std::endl;

    return 0;
  }

  catch( mc::FFGraph::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

