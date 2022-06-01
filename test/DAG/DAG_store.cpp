#include <fstream>
#include <iomanip>
#include "ffunc.hpp"

struct DagObj {
  mc::FFSubgraph subgraph;           /*!< subgraph holding list of operations in the DAG */
  std::vector< mc::FFSubgraph > storedSubgraph;
};

int main()
{
  mc::FFGraph dag;
  const unsigned int NX = 4;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &dag );

  const unsigned int NF = 2;
  mc::FFVar F[NF]
    = { X[2]*X[3]-X[0],
        X[0]*pow(exp(X[2]*X[3])+3.,4)+X[1] };
  std::cout << dag;

  DagObj pb;
  pb.subgraph = dag.subgraph( 1, F );
  pb.storedSubgraph.push_back( pb.subgraph ); // saving "subgraph"
  for( auto subgraph : pb.storedSubgraph )
    dag.output( subgraph );

  pb.subgraph = dag.subgraph( 2, F );
  pb.storedSubgraph.push_back( pb.subgraph ); // saving "subgraph"
  for( auto subgraph : pb.storedSubgraph )
    dag.output( subgraph );

  return 0;
}

