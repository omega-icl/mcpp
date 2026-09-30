#include <iostream>
#include "ocbase.hpp"
int main()
{
  mc::FFGraph DAG;
  mc::FFVar t = DAG.add_var("t"), X = DAG.add_var("x"), U = DAG.add_var("u"), S = DAG.add_var("s");
  mc::FFPartial OpP;  mc::FFIntegral OpI;  mc::FFEval OpE;
  mc::FFVar E  = OpP( X, t ) - U * X * X;
  mc::FFVar F1 = OpI( X * X, t );
  mc::FFVar F2 = OpE( X * X, t, 3.7 );
  auto show = [&]( char const* name, mc::FFVar const* f ){
    mc::FFSubgraph sg = DAG.subgraph( 1, f ); auto ex = mc::FFExpr::subgraph( &DAG, sg );
    std::cout << "    " << name << " = " << ex[0] << "\n"; };
  auto attempt = [&]( char const* what, auto fn ){
    std::cout << "  " << what << ": ";
    try{ mc::FFVar* r = fn(); std::cout << "OK\n"; show("result", r); delete[] r; }
    catch( mc::FFBase::Exceptions& e ){ std::cout << "THROWS  " << e.what() << "\n"; }
    catch( ... ){ std::cout << "THROWS (other)\n"; } };
  attempt( "DFAD  d E/dx . s     (evolution, thru OpP)", [&]{ return DAG.DFAD( 1, &E,  1, &X, &S ); } );
  attempt( "FAD   d E/du         (evolution, thru OpP)", [&]{ return DAG.FAD ( 1, &E,  1, &U ); } );
  attempt( "BAD   d E/du         (evolution, thru OpP)", [&]{ return DAG.BAD ( 1, &E,  1, &U ); } );
  attempt( "DFAD  d F1/dx . s    (integral output)    ", [&]{ return DAG.DFAD( 1, &F1, 1, &X, &S ); } );
  attempt( "FAD   d F1/dx        (integral output)    ", [&]{ return DAG.FAD ( 1, &F1, 1, &X ); } );
  attempt( "DFAD  d F2/dx . s    (point-value output) ", [&]{ return DAG.DFAD( 1, &F2, 1, &X, &S ); } );
  attempt( "FAD   d F2/dx        (point-value output) ", [&]{ return DAG.FAD ( 1, &F2, 1, &X ); } );
  mc::FFVar G = U * X * X;   // control: no external op at all
  attempt( "DFAD  d(u x^2)/dx . s  (no external op)  ", [&]{ return DAG.DFAD( 1, &G, 1, &X, &S ); } );
  return 0;
}
