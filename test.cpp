#include "interval.hpp"
typedef mc::Interval I;
 
#include "polimage.hpp"
typedef mc::PolImg<I> PI;
typedef mc::PolVar<I> PV;
int main(){
    mc::FFGraph DAG;
    mc::FFVar X[2]; X[0].set( &DAG ); X[1].set( &DAG );
    mc::FFVar F[2]; F[0] = X[0]/X[1]; F[1] = X[0]/X[1]; 
    mc::PolImg<I> Env;
    I IX[2] = { I(1,5), I(2,6) };
    mc::PolVar<I> PX[2]; PX[0].set( &Env, X[0], IX[0] ); PX[1].set( &Env, X[1], IX[1] );
    mc::PolVar<I> PF[2]; DAG.eval( 2, F, PF, 2, X, PX );
    Env.generate_cuts( 2, PF );
    std::cout << Env;
    std::cout<< DAG;
}