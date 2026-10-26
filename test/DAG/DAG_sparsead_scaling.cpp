////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DAG_sparsead_scaling -- how the cost of SFAD and SBAD grows with the size of a SPARSE Jacobian.
//
// The sparse differentiation routines return a Jacobian in coordinate form, nnz entries.  Their cost should
// scale with the work needed to form those entries -- O(nnz x depth) for a banded system -- and not with
// the number of independents times the size of the subgraph.  This driver measures it: three structures,
// doubling sizes, the fitted exponent of time against n printed per method, so that a change to SDFAD or
// SDBAD is judged by a number (2026-10-26).
//
//   TRIDIAGONAL   f_i = x_{i-1} - 2 x_i + x_{i+1} + x_i^3   -- nnz = 3n, the method-of-lines shape;
//                 cone restriction should make both routines linear in n.
//   DENSE ROW     f_0 = sum_i x_i^2, f_i = x_i for i > 0   -- one row touches every column, nnz = 2n - 1.
//                 The sum is a left-deep chain, so the forward cone of x_i holds every addition after it
//                 and yields ONE entry: forward mode stays O(n^2) here for a structural reason, while the
//                 backward cone of f_0 is the chain once, O(n).  The asymmetry the two modes exist for.
//   SHARED        f_i = x_i * s, s = sum_j x_j^2 (one shared subexpression every f_i reads) -- its Jacobian
//                 is DENSE, nnz = n^2, so a quadratic exponent is the floor here, not a defect; measured so
//                 the dense regime is on record and the per-entry cost is visible.
//
// Run with no argument for n = 100 .. 1600; the first argument sets the largest n.
//
// Before cone restriction (2026-10-26):  tridiagonal SFAD n^2.16, SBAD n^2.49 (n = 800: 298 ms, 1024 ms);
//                                        dense row SFAD n^1.97, SBAD n^2.69.
// With it:                               tridiagonal SFAD n^1.25, SBAD n^1.19 (n = 800: 5.8 ms, 5.4 ms);
//                                        dense row SFAD n^2.03 (structural, above), SBAD n^1.37.
////////////////////////////////////////////////////////////////////////////////////////////////////////////
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>
#include "ffunc.hpp"

using namespace mc;

struct Rec { size_t n; double tF, tB; unsigned nnzF, nnzB; };

static double ms( std::chrono::steady_clock::time_point const& t0 )
{ return std::chrono::duration<double,std::milli>( std::chrono::steady_clock::now() - t0 ).count(); }

template <typename BUILD>
static std::vector<Rec> run( char const* name, size_t const nmax, BUILD const& build )
{
  std::vector<Rec> out;
  // Two costs are in each time: the sweep, which the cone restriction made proportional to the cone sizes,
  // and the creation of each entry as a DAG node, inserted into a DAG that grows with n (a set insertion,
  // and the cache behaviour of a large DAG).  The second is per ENTRY and grows slowly with n; a fit over
  // the largest sizes reads mostly it.  So beside the fit: the ratio of each doubling (2 is linear), and
  // the microseconds per entry (flat is linear in nnz).
  std::cout << "\n--- " << name << "\n" << std::setw(8) << "n" << std::setw(12) << "SFAD ms" << std::setw(10) << "nnz"
            << std::setw(8) << "x2" << std::setw(9) << "us/nnz"
            << std::setw(12) << "SBAD ms" << std::setw(10) << "nnz" << std::setw(8) << "x2" << std::setw(9) << "us/nnz" << "\n";
  for( size_t n = 100; n <= nmax; n *= 2 ){
    FFGraph G;
    std::vector<FFVar> X( n ), F;
    for( size_t i = 0; i < n; ++i ) X[i].set( &G );
    build( X, F );
    Rec r{ n, 0., 0., 0, 0 };
    { auto t0 = std::chrono::steady_clock::now();
      auto J = G.SFAD( F, X );
      r.tF = ms( t0 );  r.nnzF = (unsigned)std::get<0>( J ).size(); }
    { auto t0 = std::chrono::steady_clock::now();
      auto J = G.SBAD( F, X );
      r.tB = ms( t0 );  r.nnzB = (unsigned)std::get<0>( J ).size(); }
    auto ratio = [&]( bool fwd ){ if( out.empty() ) return 0.; double const p = fwd? out.back().tF: out.back().tB; return p > 0.? ( fwd? r.tF: r.tB ) / p: 0.; };
    std::cout << std::setw(8) << n << std::setw(12) << std::fixed << std::setprecision(1) << r.tF << std::setw(10) << r.nnzF
              << std::setw(8) << std::setprecision(2) << ratio( true ) << std::setw(9) << std::setprecision(2) << 1e3*r.tF/std::max( 1u, r.nnzF )
              << std::setw(12) << std::setprecision(1) << r.tB << std::setw(10) << r.nnzB
              << std::setw(8) << std::setprecision(2) << ratio( false ) << std::setw(9) << std::setprecision(2) << 1e3*r.tB/std::max( 1u, r.nnzB ) << "\n";
    out.push_back( r );
  }
  // the fitted exponent over EVERY size from 200 on (the 100-point is in the noise): time ~ n^p
  auto slope = [&]( bool fwd ){
    if( out.size() < 3 ) return 0.;
    double sx = 0., sy = 0., sxx = 0., sxy = 0.;  size_t m = 0;
    for( size_t k = 1; k < out.size(); ++k ){
      double const x = std::log( (double)out[k].n ), y = std::log( std::max( fwd? out[k].tF: out[k].tB, 1e-3 ) );
      sx += x; sy += y; sxx += x*x; sxy += x*y; ++m;
    }
    return ( m*sxy - sx*sy ) / ( m*sxx - sx*sx ); };
  std::cout << "  exponent (fit over n >= 200):  SFAD n^" << std::setprecision(2) << slope( true )
            << "   SBAD n^" << slope( false ) << "   -- read with the x2 and us/nnz columns: a rising us/nnz is the DAG growing, not the sweep\n";
  return out;
}

int main( int argc, char** argv )
{
  size_t const nmax = argc > 1? (size_t)std::atoi( argv[1] ): 1600;
  std::cout << std::fixed;
  run( "TRIDIAGONAL  f_i = x_{i-1} - 2 x_i + x_{i+1} + x_i^3", nmax, []( std::vector<FFVar> const& X, std::vector<FFVar>& F ){
    size_t const n = X.size();
    for( size_t i = 0; i < n; ++i ){
      FFVar f = -2.*X[i] + pow( X[i], 3 );
      if( i )       f += X[i-1];
      if( i+1 < n ) f += X[i+1];
      F.push_back( f );
    } } );
  run( "DENSE ROW  f_0 = sum x_i^2,  f_i = x_i", nmax, []( std::vector<FFVar> const& X, std::vector<FFVar>& F ){
    FFVar s = 0.;
    for( auto const& x : X ) s += sqr( x );
    F.push_back( s );
    for( size_t i = 1; i < X.size(); ++i ) F.push_back( X[i] ); } );
  run( "SHARED  f_i = x_i * s,  s = sum x_j^2", nmax, []( std::vector<FFVar> const& X, std::vector<FFVar>& F ){
    FFVar s = 0.;
    for( auto const& x : X ) s += sqr( x );
    for( auto const& x : X ) F.push_back( x * s ); } );
  return 0;
}
