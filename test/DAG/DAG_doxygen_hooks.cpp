// =============================================================================
// DAG_doxygen_hooks.cpp
// The example of the documentation of ffunc.hpp on the tangent and adjoint
// hooks of external operations: the 2-norm as an external operation, with its
// directional derivative and its adjoint as two further external operations.
//
// Build (assuming MC++ headers are on the include path):
//   g++ -std=c++17 -O2 -I<path-to-mc++> -o DAG_doxygen_hooks DAG_doxygen_hooks.cpp
// =============================================================================

#include <cmath>
#include <iostream>
#include "ffunc.hpp"

namespace mc
{
// Directional derivative of the 2-norm: operands (x, t), result (x.t)/||x||
class FFnorm2Tan : public FFOp
{
public:
  FFnorm2Tan() : FFOp( (int)EXTERN ) {}
  FFVar& operator()( unsigned const nVar, FFVar const* pVar ) const
    { return *insert_external_operation( *this, 1, nVar, pVar )[0]; }
  void feval( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
              void const* vVar, unsigned const* mVar ) const
    {
      if( idU == typeid(FFVar) )
        { *static_cast<FFVar*>(vRes) = operator()( nVar, static_cast<FFVar const*>(vVar) ); return; }
      if( idU != typeid(double) ) throw std::runtime_error( "FFnorm2Tan: no evaluation in this type" );
      double const* x = static_cast<double const*>(vVar), *t = x + nVar/2;
      double xt = 0., xx = 0.;
      for( unsigned i=0; i<nVar/2; ++i ){ xt += x[i]*t[i]; xx += x[i]*x[i]; }
      *static_cast<double*>(vRes) = xt / std::sqrt(xx);
    }
  std::string name() const { return "NORM2TAN"; }
  bool commutative() const { return false; }
};

// Adjoint of the 2-norm: operands (x, w), results w x[i]/||x||
class FFnorm2Adj : public FFOp
{
public:
  FFnorm2Adj() : FFOp( (int)EXTERN ) {}
  FFVar** operator()( unsigned const nVar, FFVar const* pVar ) const
    { return insert_external_operation( *this, nVar-1, nVar, pVar ); }
  void feval( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
              void const* vVar, unsigned const* mVar ) const
    {
      if( idU == typeid(FFVar) )
        {
          FFVar** ppRes = operator()( nVar, static_cast<FFVar const*>(vVar) );
          for( unsigned i=0; i<nRes; ++i ) static_cast<FFVar*>(vRes)[i] = *ppRes[i];
          return;
        }
      if( idU != typeid(double) ) throw std::runtime_error( "FFnorm2Adj: no evaluation in this type" );
      double const* x = static_cast<double const*>(vVar), w = x[nRes];
      double xx = 0.;
      for( unsigned i=0; i<nRes; ++i ) xx += x[i]*x[i];
      for( unsigned i=0; i<nRes; ++i ) static_cast<double*>(vRes)[i] = w * x[i] / std::sqrt(xx);
    }
  std::string name() const { return "NORM2ADJ"; }
  bool commutative() const { return false; }
};

// The 2-norm, with its tangent and adjoint hooks
class FFnorm2 : public FFOp
{
public:
  FFnorm2() : FFOp( (int)EXTERN ) {}
  FFVar& operator()( unsigned const nVar, FFVar const* pVar ) const
    { return *insert_external_operation( *this, 1, nVar, pVar )[0]; }
  void feval( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
              void const* vVar, unsigned const* mVar ) const
    {
      if( idU == typeid(FFVar) )
        { *static_cast<FFVar*>(vRes) = operator()( nVar, static_cast<FFVar const*>(vVar) ); return; }
      if( idU != typeid(double) ) throw std::runtime_error( "FFnorm2: no evaluation in this type" );
      double const* x = static_cast<double const*>(vVar);
      double xx = 0.;
      for( unsigned i=0; i<nVar; ++i ) xx += x[i]*x[i];
      *static_cast<double*>(vRes) = std::sqrt(xx);
    }
  std::string name() const { return "NORM2"; }
  bool commutative() const { return true; }

  // Forward mode: the tangent of the result for given tangents of the operands
  bool has_tangent( unsigned const nSweep ) const { return true; }
  void tangent( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
                FFVar const* const* vTan, FFVar* vResTan ) const
    {
      std::vector<FFVar> vOpd( vVar, vVar+nVar );             // the operands x ...
      for( unsigned i=0; i<nVar; ++i )                        // ... then their tangents t,
        vOpd.push_back( vTan[i]? *vTan[i]: FFVar(0.) );       //     zero where the pointer is null
      vResTan[0] = FFnorm2Tan()( vOpd.size(), vOpd.data() );
    }

  // Backward mode: the contributions to the adjoints of the operands for a given adjoint of the result
  bool has_adjoint( unsigned const nSweep ) const { return true; }
  void adjoint( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
                FFVar const* const* vResAdj, FFVar* vVarAdj ) const
    {
      std::vector<FFVar> vOpd( vVar, vVar+nVar );             // the operands x ...
      vOpd.push_back( *vResAdj[0] );                          // ... then the adjoint w of the result
      FFVar** ppAdj = FFnorm2Adj()( vOpd.size(), vOpd.data() );
      for( unsigned i=0; i<nVar; ++i ) vVarAdj[i] = *ppAdj[i];
    }
};
}

int main()
{
  mc::FFGraph DAG;
  const unsigned NX = 3;
  std::vector<mc::FFVar> X = DAG.add_vars( NX ), D = DAG.add_vars( NX );
  mc::FFnorm2 norm2;
  std::vector<mc::FFVar> F = { norm2( NX, X.data() ) * X[0] };

  std::vector<mc::FFVar> dFdX_F = DAG.FAD( F, X );     // one NORM2TAN operation per variable
  std::vector<mc::FFVar> dFdX_B = DAG.BAD( F, X );     // one NORM2ADJ operation
  std::vector<mc::FFVar> dFdX_D = DAG.DFAD( F, X, D ); // one NORM2TAN operation
  std::cout << DAG;

  std::vector<double> x = { 1., -2., 2. }, d = { 0.5, 1., -1. }, dF( NX ), dir( 1 );
  DAG.eval( dFdX_F, dF, X, x );
  std::cout << "forward:     " << dF[0] << "  " << dF[1] << "  " << dF[2] << std::endl;
  DAG.eval( dFdX_B, dF, X, x );
  std::cout << "backward:    " << dF[0] << "  " << dF[1] << "  " << dF[2] << std::endl;
  std::vector<mc::FFVar> XD = X; XD.insert( XD.end(), D.begin(), D.end() );
  std::vector<double> xd = x; xd.insert( xd.end(), d.begin(), d.end() );
  DAG.eval( dFdX_D, dir, XD, xd );
  std::cout << "directional: " << dir[0] << std::endl;

  // d(||x|| x0)/dx = ( ||x|| + x0^2/||x||, x0 x1/||x||, x0 x2/||x|| ) = ( 10/3, -2/3, 2/3 ) at x = (1,-2,2)
  bool const ok = std::fabs( dF[0] - 10./3. ) < 1e-12 && std::fabs( dF[1] + 2./3. ) < 1e-12
               && std::fabs( dF[2] - 2./3. ) < 1e-12
               && std::fabs( dir[0] - ( dF[0]*d[0] + dF[1]*d[1] + dF[2]*d[2] ) ) < 1e-12;
  std::cout << "\n=== Results: " << ( ok? "passed": "failed" ) << " ===\n";
  return ok? 0: 1;
}
