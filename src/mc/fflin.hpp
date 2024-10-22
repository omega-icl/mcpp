// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FFLIN_HPP
#define MC__FFLIN_HPP

#include "ffunc.hpp"
#include "ffdep.hpp"
#include "ffinv.hpp"
#include "ffexpr.hpp"
#include "slift.hpp"
//#include "spoly.hpp"
//#include "mccormick.hpp"
#include "mcfadbad.hpp"

namespace mc
{

//! @brief C++ class defining linear combinations as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFSPoly is a C++ class for defining a linear combination as an
//! external DAG operation in MC++. The template parameter specifies the
//! type for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFLin
: public FFOp
{
private:

  // Weigth matrix
  double const* _wVar;

public:

  void set
    ( double const* wVar )
    {
      _wVar = wVar;
      data = const_cast<double*>(_wVar); // for operation comparison
    }

  // Default constructor
  FFLin
    ()
    : FFOp( EXTERN )
    {}
    
  // Copy constructor
  FFLin
    ( FFLin const& Op )
    : FFOp( Op ),
      _wVar( Op._wVar )
    {}

  // Define operation
  FFVar& operator()
    ( size_t const nVar, FFVar const* pVar, double const* wVar=nullptr )
    {
      set( wVar ); // no local copy - make sure wVar isn't going out of scope!
      return **insert_external_operation( *this, 1, nVar, pVar );
    }
    
  FFVar& operator()
    ( size_t const nVar, FFVar const*const* ppVar, double const* wVar=nullptr )
    {
      set( wVar ); // no local copy - make sure wVar isn't going out of scope!
      return **insert_external_operation( *this, 1, nVar, ppVar );
    }

  FFVar& operator()
    ( std::vector<FFVar> const& vVar, std::vector<double> const& wVar=std::vector<double>() )
    {
      set( !wVar.empty()? wVar.data(): nullptr ); // no local copy - make sure wVar isn't going out of scope!
      return **insert_external_operation( *this, 1, vVar.size(), vVar.data() );
    }
    
  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<FFVar> ) )
        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( T ) )
        return eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );

      throw std::runtime_error( "FFLin::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }
    
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      if( idU == typeid( T ) )
       return reval( nRes, static_cast<T const*>(vRes), nVar, static_cast<T*>(vVar) );
       //return true;

      throw std::runtime_error( "FFLin::reval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  template< typename U >
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::eval: generic\n"; 
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif
      if( !_wVar ){
        for( size_t i=0; i<nVar; ++i )
          if( !i ) vRes[0]  = vVar[0];
          else     vRes[0] += vVar[i];
      }
      else{
        for( size_t i=0; i<nVar; ++i )
          if( !i ) vRes[0]  = _wVar[0] * vVar[0];
          else     vRes[0] += _wVar[i] * vVar[i];     
      }
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::eval: FFVar\n"; 
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif
      vRes[0] = **insert_external_operation( *this, 1, nVar, vVar );
    }

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::eval: FFLin\n"; 
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif
      vRes[0] = 0;
      for( size_t i=0; i<nVar; ++i ) vRes[0] += vVar[i];
      vRes[0].update( FFDep::TYPE::L );
    }

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::eval: fadbad::F<T>\n";
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif

      std::vector<FFVar> vVarVal( nVar );
      for( unsigned i=0; i<nVar; ++i )
        vVarVal[i] = vVar[i].val();
      FFVar vResVal;
      eval( 1, &vResVal, nVar, vVarVal.data(), nullptr );
      vRes[0] = vResVal;
      for( unsigned i=0; i<nVar; ++i )
        vRes[0].setDepend( vVar[i] );

      for( unsigned j=0; j<vRes[0].size(); ++j ){
        vRes[0][j] = 0.;
        for( unsigned i=0; i<nVar; ++i ){
          if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
          vRes[0][j] += _wVar? _wVar[i]*vVar[i][j]: vVar[i][j];
        }
      }
    }

  void eval
    ( unsigned const nRes, fadbad::F<double>* vRes, unsigned const nVar, fadbad::F<double> const* vVar,
      unsigned const* mVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::eval: fadbad::F<double>\n";
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif

      std::vector<double> vVarVal( nVar );
      if( vVarVal.size() < nVar ) vVarVal.resize( nVar );
      for( unsigned i=0; i<nVar; ++i )
        vVarVal[i] = vVar[i].val();
      double vResVal = 0.;
      eval( 1, &vResVal, nVar, vVarVal.data(), nullptr );
      vRes[0] = vResVal;
      for( unsigned i=0; i<nVar; ++i )
        vRes[0].setDepend( vVar[i] );

      for( unsigned j=0; j<vRes[0].size(); ++j ){
        vRes[0][j] = 0.;
        for( unsigned i=0; i<nVar; ++i ){
          if( vVar[i][j] == 0. ) continue;
          vRes[0][j] += _wVar? _wVar[i]*vVar[i][j]: vVar[i][j];
        }
      }
    }

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::deriv: FFVar\n";
#endif
#ifdef MC__FFLIN_CHECK
      assert( nRes == 1 );
#endif

      for( unsigned i=0; i<nVar; ++i )
        vDer[0][i] = _wVar? _wVar[i]: 1.;
    }

  bool reval
    ( unsigned const nRes, T const* vRes, unsigned const nVar, T* vVar )
    const
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::reval: T\n";
#endif
      return true;
    }

  // Properties
  std::string name
    ()
    const
    { return "Sum"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
/*
  // Ordering
  bool lt
    ( FFOp const* op )
    const
    {
#ifdef MC__FFSPOLY_TRACE
      std::cout << "FFLin<T>::lt\n";
#endif

      return lt_SPoly<mc::lt_FFVar>()( _SPoly, dynamic_cast<FFSPoly<T> const*>(op)->_SPoly );
    }

  // Insert operation in DAG
  static FFVar insert
    ( t_SPoly const& SPoly, FFGraph* dag )
    {
      if( !dag ) return FFVar();

      std::map< FFVar const*, FFVar, lt_FFVar > x;
      for( auto const& var : SPoly.setvar() ){
#ifdef MC__FFSPOLY_DEBUG
        std::cout << "vVar[" << ivar << "] = " << *var << std::endl;
#endif
        x[var] = *var;
      }
      return SPoly.eval( x );
    }
*/
};
/*
//! @brief C++ class defining linear combinations as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFSPoly is a C++ class for defining a linear combination as an
//! external DAG operation in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFLin
: public FFOp
{

private:
  // Sparse polynomial
  mutable t_SPoly   _SPoly;

public:
  // Default constructor
  FFSPoly
    ()
    : FFOp( EXTERN ),
      _SPoly()
    {}

  // Destructor
  virtual ~FFSPoly
    ()
    {}

  // Copy constructor
  FFSPoly
    ( FFSPoly<T> const& Op )
    : FFOp( Op )
    {
#ifdef MC__FFSPOLY_TRACE
      std::cout << "FFSPoly::copy constructor\n";
#endif
      _SPoly = Op._SPoly;
    }

  // Define operation
  FFVar& operator()
    ( t_SPoly const& SPoly, bool insert=false )
    const
    {
      _SPoly = SPoly;
      data = nullptr;
      owndata = false;
#ifdef MC__FFSPOLY_TRACE
      std::cerr << "SPoly address in DAG: " << data << std::endl;
#endif
      return *insert_external_operation( *this, 1, _SPoly.setvar() )[0];
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<FFVar> ) )
        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFInv ) )
        return eval( nRes, static_cast<FFInv*>(vRes), nVar, static_cast<FFInv const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( fadbad::B<double> ) )
        return eval( nRes, static_cast<fadbad::B<double>*>(vRes), nVar, static_cast<fadbad::B<double> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( T ) )
        return eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );
      else if( idU == typeid( McCormick<T> ) )
        return eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );

      throw std::runtime_error( "FFSPoly::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  template< typename G >
  void eval
    ( unsigned const nRes, G* vRes, unsigned const nVar, G const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

  // Ordering
  bool lt
    ( FFOp const* op )
    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << this; return "SPoly[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }

  // Insert operation in DAG
  static FFVar insert
    ( t_SPoly const& SPoly, FFGraph* dag )
    {
      if( !dag ) return FFVar();

      std::map< FFVar const*, FFVar, lt_FFVar > x;
      for( auto const& var : SPoly.setvar() ){
#ifdef MC__FFSPOLY_DEBUG
        std::cout << "vVar[" << ivar << "] = " << *var << std::endl;
#endif
        x[var] = *var;
      }
      return SPoly.eval( x );
    }
};

template< typename T >
template< typename G >
inline void
FFSPoly<T>::eval
( unsigned const nRes, G* vRes, unsigned const nVar, G const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::eval: generic\n";
#endif
#ifdef MC__FFSPOLY_CHECK
  assert( nRes == 1 && nVar == _SPoly.nvar() );
#endif

  std::map< FFVar const*, G, lt_FFVar > x; // mkae thread_local static?
  unsigned ivar = 0;
  for( auto const& var : _SPoly.setvar() ){
#ifdef MC__FFSPOLY_DEBUG
    std::cout << "vVar[" << ivar << "] = " << *var << std::endl;
#endif
    x[var] = vVar[ivar++];
  }
  vRes[0] = _SPoly.eval( x );
}

template< typename T >
inline void
FFSPoly<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::eval: FFVar\n";
#endif
#ifdef MC__FFSPOLY_CHECK
  assert( nRes == 1 && nVar == _SPoly.nvar() );
#endif

  vRes[0] = operator()( _SPoly );
}

template< typename T >
inline void
FFSPoly<T>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::eval: fadbad::F<FFVar>\n";
#endif
#ifdef CRONOS__FFODE_CHECK
  assert( nRes == 1 && nVar == _SPoly.nvar() );
#endif

  FFSPoly<T> ResDer;
  vRes[0] = operator()( _SPoly );
  std::vector<FFVar> vDer( nVar );
  for( unsigned i=0; i<nVar; ++i ){
    vRes[0].setDepend( vVar[i] );
    vDer[i] = ResDer( _SPoly.diff( &vVar[i].val() ) );
  }
  for( unsigned j=0; j<vRes[0].size(); ++j ){
    vRes[0][j] = 0.;
    for( unsigned i=0; i<nVar; ++i ){
      if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
      vRes[0][j] += vDer[i] * vVar[i][j];
    }
  }
}

template< typename T >
inline void
FFSPoly<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::eval: SLiftVar\n";
#endif
#ifdef MC__FFSPOLY_CHECK
  assert( nRes == 1 && nVar == _SPoly.nvar() );
#endif

  if( vVar->env()->options.KEEPFACT || vVar->env()->options.LIFTIPOW )
    return vVar->env()->lift( nRes, vRes, nVar, vVar );
  
  std::map< FFVar const*, SLiftVar, lt_FFVar > x; // mkae thread_local static?
  unsigned ivar = 0;
  for( auto const& var : _SPoly.setvar() ){
#ifdef MC__FFSPOLY_DEBUG
    std::cout << "vVar[" << ivar << "] = " << *var << std::endl;
#endif
    x[var] = vVar[ivar++];
  }
  vRes[0] = _SPoly.eval( x );
}

template< typename T >
inline void
FFSPoly<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::deriv\n";
#endif
#ifdef MC__FFSPOLY_CHECK
  assert( nRes == 1 && nVar == _SPoly.nvar() );
#endif

  FFSPoly<T> ResDer;
  for( unsigned i=0; i<nVar; ++i ){
#ifdef MC__FFSPOLY_DEBUG
    std::cout << "&_SPoly: " << &ResDer._SPoly << std::endl;
    auto&& SPolyDer = _SPoly.diff( &vVar[i] );
    std::cout << "SPolyDer: " << SPolyDer << std::endl;
    vDer[0][i] = ResDer( SPolyDer );
#else
    vDer[0][i] = ResDer( _SPoly.diff( &vVar[i] ) );
#endif
  }
}

template< typename T >
inline bool
FFSPoly<T>::lt
( FFOp const* op )
const
{
#ifdef MC__FFSPOLY_TRACE
  std::cout << "FFSPoly<T>::lt\n";
#endif

  return lt_SPoly<mc::lt_FFVar>()( _SPoly, dynamic_cast<FFSPoly<T> const*>(op)->_SPoly );
}
*/
} // end namescape mc

#endif
