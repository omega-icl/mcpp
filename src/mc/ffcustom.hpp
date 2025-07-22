// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FFCUSTOM_HPP
#define MC__FFCUSTOM_HPP

#include <functional>

#include "ffunc.hpp"
#include "ffdep.hpp"
#include "ffinv.hpp"
#include "ffexpr.hpp"
#include "slift.hpp"
#include "spoly.hpp"
#include "mccormick.hpp"
#include "pwlu.hpp"
#include "pwcu.hpp"
#include "supmodel.hpp"
#include "polimage.hpp"

namespace mc
{

//! @brief C++ class defining operations with custom functions as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFCustom is a C++ class for defining operations with custom
//! functions as an external DAG operation in MC++. The template
//! parameter specifies the type for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template< typename T = void >
class FFCustom
: public FFOp
{
private:

  // unique identifier
  int                                                                           _uid;

  // derivative operation
  FFCustom<T>*                                                                  _Deriv;

  // real-valued evaluation
  std::function<std::vector<double>( std::vector<double> const& )>              _DEval;

  // interval-valued evaluation
  std::function<std::vector<T>( std::vector<T> const& )>                        _IEval;

  // McCormick-valued evaluation
  std::function<std::vector<McCormick<T>>( std::vector<McCormick<T>> const& )>  _MCEval;

  // default evaluation returning empty vector
  template <typename U>
  static std::vector<U> _Error
    ( std::vector<U> const& vVar )
    { return std::vector<U>(); }

public:

  // Default constructor
  FFCustom
    ()
    : FFOp     ( EXTERN ),
      _uid     ( 0 ),
      _Deriv   ( nullptr ),
      _DEval   ( _Error<double> ),
      _IEval   ( _Error<T> ),
      _MCEval   ( _Error<McCormick<T>> )
    {}

  // Copy constructor
  FFCustom
    ( FFCustom const& other )
    : FFOp     ( other ),
      _uid     ( other._uid ),
      _DEval   ( other._DEval ),
      _IEval   ( other._IEval ),
      _MCEval   ( other._MCEval )
    {
#ifdef MC__FFCUSTOM_TRACE
      std::cout << "FFCustom::copy constructor\n";
#endif
      _Deriv = (other._Deriv? new FFCustom<T>( *other._Deriv ): nullptr);
    }

  // Destructor
  virtual ~FFCustom
    ()
    {
      delete _Deriv;
    }

  // Define operation
  FFVar** operator()
    ( size_t const nDep, std::vector<FFVar> const& vVar, int const uid )
    {
#ifdef MC__FFCustom_CHECK
      assert( nDep );
#endif
      _uid = uid; // unique identifier - used for comparison
      return this->insert_external_operation( *this, nDep, vVar.size(), vVar.data() );
    }

  FFVar& operator()
    ( size_t const iDep, size_t const nDep, std::vector<FFVar> const& vVar, int const uid )
    {
#ifdef MC__FFCustom_CHECK
      assert( iDep < nDep );
#endif
      _uid = uid; // unique identifier - used for comparison
      return *(this->insert_external_operation( *this, nDep, vVar.size(), vVar.data() )[iDep]);
    }

  FFVar& operator()
    ( std::vector<FFVar> const& vVar, int const uid )
    {
      _uid = uid; // unique identifier - used for comparison
      return *(this->insert_external_operation( *this, 1, vVar.size(), vVar.data() )[0]);
    }

  // Set custom evaluation functions in double arithmetic
  void set_eval
    ( std::function<std::vector<double>( std::vector<double> const& )> const& DEval )
    {
      _DEval = DEval;
    }

  // Set custom evaluation functions in interval arithmetic
  void set_eval
    ( std::function<std::vector<T>( std::vector<T> const& )> const& IEval )
    {
      _IEval = IEval;
    }

  // Set custom evaluation functions in McCormick arithmetic
  void set_eval
    ( std::function<std::vector<McCormick<T>>( std::vector<McCormick<T>> const& )> const& MCEval )
    {
      _MCEval = MCEval;
    }

  // Set custom derivative function
  void set_deriv
    ( FFCustom<T> const& Deriv, int const uid )
    {
      _Deriv = new FFCustom<T>( Deriv );
      _Deriv->_uid = uid;
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
//      else if( idU == typeid( fadbad::F<FFVar> ) )
//        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
//      else if( idU == typeid( FFInv ) )
//        return _eval( nRes, static_cast<FFInv*>(vRes), nVar, static_cast<FFInv const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
//      else if( idU == typeid( fadbad::F<double> ) )
//        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( T ) )
        return eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );
//      else if( idU == typeid( McCormick<T> ) )
//        return _eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );
//      else if( idU == typeid( SupVar<PWCU> ) )
//        return _eval( nRes, static_cast<SupVar<PWCU>*>(vRes), nVar, static_cast<SupVar<PWCU> const*>(vVar), mVar );
//      else if( idU == typeid( SupVar<PWLU> ) )
//        return _eval( nRes, static_cast<SupVar<PWLU>*>(vRes), nVar, static_cast<SupVar<PWLU> const*>(vVar), mVar );
//      else if( idU == typeid( McCormick<SupVar<PWCU>> ) )
//        return _eval( nRes, static_cast<McCormick<SupVar<PWCU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWCU>> const*>(vVar), mVar );
//      else if( idU == typeid( McCormick<SupVar<PWLU>> ) )
//        return _eval( nRes, static_cast<McCormick<SupVar<PWLU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWLU>> const*>(vVar), mVar );
//      else if( idU == typeid( PolVar<T> ) )
//        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFCustom::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, double* vRes, size_t const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, T* vRes, size_t const nVar, T const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, McCormick<T>* vRes, size_t const nVar, McCormick<T> const* vVar, unsigned const* mVar )
    const;

//  void eval
//    ( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar, unsigned const* mVar )
//    const;

//  void eval
//    ( size_t const nRes, fadbad::F<double>* vRes, size_t const nVar, fadbad::F<double> const* vVar, unsigned const* mVar )
//    const;

  void eval
    ( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

//  void eval
//    ( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar, unsigned const* mVar )
//    const;

  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
//      if( idU == typeid( T ) )
//        return reval( nRes, static_cast<T const*>(vRes), nVar, static_cast<T*>(vVar) );
//      else if( idU == typeid( PolVar<T> ) )
//        return reval( nRes, static_cast<PolVar<T> const*>(vRes), nVar, static_cast<PolVar<T>*>(vVar) );

      throw std::runtime_error( "FFCustom::reval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

//  bool reval
//    ( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
//    const;

//  bool reval
//    ( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
//    const;

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
    { return "Custom[" + std::to_string(_uid) + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, double* pRes, size_t const nVar, double const* pVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: double\n";
#endif

  std::vector<double> const& vRes = _DEval( std::vector<double>(pVar,pVar+nVar) );
  if( vRes.empty() ) //throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    throw std::runtime_error( "FFCustom::eval ** No evaluation method for type "+std::string(typeid(*pRes).name())+"\n" );
  for( unsigned j=0; j<nRes; ++j ) pRes[j] = vRes[j];  
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, T* pRes, size_t const nVar, T const* pVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: T\n";
#endif

  std::vector<T> const& vRes = _IEval( std::vector<T>(pVar,pVar+nVar) );
  if( vRes.empty() ) //throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    throw std::runtime_error( "FFCustom::eval ** No evaluation method for type "+std::string(typeid(*pRes).name())+"\n" );
  for( unsigned j=0; j<nRes; ++j ) pRes[j] = vRes[j];  
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, McCormick<T>* pRes, size_t const nVar, McCormick<T> const* pVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: McCormick<T>\n";
#endif

  std::vector<McCormick<T>> const& vRes = _MCEval( std::vector<McCormick<T>>(pVar,pVar+nVar) );
  if( vRes.empty() ) //throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    throw std::runtime_error( "FFCustom::eval ** No evaluation method for type "+std::string(typeid(*pRes).name())+"\n" );
  for( unsigned j=0; j<nRes; ++j ) pRes[j] = vRes[j];  
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: FFDep\n";
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: FFExpr\n";
#endif

  switch( FFExpr::options.LANG ){
   case FFExpr::Options::DAG:
    for( unsigned j=0; j<nRes; ++j ){
      std::ostringstream os; os << name() << "[" << j << "]";
      vRes[j] = FFExpr::compose( os.str(), nVar, vVar );
    }
    break;
   case FFExpr::Options::GAMS:
   default:
    throw typename FFExpr::Exceptions( FFExpr::Exceptions::UNDEF );
  }
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: FFVar\n";
#endif
  FFVar** ppRes = this->insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned i=0; i<nRes; ++i )
    vRes[i] = *(ppRes[i]);
}
/*
template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, fadbad::F<double>* vRes, size_t const nVar, fadbad::F<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: fadbad::F<double>\n";
#endif
#ifdef MC__FFCUSTOM_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  std::vector<double> vVarVal( nVar );
  if( vVarVal.size() < nVar ) vVarVal.resize( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  double vResVal = 0.;
  _eval( 1, &vResVal, nVar, vVarVal.data(), nullptr );
  vRes[0] = vResVal;
  for( unsigned i=0; i<nVar; ++i )
    vRes[0].setDepend( vVar[i] );

  for( unsigned j=0; j<vRes[0].size(); ++j ){
    vRes[0][j] = 0.;
    for( unsigned i=0; i<nVar; ++i ){
      if( vVar[i][j] == 0. ) continue;
      vRes[0][j] += (_nCoef<nVar? _ptrCoef[0]*vVar[i][j]: _ptrCoef[i]*vVar[i][j]);
    }
  }
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFCUSTOM_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
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
      vRes[0][j] += (_nCoef<nVar? _ptrCoef[0]*vVar[i][j]: _ptrCoef[i]*vVar[i][j]);
    }
  }
}
*/
template< typename T >
inline void
FFCustom<T>::deriv
( unsigned const nRes, FFVar const* pRes, unsigned const nVar, FFVar const* pVar, FFVar** pDer )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::deriv\n";
#endif

  if( !_Deriv || _Deriv->type != FFOp::EXTERN ) //throw typename FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    throw std::runtime_error( "FFCustom::deriv ** No differentiation method\n" );
  FFVar** ppDer = (*_Deriv)( nVar*nRes, std::vector<FFVar>(pVar,pVar+nVar), _Deriv->_uid );

  for( unsigned j=0; j<nRes; ++j )
    for( unsigned i=0; i<nVar; ++i )
      pDer[j][i] = *(ppDer[j*nVar+i]); // [i*nRes+j];
}

template< typename T >
inline void
FFCustom<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::eval: SLiftVar\n";
#endif
  return vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline bool
FFCustom<T>::lt
( FFOp const* other )
const
{
#ifdef MC__FFCUSTOM_TRACE
  std::cout << "FFCustom::lt\n";
#endif
  auto const* oprecast = dynamic_cast<FFCustom<T> const*>(other);
  return( _uid < oprecast->_uid );
}

} // end namescape mc

#endif
