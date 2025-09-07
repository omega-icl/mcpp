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
#include "spoly.hpp"
#include "mccormick.hpp"
#include "pwlu.hpp"
#include "pwcu.hpp"
#include "supmodel.hpp"
#include "polimage.hpp"

namespace mc
{

//! @brief C++ class defining linear combinations as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFLin is a C++ class for defining a linear combination as an
//! external DAG operation in MC++. The template parameter specifies the
//! type for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template< typename T = void >
class FFLin
: public FFOp
{
private:
  // bias coefficient value
  double        _valBias;
  // pointer to the vector of weight coefficients
  size_t        _nCoef;
  // pointer to the vector of weight coefficients
  double*       _ptrCoef;
  // whether this class owns the coefficients
  bool          _ownCoef;

  // set the object and related operation in DAG
  FFVar& _set
    ( size_t const nVar, FFVar const* pVar, double const* pCoef, double const& Bias, int policy )
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::set vector\n";
      std::cout << "   Bias = " << Bias << std::endl;
      std::cout << "   Coef = [ ";
      for( unsigned i=0; i<nVar; ++i ) std::cout << pCoef[i] << " ";
      std::cout << "]\n";
#endif
      if( _ownCoef && _ptrCoef ) delete _ptrCoef;
      _ownCoef = ( policy>0? true: false );
      this->owndata = false;
      this->data = nullptr;
      _ptrCoef = const_cast<double*>(pCoef);
      _nCoef   = nVar;
      _valBias = Bias;
      
      FFVar* pRes = this->insert_external_operation( *this, 1, nVar, pVar )[0];
      
      _ownCoef = false;
      FFOp* pOp = pRes->opdef().first;
      if( policy > 0 )
        _ptrCoef = static_cast<FFLin<T>*>(pOp)->_ptrCoef; // set pointer to DAG copy
      else if( policy < 0 )
        static_cast<FFLin<T>*>(pOp)->_ownCoef = true; // transfer ownership
      // nothing to do if shallow copy requested (policy=0)
#ifdef MC__FFLIN_TRACE
      std::cerr << "FFLin operation address: " << this << std::endl;
      std::cerr << "FFLin address in DAG: "    << _ptrCoef << std::endl;
#endif
      return *pRes;
    }

  // set the object and related operation in DAG
  FFVar& _set
    ( size_t const nVar, FFVar const*const* ppVar, double const* pCoef, double const& Bias, int policy )
    {
      if( _ownCoef && _ptrCoef ) delete _ptrCoef;
      _ownCoef = ( policy>0? true: false );
      this->owndata = false;
      this->data = nullptr;
      _ptrCoef = const_cast<double*>(pCoef);
      _nCoef   = nVar;
      _valBias = Bias;

      FFVar* pRes = this->insert_external_operation( *this, 1, nVar, ppVar )[0];
      
      _ownCoef = false;
      FFOp* pOp = pRes->opdef().first;
      if( policy > 0 )
        _ptrCoef = static_cast<FFLin<T>*>(pOp)->_ptrCoef; // set pointer to DAG copy
      else if( policy < 0 )
        static_cast<FFLin<T>*>(pOp)->_ownCoef = true; // transfer ownership
      // nothing to do if shallow copy requested (policy=0)
#ifdef MC__FFLIN_TRACE
      std::cerr << "FFLin operation address: " << this << std::endl;
      std::cerr << "FFLin address in DAG: "    << _ptrCoef << std::endl;
#endif
      return *pRes;
    }

  // set the object and related operation in DAG
  FFVar& _set
    ( size_t const nVar, FFVar const* pVar, double const& Coef, double const& Bias )//, int policy )
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::set vector\n";
      std::cout << "   Bias = " << Bias << std::endl;
      std::cout << "   Coef = " << Coef << std::endl;
#endif
      if( _ownCoef && _ptrCoef ) delete _ptrCoef;
      _ownCoef = true;//( policy>0? true: false );
      this->owndata = false;
      this->data = nullptr;
      _ptrCoef = const_cast<double*>(&Coef);
      _nCoef   = 1;
      _valBias = Bias;

      FFVar* pRes = this->insert_external_operation( *this, 1, nVar, pVar )[0];
      
      _ownCoef = false;
      FFOp* pOp = pRes->opdef().first;
      //if( policy > 0 )
        _ptrCoef = static_cast<FFLin<T>*>(pOp)->_ptrCoef; // set pointer to DAG copy
      //else if( policy < 0 )
      //  static_cast<FFLin<T>*>(pOp)->_ownCoef = true; // transfer ownership
      // nothing to do if shallow copy requested (policy=0)
#ifdef MC__FFLIN_TRACE
      std::cerr << "FFLin operation address: " << this << std::endl;
      std::cerr << "FFLin address in DAG: "    << _ptrCoef << std::endl;
#endif
      return *pRes;
    }

  // set the object and related operation in DAG
  FFVar& _set
    ( size_t const nVar, FFVar const*const* ppVar, double const& Coef, double const& Bias )//, int policy )
    {
      if( _ownCoef && _ptrCoef ) delete _ptrCoef;
      _ownCoef = true;//( policy>0? true: false );
      this->owndata = false;
      this->data = nullptr;
      _ptrCoef = const_cast<double*>(&Coef);
      _nCoef   = 1;
      _valBias = Bias;

      FFVar* pRes = this->insert_external_operation( *this, 1, nVar, ppVar )[0];
      
      _ownCoef = false;
      FFOp* pOp = pRes->opdef().first;
      //if( policy > 0 )
        _ptrCoef = static_cast<FFLin<T>*>(pOp)->_ptrCoef; // set pointer to DAG copy
      //else if( policy < 0 )
      //  static_cast<FFLin<T>*>(pOp)->_ownCoef = true; // transfer ownership
      // nothing to do if shallow copy requested (policy=0)
#ifdef MC__FFLIN_TRACE
      std::cerr << "FFLin operation address: " << this << std::endl;
      std::cerr << "FFLin address in DAG: "    << _ptrCoef << std::endl;
#endif
      return *pRes;
    }

  // generic evaluation
  template< typename U >
  void _eval
    ( size_t const nRes, U* vRes, size_t const nVar, U const* vVar, unsigned const* mVar )
    const;


public:

  //! @brief Enumeration type for copy policy of DAG object
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of DAG object in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of DAG object in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of DAG object in FFGraph (with ownership transfer)
  };

  // Default constructor
  FFLin
    ()
    : FFOp     ( EXTERN ),
      _valBias ( 0. ),
      _nCoef   ( 0 ),
      _ptrCoef ( nullptr ),
      _ownCoef ( false )
    {}
    
  // Copy constructor
  FFLin
    ( FFLin const& other )
    : FFOp     ( other ),
      _valBias ( other._valBias )
    {
#ifdef MC__FFLIN_TRACE
      std::cout << "FFLin::copy constructor\n";
#endif
      if( !other._nCoef || !other._ptrCoef )
        throw std::runtime_error( "FFLin::copy constructor ** Null pointer to coefficient array\n" );

      _ownCoef = other._ownCoef;      
      if( _ownCoef ){ // deep copy
        //if( _ptrCoef ) delete[] _ptrCoef;
        _nCoef   = other._nCoef;
        _ptrCoef = new double[ _nCoef ];
        for( size_t i=0; i<_nCoef; ++i ) _ptrCoef[i] = other._ptrCoef[i];
#ifdef MC__FFLIN_TRACE
        std::cout << "FFLin copy constructor\n";
        std::cout << "   Other: " << other._ptrCoef << std::endl;
        std::cout << "   Bias = " << _valBias << std::endl;
        std::cout << "   Coef (" << _ptrCoef << ") = [ ";
        for( unsigned i=0; i<_nCoef; ++i ) std::cout << _ptrCoef[i] << " ";
        std::cout << "]\n";
#endif
      }
      else{ // shallow copy
        _nCoef   = other._nCoef;
        _ptrCoef = other._ptrCoef;
      }
    }

  // Destructor
  virtual ~FFLin
    ()
    {
      if( _ownCoef && _ptrCoef ) delete[] _ptrCoef;
    }

  // Define operation
  FFVar& operator()
    ( std::vector<FFVar> const& Var, std::vector<double> const& Coef,
      double const& Bias=0. )//, int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( Var.size() == Coef.size() );
#endif
      return _set( Var.size(), Var.data(), Coef.data(), Bias, COPY );//policy );
    }

  FFVar& operator()
    ( size_t const nVar, FFVar const* pVar, double const* pCoef,
      double const& Bias=0., int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( nVar && pVar && pCoef );
#endif
      return _set( nVar, pVar, pCoef, Bias, policy );
    }

  FFVar& operator()
    ( size_t const nVar, FFVar const*const* ppVar, double const* pCoef,
      double const& Bias=0., int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( nVar && ppVar && pCoef );
#endif
      return _set( nVar, ppVar, pCoef, Bias, policy );
    }

  FFVar& operator()
    ( std::vector<FFVar> const& Var, double const& Coef,
      double const& Bias=0. )//, int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( Var.size() );
#endif
      return _set( Var.size(), Var.data(), Coef, Bias );//, policy );
    }

  FFVar& operator()
    ( size_t const nVar, FFVar const* pVar, double const& Coef,
      double const& Bias=0. )//, int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( nVar && pVar );
#endif
      return _set( nVar, pVar, Coef, Bias );//, policy );
    }

  FFVar& operator()
    ( size_t const nVar, FFVar const*const* ppVar, double const& Coef,
      double const& Bias=0. )//, int policy=COPY )
    {
#ifdef MC__FFLIN_CHECK
      assert( nVar && ppVar );
#endif
      return _set( nVar, ppVar, Coef, Bias );//, policy );
    }

  // Retrieve bias value
  double Bias
    ()
    const
    {
      return _valBias;
    }
  // Retrieve number of coefficients
  size_t nCoef
    ()
    const
    {
      return _nCoef;
    }
  // Retrieve pointer to coefficients
  double const* Coef
    ()
    const
    {
      return _ptrCoef;
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
        return _eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return _eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFInv ) )
        return _eval( nRes, static_cast<FFInv*>(vRes), nVar, static_cast<FFInv const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return _eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return _eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( T ) )
        return _eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );
      else if( idU == typeid( McCormick<T> ) )
        return _eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWCU> ) )
        return _eval( nRes, static_cast<SupVar<PWCU>*>(vRes), nVar, static_cast<SupVar<PWCU> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWLU> ) )
        return _eval( nRes, static_cast<SupVar<PWLU>*>(vRes), nVar, static_cast<SupVar<PWLU> const*>(vVar), mVar );
      else if( idU == typeid( McCormick<SupVar<PWCU>> ) )
        return _eval( nRes, static_cast<McCormick<SupVar<PWCU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWCU>> const*>(vVar), mVar );
      else if( idU == typeid( McCormick<SupVar<PWLU>> ) )
        return _eval( nRes, static_cast<McCormick<SupVar<PWLU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWLU>> const*>(vVar), mVar );
      else if( idU == typeid( PolVar<T> ) )
        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return _eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFLin::feval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

//  void eval
//    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
//    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
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

  void eval
    ( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar, unsigned const* mVar )
    const;

  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      if( idU == typeid( T ) )
        return reval( nRes, static_cast<T const*>(vRes), nVar, static_cast<T*>(vVar) );
      else if( idU == typeid( PolVar<T> ) )
        return reval( nRes, static_cast<PolVar<T> const*>(vRes), nVar, static_cast<PolVar<T>*>(vVar) );

      throw std::runtime_error( "FFLin::reval ** No evaluation method for type"+std::string(idU.name())+"\n" );
    }

  bool reval
    ( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
    const;

  bool reval
    ( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
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
    { std::ostringstream oss; oss << this->_ptrCoef; return "Lin[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

template< typename T >
template< typename U >
inline void
FFLin<T>::_eval
( size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
  unsigned const* mVar )
const
{
  vRes[0] = _valBias;
  if( _nCoef < nVar && *_ptrCoef == 1. )
    for( size_t i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  else if( _nCoef < nVar && *_ptrCoef == -1. )
    for( size_t i=0; i<nVar; ++i ) vRes[0] -= vVar[i];
  else if( _nCoef < nVar )
    for( size_t i=0; i<nVar; ++i ){
      if( *_ptrCoef == 0. ) continue;
      vRes[0] += *_ptrCoef * vVar[i];
    }
  else
    for( size_t i=0; i<nVar; ++i ){
      if( _ptrCoef[i] == 0. ) continue;
      vRes[0] += _ptrCoef[i] * vVar[i];
    }
}
/*
template< typename T >
inline void
FFLin<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: FFDep\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif
  vRes[0] = 0;
  for( size_t i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::L );
}
*/
template< typename T >
inline void
FFLin<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: FFVar\n";
  std::cerr << "DAG operation address: " << this     << std::endl;
  std::cerr << "DAG address in DAG: "    << _ptrCoef << std::endl;
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  vRes[0] = *(this->insert_external_operation( *this, 1, nVar, vVar )[0]);
}
/*
template< typename T >
inline void
FFLin<T>::eval
( size_t const nRes, fadbad::F<double>* vRes, size_t const nVar, fadbad::F<double> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: fadbad::F<double>\n";
#endif
#ifdef MC__FFLIN_CHECK
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
FFLin<T>::eval
( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFLIN_CHECK
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
FFLin<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::deriv: FFVar\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  for( unsigned i=0; i<nVar; ++i )
    vDer[0][i] = (_nCoef<nVar? _ptrCoef[0]: _ptrCoef[i]);
}

template< typename T >
inline void
FFLin<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: SLiftVar\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  if( vVar->env()->options.KEEPFACT )// || vVar->env()->options.LIFTIPOW )
    return vVar->env()->lift( nRes, vRes, nVar, vVar );

  _eval( 1, vRes, nVar, vVar, nullptr ); // Expand operation
}

template< typename T >
inline void
FFLin<T>::eval
( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::eval: PolVar\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__FFLIN_CHECK
  assert( img && dag );
#endif
  FFVar* pRes = dag->curOp()->varout[0];
#ifdef MC__FFLIN_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  std::vector<T> vTVar( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vTVar[i] = vVar[i].range();
  T TRes;
  _eval( 1, &TRes, nVar, vTVar.data(), nullptr ); // Expand operation
  vRes[0].set( img, *pRes, TRes );
}

template< typename T >
inline bool
FFLin<T>::reval
( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::reval: PolVar\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__FFLIN_CHECK
  assert( img && pop );
#endif

  //if( _nCoef < nVar )
  //  img->add_cut( pop, PolCut<T>::EQ, -_valBias, nVar, vVar, _ptrCoef[0], vRes[0], -1. );
  //else
  //  img->add_cut( pop, PolCut<T>::EQ, -_valBias, nVar, vVar, _ptrCoef, vRes[0], -1. );

  // Exclude zero coefficients in cut
  auto cut = *img->add_cut( pop, PolCut<T>::EQ, -_valBias, vRes[0], -1. );
  if( _nCoef < nVar )
    for( size_t i=0; i<nVar; ++i ){
      if( *_ptrCoef == 0. ) continue;
      cut->append( vVar[i], _ptrCoef[0] );
    }
  else
    for( size_t i=0; i<nVar; ++i ){
      if( _ptrCoef[i] == 0. ) continue;
      cut->append( vVar[i], _ptrCoef[i] );
    }
  
  return true;
}

template< typename T >
inline bool
FFLin<T>::reval
( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
const
{
#ifdef MC__FFLIN_TRACE
  std::cout << "FFLin::reval: T\n";
#endif
#ifdef MC__FFLIN_CHECK
  assert( _nCoef && _ptrCoef && nRes == 1 );
#endif

  for( unsigned i=0; i<nVar; ++i ){
    if( isequal( _nCoef<nVar? _ptrCoef[0]: _ptrCoef[i], 0. ) ) continue;
    T Var_i( vRes[0] - _valBias ); 
    for( unsigned j=0; j<nVar; ++j ){
      if( i == j ) continue;
      Var_i -= (_nCoef<nVar? _ptrCoef[0]: _ptrCoef[j]) * vVar[j];
    }
    Var_i /= ( _nCoef<nVar? _ptrCoef[0]: _ptrCoef[i] );
    if( !Op<T>::inter( vVar[i], Var_i, vVar[i] ) ) return false;
  }
  return true;
}

template< typename T >
inline bool
FFLin<T>::lt
( FFOp const* other )
const
{
  auto const* oprecast = dynamic_cast<FFLin<T> const*>(other);
#ifdef MC__FFLin_TRACE
  std::cout << "FFLin<T>::lt\n";
  std::cout << _ptrCoef << " (" << data << ") <=> " << oprecast->_ptrCoef << " (" << oprecast->data << ")" << std::endl;
#endif
//  if( _ptrCoef < oprecast->_ptrCoef ) return true;
//  if( _ptrCoef > oprecast->_ptrCoef ) return false;
  // Compare linear coefficients one by one, then bias
  for( unsigned i=0; i<_nCoef || i<oprecast->_nCoef; ++i ){
    if( _ptrCoef[i<_nCoef?i:0] < oprecast->_ptrCoef[i<oprecast->_nCoef?i:0] ) return true;
    if( _ptrCoef[i<_nCoef?i:0] > oprecast->_ptrCoef[i<oprecast->_nCoef?i:0] ) return false;
  }
  return( _valBias < oprecast->_valBias );
}

} // end namescape mc

#endif
