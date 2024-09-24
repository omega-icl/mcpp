#ifndef MC__FFVECT_HPP
#define MC__FFVECT_HPP

#undef  MC__FFVECT_DEBUG
#define MC__FFVECT_CHECK

#include "ffunc.hpp"
#include "ffdep.hpp"
#include "slift.hpp"

namespace mc
{

//! @brief C++ class for vectorizing DAG evaluation
////////////////////////////////////////////////////////////////////////
//! mc::Vect is a C++ class for vectorizing DAG evaluation in MC++.
////////////////////////////////////////////////////////////////////////
class Vect
////////////////////////////////////////////////////////////////////////
{
private:

  //! @brief DAG duplication for thread evaluation
  struct Worker
  {
    //! @brief Default constructor
    Worker
      ()
      : pDAG( nullptr )
      {}
      
    //! @brief Destructor
    ~Worker
      ()
      {
        if( pDAG ) delete pDAG;
      }

    //! @brief local copy of DAG
    FFGraph*                pDAG;

    //! @brief DAG subgraph
    FFSubgraph              sgDep;

    //! @brief vector of dependent variables
    std::vector<FFVar>      vDep;

    //! @brief vector of independent variables
    std::vector<FFVar>      vVar;

    //! @brief vector of constants
    std::vector<FFVar>      vCst;
 
    //! @brief copy worker
    bool copy
      ( Worker const& wkr )
      {
        return set( wkr.pDAG, wkr.vDep, wkr.vVar, wkr.vCst );
      }
 
    //! @brief set DAG copy
    bool set
      ( FFGraph* _pDAG, std::vector<FFVar> const& _vDep, std::vector<FFVar> const& _vVar,
        std::vector<FFVar> const& _vCst )
      {
        try{
          if( pDAG ) delete pDAG;
          pDAG = new FFGraph;
          pDAG->options = _pDAG->options;
          pDAG->insert( _pDAG, _vVar, vVar );
          if( _vCst.empty() )
            vCst.clear();
          else
            pDAG->insert( _pDAG, _vCst, vCst );
          pDAG->insert( _pDAG, _vDep, vDep );
          sgDep = pDAG->subgraph( vDep );
        }
  
        catch( mc::FFBase::Exceptions &eObj ){
          std::cerr << "Error " << eObj.ierr()
                    << " in DAG copy for parallel evaluation:" << std::endl
                    << eObj.what() << std::endl;
          return false;
        }

        return true;
      }

    //! @brief evaluate DAG
    template <typename U>
    void eval
      ( std::vector<U>& wkU, U* uDep, U const* uVar, U const* uCst )
      {
#ifdef MC__FFVECT_DEBUG
        pDAG->output( sgDep );
#endif
        if( uCst )
          pDAG->eval( sgDep, wkU, vDep.size(), vDep.data(), uDep, vVar.size(), vVar.data(), uVar,
                      vCst.size(), vCst.data(), uCst );
        else
          pDAG->eval( sgDep, wkU, vDep.size(), vDep.data(), uDep, vVar.size(), vVar.data(), uVar );
      }
#ifdef MC__FFVECT_DEBUG
    void eval
      ( std::vector<double>& wkU, double* uDep, double const* uVar, double const* uCst )
      {
        pDAG->output( sgDep );
        for( size_t i=0; i<vVar.size(); ++i ) std::cout << "Var: " << vVar[i] << " = " << uVar[i] << std::endl;
        for( size_t i=0; i<vCst.size(); ++i ) std::cout << "Cst: " << vCst[i] << " = " << uCst[i] << std::endl;
        if( uCst )
          pDAG->eval( sgDep, wkU, vDep.size(), vDep.data(), uDep, vVar.size(), vVar.data(), uVar,
                      vCst.size(), vCst.data(), uCst );
        else
          pDAG->eval( sgDep, wkU, vDep.size(), vDep.data(), uDep, vVar.size(), vVar.data(), uVar );
        for( size_t i=0; i<vDep.size(); ++i ) std::cout << "Dep:" << vDep[i] << " = " << uDep[i] << std::endl;
      }
#endif
  };

  //! @brief Underlying DAG
  FFGraph*                                      _pDAG;
  //! @brief Participating variables
  std::vector<FFVar>                            _vVar;
  //! @brief Number of variables
  size_t                                        _nVar;
  //! @brief Participating variables
  std::vector<FFVar>                            _vCst;
  //! @brief Number of variables
  size_t                                        _nCst;
  //! @brief Function structure
  std::vector<std::vector<FFVar>>               _vFun;
  //! @brief Number of functions
  size_t                                        _nFun;
  //! @brief Workers
  std::vector<Worker>                           _vWkr;
  //! @brief Threads
#ifdef MC__USE_THREAD
  std::vector<std::thread>                      _vThd;
#endif

public:

  //! @brief Default constructor
  Vect
    ()
    : _pDAG( nullptr ),
      _nVar( 0 ), _nCst( 0 ), _nFun( 0 )
    {}

  //! @brief Default constructor
  Vect
    ( FFGraph* pDAG, std::vector<FFVar> const& vVar,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      set( pDAG, vVar, std::vector<FFVar>(), vFun );
    }

  Vect
    ( FFGraph* pDAG, std::vector<FFVar> const& vVar, std::vector<FFVar> const& vCst,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      set( pDAG, vVar, vCst, vFun );
    }

  //! @brief Copy constructor
  Vect
    ( Vect const& v )
    : _pDAG( v._pDAG ),
      _vVar( v._vVar ), _nVar( v._nVar ),
      _vCst( v._vCst ), _nCst( v._nCst ),
      _vFun( v._vFun ), _nFun( v._nFun ),
      options( v.options )
    {
      _vWkr.resize( v._vWkr.size() );
      for( size_t iwkr=0; iwkr<v._vWkr.size(); ++iwkr ){
#ifdef MC__FFVECT_DEBUG
        std::cout << "Copying worker #" << iwkr << std::endl;
#endif
        _vWkr[iwkr].copy( v._vWkr[iwkr] );

      }
    }

  ~Vect
    () 
    {}

  //! @brief Set function structure
  bool set
    ( FFGraph* pDAG, size_t nVar, FFVar const* pVar,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      return set( pDAG, std::vector<FFVar>( pVar, pVar+nVar ), std::vector<FFVar>(), vFun );
    }

  bool set
    ( FFGraph* pDAG, size_t nVar, FFVar const* pVar, size_t nCst, FFVar const* pCst,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      return set( pDAG, std::vector<FFVar>( pVar, pVar+nVar ), std::vector<FFVar>( pCst, pCst+nCst ), vFun );
    }

  bool set
    ( FFGraph* pDAG, std::vector<FFVar> const& vVar,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      return set( pDAG, vVar, std::vector<FFVar>(), vFun );
    }
    
  bool set
    ( FFGraph* pDAG, std::vector<FFVar> const& vVar, std::vector<FFVar> const& vCst,
      std::vector<std::vector<FFVar>> const& vFun=std::vector<std::vector<FFVar>>() )
    {
      _pDAG = pDAG;

      _vVar = vVar;
      _nVar = vVar.size();

      if( vCst.empty() ){
        _vCst.clear();
        _nCst = 0;
      }
      else{
        _vCst = vCst;
        _nCst = vCst.size();
      }
      
      _nFun = 0;
      _vFun = vFun;
      _vWkr.clear();
      _vWkr.reserve( vFun.size() );
      for( size_t i=0; i<vFun.size(); ++i ){
        _nFun += vFun[i].size();
        _vWkr.push_back( Worker() );
        if( !_vWkr.back().set( pDAG, vFun[i], vVar, vCst ) )
          return false;
      }
      return true;
    }

//  //! @brief Add function to vector structure
//  bool add
//    ( std::vector<FFVar> const& vFun )
//    {
//      _nFun += vFun.size();
//      _vFun.push_back( vFun );
//      _vWkr.push_back( Worker() );
//      return _vWkr.back().set( _pDAG, vFun, _vVar, _vCst );
//    }

  //! @brief Original DAG
  FFGraph* pDAG
    ()
    const
    { return _pDAG; }

  //! @brief Number of constants
  size_t nCst
    ()
    const
    { return _nCst; }

  //! @brief Participating constats
  std::vector<FFVar> const& vCst
    ()
    const
    { return _vCst; }

  //! @brief Number of variables
  size_t nVar
    ()
    const
    { return _nVar; }

  //! @brief Participating variables
  std::vector<FFVar> const& vVar
    ()
    const
    { return _vVar; }

  //! @brief Number of functions
  size_t nFun
    ()
    const
    { return _nFun; }

  //! @brief Vectorized function
  std::vector<std::vector<FFVar>> const& vFun
    ()
    const
    { return _vFun; }

  //! @brief MLP options
  struct Options
  {
    //! @brief Default constructor
    Options():
      AUTODIFF(F)
      {}

    //! @brief Assignment operator
    Options& operator=
      ( Options const& opt ){
        AUTODIFF  = opt.AUTODIFF;
        return *this;
      }

    //! @brief Enumeration type for AD strategy
    enum AD{
      F=0,	//!< Forward differentiation
      B		//!< Backward differentiation
    };

    //! @brief Whether to apply forward or reverse automatic differentiation
    int       AUTODIFF;
  } options;

  //! @brief Evaluate function
  template <typename U>
  void eval_serial
    ( U* y, U const* x, U const* c );

  //! @brief Evaluate function
  template <typename U>
  void eval_parallel
    ( U* y, U const* x, U const* c );
};

template< typename U >
inline void
Vect::eval_serial
( U* y, U const* x, U const* c )
{
  static thread_local std::vector<U> wkU;
  U* yi = y;
  for( size_t i=0; i<_vFun.size(); yi+=_vFun[i].size(), ++i ) // increment i last
    _vWkr[i].eval( wkU, yi, x, c ); 
}


template< typename U >
inline void
Vect::eval_parallel
( U* y, U const* x, U const* c )
{
  static thread_local std::vector<std::vector<U>> vwkU( _vFun.size() );

  // Run evaluations on current and auxiliary threads
  U* yi = y;
#ifdef MC__USE_THREAD
  if( _vThd.size() < _vWkr.size()-1 ) _vThd.resize( _vWkr.size()-1 ); 
  for( size_t i=0; i<_vFun.size(); yi+=_vFun[i].size(), ++i ){ // increment i last
    // Final element to run on current thread
    if( i == _vFun.size()-1 ){
#ifdef MC__FFVECT_DEBUG
      std::cout << "Running on current thread" << std::endl;
#endif
      _vWkr[i].eval( vwkU[i], yi, x, c );
    }
    // Other elements to run on auxiliary threads
    else{
#ifdef MC__FFVECT_DEBUG
      std::cout << "Starting thread #" << i << std::endl;
#endif
      //_vWkr[i].eval( vwkU[i], yi, x, c );
      _vThd[i] = std::thread( &Vect::Worker::eval<U>, &_vWkr[i], std::ref(vwkU[i]), yi, x, c );
    }
  }
#else
  for( size_t i=0; i<_vFun.size(); yi+=_vFun[i].size(), ++i ){ // increment i last
    _vWkr[i].eval( vwkU[i], yi, x, c );
  }
#endif

  // Join all the threads to the main one
  for( size_t i=0; i<_vFun.size()-1; ++i )
    _vThd[i].join();
}

//! @brief C++ class defining vectorized functions as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFVect is a C++ class for defining vectorized functions as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFVect
////////////////////////////////////////////////////////////////////////
: public FFOp
{

private:
  // pointer to vectorized function
  Vect*   _pFun;
  // Whether this class owns _pFun
  bool    _ownFun;

  FFVar** _set
    ( Vect* pFun )
    {
      if( _ownFun && _pFun ) delete _pFun;
      _ownFun = true;
      _pFun = pFun;
      owndata = false;
      data = _pFun;

      FFVar** ppRes = ( _pFun->nCst()?
                        insert_external_operation( *this, _pFun->nFun(), _pFun->nVar(), _pFun->vVar().data(),
                                                   _pFun->nCst(), _pFun->vCst().data() ):
                        insert_external_operation( *this, _pFun->nFun(), _pFun->nVar(), _pFun->vVar().data() ) );

      _ownFun = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      
      _pFun = dynamic_cast<FFVect<T>*>(pOp)->_pFun;
#ifdef MC__FFVECT_TRACE
      std::cerr << "Vect address in DAG: " << _pFun << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Default Constructor
  FFVect
    ()
    : FFOp( EXTERN ),
      _pFun( nullptr ),
      _ownFun( false )
    {}

  // Destructor
  virtual ~FFVect
    ()
    {
      if( _ownFun ) delete _pFun;
    }

  // Copy constructor
  FFVect
    ( FFVect<T> const& Op )
    : FFOp( Op )
    {
#ifdef MC__FFVECT_TRACE
      std::cout << "FFVect::copy constructor\n";
#endif
      if( !Op._pFun )
        throw std::runtime_error( "FFVect::copy constructor ** Undefined function\n" );

      _ownFun = Op._ownFun;
      if( _ownFun )
        _pFun = new Vect( *Op._pFun );
      else
        _pFun = Op._pFun;
    }

  // Define operation
  FFVar& operator()
    ( size_t const iFun, Vect* pFun )
    {
#ifdef MC__FFVECT_CHECK
      assert( iFun < pFun->nFun() );
#endif
      return *(_set( pFun )[iFun]);
    }

  FFVar** operator()
    ( Vect* pFun )
    {
      return _set( pFun );
    }

  // Forward evaluation overloads
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
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
//      else if( idU == typeid( T ) )
//        return eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );
//      else if( idU == typeid( McCormick<T> ) )
//        return eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );
//      else if( idU == typeid( ISVar<T> ) )
//        return eval( nRes, static_cast<ISVar<T>*>(vRes), nVar, static_cast<ISVar<T> const*>(vVar), mVar );
//      else if( idU == typeid( ASVar<T> ) )
//        return eval( nRes, static_cast<ASVar<T>*>(vRes), nVar, static_cast<ASVar<T> const*>(vVar), mVar );
//      else if( idU == typeid( PolVar<T> ) )
//        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );

      throw std::runtime_error( "FFVect::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  template< typename U >
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const;

#ifdef MC__FFVECT_DEBUG
  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;
#endif

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

//  void eval
//    ( unsigned const nRes, PolVar<T>* vRes, unsigned const nVar, PolVar<T> const* vVar, unsigned const* mVar )
//    const;


  // Backward evaluation overloads
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
//      if( idU == typeid( PolVar<T> ) )
//        return reval( nRes, static_cast<PolVar<T> const*>(vRes), nVar, static_cast<PolVar<T>*>(vVar) );

      throw std::runtime_error( "FFVect::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

//  bool reval
//    ( unsigned const nRes, PolVar<T> const* vRes, unsigned const nVar, PolVar<T>* vVar )
//    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;
    
  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << data; return "Vect[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

template< typename T >
template< typename U >
inline void
FFVect<T>::eval
( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() );
#endif

  //_pFun->eval_serial( vRes, vVar );
  _pFun->eval_parallel( vRes, vVar, vVar+_pFun->nVar() );
}

#ifdef MC__FFVECT_DEBUG
template< typename T >
inline void
FFVect<T>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::eval: double\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() );
#endif

  //_pFun->eval_serial( vRes, vVar );
  _pFun->eval_parallel( vRes, vVar, vVar+_pFun->nVar() );
  for( size_t i=0; i<nRes; ++i ) std::cout << "vRes[" << i << "] = " << vRes[i] << std::endl;
}
#endif

template< typename T >
inline void
FFVect<T>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::eval: FFDep\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() );
#endif

  _pFun->eval_serial( vRes, vVar, vVar+_pFun->nVar() );

//  vRes[0] = 0;
//  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
//  vRes[0].update( FFDep::TYPE::N );
//  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFVect<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef NC__FFVECT_TRACE
  std::cout << "FFVect::eval: FFVar\n";
  std::cerr << "FFVect operation address: " << this << std::endl;
  std::cerr << "Vect address in DAG: " << _pFun << std::endl;
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() );
  #endif

  FFVar** ppRes = insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
}
/*
template< typename T >
inline void
FFVect<T>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar() && nRes == _pFun->nFun() && _pFun->pDAG() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = insert_external_operation( *this, nRes, nVar, vVarVal.data() );

  std::vector<std::vector<FFVar>> vFun;
  vFun.reserve( _pFun->vFun().size() );
//  Vect FunDer( _pFun->pDAG(), vVarVal );
  for( auto const& vFi : _pFun->vFun() )
    switch( _pFun->options.AUTODIFF ){
      default:
      case Vect::Options::AD::F:
        vFun.push_back( _pFun->pDAG()->FAD( vFi, vVarVal ) );
//        FunDer.add( _pFun->pDAG()->FAD( vFi, vVarVal ) );
        break;
      case Vect::Options::AD::B:
        vFun.push_back( _pFun->pDAG()->BAD( vFi, vVarVal ) );
//        FunDer.add( _pFun->pDAG()->BAD( vFi, vVarVal ) );
        break;
    }

  Vect FunDer( _pFun->pDAG(), vVarVal, vFun );
  FFVect<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( &FunDer ); 

  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = *vResVal[k];
    for( unsigned i=0; i<nVar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j ){
      vRes[k][j] = 0.;
      for( unsigned i=0; i<nVar; ++i ){
        if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
        vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
      }
    }
  }
}
*/
template< typename T >
inline void
FFVect<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* pVar, FFVar** vDer )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::deriv: FFVar\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() && _pFun->pDAG() );
#endif

  auto&& vVar = std::vector<FFVar>( pVar, pVar+_pFun->nVar() );
  auto&& vCst = ( _pFun->nCst()? std::vector<FFVar>( pVar+_pFun->nVar(), pVar+nVar ): std::vector<FFVar>() );
  std::vector<std::vector<FFVar>> vFun;
  vFun.reserve( _pFun->vFun().size() );
//  Vect FunDer( _pFun->pDAG(), vVar, vCst );
  for( auto const& vFi : _pFun->vFun() )
    switch( _pFun->options.AUTODIFF ){
      default:
      case Vect::Options::AD::F:
        vFun.push_back( _pFun->pDAG()->FAD( vFi, vVar ) );
        //for( auto const& vFuni : vFun.back() )
          //std::cout << vFuni << std::endl;
//        FunDer.add( _pFun->pDAG()->FAD( vFi, vVar ) );
        break;
      case Vect::Options::AD::B:
        vFun.push_back( _pFun->pDAG()->BAD( vFi, vVar ) );
//        FunDer.add( _pFun->pDAG()->BAD( vFi, vVar ) );
        break;
    }

  Vect FunDer( _pFun->pDAG(), vVar, vCst, vFun );
  FFVect<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( &FunDer ); 

  for( unsigned b=0, bik=0, boff=0; b<_pFun->vFun().size(); boff+=_pFun->vFun()[b].size(), ++b )
    for( unsigned k=0; k<_pFun->vFun()[b].size(); ++k ){
      for( unsigned i=0; i<_pFun->nVar(); ++i, ++bik ){
        //std::cout << boff+k << "," << i << " ==> " << bik << std::endl;
        vDer[boff+k][i] = *vResDer[bik];
      }
      for( unsigned i=_pFun->nVar(); i<nVar; ++i ){
        //std::cout << boff+k << "," << i << " ==> --" << std::endl;
        vDer[boff+k][i] = 0;
      }
    }
}

//F1[1]
//...
//F1[N1]
//...
//FM[1]
//...
//FM[NM]


//dF1[1]dP[1]
//...
//dF1[1]dP[NP]
//...
//dF1[N1]dP[1]
//...
//dF1[N1]dP[NP]
//...
//...
//dFM[1]dP[1]
//...
//dFM[1]dP[NP]
//...
//dFM[NM]dP[1]
//...
//dFM[NM]dP[NP]


template< typename T >
inline void
FFVect<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFVECT_TRACE
  std::cout << "FFVect::eval: SLiftVar\n";
#endif
#ifdef MC__FFVECT_CHECK
  assert( _pFun && nVar == _pFun->nVar()+_pFun->nCst() && nRes == _pFun->nFun() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

} // end namespace mc

#endif
