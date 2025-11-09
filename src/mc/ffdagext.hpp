#ifndef MC__FFDAGEXT_HPP
#define MC__FFDAGEXT_HPP

//#define MC__FFDAGEXT_DEBUG
#define MC__FFDAGEXT_CHECK

#include "ffextern.hpp"

namespace mc
{

//! @brief C++ class for evaluation and relaxation of an operation defined by an expression tree
////////////////////////////////////////////////////////////////////////
//! mc::DAGEXT is a C++ class for evaluation and relaxation of an
//! operation defined by an expression tree using MC++
////////////////////////////////////////////////////////////////////////
template <typename T> 
class DAGEXT
////////////////////////////////////////////////////////////////////////
{
private:

  //! @brief Expression tree
  FFGraph*                               _dag;
  //! @brief Variables
  std::vector<FFVar>                     _varin;
  //! @brief Dependents
  std::vector<FFVar>                     _varout;
  //! @brief Codelist
  FFSubgraph                             _codelist;

  //! @brief Intermediate storage for DAG evaluation
  std::vector<double>                   _wkD;
  std::vector<fadbad::F<double>>        _wkFD;
  std::vector<T>                        _wkI;
  std::vector<T>                        _wkCPI;
  std::vector<McCormick<T>>             _wkMC;
  std::vector<SupVar<PWCU>>             _wkPWCS;
  std::vector<McCormick<SupVar<PWCU>>>  _wkMCPWCS;
  std::vector<SupVar<PWLU>>             _wkPWLS;
  std::vector<McCormick<SupVar<PWLU>>>  _wkMCPWLS;
  std::vector<PolVar<T>>                _wkPOL;
  std::vector<SCVar<T>>                 _wkSC;

  //! @brief Set expression tree
  void _set
    ( FFGraph* dag, std::vector<FFVar> const& varin, std::vector<FFVar> const& varout )
    {
      delete _dag; _dag = new FFGraph;
#ifdef MC__DAGEXT_DEBUG
      std::cout << "DAGEXT:: Original DAG: " <<  dag << std::endl;
      std::cout << "DAGEXT:: Copied DAG:   " << _dag << std::endl;
#endif
      _dag->insert( dag, varin,  _varin  );
      _dag->insert( dag, varout, _varout );
#ifdef MC__DAGEXT_DEBUG
      _codelist = _dag->subgraph( _varout );
      std::vector<FFExpr> strout = FFExpr::subgraph( _dag, _codelist );
      for( size_t i=0; i<strout.size(); ++i )
        std::cout << "F" << i << ": " << strout[i] << std::endl;
#else
      _codelist.clear();
#endif
    }

  //! @brief Evaluate expression tree
  template <typename U>
  void _eval
    ( U const* valin, U* valout, std::vector<U>& wk )
//    const
    {
      _dag->eval( _codelist, wk, _varout.size(), _varout.data(), valout, _varin.size(), _varin.data(), valin );
    }

  //! @brief Evaluate expression tree
  template <typename U>
  bool _reval
    ( U* valin, U const* valout, std::vector<U>& wk, U const& inf )
//    const
    {
      int flag = _dag->reval( _codelist, wk, _varout.size(), _varout.data(), const_cast<U*>( valout ),
                              _varin.size(), _varin.data(), valin, inf, options.CPMAX, options.CPTHRES );
#ifdef MC__DAGEXT_DEBUG
      std::cout << "DAGEXT:: Work array: " << flag << " passes\n";
      for( unsigned i=0; i<_codelist.len_tap-_codelist.len_wrk; i++ )
        std::cout << "wk[" << i << "] = " << wk[i] << std::endl;
      for( unsigned i=0; i<_varin.size(); i++ )
        std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif
      return( flag<0? false: true );
    }

public:

  //! @brief DAGEXT options
  struct Options
  {
    //! @brief Default constructor
    Options()
      {
        reset();
      }

    //! @brief Assignment operator
    Options& operator=
      ( Options const& other ){
        AUTODIFF = other.AUTODIFF;
        CPMAX    = other.CPMAX;
        CPTHRES  = other.CPTHRES;
        CPINF    = other.CPINF;

        return *this;
      }

    //! @brief Reset options
    void reset
      ()
      {
        AUTODIFF = F;
        CPMAX    = 1; // Leave it to the outer DAG to iterate as necessary by default
        CPTHRES  = 1e4*DBL_EPSILON;
        CPINF    = 1e30;
      }

    //! @brief Enumeration type for AD strategy
    enum AD{
      F=0,	//!< Forward differentiation
      B		//!< Backward differentiation
    };

    //! @brief Whether to apply forward or reverse automatic differentiation
    int      AUTODIFF;
    //! @brief Maximum rounds of constraint propagation
    size_t   CPMAX;
    //! @brief Threshold for repeating constraint propagation (minimum relative reduction in any variable)
    double   CPTHRES;
    //! @brief Infinite value for unbounded variables in constraint propagation
    double   CPINF;
  } options;

  //! @brief Default constructor
  DAGEXT
    ():
    _dag( nullptr )
    {}

  //! @brief Data constructor
  DAGEXT
    ( FFGraph* dag, std::vector<FFVar> const& varin, std::vector<FFVar> const& varout ):
    _dag( nullptr )
    {
      _set( dag, varin, varout );
    }

  //! @brief Copy constructor
  DAGEXT
    ( DAGEXT const& other ):
    _dag( nullptr )
    {
      _set( other._dag, other._varin, other._varout );
      options = other.options;
    }

  virtual ~DAGEXT() 
    {
      delete _dag;
    }

  //! @brief Set DAGEXT data
  void set
    ( FFGraph* dag, std::vector<FFVar> const& varin, std::vector<FFVar> const& varout )
    {
      _set( dag, varin, varout );
    }

  //! @brief DAG evaluation
  void eval
    ( double const* valin, double* valout )
    {
      _eval( valin, valout, _wkD );
    }
  void eval
    ( fadbad::F<double> const* valin, fadbad::F<double>* valout )
    {
      _eval( valin, valout, _wkFD );
    }
  void eval
    ( T const* valin, T* valout )
    {
      _eval( valin, valout, _wkI );
    }
  void eval
    ( McCormick<T> const* valin, McCormick<T>* valout )
    {
      _eval( valin, valout, _wkMC );
    }
  void eval
    ( SCVar<T> const* valin, SCVar<T>* valout )
    {
      _eval( valin, valout, _wkSC );
    }
  void eval
    ( SupVar<PWCU> const* valin, SupVar<PWCU>* valout )
    {
      _eval( valin, valout, _wkPWCS );
    }
  void eval
    ( McCormick<SupVar<PWCU>> const* valin, McCormick<SupVar<PWCU>>* valout )
    {
      _eval( valin, valout, _wkMCPWCS );
    }
  void eval
    ( SupVar<PWLU> const* valin, SupVar<PWLU>* valout )
    {
      _eval( valin, valout, _wkPWLS );
    }
  void eval
    ( McCormick<SupVar<PWLU>> const* valin, McCormick<SupVar<PWLU>>* valout )
    {
      _eval( valin, valout, _wkMCPWLS );
    }
  void eval
    ( PolVar<T> const* valin, PolVar<T>* valout )
    {
      _eval( valin, valout, _wkPOL );
    }
  template <typename U>
  void eval
    ( U const* valin, U* valout )
    { 
      static thread_local std::vector<U> wkU;
      _eval( valin, valout, wkU );
    }

  //! @brief DAG reverse evaluation
  bool reval
    ( T* valin, T const* valout )
    {
      return _reval( valin, valout, _wkCPI, options.CPINF*T(-1,1) );
    }

  //! @brief DAG differentiation
  std::tuple< std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar> > diff
    ()
    const
    {
      switch( options.AUTODIFF ){
      case Options::F:
      default:
        return _dag->SFAD( _varout, _varin );
      case Options::B:
        return _dag->SBAD( _varout, _varin );
      }
    }

  //! @brief Query members
  size_t nin
    ()
    const
    {
      return _varin.size();
    }
  size_t nout
    ()
    const
    {
      return _varout.size();
    }
  std::vector<FFVar> const& varin
    ()
    const
    {
      return _varin;
    }
  std::vector<FFVar> const& varout
    ()
    const
    {
      return _varout;
    }
  FFGraph* dag
    ()
    const
    {
      return _dag;
    }
};

//! @brief C++ class defining expression tree as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFDAGEXT is a C++ class for defining expression trees as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFDAGEXT
////////////////////////////////////////////////////////////////////////
: public FFEXTERN<T,DAGEXT<T>>
{

protected:

  using FFEXTERN<T,DAGEXT<T>>::_ptrObj;
  using FFEXTERN<T,DAGEXT<T>>::_ownObj;
    
  // set the object and related operation in DAG
  FFVar** _set
    ( size_t const nVar, FFVar const* pVar, DAGEXT<T>* pDAG, int policy )
    {
#ifdef MC__FFDAGEXT_CHECK
      assert( nVar == pDAG->nin() );
#endif
      if( _ownObj && _ptrObj ) delete _ptrObj;
      _ownObj = ( policy>0? true: false ); //copy;
      //_ownObj = true;
      this->owndata = false;
      this->data = _ptrObj = pDAG;
      //this->sparse = sparse;

      FFVar** ppRes = this->insert_external_operation( *this, pDAG->nout(), nVar, pVar );

      _ownObj = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _ptrObj = static_cast<FFDAGEXT<T>*>(pOp)->_ptrObj; // set pointer to DAGEXT copy
      else if( policy < 0 )
        static_cast<FFDAGEXT<T>*>(pOp)->_ownObj = true; // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFDAGEXT_DEBUG
      std::cerr << "DAGEXT operation address: " << this << std::endl;
      std::cerr << "DAGEXT address in DAG: " << _ptrObj << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for copy policy of DAGEXT object
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of DAGEXT object in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of DAGEXT object in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of DAGEXT object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFDAGEXT
    ( bool const sparse=true )
    : FFEXTERN<T,DAGEXT<T>>()
    {
      this->sparse = sparse;
    }

  // Destructor
  virtual ~FFDAGEXT
    ()
    {}

  // Copy constructor
  FFDAGEXT
    ( FFDAGEXT<T> const& Other )
    : FFEXTERN<T,DAGEXT<T>>( Other )
    {}

  // Define operation
  //FFVar** operator()
  std::vector<FFVar> operator()
    ( std::vector<FFVar> const& vVar, DAGEXT<T>* pDAG, int policy=COPY )
    {
      //return _set( vVar.size(), vVar.data(), pDAG, policy );
      FFVar** ppDer = _set( vVar.size(), vVar.data(), pDAG, policy );
      std::vector<FFVar> vDer( vVar.size() );
      for( size_t i=0; i<vVar.size(); ++i ) vDer[i] = *ppDer[i];
      return std::move( vDer );
    }

  FFVar& operator()
    ( size_t const idep, std::vector<FFVar> const& vVar, DAGEXT<T>* pDAG, int policy=COPY )
    {
#ifdef MC__FFDAGEXT_CHECK
      assert( idep < pDAG->nout() );
#endif
      return *(_set( vVar.size(), vVar.data(), pDAG, policy )[idep]);
    }

  FFVar** operator()
    ( size_t const nVar, FFVar const* pVar, DAGEXT<T>* pDAG, int policy=COPY )
    {
      return _set( nVar, pVar, pDAG, policy );
    }

  FFVar& operator()
    ( size_t const idep, size_t const nVar, FFVar const* pVar, DAGEXT<T>* pDAG, int policy=COPY )
    {
#ifdef MC__FFDAGEXT_CHECK
      assert( idep < pDAG->nout() );
#endif
      return *(_set( nVar, pVar, pDAG, policy )[idep]);
    }

  DAGEXT<T>* pDAGEXT
    ()
    const
    {
      //std::cerr << "DAGEXT address retreived: " << _ptrObj << std::endl;
      return _ptrObj;
    }

  // Forward evaluation overloads
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
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<double> ) )
        return eval( nRes, static_cast<fadbad::F<double>*>(vRes), nVar, static_cast<fadbad::F<double> const*>(vVar), mVar );
      else if( idU == typeid( T ) )
        return eval( nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar );
      else if( idU == typeid( PolVar<T> ) )
        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );
      else if( idU == typeid( McCormick<T> ) )
        return eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWCU> ) )
        return eval( nRes, static_cast<SupVar<PWCU>*>(vRes), nVar, static_cast<SupVar<PWCU> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWLU> ) )
        return eval( nRes, static_cast<SupVar<PWLU>*>(vRes), nVar, static_cast<SupVar<PWLU> const*>(vVar), mVar );
      else if( idU == typeid( McCormick<SupVar<PWCU>> ) )
        return eval( nRes, static_cast<McCormick<SupVar<PWCU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWCU>> const*>(vVar), mVar );
      else if( idU == typeid( McCormick<SupVar<PWLU>> ) )
        return eval( nRes, static_cast<McCormick<SupVar<PWLU>>*>(vRes), nVar, static_cast<McCormick<SupVar<PWLU>> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFDAGEXT::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  template< typename U >
  void eval
    ( size_t const nRes, U* vRes, size_t const nVar, U const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar, unsigned const* mVar )
    const;

  // Backward evaluation overloads
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      if( idU == typeid( T ) )
        return reval( nRes, static_cast<T const*>(vRes), nVar, static_cast<T*>(vVar) );
      else if( idU == typeid( PolVar<T> ) )
        return reval( nRes, static_cast<PolVar<T> const*>(vRes), nVar, static_cast<PolVar<T>*>(vVar) );

      throw std::runtime_error( "FFDAGEXT::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  bool reval
    ( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
    const;

  bool reval
    ( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
      FFVar** vDer )
    const;

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
      FFVar** vDer, size_t* nnz, size_t** colnz )
    const;

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << this->data; return "DAGEXT[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};
/*
template< typename T >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: FFDep\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  // COULD THIS LEVERAGE _ptrObj->eval in FFDep?
  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}
*/
template< typename T >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: FFVar\n";
  std::cerr << "DAGEXT operation address: " << this    << std::endl;
  std::cerr << "DAGEXT address in DAG: "    << _ptrObj << std::endl;
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  FFVar** ppRes = this->insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template< typename T >
template< typename U >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  _ptrObj->eval( vVar, vRes );
}

template< typename T >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = this->insert_external_operation( *this, nRes, nVar, vVarVal.data() );

  auto&& sdervarout = _ptrObj->diff();
  DAGEXT<T>* pDAGDer = new DAGEXT<T>( _ptrObj->dag(), _ptrObj->varin(), std::get<2>(sdervarout) );
  FFDAGEXT<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVarVal.data(), pDAGDer, TRANSFER );

  // initialize forward AD and set all derivatives to zero
  for( unsigned k=0; k<nRes; ++k ){
    vRes[k] = *vResVal[k];
    for( unsigned i=0; i<nVar; ++i )
      vRes[k].setDepend( vVar[i] );
    for( unsigned j=0; j<vRes[k].size(); ++j )
      vRes[k][j] = 0.;
  }

  // accumulate nonzero derivatives from sparse Jacobian
  for( unsigned int ie=0; ie<std::get<0>(sdervarout).size(); ++ie ){
    for( unsigned j=0; j<vRes[std::get<0>(sdervarout)[ie]].size(); ++j ){
      auto const& k = std::get<0>(sdervarout)[ie];
      auto const& i = std::get<1>(sdervarout)[ie];
      if( vVar[i][j].cst() && vVar[i][j].num().val() == 0. ) continue;
      vRes[k][j] += *vResDer[ie] * vVar[i][j];
      //vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
    }
  }
}

template< typename T >
inline void
FFDAGEXT<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer, size_t* nnz, size_t** colnz )
const
{
//#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::deriv (sparse): FFVar\n";
//#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
  assert( this->sparse && nnz && colnz );
#endif

  auto&& sdervarout = _ptrObj->diff();
  DAGEXT<T>* pDAGDer = new DAGEXT<T>( _ptrObj->dag(), _ptrObj->varin(), std::get<2>(sdervarout) );
  FFDAGEXT<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVar, pDAGDer, TRANSFER );

  // count non-zero dependencies for each output
  for( size_t k=0; k<nRes; ++k ) nnz[k] = 0;
  for( size_t ie=0; ie<std::get<0>(sdervarout).size(); ++ie )
    ++nnz[std::get<0>(sdervarout)[ie]];

  // resize each derivative and component arrays
  for( size_t k=0; k<nRes; ++k ){
    assert( !vDer[k] && !colnz[k] );
    if( !nnz[k] ) continue;
    vDer[k]  = new FFVar[nnz[k]];
    colnz[k] = new size_t[nnz[k]];
//#ifdef MC__FFDAGEXT_DEBUG
    std::cout << "NNZ[" << k << "] = " << nnz[k] << std::endl;
//#endif
  }

  // copy nonzero derivatives from sparse Jacobian
  for( size_t k=0; k<nRes; ++k ) nnz[k] = 0;
  for( size_t ie=0; ie<std::get<0>(sdervarout).size(); ++ie ){
    auto const& k = std::get<0>(sdervarout)[ie];
    auto const& i = std::get<1>(sdervarout)[ie];
    colnz[k][nnz[k]] = i;
    vDer[k][nnz[k]]  = *vResDer[ie];
    ++nnz[k];
  }
}

template< typename T >
inline void
FFDAGEXT<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
  FFVar** vDer )
const
{
//#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::deriv (dense): FFVar\n";
//#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( !this->sparse && _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  auto&& sdervarout = _ptrObj->diff();
  DAGEXT<T>* pDAGDer = new DAGEXT<T>( _ptrObj->dag(), _ptrObj->varin(), std::get<2>(sdervarout) );
  FFDAGEXT<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVar, pDAGDer, TRANSFER );

  // set all derivatives to zero
  for( size_t k=0; k<nRes; ++k )
    for( size_t i=0; i<nVar; ++i )
        vDer[k][i] = 0;

  // copy nonzero derivatives from sparse Jacobian
  for( size_t ie=0; ie<std::get<0>(sdervarout).size(); ++ie ){
    auto const& k = std::get<0>(sdervarout)[ie];
    auto const& i = std::get<1>(sdervarout)[ie];
    vDer[k][i] = *vResDer[ie]; //[k+nRes*i];
  }
}

template< typename T >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: SLiftVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFDAGEXT<T>::eval
( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: PolVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__FFDAGEXT_CHECK
  assert( img && dag );
#endif
  FFVar** ppRes = dag->curOp()->varout.data(); // ACCOUNT FOR MULTIPLE OUTPUTS
#ifdef MC__FFDAGEXT_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  this->_resize_relax( img );
  this->_propagate_relax( img, ppRes, vRes, vVar );
}

template< typename T >
inline bool
FFDAGEXT<T>::reval
( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::reval: PolVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__FFDAGEXT_CHECK
  assert( img && pop );
#endif

  this->_backpropagate_relax( img, pop, vRes, vVar );
  return true;
}

template< typename T >
inline bool
FFDAGEXT<T>::reval
( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::reval: T\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  return _ptrObj->reval( vVar, vRes );
}

} // end namespace mc

#endif
