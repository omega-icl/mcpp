#ifndef MC__FFMLP_HPP
#define MC__FFMLP_HPP

//#define MC__FFMLP_DEBUG
#define MC__FFMLP_CHECK

#include "fflin.hpp"
#include "ffextern.hpp"

namespace mc
{

//! @brief C++ class for evaluation and relaxation of multilayer perceptrons
////////////////////////////////////////////////////////////////////////
//! mc::MLP is a C++ class for evaluation and relaxation of (trained)
//! multilayer perceptrons (MLP) that leverages expression trees and
//! arithmetics available through MC++
////////////////////////////////////////////////////////////////////////
template <typename T> 
class MLP
////////////////////////////////////////////////////////////////////////
{
public:

  //! @brief Enumeration type for activation function
  enum ACTIVTYPE{
    LINEAR=0, //!< Linear activation function
    RELU,     //!< ReLU activation function
    TANH,     //!< tanh activation function
    SIGMOID   //!< Sigmoid activation function
  };

  //! @brief MLP data
  std::vector<std::pair<std::vector<std::vector<double>>,int>> data;

  //! @brief Reset MLP data
  void reset_data
    ();

  //! @brief Set MLP data
  bool set_data
    ( std::vector<std::pair<std::vector<std::vector<double>>,int>> const& data );

  //! @brief Append multi-neuron layer to MLP data: outer vector size is #neurons in hidden layer; inner vector has scalar weights multiplying each neuron in previous (input or hidden) layer and a bias term and size 1+#neuron in previous layer 
  bool append_data
    ( std::vector<std::vector<double>> const& data, int const activ=LINEAR, bool const reset=false );

  //! @brief Append single-neuron layer to MLP data: vector has scalar weights multiplying each neuron in previous (input or hidden) layer and a bias term and size 1+#neuron in previous layer 
  bool append_data
    ( std::vector<double> const& data, int const activ=LINEAR, bool const reset=false );

  //! @brief MLP options
  struct Options
  {
    //! @brief Default constructor
    Options()
      {
        reset();
      }

    //! @brief Assignment operator
    Options& operator=
      ( Options const& other )
      {
        ZEROTOL   = other.ZEROTOL;
        RELU2ABS  = other.RELU2ABS;
        SIG2EXP   = other.SIG2EXP;
        AUTODIFF  = other.AUTODIFF;
        CPMAX     = other.CPMAX;
        CPTHRES   = other.CPTHRES;
        CPINF     = other.CPINF;

        return *this;
      }

    //! @brief Reset options
    void reset
      ()
      {
        ZEROTOL  = DBL_EPSILON;
        RELU2ABS = false;
        SIG2EXP  = false;
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

    //! @brief Threshold for zero coefficient in neural network
    double   ZEROTOL;
    //! @brief Whether to convert ReLU to abs (true) or max (false) function
    bool     RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false)
    bool     SIG2EXP;
    //! @brief Whether to apply forward or reverse automatic differentiation
    int      AUTODIFF;
    //! @brief Maximum rounds of constraint propagation
    size_t   CPMAX;
    //! @brief Threshold for repeating constraint propagation (minimum relative reduction in any variable)
    double   CPTHRES;
    //! @brief Infinite value for unbounded variables in constraint propagation
    double   CPINF;
  } options;

private:

  //! @brief Expression tree
  FFGraph*                               _dag;
  //! @brief Expression tree needs updating
  bool                                   _update;
  //! @brief Number of variables
  size_t                                 _nin;
  //! @brief Number of dependents
  size_t                                 _nout;
  //! @brief Number of hidden layers
  size_t                                 _nhid;
  //! @brief Variables
  std::vector<FFVar>                     _varin;
  //! @brief Dependents
  std::vector<FFVar>                     _varout;
  //! @brief Codelist
  FFSubgraph                             _codelist;

  //! @brief Intermediate storage for DAG evaluation
  std::vector<std::vector<FFVar>>                    _wkFF;
  std::vector<std::vector<double>>                   _wkD;
  std::vector<std::vector<fadbad::F<double>>>        _wkFD;
  std::vector<std::vector<fadbad::B<double>>>        _wkBD;
  std::vector<std::vector<T>>                        _wkI;
  std::vector<std::vector<McCormick<T>>>             _wkMC;
  std::vector<std::vector<SupVar<PWCU>>>             _wkPWCS;
  std::vector<std::vector<McCormick<SupVar<PWCU>>>>  _wkMCPWCS;
  std::vector<std::vector<SupVar<PWLU>>>             _wkPWLS;
  std::vector<std::vector<McCormick<SupVar<PWLU>>>>  _wkMCPWLS;
  std::vector<std::vector<SCVar<T>>>                 _wkSC;
  std::vector<PolVar<T>>                             _wkPOL;
  std::vector<T>                                     _wkCPI;

  //! @brief ReLU activation
  template <typename U>
  U ReLU
    ( U const& x )
    const
    {
      return options.RELU2ABS?
             ( x + Op<U>::fabs(x) ) * 0.5:
             Op<U>::max( x, 0. );
    }
  template <typename U>
  fadbad::F<U> ReLU
    ( fadbad::F<U> const& x )
    const
    {
      fadbad::F<U> z = ReLU( x.val() );
      z.setDepend( x );
      for( unsigned j=0; j<z.size(); ++j )
        z[j] = Op<U>::fstep( x.val() ) * x[j];
      return z;
    }

  //! @brief Set expression tree
  void _set_dag
    ( std::vector<std::vector<FFVar>>& wkhid );

  //! @brief Evaluate MLP
  template <typename U>
  void _eval
    ( U const* valin, U* valout, std::vector<std::vector<U>>& wkhid )
    const;

  //! @brief Evaluate MLP through expression tree
  template <typename U>
  void _eval
    ( U const* valin, U* valout, std::vector<U>& wk );
//    const;

  //! @brief Reverse evaluate MLP through expression tree
  template <typename U>
  bool _reval
    ( U* valin, U const* valout, std::vector<U>& wk, U const& inf );
//    const;

public:

  //! @brief Default constructor
  MLP
  ()
  : _dag    ( nullptr ),
    _update ( false ),
    _nin     ( 0 ),
    _nout    ( 0 ),
    _nhid    ( 0 )
  {}
/*
  //! @brief Data constructor
  DAG
    ( FFGraph* dag, std::vector<FFVar> const& varin, std::vector<FFVar> const& varout ):
    _dag( nullptr )
    {
      _set( dag, varin, varout );
    }
*/
  //! @brief Copy constructor
  MLP
    ( MLP const& other )
    : data    ( other.data ),
      options ( other.options ),
      _dag    ( nullptr ),
      _update ( false ),
      _nin    ( other._nin ),
      _nout   ( other._nout ),
      _nhid   ( other._nhid )
    {}

  virtual ~MLP() 
    {
      delete _dag;
    }

  //! @brief Set MLP DAG
  void set_dag
    ()
    {
      _set_dag( _wkFF );
    }

  //! @brief MLP evaluation
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
    ( fadbad::B<double> const* valin, fadbad::B<double>* valout )
    {
      _eval( valin, valout, _wkBD );
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
      // Update MLP DAG first
      _set_dag( _wkFF );
      _eval( valin, valout, _wkPOL );
    }
  template <typename U>
  void eval
    ( U const* valin, U* valout )
    { 
      static thread_local std::vector<std::vector<U>> wkU;
      _eval( valin, valout, wkU );
    }

  //! @brief MLP reverse evaluation
  bool reval
    ( T* valin, T const* valout )
    {
      // Update MLP DAG first
      _set_dag( _wkFF );
      return _reval( valin, valout, _wkCPI, options.CPINF*T(-1,1) );
    }

  //! @brief Query members
  size_t nin
    ()
    const
    {
      return _nin;
    }
  size_t nout
    ()
    const
    {
      return _nout;
    }
  size_t nhid
    ()
    const
    {
      return _nhid;
    }
  std::vector<FFVar> const& varin
    ()
    //const
    {
      _set_dag( _wkFF );
      return _varin;
    }
  std::vector<FFVar> const& varout
    ()
    //const
    {
      _set_dag( _wkFF );
      return _varout;
    }
  FFGraph* dag
    ()
    //const
    {
      _set_dag( _wkFF );
      return _dag;
    }
};

template< typename T >
inline bool
MLP<T>::set_data
( std::vector<std::pair<std::vector<std::vector<double>>,int>> const& mlp )
{
  if( !mlp.size()
   || !mlp[0].first.size()
   || (mlp[0].first)[0].size() <= 1 )
    return false;

  // Set data
  data  = mlp;
  _nin  = data.front().first.front().size() - 1;
  _nout = data.back().first.size();
  _nhid = data.size() - 1;
  _update = true;

  // Cleanse data
  for( auto& [layer,activ] : data )
    for( auto& neuron : layer )
      for( auto& weight : neuron )
        if( std::fabs(weight) < options.ZEROTOL )
          weight = 0.;

  return true;
}

template< typename T >
inline void
MLP<T>::reset_data
()
{
  // Reset data
  _nin = _nout = _nhid = 0;
  data.clear();
  _update = true;
}

template< typename T >
inline bool
MLP<T>::append_data
( std::vector<std::vector<double>> const& layer, int const activ, bool const reset )
{
  if( reset ){
    reset_data();
    _update = true;
  }
  if( !layer.size()
   || layer[0].size() <= 1
   || ( !data.empty() && data.back().first.size() != layer[0].size()-1 ) )
    return false;

  // Set data
  data.push_back( std::make_pair( layer, activ ) );
  _nin  = data.front().first.front().size()-1;
  _nout = data.back().first.size();
  _nhid = data.size()-1;
  _update = true;

  // Cleanse data
  for( auto& neuron : data.back().first )
    for( auto& weight : neuron )
      if( std::fabs(weight) < options.ZEROTOL )
        weight = 0.;

  return true;
}

template< typename T >
inline bool
MLP<T>::append_data
( std::vector<double> const& layer, int const activ, bool const reset )
{
  return append_data( std::vector<std::vector<double>>({layer}), activ, reset );
}

template< typename T >
inline void
MLP<T>::_set_dag
( std::vector<std::vector<FFVar>>& wkhid )
{
  // Nothing to be updated
  if( _dag && !_update ) return;

  delete _dag;
  _dag = nullptr;
  if( !_nin || !_nout ) return;

  // Created new DAG for MLP
  _dag = new FFGraph;
  _varin  = _dag->add_vars( _nin,  "X" );
  _varout = _dag->add_vars( _nout, "Y" );

  // Propagate DAG through hidden layers
  wkhid.resize( _nhid );
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << _nhid << std::endl;
#endif
  for( unsigned l=0; l<_nhid; ++l ){
    size_t const nneu = data[l].first.size();
#ifdef MC__FFMLP_CHECK
    assert( nneu ); // number of neurons in hidden layer l+1
#endif
    wkhid[l].resize( nneu );
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i ){
      FFLin<I> sum;
      // Need to clean data beforehand
      wkhid[l][i] = sum( (data[l].first)[i].size()-1, l? wkhid[l-1].data(): _varin.data(),
                         (data[l].first)[i].data()+1, (data[l].first)[i][0], FFLin<I>::SHALLOW );
#ifdef MC__FFMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": " << (data[l].first)[i].size()-1 << std::endl;
#endif
      switch( data[l].second ){
        case LINEAR:  default:                                             break;
        case RELU:    wkhid[l][i] = ReLU( wkhid[l][i] );                   break;
        case TANH:    wkhid[l][i] = tanh( wkhid[l][i] );                   break;
        case SIGMOID: wkhid[l][i] = (options.SIG2EXP?
                                     1. / ( exp( -wkhid[l][i] ) + 1. ):
                                     tanh( wkhid[l][i]*0.5 ) * 0.5 + 0.5); break;
      }
    }
  }

  // Propagate DAG through output layer
#ifdef MC__FFMLP_CHECK
  assert( data.back().first.size() == _nout ); // number of neurons in output layer
#endif
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << _nhid << ": " << _nout << std::endl;
#endif
  for( unsigned i=0; i<_nout; ++i ){
    FFLin<I> sum;
    // Need to clean data beforehand
    _varout[i] = sum( (data.back().first)[i].size()-1, _nhid? wkhid[_nhid-1].data(): _varin.data(),
                      (data.back().first)[i].data()+1, (data.back().first)[i][0], FFLin<I>::SHALLOW );
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << _nhid << ": " << (data.back().first)[i].size()-1 << std::endl;
#endif
    switch( data.back().second ){
      case LINEAR:  default:                                             break;
      case RELU:    _varout[i] = ReLU( _varout[i] );                     break;
      case TANH:    _varout[i] = tanh( _varout[i] );                     break;
      case SIGMOID: _varout[i] = (options.SIG2EXP?
                                  1. / ( exp( -_varout[i] ) + 1. ):
                                  tanh( _varout[i]*0.5 ) * 0.5 + 0.5);   break;
    }
  }

#ifdef MC__FFMLP_DEBUG
  _codelist = _dag->subgraph( _varout );
//  _dag->output( _codelist );
  std::vector<FFExpr> strout = FFExpr::subgraph( _dag, _codelist );
  for( size_t i=0; i<_nout; ++i )
    std::cout << "F" << i << ": " << strout[i] << std::endl;
//  std::ofstream oFile( "MLP.dot", std::ios_base::out );
//  _dag->dot_script( _varout, oFile );
//  oFile.close();
#else
  _codelist.clear();
#endif
  _update = false;
}

template< typename T >
template< typename U >
inline void
MLP<T>::_eval
( U const* valin, U* valout, std::vector<std::vector<U>>& wkhid )
const
{
  // Propagate through hidden layers
  wkhid.resize( _nhid );
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << _nhid << std::endl;
#endif
  for( unsigned l=0; l<_nhid; ++l ){
#ifdef MC__FFMLP_CHECK
    assert( data[l].first.size() ); // number of neurons in layer l+1
#endif
    size_t const nneu = data[l].first.size();
    wkhid[l].resize( nneu );
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i ){
      wkhid[l][i] = (data[l].first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": " << (data[l].first)[i].size()-1 << std::endl;
#endif
      for( unsigned j=0; j<(data[l].first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
        std::cout << "layer:" << l << " neuron:" << i << " input:" << j << std::endl;
#endif
        if( std::fabs((data[l].first)[i][1+j]) < options.ZEROTOL ) continue;
        wkhid[l][i] += (l? wkhid[l-1][j]: valin[j]) * (data[l].first)[i][1+j];
      }
      switch( data[l].second ){
        case LINEAR:  default:                                                    break;
        case RELU:    wkhid[l][i] = ReLU( wkhid[l][i] );                          break;
        case TANH:    wkhid[l][i] = Op<U>::tanh( wkhid[l][i] );                   break;
        case SIGMOID: wkhid[l][i] = (options.SIG2EXP?
                                     1. / ( Op<U>::exp( -wkhid[l][i] ) + 1. ):
                                     Op<U>::tanh( wkhid[l][i]*0.5 ) * 0.5 + 0.5); break;
      }
    }
  }

  // Propagate through output layers
#ifdef MC__FFMLP_CHECK
  assert( data.back().first.size() ); // number of neurons in layer l+1
#endif
  size_t const nneu = data.back().first.size();
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << _nhid << ": " << nneu << std::endl;
#endif
  for( unsigned i=0; i<nneu; ++i ){
    valout[i] = (data.back().first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << _nhid << ": " << (data.back().first)[i].size()-1 << std::endl;
#endif
    for( unsigned j=0; j<(data.back().first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
      std::cout << "layer:" << _nhid << " neuron:" << i << " input:" << j << std::endl;
#endif
      if( std::fabs((data.back().first)[i][1+j]) < options.ZEROTOL ) continue;
      valout[i] += (_nhid? wkhid[_nhid-1][j]: valin[j]) * (data.back().first)[i][1+j];
    }
    switch( data.back().second ){
      case LINEAR:  default:                                                break;
      case RELU:    valout[i] = ReLU( valout[i] );                          break;
      case TANH:    valout[i] = Op<U>::tanh( valout[i] );                   break;
      case SIGMOID: valout[i] = (options.SIG2EXP?
                                 Op<U>::inv( Op<U>::exp( -valout[i] ) + 1. ):
                                 Op<U>::tanh( valout[i]*0.5 ) * 0.5 + 0.5); break;
    }
  }
}

template< typename T >
template< typename U >
inline void
MLP<T>::_eval
( U const* valin, U* valout, std::vector<U>& wk )
//const
{
  // Run eval on MLP DAG
  _dag->eval( _codelist, wk, _varout.size(), _varout.data(), valout, _varin.size(), _varin.data(), valin );
}

template< typename T >
template< typename U >
inline bool
MLP<T>::_reval
( U* valin, U const* valout, std::vector<U>& wk, U const& inf )
//const
{
  // Run reval on MLP DAG
  int flag = _dag->reval( _codelist, wk, _varout.size(), _varout.data(), const_cast<U*>( valout ),
                          _varin.size(), _varin.data(), valin, inf, options.CPMAX, options.CPTHRES );

#ifdef MC__FFMLP_DEBUG
  std::cout << "MLP:: Work array: " << flag << " passes\n";
  for( unsigned i=0; i<_codelist.len_tap-_codelist.len_wrk; i++ )
    std::cout << "wk[" << i << "] = " << wk[i] << std::endl;
  for( unsigned i=0; i<_varin.size(); i++ )
    std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif

  return( flag<0? false: true );
}

//! @brief C++ class defining neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLP is a C++ class for defining neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFMLP
////////////////////////////////////////////////////////////////////////
: public FFEXTERN<T,MLP<T>>
{

protected:

  using FFEXTERN<T,MLP<T>>::_ptrObj;
  using FFEXTERN<T,MLP<T>>::_ownObj;

  // set the object and related operation in DAG
  FFVar** _set
    ( size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() );
#endif
      if( _ownObj && _ptrObj ) delete _ptrObj;
      _ownObj = ( policy>0? true: false ); //copy;
      //_ownObj = true;
      this->owndata = false;
      this->data = _ptrObj = pMLP;

      FFVar** ppRes = this->insert_external_operation( *this, pMLP->nout(), nVar, pVar );

      _ownObj = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _ptrObj = static_cast<FFMLP<T>*>(pOp)->_ptrObj; // set MLP pointer to DAG copy
      else if( policy < 0 )
        static_cast<FFMLP<T>*>(pOp)->_ownObj = true; // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFMLP_DEBUG
      std::cerr << "MLP operation address: " << this << std::endl;
      std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFMLP
    ()
    : FFEXTERN<T,MLP<T>>()
    {}

  // Destructor
  virtual ~FFMLP
    ()
    {}

  // Copy constructor
  FFMLP
    ( FFMLP<T> const& Other )
    : FFEXTERN<T,MLP<T>>( Other )
    {}

  // Define operation
  FFVar** operator()
    ( std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( vVar.size() == pMLP->nin() );
#endif
      return _set( vVar.size(), vVar.data(), pMLP, policy );
    }

  FFVar& operator()
    ( size_t const idep, std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( vVar.size() == pMLP->nin() && idep < pMLP->nout() );
#endif
      return *(_set( vVar.size(), vVar.data(), pMLP, policy )[idep]);
    }

  FFVar** operator()
    ( size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() );
#endif
      return _set( nVar, pVar, pMLP, policy );
    }

  FFVar& operator()
    ( size_t const idep, size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() && idep < pMLP->nout() );
#endif
      return *(_set( nVar, pVar, pMLP, policy )[idep]);
    }

  // MLP pointer
  MLP<T>* pMLP
    ()
    const
    {
      //std::cerr << "MLP address retreived: " << _ptrObj << std::endl;
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

      throw std::runtime_error( "FFMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  template< typename U >
  void eval
    ( size_t const nRes, U* vRes, size_t const nVar, U const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar, unsigned const* mVar )
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

      throw std::runtime_error( "FFMLP::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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
    
  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << this->data; return "MLP[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

//! @brief C++ class defining gradient of neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradMLP is a C++ class for defining gradient of neural 
//! networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFGradMLP
////////////////////////////////////////////////////////////////////////
: public FFEXTERN<T,MLP<T>>
{
  friend class FFMLP<T>;
  
protected:

  using FFEXTERN<T,MLP<T>>::_ptrObj;
  using FFEXTERN<T,MLP<T>>::_ownObj;

  // set the object and related operation in DAG
  FFVar** _set
    ( size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() );
#endif
      if( _ownObj && _ptrObj ) delete _ptrObj;
      _ownObj = ( policy>0? true: false ); //copy;
      //_ownObj = true;
      this->owndata = false;
      this->data = _ptrObj = pMLP;

      FFVar** ppRes = this->insert_external_operation( *this, nVar*pMLP->nout(), nVar, pVar );

      _ownObj = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _ptrObj = static_cast<FFGradMLP<T>*>(pOp)->_ptrObj; // set MLP pointer to DAG copy
      else if( policy < 0 )
        static_cast<FFGradMLP<T>*>(pOp)->_ownObj = true; // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFMLP_DEBUG
      std::cerr << "GradMLP operation address: " << this << std::endl;
      std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFGradMLP
    ()
    : FFEXTERN<T,MLP<T>>()
    {}

  // Destructor
  virtual ~FFGradMLP
    ()
    {}

  // Copy constructor
  FFGradMLP
    ( FFGradMLP<T> const& Other )
    : FFEXTERN<T,MLP<T>>( Other )
    {}

  // Define operation
  FFVar** operator()
    ( std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( vVar.size() == pMLP->nin() );
#endif
      return _set( vVar.size(), vVar.data(), pMLP, policy );
    }

  FFVar& operator()
    ( size_t const idep, std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( vVar.size() == pMLP->nin() && idep < pMLP->nout()*pMLP->nin() );
#endif
      return *(_set( vVar.size(), vVar.data(), pMLP, policy )[idep]);
    }

  FFVar** operator()
    ( size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() );
#endif
      return _set( nVar, pVar, pMLP, policy );
    }

  FFVar& operator()
    ( size_t const idep, size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy=COPY )
    {
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin() && idep < pMLP->nout()*pMLP->nin() );
#endif
      return *(_set( nVar, pVar, pMLP, policy )[idep]);
    }

  // MLP pointer
  MLP<T>* pMLP
    ()
    const
    {
      //std::cerr << "MLP address retreived: " << _ptrObj << std::endl;
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
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );

      throw std::runtime_error( "FFGradMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( size_t const nRes, double* vRes, size_t const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  // Backward evaluation overloads
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      throw std::runtime_error( "FFGradMLP::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << this->data; return "GradMLP[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

template< typename T >
inline void
FFMLP<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFDep\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFGradMLP<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFDep\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout()*_ptrObj->nin() && nVar == _ptrObj->nin() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFMLP<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFVar\n";
  std::cerr << "FFMLP operation address: " << this    << std::endl;
  std::cerr << "MLP address in DAG: "    << _ptrObj << std::endl;
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  FFVar** ppRes = this->insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template< typename T >
inline void
FFGradMLP<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFVar\n";
  std::cerr << "FFGradMLP operation address: " << this    << std::endl;
  std::cerr << "MLP address in DAG: "    << _ptrObj << std::endl;
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout()*_ptrObj->nin() && nVar == _ptrObj->nin() );
#endif

  FFVar** ppRes = this->insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template< typename T >
template< typename U >
inline void
FFMLP<T>::eval
( size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  _ptrObj->eval( vVar, vRes );
}

template< typename T >
inline void
FFGradMLP<T>::eval
( size_t const nRes, double* vRes, size_t const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: double\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout()*_ptrObj->nin() && nVar == _ptrObj->nin() );
#endif

  switch( _ptrObj->options.AUTODIFF ){
   default:
   case MLP<T>::Options::AD::F:{
    std::vector<fadbad::F<double>> vFVar( _ptrObj->nin() );
    for( unsigned i=0; i<_ptrObj->nin(); ++i ){
      vFVar[i] = vVar[i];
      vFVar[i].diff(i,_ptrObj->nin());
    }
    std::vector<fadbad::F<double>> vFRes( _ptrObj->nout() ); 
    _ptrObj->eval( vFVar.data(), vFRes.data() );
    for( unsigned k=0; k<_ptrObj->nout(); ++k )
      for( unsigned i=0; i<_ptrObj->nin(); ++i )
        vRes[k*_ptrObj->nin()+i] = vFRes[k].d(i);
    break;
   }
   
   case MLP<T>::Options::AD::B:{
    std::vector<fadbad::B<double>> vBVar( _ptrObj->nin() );
    for( unsigned i=0; i<_ptrObj->nin(); ++i )
      vBVar[i] = vVar[i];
    std::vector<fadbad::B<double>> vBRes( _ptrObj->nout() ); 
    _ptrObj->eval( vBVar.data(), vBRes.data() );
    for( unsigned k=0; k<_ptrObj->nout(); ++k )
      vBRes[k].diff( k, _ptrObj->nout() );
    for( unsigned k=0; k<_ptrObj->nout(); ++k )
      for( unsigned i=0; i<_ptrObj->nin(); ++i )
        vRes[k*_ptrObj->nin()+i] = vBVar[i].d(k);
    break;
   }
  }
}

template< typename T >
inline void
FFMLP<T>::eval
( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = this->insert_external_operation( *this, nRes, nVar, vVarVal.data() );

  FFGradMLP<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVarVal.data(), _ptrObj, COPY );

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

template< typename T >
inline void
FFMLP<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::deriv\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  FFGradMLP<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVar, _ptrObj, COPY );
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
      vDer[k][i] = *vResDer[k+nRes*i];
}

template< typename T >
inline void
FFMLP<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFGradMLP<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout()*_ptrObj->nin() && nVar == _ptrObj->nin() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFMLP<T>::eval
( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFExpr\n";
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
FFGradMLP<T>::eval
( size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFExpr\n";
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
FFMLP<T>::eval
( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: PolVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__FFMLP_CHECK
  assert( img && dag );
#endif
  FFVar** ppRes = dag->curOp()->varout.data();
#ifdef MC__FFMLP_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  this->_resize_relax( img );
  this->_propagate_relax( img, ppRes, vRes, vVar );
}

template< typename T >
inline bool
FFMLP<T>::reval
( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: PolVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__FFMLP_CHECK
  assert( img && pop );
#endif

  this->_backpropagate_relax( img, pop, vRes, vVar );
  return true;
}

template< typename T >
inline bool
FFMLP<T>::reval
( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: T\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  return _ptrObj->reval( vVar, vRes );
}

} // end namespace mc

#endif
