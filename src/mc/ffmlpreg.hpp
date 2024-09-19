#ifndef MC__FFMLPREG_HPP
#define MC__FFMLPREG_HPP

//#define MC__FFMLPREG_DEBUG
#define MC__FFMLPREG_CHECK

#include "mcfunc.hpp"
#include "ffunc.hpp"
#include "slift.hpp"

namespace mc
{

//! @brief C++ class for training of multilayer perceptrons
////////////////////////////////////////////////////////////////////////
//! mc::MLPREG is a C++ class for training of multilayer perceptrons
//! (MLP) that leverages expression trees and arithmetics available
//! through MC++.
////////////////////////////////////////////////////////////////////////
template <typename T> 
class MLPREG
////////////////////////////////////////////////////////////////////////
{
private:

  //! @brief Number of inputs
  size_t                                                _nin;
  //! @brief Number of outputs
  size_t                                                _nout;
  //! @brief Number of hidden layers
  size_t                                                _nhid;
  //! @brief Number of weights
  size_t                                                _nwei;

  //! @brief MLP structure
  std::pair<size_t,std::vector<std::pair<size_t,int>>>  _struct;

  //! @brief Intermediate storage
  std::vector<std::vector<double>>                      _Dhid;
  //! @brief Intermediate storage
  std::vector<std::vector<fadbad::F<double>>>           _FDhid;
  //! @brief Intermediate storage
  std::vector<std::vector<fadbad::B<double>>>           _BDhid;
    
public:

  //! @brief Default constructor
  MLPREG
    ()
    : _nin(0), _nout(0), _nhid(0), _nwei(0)
    {}

  //! @brief Default constructor
  MLPREG
    ( std::pair<size_t,std::vector<std::pair<size_t,int>>> const& structure )
    {
      set_structure( structure );
    }

  //! @brief Copy constructor
  MLPREG
    ( MLPREG<T> const& MLP )
    : _nin(MLP._nin), _nout(MLP._nout), _nhid(MLP._nhid), _nwei(MLP._nwei),
      _struct(MLP._struct), options(MLP.options)
    {}

  ~MLPREG
    () 
    {}

  //! @brief MLP structure
  void set_structure
    ( std::pair<size_t,std::vector<std::pair<size_t,int>>> const& structure )
    {
      _struct = structure;
      _nin    = _struct.first;
      _nout   = _struct.second.back().first;
      _nhid   = _struct.second.size()-1;

      _nwei = 0;
      size_t nlayin = _nin;
      for( auto const& [nneu,activ] : _struct.second ){
        _nwei += (nlayin+1)*nneu; // +1 accounts for bias term
        nlayin = nneu;
      }
    }

  //! @brief Number of inputs
  size_t nin
    ()
    const
    { return _nin; }

  //! @brief Number of outputs
  size_t nout
    ()
    const
    { return _nout; }

  //! @brief Number of hidden layers
  size_t nhid
    ()
    const
    { return _nhid; }

  //! @brief MLP structure
  size_t nwei
    ()
    const
    { return _nwei; }

  //! @brief MLP options
  struct Options
  {
    //! @brief Activation type
    enum ACTIVTYPE{
      LINEAR=0, //!< Linear activation function
      RELU,     //!< ReLU activation function
      TANH,     //!< tanh activation function
      SIGMOID   //!< Sigmoid activation function
    };

    //! @brief Default constructor
    Options():
      RELU2ABS(false), SIG2EXP(false), AUTODIFF(F)
      {}

    //! @brief Assignment operator
    Options& operator=
      ( Options const& opt ){
        RELU2ABS  = opt.RELU2ABS;
        SIG2EXP   = opt.SIG2EXP;
        AUTODIFF  = opt.AUTODIFF;
        return *this;
      }

    //! @brief Enumeration type for AD strategy
    enum AD{
      F=0,	//!< Forward differentiation
      B		//!< Backward differentiation
    };

    //! @brief Whether to convert ReLU to abs (true) or max (false) function
    bool     RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false)
    bool     SIG2EXP;
    //! @brief Whether to apply forward or reverse automatic differentiation
    int      AUTODIFF;
  } options;

  //! @brief Evaluate neural network
  template <typename U>
  void eval
    ( U* y, U const* x, U const* w )
    const
    { static thread_local std::vector<std::vector<U>> Uhid;
      _eval( y, x, w, Uhid ); }
  void eval
    ( double* y, double const* x, double const* w )
    { _eval( y, x, w, _Dhid ); }
  void eval
    ( fadbad::F<double>* y, fadbad::F<double> const* x, fadbad::F<double> const* w )
    { _eval( y, x, w, _FDhid ); }
  void eval
    ( fadbad::B<double>* y, fadbad::B<double> const* x, fadbad::B<double> const* w )
    { _eval( y, x, w, _BDhid ); _BDhid.clear(); }

private:

  template< typename U >
  void _eval
    ( U* y, U const* x, U const* w, std::vector<std::vector<U>>& vhid )
    const;

  template <typename U>
  U ReLU
    ( U const& x )
    const
    {
      return options.RELU2ABS?
             ( x + Op<U>::fabs(x) ) * 0.5:
             Op<U>::max( x, U(0.) );
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
};


template< typename T >
template< typename U >
inline void
MLPREG<T>::_eval
( U* y, U const* x, U const* w, std::vector<std::vector<U>>& vhid )
const
{
  // Propagate through hidden layers
#ifdef MC__FFMLPREG_DEBUG
  std::cerr << std::scientific << std::setprecision(8);
  std::cerr << "No hidden layers: " << _nhid << std::endl;
#endif
  vhid.resize( _nhid );
  U const* pwei = w;
  size_t nneu = _nin;
  for( unsigned l=0; l<_nhid; ++l ){
    size_t const ncoef = nneu+1;
    nneu = _struct.second[l].first;
#ifdef MC__FFMLPREG_CHECK
    assert( nneu ); // number of neurons in layer l+1
#endif
    vhid[l].resize( nneu );
#ifdef MC__FFMLPREG_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i, pwei+=ncoef ){
      vhid[l][i] = pwei[0]; // bias term
#ifdef MC__FFMLPREG_DEBUG
      std::cerr << "No coefficients in neuron " << i << " in layer " << l << ": " << ncoef << std::endl;
#endif
      for( unsigned j=0; j<ncoef-1; ++j ){
#ifdef MC__FFMLPREG_DEBUG
        std::cerr << "layer:" << l << " neuron:" << i << " input:" << j << std::endl;
#endif
        //if( std::fabs((pwei[1+j]) < options.ZEROTOL ) continue;
        vhid[l][i] += (l? vhid[l-1][j]: x[j]) * pwei[1+j];
      }
      switch( _struct.second[l].second ){
        case Options::LINEAR:  default:                                                  break;
        case Options::RELU:    vhid[l][i] = ReLU( vhid[l][i] );                          break;
        case Options::TANH:    vhid[l][i] = Op<U>::tanh( vhid[l][i] );                   break;
        case Options::SIGMOID: vhid[l][i] = (options.SIG2EXP?
                                             1. / ( Op<U>::exp( -vhid[l][i] ) + 1. ):
                                             Op<U>::tanh( vhid[l][i]*0.5 ) * 0.5 + 0.5); break;
      }
#ifdef MC__FFMLPREG_DEBUG
      //std::cerr << "layer:" << l << " neuron:" << i << " value: " << vhid[l][i] << std::endl;
#endif
    }
  }

  // Propagate through output layers
  size_t const ncoef = nneu+1;
  nneu = _struct.second.back().first;
#ifdef MC__FFMLPREG_CHECK
  assert( nneu ); // number of neurons in output layer
#endif
#ifdef MC__FFMLPREG_DEBUG
  std::cerr << "No neurons in layer " << _nhid << ": " << nneu << std::endl;
#endif
  for( unsigned i=0; i<nneu; ++i, pwei+=ncoef ){
    y[i] = pwei[0]; // bias term
#ifdef MC__FFMLPREG_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << _nhid << ": " << ncoef << std::endl;
#endif
    for( unsigned j=0; j<ncoef-1; ++j ){
#ifdef MC__FFMLPREG_DEBUG
      std::cerr << "layer:" << _nhid << " neuron:" << i << " input:" << j << std::endl;
#endif
      //if( std::fabs((pwei[1+j]) < options.ZEROTOL ) continue;
      y[i] += (_nhid? vhid[_nhid-1][j]: x[j]) * pwei[1+j];
    }
    switch( _struct.second.back().second ){
      case Options::LINEAR:  default:                                      break;
      case Options::RELU:    y[i] = ReLU( y[i] );                          break;
      case Options::TANH:    y[i] = Op<U>::tanh( y[i] );                   break;
      case Options::SIGMOID: y[i] = (options.SIG2EXP?
                                     Op<U>::inv( Op<U>::exp( -y[i] ) + 1. ):
                                     Op<U>::tanh( y[i]*0.5 ) * 0.5 + 0.5); break;
    }
#ifdef MC__FFMLPREG_DEBUG
    //std::cerr << "layer:" << _nhid << " neuron:" << i << " value: " << y[i] << std::endl;
#endif
  }
}
/*
template< typename T >
template< typename U >
inline void
MLPREG<T>::_eval
( U* y, U const* x, std::vector<std::pair<std::vector<std::vector<U>>,int>> vcoef,
  std::vector<std::vector<U>>& vhid )
const
{
  // Propagate through hidden layers
  vhid.resize( nhid );
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << nhid << std::endl;
#endif
  for( unsigned l=0; l<nhid; ++l ){
#ifdef MC__FFMLP_CHECK
    assert( vcoef[l].first.size() ); // number of neurons in layer l+1
#endif
    size_t const nneu = vcoef[l].first.size();
    vhid[l].resize( nneu );
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i ){
      vhid[l][i] = (vcoef[l].first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": " << vcoef[l][i].size()-1 << std::endl;
#endif
      for( unsigned j=0; j<(vcoef[l].first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
        std::cout << "layer:" << l << " neuron:" << i << " input:" << j << std::endl;
#endif
        if( std::fabs((vcoef[l].first)[i][1+j]) < options.ZEROTOL ) continue;
        vhid[l][i] += (l? vhid[l-1][j]: x[j]) * (vcoef[l].first)[i][1+j];
      }
      switch( vcoef[l].second ){
        case Options::LINEAR:  default:                                                  break;
        case Options::RELU:    vhid[l][i] = ReLU( vhid[l][i] );                          break;
        case Options::TANH:    vhid[l][i] = Op<U>::tanh( vhid[l][i] );                   break;
        case Options::SIGMOID: vhid[l][i] = (options.SIG2EXP?
                                             1. / ( Op<U>::exp( -vhid[l][i] ) + 1. ):
                                             Op<U>::tanh( vhid[l][i]*0.5 ) * 0.5 + 0.5); break;
      }
    }
  }

  // Propagate through output layers
#ifdef MC__FFMLP_CHECK
  assert( vcoef.back().first.size() ); // number of neurons in layer l+1
#endif
  size_t const nneu = vcoef.back().first.size();
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << nhid << ": " << nneu << std::endl;
#endif
  for( unsigned i=0; i<nneu; ++i ){
    y[i] = (vcoef.back().first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << nhid << ": " << (vcoef.back().first)[i].size()-1 << std::endl;
#endif
    for( unsigned j=0; j<(vcoef.back().first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
      std::cout << "layer:" << nhid << " neuron:" << i << " input:" << j << std::endl;
#endif
      if( std::fabs((vcoef.back().first)[i][1+j]) < options.ZEROTOL ) continue;
      y[i] += (nhid? vhid[nhid-1][j]: x[j]) * (vcoef.back().first)[i][1+j];
    }
    switch( vcoef.back().second ){
      case Options::LINEAR:  default:                                      break;
      case Options::RELU:    y[i] = ReLU( y[i] );                          break;
      case Options::TANH:    y[i] = Op<U>::tanh( y[i] );                   break;
      case Options::SIGMOID: y[i] = (options.SIG2EXP?
                                     Op<U>::inv( Op<U>::exp( -y[i] ) + 1. ):
                                     Op<U>::tanh( y[i]*0.5 ) * 0.5 + 0.5); break;
    }
  }
}
*/
//! @brief C++ class defining neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLPREG is a C++ class for defining neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFMLPREG
////////////////////////////////////////////////////////////////////////
: public FFOp
{

private:
  // pointer to MLP
  MLPREG<T>*   _pMLP;
  // Whether this class owns _pMLP
  bool         _ownMLP;
  
  FFVar** _set
    ( unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei, MLPREG<T>* pMLP )
    {
#ifdef MC__FFMLPREG_CHECK
      assert( nIn == pMLP->nin() && nWei == pMLP->nwei() );
#endif
      if( _ownMLP && _pMLP ) delete _pMLP;
      _ownMLP = true;
      _pMLP = pMLP;
      owndata = false;
      data = _pMLP;

      FFVar** ppRes = insert_external_operation( *this, _pMLP->nout(), nIn, pIn, nWei, pWei );

      _ownMLP = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      _pMLP = dynamic_cast<FFMLPREG<T>*>(pOp)->_pMLP;
#ifdef MC__FFMLPREG_TRACE
      std::cerr << "MLPREG address in DAG: " << _pMLP << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Default Constructor
  FFMLPREG
    ()
    : FFOp( EXTERN ),
      _pMLP( nullptr ),
      _ownMLP( false )
    {}

  // Destructor
  virtual ~FFMLPREG
    ()
    {
      if( _ownMLP ) delete _pMLP;
    }

  // Copy constructor
  FFMLPREG
    ( FFMLPREG<T> const& Op )
    : FFOp( Op )
    {
#ifdef MC__FFMLPREG_TRACE
      std::cout << "FFMLPREG::copy constructor\n";
#endif
      if( !Op._pMLP )
        throw std::runtime_error( "FFMLPREG::copy constructor ** Undefined MLP\n" );

      _ownMLP = Op._ownMLP;      
      if( _ownMLP )
        _pMLP = new MLPREG<T>( *Op._pMLP );
      else
        _pMLP = Op._pMLP;
    }

  // Define operation
  FFVar& operator()
    ( unsigned const iOut, unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei,
      MLPREG<T>* pMLP )
    {
#ifdef MC__FFMLPREG_CHECK
      assert( iOut < pMLP->nout() );
#endif
      return *(_set( nIn, pIn, nWei, pWei, pMLP )[iOut]);
    }

  FFVar** operator()
    ( unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei, MLPREG<T>* pMLP )
    {
      return _set( nIn, pIn, nWei, pWei, pMLP );
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

      throw std::runtime_error( "FFMLPREG::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  template< typename U >
  void eval
    ( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar, unsigned const* mVar )
    const;

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

      throw std::runtime_error( "FFMLPREG::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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
    { std::ostringstream oss; oss << data; return "MLPREG[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

//! @brief C++ class defining gradient of neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGRADMLPREG is a C++ class for defining gradient of neural
//! networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFGRADMLPREG
////////////////////////////////////////////////////////////////////////
: public FFOp
{

private:

  // pointer to MLP
  MLPREG<T>*   _pMLP;
  // Whether this class owns _pMLP
  bool         _ownMLP;
  
  FFVar** _set
    ( unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei, MLPREG<T>* pMLP )
    {
#ifdef MC__FFMLPREG_CHECK
      assert( nIn == pMLP->nin() && nWei == pMLP->nwei() );
#endif
      if( _ownMLP && _pMLP ) delete _pMLP;
      _ownMLP = true;
      _pMLP = pMLP;
      owndata = false;
      data = _pMLP;

      FFVar** ppRes = insert_external_operation( *this, (nIn+nWei)*_pMLP->nout(), nIn, pIn, nWei, pWei );
      
      _ownMLP = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      _pMLP = dynamic_cast<FFGRADMLPREG<T>*>(pOp)->_pMLP;
#ifdef MC__FFMLPREG_TRACE
      std::cerr << "GRADMLPREG address in DAG: " << _pMLP << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Default Constructor
  FFGRADMLPREG
    ()
    : FFOp( EXTERN ),
      _pMLP( nullptr ),
      _ownMLP( false )
    {}

  // Destructor
  virtual ~FFGRADMLPREG
    ()
    {
      if( _ownMLP ) delete _pMLP;
    }

  // Copy constructor
  FFGRADMLPREG
    ( FFGRADMLPREG<T> const& Op )
    : FFOp( Op )
    {
#ifdef MC__FFMLPREG_TRACE
      std::cout << "FFGRADMLPREG::copy constructor\n";
#endif
      if( !Op._pMLP )
        throw std::runtime_error( "FFGRADMLPREG::copy constructor ** Undefined MLP\n" );

      _ownMLP = Op._ownMLP;      
      if( _ownMLP )
        _pMLP = new MLPREG<T>( *Op._pMLP );
      else
        _pMLP = Op._pMLP;
    }

  // Define operation
  FFVar& operator()
    ( unsigned const iOut, unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei,
      MLPREG<T>* pMLP )
    {
#ifdef MC__FFMLPREG_CHECK
      assert( iOut < pMLP->nout() );
#endif
      return *(_set( nIn, pIn, nWei, pWei, pMLP )[iOut]);
    }

  FFVar** operator()
    ( unsigned const nIn, FFVar const* pIn, unsigned const nWei, FFVar const* pWei,
      MLPREG<T>* pMLP )
    {
      return _set( nIn, pIn, nWei, pWei, pMLP );
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

      throw std::runtime_error( "FFGRADMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const;

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const;

  // Backward evaluation overloads
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      throw std::runtime_error( "FFGRADMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << data; return "GRADMLPREG[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};


template< typename T >
template< typename U >
inline void
FFMLPREG<T>::eval
( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  _pMLP->eval( vRes, vVar, vVar+_pMLP->nin() );
}

template< typename T >
inline void
FFMLPREG<T>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::eval: FFDep\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFMLPREG<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::eval: FFVar\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  FFVar** ppRes = insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
}

template< typename T >
inline void
FFMLPREG<T>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = insert_external_operation( *this, nRes, nVar, vVarVal.data() );
  FFGRADMLPREG<T> ResDer;
  FFVar const*const* vResDer = ResDer( _pMLP->nin(), vVarVal.data(), _pMLP->nwei(), vVarVal.data()+_pMLP->nin(), _pMLP );

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
FFMLPREG<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::deriv: FFVar\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  FFGRADMLPREG<T> ResDer;
  FFVar const*const* vResDer = ResDer( _pMLP->nin(), vVar, _pMLP->nwei(), vVar+_pMLP->nin(), _pMLP );
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
      vDer[k][i] = *vResDer[k+nRes*i];
}

template< typename T >
inline void
FFMLPREG<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFMLPREG::eval: SLiftVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == _pMLP->nout() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFGRADMLPREG<T>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADMLPREG::eval: double\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == nVar * _pMLP->nout() );
#endif

  size_t nIndep = _pMLP->nin() + _pMLP->nwei();
  switch( _pMLP->options.AUTODIFF ){
   default:
   case MLPREG<T>::Options::AD::F:{
    std::vector<fadbad::F<double>> vFVar( nIndep );
    for( unsigned i=0; i<nIndep; ++i ){
      vFVar[i] = vVar[i];
      vFVar[i].diff(i,nIndep);
    }
    std::vector<fadbad::F<double>> vFRes( _pMLP->nout() ); 
    _pMLP->eval( vFRes.data(), vFVar.data(), vFVar.data()+_pMLP->nin() );
    for( unsigned k=0; k<_pMLP->nout(); ++k )
      for( unsigned i=0; i<nIndep; ++i )
        vRes[k*nIndep+i] = vFRes[k].d(i);
    break;
   }
   
   case MLPREG<T>::Options::AD::B:{
    std::vector<fadbad::B<double>> vBVar( nIndep );
    for( unsigned i=0; i<nIndep; ++i )
      vBVar[i] = vVar[i];
    std::vector<fadbad::B<double>> vBRes( _pMLP->nout() ); 
    _pMLP->eval( vBRes.data(), vBVar.data(), vBVar.data()+_pMLP->nin() );
    for( unsigned k=0; k<_pMLP->nout(); ++k )
      vBRes[k].diff( k, _pMLP->nout() );
    for( unsigned k=0; k<_pMLP->nout(); ++k )
      for( unsigned i=0; i<nIndep; ++i )
        vRes[k*nIndep+i] = vBVar[i].d(k);
    break;
   }
  }
}

template< typename T >
inline void
FFGRADMLPREG<T>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFGRADMLPREG::eval: FFDep\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == nVar * _pMLP->nout() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFGRADMLPREG<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFGRADMLPREG::eval: FFVar\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == nVar * _pMLP->nout() );
#endif

  FFVar** ppRes = insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
}

template< typename T >
inline void
FFGRADMLPREG<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLPREG_TRACE
  std::cout << "FFGRADMLPREG::eval: SLiftFVar\n";
#endif
#ifdef MC__FFMLPREG_CHECK
  assert( _pMLP && nVar == _pMLP->nin() + _pMLP->nwei() && nRes == nVar * _pMLP->nout() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

} // end namespace mc

#endif
