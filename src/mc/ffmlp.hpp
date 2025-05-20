#ifndef MC__FFMLP_HPP
#define MC__FFMLP_HPP

//#define MC__FFMLP_DEBUG
#define MC__FFMLP_CHECK

#include "mccormick.hpp"
#include "supmodel.hpp"
#include "pwcu.hpp"
#include "pwlu.hpp"
#include "mcfunc.hpp"
#include "ffunc.hpp"
#include "ffdep.hpp"
#include "polimage.hpp"
#include "slift.hpp"

namespace mc
{

// IDEAS:
// - read data from {structure, weights} from MLPREG
// - pass data from pytorch trained ANN

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

  //! @brief Number of inputs
  size_t                                 nin;
  //! @brief Number of outputs
  size_t                                 nout;
  //! @brief Number of hidden layers
  size_t                                 nhid;
  //! @brief MLP data
  std::vector<std::pair<std::vector<std::vector<double>>,int>> data;
  //! @brief Whether MLP data changed since DAG setup
  bool                                   DAGupdt;

  //! @brief Storage for Polyhedral relaxation
  FFGraph*                               DAG;
  FFSubgraph                             DAGOps;
  std::vector<FFVar>                     DAGVar;
  std::vector<FFVar>                     DAGRes;
  PolImg<T>*                             POLEnv;
  std::vector<PolVar<T>>                 POLVar;
  std::vector<PolVar<T>>                 POLRes;
  std::map<PolVar<T> const*,PolVar<T>,lt_PolVar<T>> POLMap;

  //! @brief Storage for Interval bounds
  std::vector<T>                         IVar;
  std::vector<T>                         IRes;

  //! @brief Storage for McCormick relaxation
  std::vector<McCormick<T>>              MCVar;
  std::vector<McCormick<T>>              MCRes;

  //! @brief Storage for interval superposition models
  SupModel<PWCU>*                        PWCSEnv;
  std::vector<SupVar<PWCU>>              PWCSVar;
  std::vector<SupVar<PWCU>>              PWCSRes;
  std::vector<std::vector<PolVar<T>>>    POLPWCSAux;
  std::vector<double>                    DLPWCSAux;
  std::vector<double>                    DUPWCSAux;
  std::vector<McCormick<SupVar<PWCU>>>   MCPWCSVar;
  std::vector<McCormick<SupVar<PWCU>>>   MCPWCSRes;

  //! @brief Storage for affine superposition models
  SupModel<PWLU>*                        PWLSEnv;
  std::vector<SupVar<PWLU>>              PWLSVar;
  std::vector<SupVar<PWLU>>              PWLSRes;
  std::vector<PolVar<T>>                 POLPWLSAux;
  std::vector<double>                    DXPWLSAux;
  std::vector<double>                    DYPWLSAux;
  std::vector<McCormick<SupVar<PWLU>>>   MCPWLSVar;
  std::vector<McCormick<SupVar<PWLU>>>   MCPWLSRes;

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

    //! @brief Relaxation type
    enum RELAXTYPE{
      AUX=0,   //!< Auxiliary variable polyhedral relaxation
      INT,     //!< Interval bounds
      MC,      //!< McCormick relaxation with interval bounds
      PWCS,    //!< Piecewise-constant superposition model
      MCPWCS,  //!< McCormick relaxation with piecewise-constant superposition bounds
      PWLS,    //!< Piecewise-linear superposition model
      MCPWLS   //!< McCormick relaxation with piecewise-linear superposition bounds
    };

    //! @brief Default constructor
    Options():
      RELAX(MC), SUPDIV(64), INIPWL(1), SUPCONT(true), SUPSHADOW(true), CUTSHADOW(false), 
      ZEROTOL(DBL_EPSILON), RELU2ABS(false), SIG2EXP(false), AUTODIFF(F)
      {}

    //! @brief Assignment operator
    Options& operator=
      ( Options const& opt ){
        RELAX     = opt.RELAX;
        SUPDIV    = opt.SUPDIV;
        INIPWL    = opt.INIPWL;
        SUPCONT   = opt.SUPCONT;
        SUPSHADOW = opt.SUPSHADOW;
        CUTSHADOW = opt.CUTSHADOW;
        ZEROTOL   = opt.ZEROTOL;
        RELU2ABS  = opt.RELU2ABS;
        SIG2EXP   = opt.SIG2EXP;
        return *this;
      }

    //! @brief Enumeration type for AD strategy
    enum AD{
      F=0,	//!< Forward differentiation
      B		//!< Backward differentiation
    };

    //! @brief Type of relaxation
    RELAXTYPE RELAX;
    //! @brief Partition size in superposition model
    size_t SUPDIV;
    //! @brief Initial partition size in piecewise-linear superposition model
    size_t INIPWL;
    //! @brief Whether to construct continuous or binary relaxation of superposition model
    bool     SUPCONT;
    //! @brief Whether to propagate shadow estimators in superposition model
    bool     SUPSHADOW;
    //! @brief Whether to append cuts from shadow estimators in superposition model relaxation
    bool     CUTSHADOW;
    //! @brief Threshold for zero coefficient in neural network
    double   ZEROTOL;
    //! @brief Whether to convert ReLU to abs (true) or max (false) function
    bool     RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false)
    bool     SIG2EXP;
    //! @brief Whether to apply forward or reverse automatic differentiation
    int      AUTODIFF;
  } options;

  //! @brief Default constructor
  MLP
  ()
  : nin(0), nout(0), DAGupdt(false),
    DAG(nullptr), POLEnv(nullptr), PWCSEnv(nullptr), PWLSEnv(nullptr)
  {}

  ~MLP() 
  {
    delete POLEnv;
    delete DAG;
    delete PWCSEnv;
    delete PWLSEnv;
  }

  //! @brief Reset MLP data
  void reset_data
    ();

  //! @brief Set MLP data
  bool set_data
    ( std::vector<std::pair<std::vector<std::vector<double>>,int>> const& data );

  //! @brief Append multi-neuron layer to MLP data: outer vector size is #neurons in hidden layer; inner vector has scalar weights multiplying each neuron in previous (input or hidden) layer and a bias term and size 1+#neuron in previous layer 
  bool append_data
    ( std::vector<std::vector<double>> const& data, int const activ=Options::LINEAR, bool const reset=false );

  //! @brief Append single-neuron layer to MLP data: vector has scalar weights multiplying each neuron in previous (input or hidden) layer and a bias term and size 1+#neuron in previous layer 
  bool append_data
    ( std::vector<double> const& data, int const activ=Options::LINEAR, bool const reset=false );

  //! @brief Resize relaxation containers
  void resize_relax
    ();

  //! @brief Evaluate neural network
  template <typename U>
  void evaluate
    ( U* y, U const* x )
    const
    { static thread_local std::vector<std::vector<U>> Uhid;
      evaluate( y, x, Uhid ); }
  void evaluate
    ( double* y, double const* x )
    { evaluate( y, x, Dhid ); }
  void evaluate
    ( fadbad::F<double>* y, fadbad::F<double> const* x )
    { evaluate( y, x, FDhid ); }
  void evaluate
    ( T* y, T const* x )
    { evaluate( y, x, Ihid ); }
  void evaluate
    ( McCormick<T>* y, McCormick<T> const* x )
    { evaluate( y, x, MCIhid ); }
  void evaluate
    ( SupVar<PWCU>* y, SupVar<PWCU> const* x )
    { evaluate( y, x, PWCShid ); }
  void evaluate
    ( McCormick<SupVar<PWCU>>* y, McCormick<SupVar<PWCU>> const* x )
    { evaluate( y, x, MCPWCShid ); }
  void evaluate
    ( SupVar<PWLU>* y, SupVar<PWLU> const* x )
    { evaluate( y, x, PWLShid ); }
  void evaluate
    ( McCormick<SupVar<PWLU>>* y, McCormick<SupVar<PWLU>> const* x )
    { evaluate( y, x, MCPWLShid ); }
  void evaluate
    ( PolVar<T>* y, PolVar<T> const* x )
    { evaluate( y, x, POLhid ); }

  //! @brief Propagate polyhedral image through neural network
  void propagate_relax
    ( PolImg<T>* img, FFVar** pRes, PolVar<T>* vRes, PolVar<T> const* vVar );

  //! @brief Append polyhedral image cuts for neural network
  void backpropagate_relax
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const* vRes, PolVar<T>* vVar );

private:

  template <typename U>
  void evaluate
    ( U* y, U const* x, std::vector<std::vector<U>>& vhid )
    const;

  void append_PWLScuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const& vRes, PolVar<T>* vVar,
      std::vector<PWLU> const& est, std::set<unsigned int> const& dep, bool const under );

  template <typename U>
  U ReLU
    ( U const& x )
    const
    {
      return options.RELU2ABS?
             ( x + Op<U>::fabs(x) ) * 0.5:
             Op<U>::max( x, 0. );
    }
/*
  SupVar<PWCU> ReLU
    ( SupVar<PWCU> const& x )
    const
    {
      return relu( x );
    }

  template <typename U>
  SupVar<PWLU> ReLU
    ( SupVar<PWCU> const& x )
    const
    {
      return relu( x );
    }
*/
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

  //! @brief Intermediate storage for MLP evaluation
  static thread_local std::vector<std::vector<double>>                   Dhid;
  static thread_local std::vector<std::vector<fadbad::F<double>>>        FDhid;
  static thread_local std::vector<std::vector<T>>                        Ihid;
  static thread_local std::vector<std::vector<McCormick<T>>>             MCIhid;
  static thread_local std::vector<std::vector<SupVar<PWCU>>>             PWCShid;
  static thread_local std::vector<std::vector<McCormick<SupVar<PWCU>>>>  MCPWCShid;
  static thread_local std::vector<std::vector<SupVar<PWLU>>>             PWLShid;
  static thread_local std::vector<std::vector<McCormick<SupVar<PWLU>>>>  MCPWLShid;
  static thread_local std::vector<std::vector<PolVar<T>>>                POLhid;
};

template <typename T> inline thread_local std::vector<std::vector<double>>                  MLP<T>::Dhid      = std::vector<std::vector<double>>();
template <typename T> inline thread_local std::vector<std::vector<fadbad::F<double>>>       MLP<T>::FDhid     = std::vector<std::vector<fadbad::F<double>>>();
template <typename T> inline thread_local std::vector<std::vector<T>>                       MLP<T>::Ihid      = std::vector<std::vector<T>>();
template <typename T> inline thread_local std::vector<std::vector<McCormick<T>>>            MLP<T>::MCIhid    = std::vector<std::vector<McCormick<T>>>();
template <typename T> inline thread_local std::vector<std::vector<SupVar<PWCU>>>            MLP<T>::PWCShid   = std::vector<std::vector<SupVar<PWCU>>>();
template <typename T> inline thread_local std::vector<std::vector<McCormick<SupVar<PWCU>>>> MLP<T>::MCPWCShid = std::vector<std::vector<McCormick<SupVar<PWCU>>>>();
template <typename T> inline thread_local std::vector<std::vector<SupVar<PWLU>>>            MLP<T>::PWLShid   = std::vector<std::vector<SupVar<PWLU>>>();
template <typename T> inline thread_local std::vector<std::vector<McCormick<SupVar<PWLU>>>> MLP<T>::MCPWLShid = std::vector<std::vector<McCormick<SupVar<PWLU>>>>();
template <typename T> inline thread_local std::vector<std::vector<PolVar<T>>>               MLP<T>::POLhid    = std::vector<std::vector<PolVar<T>>>();

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
  data = mlp;
  nin  = data.front().first.front().size() - 1;
  nout = data.back().first.size();
  nhid = data.size() - 1;
  DAGupdt = true;

  return true;
}

template< typename T >
inline void
MLP<T>::reset_data
()
{
  // Reset data
  nin  = nout = nhid = 0;
  data.clear();
  DAGupdt = true;
}

template< typename T >
inline bool
MLP<T>::append_data
( std::vector<std::vector<double>> const& layer, int const activ, bool const reset )
{
  if( reset ){
    reset_data();
    DAGupdt = true;
  }
  if( !layer.size()
   || layer[0].size() <= 1
   || ( !data.empty() && data.back().first.size() != layer[0].size() - 1 ) )
    return false;

  // Set data
  data.push_back( std::make_pair( layer, activ ) );
  nin  = data.front().first.front().size() - 1;
  nout = data.back().first.size();
  nhid = data.size() - 1;
  DAGupdt = true;

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
template< typename U >
inline void
MLP<T>::evaluate
( U* y, U const* x, std::vector<std::vector<U>>& vhid )
const
{
  // Propagate through hidden layers
  vhid.resize( nhid );
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << nhid << std::endl;
#endif
  for( unsigned l=0; l<nhid; ++l ){
#ifdef MC__FFMLP_CHECK
    assert( data[l].first.size() ); // number of neurons in layer l+1
#endif
    size_t const nneu = data[l].first.size();
    vhid[l].resize( nneu );
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i ){
      vhid[l][i] = (data[l].first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": " << (data[l].first)[i].size()-1 << std::endl;
#endif
      for( unsigned j=0; j<(data[l].first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
        std::cout << "layer:" << l << " neuron:" << i << " input:" << j << std::endl;
#endif
        if( std::fabs((data[l].first)[i][1+j]) < options.ZEROTOL ) continue;
        vhid[l][i] += (l? vhid[l-1][j]: x[j]) * (data[l].first)[i][1+j];
      }
      switch( data[l].second ){//options.ACTIV ){
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
  assert( data.back().first.size() ); // number of neurons in layer l+1
#endif
  size_t const nneu = data.back().first.size();
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << nhid << ": " << nneu << std::endl;
#endif
  for( unsigned i=0; i<nneu; ++i ){
    y[i] = (data.back().first)[i][0]; // bias term
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << nhid << ": " << (data.back().first)[i].size()-1 << std::endl;
#endif
    for( unsigned j=0; j<(data.back().first)[i].size()-1; ++j ){
#ifdef MC__FFMLP_DEBUG
      std::cout << "layer:" << nhid << " neuron:" << i << " input:" << j << std::endl;
#endif
      if( std::fabs((data.back().first)[i][1+j]) < options.ZEROTOL ) continue;
      y[i] += (nhid? vhid[nhid-1][j]: x[j]) * (data.back().first)[i][1+j];
    }
    switch( data.back().second ){
      case Options::LINEAR:  default:                                      break;
      case Options::RELU:    y[i] = ReLU( y[i] );                          break;
      case Options::TANH:    y[i] = Op<U>::tanh( y[i] );                   break;
      case Options::SIGMOID: y[i] = (options.SIG2EXP?
                                     Op<U>::inv( Op<U>::exp( -y[i] ) + 1. ):
                                     Op<U>::tanh( y[i]*0.5 ) * 0.5 + 0.5); break;
    }
  }
}

template< typename T >
inline void
MLP<T>::resize_relax
()
{
  // Set relaxation environment and containers
  switch( options.RELAX ){
    case Options::AUX:
    default:
      if( !DAG || DAGupdt ){
        delete DAG;
        DAG = new FFGraph;
        DAGVar.resize( nin );
        DAGRes.resize( nout );
        for( unsigned i=0; i<nin; ++i )
          DAGVar[i].set( DAG );
        evaluate( DAGRes.data(), DAGVar.data() );
        DAGOps = DAG->subgraph( nout, DAGRes.data() );
#ifdef MC__FFMLP_DEBUG
        DAG->output( DAGOps, " MLP", std::cerr );
        std::ofstream oFile( "MLP.dot", std::ios_base::out );
        DAG->dot_script( nout, DAGRes.data(), oFile );
        oFile.close();
#endif
        DAGupdt = false;
      }
      if( !POLEnv ) POLEnv = new PolImg<T>;
      POLVar.resize( nin );
      POLRes.resize( nout );
      break;

    case Options::INT:
      IVar.resize( nin );
      IRes.resize( nout );
      break;

    case Options::MC:
      MCVar.resize( nin );
      MCRes.resize( nout );
      break;

    case Options::MCPWCS:
      MCPWCSVar.resize( nin );
      MCPWCSRes.resize( nout );
      // no break

    case Options::PWCS:
      if( !PWCSEnv || PWCSEnv->nvar() != nin ){
        delete PWCSEnv;
        PWCSEnv = new SupModel<PWCU>( nin );
      }
      PWCSEnv->options.USE_SHADOW = options.SUPSHADOW;
      PWCSVar.resize( nin );
      PWCSRes.resize( nout );
      POLPWCSAux.resize( nin );
      DLPWCSAux.resize( options.SUPDIV );
      DUPWCSAux.resize( options.SUPDIV );
      break;

    case Options::MCPWLS:
      MCPWLSVar.resize( nin );
      MCPWLSRes.resize( nout );
      // no break

    case Options::PWLS:
      if( !PWLSEnv || PWLSEnv->nvar() != nin ){
        delete PWLSEnv;
        PWLSEnv = new SupModel<PWLU>( nin );
      }
      PWLSEnv->options.USE_SHADOW = options.SUPSHADOW;
      PWLSVar.resize( nin );
      PWLSRes.resize( nout );
      POLPWLSAux.resize( nin );
      break;
  }
}

template< typename T >
inline void
MLP<T>::propagate_relax
( PolImg<T>* img, FFVar ** pRes, PolVar<T>* vRes, PolVar<T> const* vVar )
{
  switch( options.RELAX ){
    // Polyhedral relaxation with auxiliary variables
    case Options::AUX:
    default:
      POLEnv->options = img->options;
      POLEnv->reset();
      POLMap.clear();
      for( unsigned i=0; i<nin; ++i ){
        POLVar[i].set( POLEnv, DAGVar[i], vVar[i].range() );
        POLMap[&POLVar[i]] = vVar[i];
      }
#ifdef MC__FFMLP_DEBUG
      { std::vector<PolVar<T>> POLwk;
        DAG->eval( DAGOps, POLwk, nout, DAGRes.data(), POLRes.data(), nin, DAGVar.data(), POLVar.data() );
        for( auto polv : POLwk ) std::cout << polv.name() << ": " << polv.range() << std::endl; }
#else
      DAG->eval( DAGOps, nout, DAGRes.data(), POLRes.data(), nin, DAGVar.data(), POLVar.data() );
#endif
      // COULD DO CONSTRAINT PROPAGATION BASED ON vRes.range()
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], POLRes[j].range() );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "vRes[" << j << "] in " << vRes[j] << std::endl;
#endif
        POLMap[&POLRes[j]] = vRes[j];
      }
      break;

    // Interval bounds
    case MLP<T>::Options::RELAXTYPE::INT:
      for( unsigned i=0; i<nin; ++i )
        IVar[i] = vVar[i].range();
      evaluate( IRes.data(), IVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], IRes[j] );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "vRes[" << j << "] in " << vRes[j] << std::endl;
#endif
      }
      break;
      
    // McCormick relaxations at mid-point with subgradient in each direction
    case MLP<T>::Options::RELAXTYPE::MC:
      for( unsigned i=0; i<nin; ++i )
      MCVar[i] = McCormick<T>( vVar[i].range(), Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
      evaluate( MCRes.data(), MCVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], MCRes[j].I() );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "vRes[" << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      break;

    // Piecewise-constant superposition models
    case MLP<T>::Options::RELAXTYPE::PWCS:
      for( unsigned i=0; i<nin; ++i )
        PWCSVar[i].set( *PWCSEnv, i, vVar[i].range(), options.SUPDIV );
      evaluate( PWCSRes.data(), PWCSVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], T(PWCSRes[j].l(),PWCSRes[j].u()) );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWCSRes[" << j << "] in " << PWCSRes[j]      << std::endl;
        std::cerr << "vRes["    << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      break;
 
    // McCormick relaxation with piecewise-constant superposition bounds and subgradient calculated at mid-point
    case MLP<T>::Options::RELAXTYPE::MCPWCS:
      for( unsigned i=0; i<nin; ++i )
        MCPWCSVar[i] = McCormick<SupVar<PWCU>>( SupVar<PWCU>( *PWCSEnv, i, vVar[i].range(), options.SUPDIV ),
                                                Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
      evaluate( MCPWCSRes.data(), MCPWCSVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], T(MCPWCSRes[j].I().l(),MCPWCSRes[j].I().u()) );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWCSRes["   << j << "] in " << MCPWCSRes[j].I() << std::endl;
        std::cerr << "MCPWCSRes[" << j << "] in " << MCPWCSRes[j]     << std::endl;
        std::cerr << "vRes["      << j << "] in " << vRes[j].range()  << std::endl;
#endif
      }
      break;

    // Piecewise-linear superposition models
    case MLP<T>::Options::RELAXTYPE::PWLS:
      //UnivarPWLE<double>::nbpsMax = options.ASMBPS;       
      for( unsigned i=0; i<nin; ++i )
        PWLSVar[i].set( *PWLSEnv, i, vVar[i].range(), options.INIPWL );
      evaluate( PWLSRes.data(), PWLSVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], T(PWLSRes[j].l(),PWLSRes[j].u()) );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWLSRes[" << j << "] in " << PWLSRes[j]      << std::endl;
        std::cerr << "vRes["    << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      break;

    // McCormick relaxation with piecewise-linear superposition bounds and subgradient calculated at mid-point
    case MLP<T>::Options::RELAXTYPE::MCPWLS:
      for( unsigned i=0; i<nin; ++i )
        MCPWLSVar[i] = McCormick<SupVar<PWLU>>( SupVar<PWLU>( *PWLSEnv, i, vVar[i].range(), options.INIPWL ),
                                                Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
      evaluate( MCPWLSRes.data(), MCPWLSVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], T(MCPWLSRes[j].I().l(),MCPWLSRes[j].I().u()) );
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWLSRes["   << j << "] in " << MCPWLSRes[j].I() << std::endl;
        std::cerr << "MCPWLSRes[" << j << "] in " << MCPWLSRes[j]     << std::endl;
        std::cerr << "vRes["      << j << "] in " << vRes[j].range()  << std::endl;
#endif
      }
      break;
  }
}

template< typename T >
inline void
MLP<T>::backpropagate_relax
( PolImg<T>* img, FFOp* pOp, PolVar<T> const* vRes, PolVar<T>* vVar )
{
  switch( options.RELAX ){
    // Polyhedral relaxation with auxiliary variables
    case Options::AUX:
    default:
      POLEnv->generate_cuts( nout, POLRes.data() );
#ifdef MC__FFMLP_DEBUG
      std::cerr << "POLEnv:" << *POLEnv << std::endl;
#endif
      img->insert_cuts( POLEnv, POLMap );
      break;

    // Interval bounds
    case MLP<T>::Options::RELAXTYPE::INT:
      break;
      
    // McCormick relaxations at mid-point with subgradient in each direction
    case MLP<T>::Options::RELAXTYPE::MC:
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__FFMLP_DEBUG
        std::cerr << "MCRes[" << j << "] in " << MCRes[j] << std::endl;
#endif
        double rhs1 = -MCRes[j].cv(),
               rhs2 = -MCRes[j].cc();
        for( unsigned i=0; i<nin; ++i ){
          rhs1 += MCRes[j].cvsub(i)*MCVar[i].cv();
          rhs2 += MCRes[j].ccsub(i)*MCVar[i].cc();
        }
        img->add_cut( pOp, PolCut<T>::LE, rhs1, nin, vVar, MCRes[j].cvsub(), vRes[j], -1. );
        img->add_cut( pOp, PolCut<T>::GE, rhs2, nin, vVar, MCRes[j].ccsub(), vRes[j], -1. );
      }
      break;

    // Piecewise-constant superposition models
    case MLP<T>::Options::RELAXTYPE::PWCS:
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        auto& PWCSVar_i = PWCSVar[i];
#ifdef MC__FFMLP_CHECK
        assert( PWCSVar_i.uest()[i].yL().size() == options.SUPDIV 
             && PWCSVar_i.uest()[i].yL().size() == PWCSVar_i.oest()[i].yU().size() );
#endif
        POLPWCSAux[i].resize( options.SUPDIV );
        for( unsigned k=0; k<options.SUPDIV; ++k )
          POLPWCSAux[i][k].set( img, Op<T>::zeroone(), options.SUPCONT );

        // sum of auxiliaries equal to 1
        img->add_cut( pOp, PolCut<T>::EQ, 1., options.SUPDIV, POLPWCSAux[i].data(), 1. );

        // link auxiliaries to independent variables
        img->add_cut( pOp, PolCut<T>::LE, 0., options.SUPDIV, POLPWCSAux[i].data(), PWCSVar_i.uest()[i].yL().data(), vVar[i], -1. );
        img->add_cut( pOp, PolCut<T>::GE, 0., options.SUPDIV, POLPWCSAux[i].data(), PWCSVar_i.oest()[i].yU().data(), vVar[i], -1. );
      }

      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
        auto& PWCSRes_j = PWCSRes[j];
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWCSRes["   << j << "] in " << PWCSRes_j    << std::endl;
#endif
        // constant superposition relaxation
        if( PWCSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWCSRes_j.cst(), vRes[j], 1. );
          continue;
        }

        // polyhedral cuts for superposition relaxation
        auto cutF1 = *img->add_cut( pOp, PolCut<T>::LE, 0., vRes[j], -1. );
        auto cutF2 = *img->add_cut( pOp, PolCut<T>::GE, 0., vRes[j], -1. );
        for( auto const& i : PWCSRes_j.sdep() ){
#ifdef MC__FFMLP_CHECK
          assert( PWCSRes_j.uest()[i].yL().size() == options.SUPDIV 
               && PWCSRes_j.uest()[i].yL().size() == PWCSRes_j.oest()[i].yU().size() );
#endif
          cutF1->append( options.SUPDIV, POLPWCSAux[i].data(), PWCSRes_j.uest()[i].yL().data() );
          cutF2->append( options.SUPDIV, POLPWCSAux[i].data(), PWCSRes_j.oest()[i].yU().data() );
        }
      }
      break;

    // McCormick relaxation with piecewise-constant superposition bounds and subgradient calculated at mid-point
    case MLP<T>::Options::RELAXTYPE::MCPWCS:
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        auto& PWCSVar_i = MCPWCSVar[i].I();
#ifdef MC__FFMLP_CHECK
        assert( PWCSVar_i.uest()[i].yL().size() == options.SUPDIV 
             && PWCSVar_i.uest()[i].yL().size() == PWCSVar_i.oest()[i].yU().size() );
#endif
        POLPWCSAux[i].resize( options.SUPDIV );
        for( unsigned k=0; k<options.SUPDIV; ++k )
          POLPWCSAux[i][k].set( img, Op<T>::zeroone(), options.SUPCONT );

        // sum of auxiliaries equal to 1
        img->add_cut( pOp, PolCut<T>::EQ, 1., options.SUPDIV, POLPWCSAux[i].data(), 1. );

        // link auxiliaries to independent variables
        img->add_cut( pOp, PolCut<T>::LE, 0., options.SUPDIV, POLPWCSAux[i].data(), PWCSVar_i.uest()[i].yL().data(), vVar[i], -1. );
        img->add_cut( pOp, PolCut<T>::GE, 0., options.SUPDIV, POLPWCSAux[i].data(), PWCSVar_i.oest()[i].yU().data(), vVar[i], -1. );
      }

      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
        auto& PWCSRes_j = MCPWCSRes[j].I();
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWCSRes["   << j << "] in " << PWCSRes_j    << std::endl;
#endif
        // constant superposition relaxation
        if( PWCSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWCSRes_j.cst(), vRes[j], 1. );
          continue;
        }

        // polyhedral cuts for superposition relaxation
        auto cutF1 = *img->add_cut( pOp, PolCut<T>::LE, 0., vRes[j], -1. );
        auto cutF2 = *img->add_cut( pOp, PolCut<T>::GE, 0., vRes[j], -1. );
        for( auto const& i : PWCSRes_j.sdep() ){
#ifdef MC__FFMLP_CHECK
          assert( PWCSRes_j.uest()[i].yL().size() == options.SUPDIV 
               && PWCSRes_j.uest()[i].yL().size() == PWCSRes_j.oest()[i].yU().size() );
#endif
          cutF1->append( options.SUPDIV, POLPWCSAux[i].data(), PWCSRes_j.uest()[i].yL().data() );
          cutF2->append( options.SUPDIV, POLPWCSAux[i].data(), PWCSRes_j.oest()[i].yU().data() );
        }

        // polyhedral cut generation for MC  
#ifdef MC__FFMLP_DEBUG
        std::cerr << "MCPWCSRes[" << j << "] in " << MCPWCSRes[j] << std::endl;
#endif
        double rhs1 = -MCPWCSRes[j].cv(),  
               rhs2 = -MCPWCSRes[j].cc();
        for( unsigned i=0; i<nin; ++i ){ 
          rhs1 += MCPWCSRes[j].cvsub(i) * MCPWCSVar[i].cv();
          rhs2 += MCPWCSRes[j].ccsub(i) * MCPWCSVar[i].cc();
        }
        img->add_cut( pOp, PolCut<T>::LE, rhs1, nin, vVar, MCPWCSRes[j].cvsub(), vRes[j], -1. );
        img->add_cut( pOp, PolCut<T>::GE, rhs2, nin, vVar, MCPWCSRes[j].ccsub(), vRes[j], -1. );
      }
      break;
      
    // Piecewise-linear superposition models
    case MLP<T>::Options::RELAXTYPE::PWLS:
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
        auto& PWLSRes_j = PWLSRes[j];
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWLSRes["   << j << "] in " << PWLSRes_j    << std::endl;
#endif
        // constant superposition relaxation
        if( PWLSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWLSRes_j.cst(), vRes[j], 1. );
          continue;
        }

        // polyhedral cuts for superposition relaxation
        append_PWLScuts( img, pOp, vRes[j], vVar, PWLSRes_j.uest(), PWLSRes_j.sdep(), true  );          
        append_PWLScuts( img, pOp, vRes[j], vVar, PWLSRes_j.oest(), PWLSRes_j.sdep(), false );          
      }
      break;

    // McCormick relaxation with piecewise-linear superposition bounds and subgradient calculated at mid-point
    case MLP<T>::Options::RELAXTYPE::MCPWLS:
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
        auto& PWLSRes_j = MCPWLSRes[j].I();
#ifdef MC__FFMLP_DEBUG
        std::cerr << "PWLSRes["   << j << "] in " << PWLSRes_j    << std::endl;
#endif
        // constant superposition relaxation
        if( PWLSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWLSRes_j.cst(), vRes[j], 1. );
          continue;
        }

        // polyhedral cuts for superposition relaxation
        append_PWLScuts( img, pOp, vRes[j], vVar, PWLSRes_j.uest(), PWLSRes_j.sdep(), true  );          
        append_PWLScuts( img, pOp, vRes[j], vVar, PWLSRes_j.oest(), PWLSRes_j.sdep(), false );          

        // polyhedral cut generation for MC  
#ifdef MC__FFMLP_DEBUG
        std::cerr << "MCPWLSRes[" << j << "] in " << MCPWLSRes[j] << std::endl;
#endif
        double rhs1 = -MCPWLSRes[j].cv(),  
               rhs2 = -MCPWLSRes[j].cc();
        for( unsigned i=0; i<nin; ++i ){ 
          rhs1 += MCPWLSRes[j].cvsub(i) * MCPWLSVar[i].cv();
          rhs2 += MCPWLSRes[j].ccsub(i) * MCPWLSVar[i].cc();
        }
        img->add_cut( pOp, PolCut<T>::LE, rhs1, nin, vVar, MCPWLSRes[j].cvsub(), vRes[j], -1. );
        img->add_cut( pOp, PolCut<T>::GE, rhs2, nin, vVar, MCPWLSRes[j].ccsub(), vRes[j], -1. );
      }
      break;
  } // end switch

#ifdef MC__FFMLP_DEBUG
  std::cerr << *img;
  {int dum; std::cout << "PAUSED, ENTER 1"; std::cin >> dum;}
#endif
}


template< typename T >
inline void
MLP<T>::append_PWLScuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& vRes, PolVar<T>* vVar,
  std::vector<PWLU> const& est, std::set<unsigned int> const& dep, bool const under )
{
  for( auto const& i : dep ){
    POLPWLSAux[i].set( img, T(est.at(i).xL(),est.at(i).xU()), true );

    // Generate breakpoints from segment-based representation in mc::PWLU
    auto const& est_i = est.at(i);
    size_t const NK = est_i.dx().size()+1;
    DXPWLSAux.resize( NK ); DXPWLSAux[0] = est_i.xL();
    DYPWLSAux.resize( NK ); DYPWLSAux[0] = est_i.yL();
    auto idx = est_i.dx().cbegin(), idy = est_i.dy().cbegin();
    for( unsigned k=1; k<NK; ++k, ++idx, ++idy ){
      DXPWLSAux[k] = DXPWLSAux[k-1] + *idx;
      DYPWLSAux[k] = DYPWLSAux[k-1] + *idx * *idy;
    }

    img->add_semilinear_cuts( pOp, NK, vVar[i], DXPWLSAux.data(), POLPWLSAux[i], DYPWLSAux.data(), mc::PolCut<T>::EQ );
  }
     
  img->add_cut( pOp, (under? PolCut<T>::LE: PolCut<T>::GE), 0., dep, POLPWLSAux.data(), 1., vRes, -1. );
} 

//! @brief C++ class defining neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLP is a C++ class for defining neural networks as external
//! DAG operations in MC++. The template parameter specifies the type
//! for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFMLP
////////////////////////////////////////////////////////////////////////
: public FFOp
{
public:

  //! @brief Constructors
  FFMLP
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      owndata = false;
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin && idep < pMLP->nout );
#endif
      return *(insert_external_operation( *this, pMLP->nout, nVar, pVar )[idep]);
    }

  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      owndata = false;
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin );
#endif
      return insert_external_operation( *this, pMLP->nout, nVar, pVar );
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
      else if( idU == typeid( McCormick<T> ) )
        return eval( nRes, static_cast<McCormick<T>*>(vRes), nVar, static_cast<McCormick<T> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWCU> ) )
        return eval( nRes, static_cast<SupVar<PWCU>*>(vRes), nVar, static_cast<SupVar<PWCU> const*>(vVar), mVar );
      else if( idU == typeid( SupVar<PWLU> ) )
        return eval( nRes, static_cast<SupVar<PWLU>*>(vRes), nVar, static_cast<SupVar<PWLU> const*>(vVar), mVar );
      else if( idU == typeid( PolVar<T> ) )
        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );

      throw std::runtime_error( "FFMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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

  void eval
    ( unsigned const nRes, PolVar<T>* vRes, unsigned const nVar, PolVar<T> const* vVar, unsigned const* mVar )
    const;


  // Backward evaluation overloads
  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      if( idU == typeid( PolVar<T> ) )
        return reval( nRes, static_cast<PolVar<T> const*>(vRes), nVar, static_cast<PolVar<T>*>(vVar) );

      throw std::runtime_error( "FFMLP::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  bool reval
    ( unsigned const nRes, PolVar<T> const* vRes, unsigned const nVar, PolVar<T>* vVar )
    const;

  // Derivatives
  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;
    
  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << data; return "MLP[" + oss.str() + "]"; }

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
: public FFOp
{
public:

  //! @brief Constructors
  FFGradMLP
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      owndata = false;
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin && idep < pMLP->nout );
#endif
      return *(insert_external_operation( *this, nVar*pMLP->nout, nVar, pVar )[idep]);
    }

  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      owndata = false;
#ifdef MC__FFMLP_CHECK
      assert( nVar == pMLP->nin );
#endif
      return insert_external_operation( *this, nVar*pMLP->nout, nVar, pVar );
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

      throw std::runtime_error( "FFGradMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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
      throw std::runtime_error( "FFGradMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  // Properties
  std::string name
    ()
    const
    { std::ostringstream oss; oss << data; return "GRADMLP[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};


template< typename T >
template< typename U >
inline void
FFMLP<T>::eval
( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__MLPOP_TRACE
  std::cout << "FFMLP::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  pMLP->evaluate( vRes, vVar );
}

template< typename T >
inline void
FFMLP<T>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFDep\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFMLP<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  FFVar** pRes = operator()( nVar, vVar, pMLP );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template< typename T >
inline void
FFMLP<T>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: fadbad::F<FFVar>\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = operator()( nVar, vVarVal.data(), pMLP );
  FFGradMLP<T> ResDer;
  FFVar const*const* vResDer = ResDer( nVar, vVarVal.data(), pMLP );

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
  std::cout << "FFMLP::deriv: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  FFGradMLP<T> ResDer;
  FFVar const*const* vResDer = ResDer( nVar, vVar, pMLP );
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
      vDer[k][i] = *vResDer[k+nRes*i];
}

template< typename T >
inline void
FFMLP<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: SLiftVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFMLP<T>::eval
( unsigned const nRes, PolVar<T>* vRes, unsigned const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: PolVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__FFMLP_CHECK
  assert( img && dag );
#endif
  FFVar** pRes = dag->curOp()->varout.data(); // ACCOUNT FOR MULTIPLE OUTPUTS
#ifdef MC__FFMLP_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  pMLP->resize_relax();
  pMLP->propagate_relax( img, pRes, vRes, vVar );
}

template< typename T >
inline bool
FFMLP<T>::reval
( unsigned const nRes, PolVar<T> const* vRes, unsigned const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: PolVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__FFMLP_CHECK
  assert( img && pop );
#endif

  pMLP->backpropagate_relax( img, pop, vRes, vVar );
  return true;
}

template< typename T >
inline void
FFGradMLP<T>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGradMLP::eval: double\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  switch( pMLP->options.AUTODIFF ){
   default:
   case MLP<T>::Options::AD::F:{
    std::vector<fadbad::F<double>> vFVar( pMLP->nin );
    for( unsigned i=0; i<pMLP->nin; ++i ){
      vFVar[i] = vVar[i];
      vFVar[i].diff(i,pMLP->nin);
    }
    std::vector<fadbad::F<double>> vFRes( pMLP->nout ); 
    pMLP->evaluate( vFRes.data(), vFVar.data() );
    for( unsigned k=0; k<pMLP->nout; ++k )
      for( unsigned i=0; i<pMLP->nin; ++i )
        vRes[k*pMLP->nin+i] = vFRes[k].d(i);
    break;
   }
   
   case MLP<T>::Options::AD::B:{
    std::vector<fadbad::B<double>> vBVar( pMLP->nin );
    for( unsigned i=0; i<pMLP->nin; ++i )
      vBVar[i] = vVar[i];
    std::vector<fadbad::B<double>> vBRes( pMLP->nout ); 
    pMLP->evaluate( vBRes.data(), vBVar.data() );
    for( unsigned k=0; k<pMLP->nout; ++k )
      vBRes[k].diff( k, pMLP->nout );
    for( unsigned k=0; k<pMLP->nout; ++k )
      for( unsigned i=0; i<pMLP->nin; ++i )
        vRes[k*pMLP->nin+i] = vBVar[i].d(k);
    break;
   }
  }
}

template< typename T >
inline void
FFGradMLP<T>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFDep\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFGradMLP<T>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGradMLP::eval: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  FFVar** pRes = operator()( nVar, vVar, pMLP );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template< typename T >
inline void
FFGradMLP<T>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGradMLP::eval: SLiftFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__FFMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

} // end namespace mc

#endif
