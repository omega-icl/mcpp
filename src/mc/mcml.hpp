#ifndef MC__MCML_HPP
#define MC__MCML_HPP

//#define MC__MCMLP_DEBUG
#define MC__MCMLP_CHECK

#include "mccormick.hpp"
#include "ismodel.hpp"
#include "asmodel.hpp"
#include "mcfunc.hpp"
#include "ffunc.hpp"
#include "polimage.hpp"
#include "slift.hpp"

namespace mc
{

//! @brief C++ class for evaluation and relaxation of multilayer perceptrons
////////////////////////////////////////////////////////////////////////
//! mc::MLP is a C++ class for evaluation and relaxation of multilayer
//! perceptrons (MLP) that leverages expression trees and arithmetics
//! available through MC++
////////////////////////////////////////////////////////////////////////
template <typename T> 
class MLP
////////////////////////////////////////////////////////////////////////
{
public:

  //! @brief Number of inputs
  size_t                               nin;
  //! @brief Number of outputs
  size_t                               nout;
  //! @brief Number of hidden layers
  size_t                               nhid;
  //! @brief MLP data
  std::vector<std::pair<std::vector<std::vector<double>>,int>> data;

  //! @brief Storage for Polhedral relaxation
  std::pair<FFGraph*,void*>            DAG;
  FFSubgraph                           DAGOps;
  std::vector<FFVar>                   DAGVar;
  std::vector<FFVar>                   DAGRes;
  PolImg<T>*                           POLEnv;
  std::vector<PolVar<T>>               POLVar;
  std::vector<PolVar<T>>               POLRes;
  std::map<PolVar<T> const*,PolVar<T>,lt_PolVar<T>> POLMap;

  //! @brief Storage for Interval bounds
  std::vector<T>                       IVar;
  std::vector<T>                       IRes;

  //! @brief Storage for McCormick relaxation
  std::vector<McCormick<T>>            MCVar;
  std::vector<McCormick<T>>            MCRes;

  //! @brief Storage for interval superposition models
  ISModel<T>*                          ISMEnv;
  std::vector<ISVar<T>>                ISMVar;
  std::vector<ISVar<T>>                ISMRes;
  std::vector<std::vector<PolVar<T>>>  POLISMAux;
  std::vector<double>                  DLISMAux;
  std::vector<double>                  DUISMAux;
  std::vector<McCormick<ISVar<T>>>     MCISMVar;
  std::vector<McCormick<ISVar<T>>>     MCISMRes;

  //! @brief Storage for affine superposition models
  ASModel<T>*                          ASMEnv;
  std::vector<ASVar<T>>                ASMVar;
  std::vector<ASVar<T>>                ASMRes;
  std::vector<PolVar<T>>               POLLASMAux;
  std::vector<PolVar<T>>               POLUASMAux;
  std::vector<double>                  DXASMAux;
  std::vector<double>                  DYASMAux;

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
      AUX=0, //!< Auxiliary variable polyhedral relaxation
      INT,   //!< Interval bounds
      MC,    //!< McCormick relaxation with interval bounds
      ISM,   //!< Interval superposition model
      MCISM, //!< McCormick relaxation with interval superposition bounds
      ASM    //!< Affine superposition model
    };

    //! @brief Default constructor
    Options():
      RELAX(MC), ISMDIV(64), ASMBPS(8), ISMCONT(true), ISMSLOPE(true), ISMSHADOW(true), CUTSHADOW(false), 
      ZEROTOL(machprec()), RELU2ABS(false), SIG2EXP(false), AUTODIFF(F)
      {}

    //! @brief Assignment operator
    Options& operator=
      ( Options const& opt ){
        RELAX     = opt.RELAX;
        ISMDIV    = opt.ISMDIV;
        ASMBPS    = opt.ASMBPS;
        ISMCONT   = opt.ISMCONT;
        ISMSLOPE  = opt.ISMSLOPE;
        ISMSHADOW = opt.ISMSHADOW;
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
    //! @brief Number of subdivisions in superposition model
    unsigned ISMDIV;
    //! @brief Number of ??? in superposition model   
    unsigned ASMBPS;
    //! @brief Whether to construct continuous or binary relaxation of superposition model
    bool     ISMCONT;
    //! @brief Whether to propagate slopes in superposition model
    bool     ISMSLOPE;
    //! @brief Whether to propagate shadow remainders in superposition model
    bool     ISMSHADOW;
    //! @brief Whether to append cuts from shadow remainders in superposition model relaxation
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
  : nin(0), nout(0), DAG(nullptr,nullptr), POLEnv(nullptr), ISMEnv(nullptr), ASMEnv(nullptr)
  {}

  ~MLP() 
  {
    delete POLEnv;
    delete DAG.first;
    delete ISMEnv;
    delete ASMEnv;
  }

  //! @brief Clear MLP data
  void clear_data
    ();

  //! @brief Set MLP data
  bool set_data
    ( std::vector<std::pair<std::vector<std::vector<double>>,int>> const& data );

  //! @brief Append multi-neuron layer to MLP data
  bool append_data
    ( std::vector<std::vector<double>> const& data, int const activ=Options::LINEAR, bool const reset=false );

  //! @brief Append single-neuron layer to MLP data
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
    ( ISVar<T>* y, ISVar<T> const* x )
    { evaluate( y, x, ISMhid ); }
  void evaluate
    ( McCormick<ISVar<T>>* y, McCormick<ISVar<T>> const* x )
    { evaluate( y, x, MCISMhid ); }
  void evaluate
    ( ASVar<T>* y, ASVar<T> const* x )
    { evaluate( y, x, ASMhid ); }
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

  void append_ASMcuts
    ( PolImg<T>* img, FFOp* pop, PolVar<T> const& vRes, PolVar<T>* vVar,
      std::vector<UnivarPWL<T>> const& pwlEst );

  void append_ASMcuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const& vRes, PolVar<T>* vVar,
      double const& rhsEst, std::vector<double> const& lnrEst=std::vector<double>() );
  
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
  ISVar<U> ReLU
    ( ISVar<U> const& x )
    const
    {
      return relu( x );
    }

  template <typename U>
  ASVar<U> ReLU
    ( ASVar<U> const& x )
    const
    {
      return relu( x );
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

  //! @brief Intermediate storage for MLP evaluation
  static thread_local std::vector<std::vector<double>>              Dhid;
  static thread_local std::vector<std::vector<fadbad::F<double>>>   FDhid;
  static thread_local std::vector<std::vector<T>>                   Ihid;
  static thread_local std::vector<std::vector<McCormick<T>>>        MCIhid;
  static thread_local std::vector<std::vector<ISVar<T>>>            ISMhid;
  static thread_local std::vector<std::vector<McCormick<ISVar<T>>>> MCISMhid;
  static thread_local std::vector<std::vector<ASVar<T>>>            ASMhid;
  static thread_local std::vector<std::vector<PolVar<T>>>           POLhid;
};

template <typename T> inline thread_local std::vector<std::vector<double>>              MLP<T>::Dhid     = std::vector<std::vector<double>>();
template <typename T> inline thread_local std::vector<std::vector<fadbad::F<double>>>   MLP<T>::FDhid    = std::vector<std::vector<fadbad::F<double>>>();
template <typename T> inline thread_local std::vector<std::vector<T>>                   MLP<T>::Ihid     = std::vector<std::vector<T>>();
template <typename T> inline thread_local std::vector<std::vector<McCormick<T>>>        MLP<T>::MCIhid   = std::vector<std::vector<McCormick<T>>>();
template <typename T> inline thread_local std::vector<std::vector<ISVar<T>>>            MLP<T>::ISMhid   = std::vector<std::vector<ISVar<T>>>();
template <typename T> inline thread_local std::vector<std::vector<McCormick<ISVar<T>>>> MLP<T>::MCISMhid = std::vector<std::vector<McCormick<ISVar<T>>>>();
template <typename T> inline thread_local std::vector<std::vector<ASVar<T>>>            MLP<T>::ASMhid   = std::vector<std::vector<ASVar<T>>>();
template <typename T> inline thread_local std::vector<std::vector<PolVar<T>>>           MLP<T>::POLhid   = std::vector<std::vector<PolVar<T>>>();

template<typename T>
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
  return true;
}

template<typename T>
inline void
MLP<T>::clear_data
()
{
  // Reset data
  nin  = nout = nhid = 0;
  data.clear();
}

template<typename T>
inline bool
MLP<T>::append_data
( std::vector<std::vector<double>> const& layer, int const activ, bool const reset )
{
  if( reset ) clear_data();
  if( !layer.size()
   || layer[0].size() <= 1
   || ( !data.empty() && data.back().first.size() != layer[0].size() - 1 ) )
    return false;

  // Set data
  data.push_back( std::make_pair( layer, activ ) );
  nin  = data.front().first.front().size() - 1;
  nout = data.back().first.size();
  nhid = data.size() - 1;
  return true;
}

template<typename T>
inline bool
MLP<T>::append_data
( std::vector<double> const& layer, int const activ, bool const reset )
{
  return append_data( std::vector<std::vector<double>>({layer}), activ, reset );
}

template<typename T>
template<typename U>
inline void
MLP<T>::evaluate
( U* y, U const* x, std::vector<std::vector<U>>& vhid )
const
{
  // Propagate through hidden layers
  vhid.resize( nhid );
#ifdef MC__MCMLP_DEBUG
  std::cerr << "No hidden layers: " << nhid << std::endl;
#endif
  for( unsigned l=0; l<nhid; ++l ){
#ifdef MC__MCMLP_CHECK
    assert( data[l].first.size() ); // number of neurons in layer l+1
#endif
    size_t const nneu = data[l].first.size();
    vhid[l].resize( nneu );
#ifdef MC__MCMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for( unsigned i=0; i<nneu; ++i ){
      vhid[l][i] = (data[l].first)[i][0]; // bias term
#ifdef MC__MCMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": " << data[l][i].size()-1 << std::endl;
#endif
      for( unsigned j=0; j<(data[l].first)[i].size()-1; ++j ){
#ifdef MC__MCMLP_DEBUG
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
#ifdef MC__MCMLP_CHECK
  assert( data.back().first.size() ); // number of neurons in layer l+1
#endif
  size_t const nneu = data.back().first.size();
#ifdef MC__MCMLP_DEBUG
  std::cerr << "No neurons in layer " << nhid << ": " << nneu << std::endl;
#endif
  for( unsigned i=0; i<nneu; ++i ){
    y[i] = (data.back().first)[i][0]; // bias term
#ifdef MC__MCMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << nhid << ": " << (data.back().first)[i].size()-1 << std::endl;
#endif
    for( unsigned j=0; j<(data.back().first)[i].size()-1; ++j ){
#ifdef MC__MCMLP_DEBUG
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

template<typename T>
inline void
MLP<T>::resize_relax
()
{
  // Set relaxation environment and containers
  switch( options.RELAX ){
    case Options::AUX:
    default:
      if( !DAG.first || DAG.second != &data ){
        delete DAG.first;
        DAG = std::make_pair( new FFGraph, &data );
        DAGVar.resize( nin );
        DAGRes.resize( nout );
        for( unsigned i=0; i<nin; ++i )
          DAGVar[i].set( DAG.first );
        evaluate( DAGRes.data(), DAGVar.data() );
        DAGOps = DAG.first->subgraph( nout, DAGRes.data() );
#ifdef MC__MCMLP_DEBUG
        DAG.first->output( DAGOps, " MLP", std::cerr );
        std::ofstream oFile( "MLP.dot", std::ios_base::out );
        DAG.first->dot_script( nout, DAGRes.data(), oFile );
        oFile.close();
#endif
      }
      if( !POLEnv ) POLEnv = new PolImg<T>;
      POLVar.resize( nin );
      POLRes.resize( nout );
      break;

    case Options::INT:
      IVar.resize( nin );
      IRes.resize( nout );
      return;

    case Options::MC:
      MCVar.resize( nin );
      MCRes.resize( nout );
      return;

    case Options::MCISM:
      MCISMVar.resize( nin );
      MCISMRes.resize( nout );
      // no break

    case Options::ISM:
      if( !ISMEnv || ISMEnv->nvar() != nin || ISMEnv->ndiv() != options.ISMDIV ){
        delete ISMEnv;
        ISMEnv = new ISModel<T>( nin, options.ISMDIV );
      }
      ISMEnv->options.SLOPE_USE  = options.ISMSLOPE;  
      ISMEnv->options.SHADOW_USE = options.ISMSHADOW;
      ISMVar.resize( nin );
      ISMRes.resize( nout );
      POLISMAux.resize( nin );
      DLISMAux.resize( options.ISMDIV );
      DUISMAux.resize( options.ISMDIV );
      return;

    case Options::ASM:
      if( !ASMEnv || ASMEnv->nvar() != nin || ASMEnv->ndiv() != options.ISMDIV ){
        delete ASMEnv;
        ASMEnv = new ASModel<T>( nin, options.ISMDIV );
      }
      ASMEnv->options.SLOPE_USE  = options.ISMSLOPE;
      ASMEnv->options.SHADOW_USE = options.ISMSHADOW;
      ASMVar.resize( nin );
      ASMRes.resize( nout );
      POLLASMAux.resize( nin );
      POLUASMAux.resize( nin );
      return;
  }
}

template<typename T>
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
#ifdef MC__MCMLP_DEBUG
      std::vector<PolVar<T>> POLwk;
      DAG.first->eval( DAGOps, POLwk, nout, DAGRes.data(), POLRes.data(), nin, DAGVar.data(), POLVar.data() );
      for( auto polv : POLwk ) std::cout << polv.name() << ": " << polv.range() << std::endl;
#else
      DAG.first->eval( DAGOps, nout, DAGRes.data(), POLRes.data(), nin, DAGVar.data(), POLVar.data() );
#endif
      // COULD DO CONSTRAINT PROPAGATION BASED ON vRes.range()
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], POLRes[j].range() );
#ifdef MC__MCMLP_DEBUG
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
#ifdef MC__MCMLP_DEBUG
        std::cerr << "vRes[" << j << "] in " << vRes[j] << std::endl;
#endif
      }
      return;
      
    // McCormick relaxations at mid-point with subgradient in each direction
    case MLP<T>::Options::RELAXTYPE::MC:
      for( unsigned i=0; i<nin; ++i )
      MCVar[i] = McCormick<T>( vVar[i].range(), Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
      evaluate( MCRes.data(), MCVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], MCRes[j].I() );
#ifdef MC__MCMLP_DEBUG
        std::cerr << "vRes[" << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      return;

    // Interval superposition models
    case MLP<T>::Options::RELAXTYPE::ISM:
      for( unsigned i=0; i<nin; ++i )
        ISMVar[i].set( ISMEnv, i, vVar[i].range() );
      evaluate( ISMRes.data(), ISMVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], ISMRes[j].B() );
#ifdef MC__MCMLP_DEBUG
        std::cerr << "ISMRes[" << j << "] in " << ISMRes[j] << std::endl;
        std::cerr << "vRes[" << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      return;

 
    // McCormick relaxation with ISM bounds at mid-point with subgradient in each direction
    case MLP<T>::Options::RELAXTYPE::MCISM:
      for( unsigned i=0; i<nin; ++i )
        MCISMVar[i] = McCormick<ISVar<T>>( ISVar<T>( ISMEnv, i, vVar[i].range() ), Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
      evaluate( MCISMRes.data(), MCISMVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], MCISMRes[j].I().B() );
#ifdef MC__MCMLP_DEBUG
        std::cerr << "ISMRes[" << j << "] in " << MCISMRes[j].I() << std::endl;
        std::cerr << "MCISMRes[" << j << "] in " << MCISMRes[j] << std::endl;
        std::cerr << "vRes[" << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      return;

    // Affine superposition models
    case MLP<T>::Options::RELAXTYPE::ASM:
      UnivarPWLE<double>::nbpsMax = options.ASMBPS;       
      for( unsigned i=0; i<nin; ++i )
        ASMVar[i].set( ASMEnv, i, vVar[i].range() );
      evaluate( ASMRes.data(), ASMVar.data() );
      for( unsigned j=0; j<nout; ++j ){
        vRes[j].set( img, *pRes[j], ASMRes[j].B() );
#ifdef MC__MCMLP_DEBUG
        std::cerr << "ASMRes[" << j << "] in " << ASMRes[j] << std::endl;
        std::cerr << "vRes[" << j << "] in " << vRes[j].range() << std::endl;
#endif
      }
      return;
  }
}

template<typename T>
inline void
MLP<T>::backpropagate_relax
( PolImg<T>* img, FFOp* pOp, PolVar<T> const* vRes, PolVar<T>* vVar )
{
  switch( options.RELAX ){
    // Polyhedral relaxation with auxiliary variables
    case Options::AUX:
    default:
      POLEnv->generate_cuts( nout, POLRes.data() );
#ifdef MC__MCMLP_DEBUG
      std::cerr << "POLEnv:" << *POLEnv << std::endl;
#endif
      img->insert_cuts( POLEnv, POLMap );
      break;

    // Interval bounds
    case MLP<T>::Options::RELAXTYPE::INT:
      return;
      
    // McCormick relaxations at mid-point with subgradient in each direction
    case MLP<T>::Options::RELAXTYPE::MC:
      // polyhedral cut generation
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__MCMLP_DEBUG
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
      return;

    // Interval superposition models
    case MLP<T>::Options::RELAXTYPE::ISM:
#ifdef MC__MCMLP_CHECK
      assert( ISMEnv->ndiv() == options.ISMDIV );
#endif
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        POLISMAux[i].resize( ISMEnv->ndiv() );
        for( unsigned k=0; k<ISMEnv->ndiv(); ++k )
          POLISMAux[i][k].set( img, Op<T>::zeroone(), options.ISMCONT );
      }

      // polyhedral cut generation
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__MCMLP_DEBUG
        std::cerr << "MCRes[" << j << "] in " << MCRes[j] << std::endl;
#endif
        auto cutF1 = *img->add_cut( pOp, PolCut<T>::LE, 0., vRes[j], -1. );
        auto cutF2 = *img->add_cut( pOp, PolCut<T>::GE, 0., vRes[j], -1. );
        for( unsigned i=0; i<nin; ++i ){
          auto const& ISMi = ISMRes[j].C()[i];
          if( ISMi.empty() ) continue;
          for( unsigned k=0; k<ISMEnv->ndiv(); ++k ){
            DLISMAux[k] = Op<T>::l( ISMi[k] );
            DUISMAux[k] = Op<T>::u( ISMi[k] );
          }
          cutF1->append( ISMEnv->ndiv(), POLISMAux[i].data(), DLISMAux.data() );
          cutF2->append( ISMEnv->ndiv(), POLISMAux[i].data(), DUISMAux.data() );
        }

        // add polyhedral cuts for ISM-participating variables
        for( unsigned i=0; i<nin; i++ ){
          if( POLISMAux[i].empty() ) continue;
          // auxiliaries add up to 1
          img->add_cut( pOp, PolCut<T>::EQ, 1., ISMEnv->ndiv(), POLISMAux[i].data(), 1. );

          // link auxiliaries to model variables
          PolVar<T> POLvarL( 0. ), POLvarU( 0. );
          auto&& ISMi = ISMVar[i].C()[i];
          assert( !ISMi.empty() );
          for( unsigned k=0; k<ISMEnv->ndiv(); k++ ){
            DLISMAux[k] = Op<T>::l(ISMi[k]);
            DUISMAux[k] = Op<T>::u(ISMi[k]);
          } 
          img->add_cut( pOp, PolCut<T>::LE, 0., ISMEnv->ndiv(), POLISMAux[i].data(), DLISMAux.data(), vVar[i], -1. );
          img->add_cut( pOp, PolCut<T>::GE, 0., ISMEnv->ndiv(), POLISMAux[i].data(), DUISMAux.data(), vVar[i], -1. );
        }
      }
      return;

    // McCormick relaxation with ISM bounds at mid-point with subgradient in each direction    
    case MLP<T>::Options::RELAXTYPE::MCISM:
#ifdef MC__MCMLP_CHECK
      assert( ISMEnv->ndiv() == options.ISMDIV );
#endif
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        POLISMAux[i].resize( ISMEnv->ndiv() );
        for( unsigned k=0; k<ISMEnv->ndiv(); ++k )
          POLISMAux[i][k].set( img, Op<T>::zeroone(), options.ISMCONT );  
      }
  
      // polyhedral cut generation
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__MCMLP_DEBUG
        std::cerr << "MCISMRes[" << j << "] in " << MCISMRes[j] << std::endl;
#endif
        auto cutF1 = *img->add_cut( pOp, PolCut<T>::LE, 0., vRes[j], -1. );
        auto cutF2 = *img->add_cut( pOp, PolCut<T>::GE, 0., vRes[j], -1. );
        for( unsigned i=0; i<nin; ++i ){
          auto const& ISMi = MCISMRes[j].I().C()[i];
          if( ISMi.empty() ) continue;
          for( unsigned k=0; k<ISMEnv->ndiv(); ++k ){
            DLISMAux[k] = Op<T>::l( ISMi[k] );
            DUISMAux[k] = Op<T>::u( ISMi[k] );
          }
          cutF1->append( ISMEnv->ndiv(), POLISMAux[i].data(), DLISMAux.data() );
          cutF2->append( ISMEnv->ndiv(), POLISMAux[i].data(), DUISMAux.data() );
        }

        // add polyhedral cuts for ISM-participating variables
        for( unsigned i=0; i<nin; i++ ){
          if( POLISMAux[i].empty() ) continue;
          // auxiliaries add up to 1
          img->add_cut( pOp, PolCut<T>::EQ, 1., ISMEnv->ndiv(), POLISMAux[i].data(), 1. );

          // link ISM auxiliaries to model variables
          PolVar<T> POLvarL( 0. ), POLvarU( 0. );
          auto&& ISMi = MCISMVar[i].I().C()[i];
          assert( !ISMi.empty() );
          for( unsigned k=0; k<ISMEnv->ndiv(); k++ ){
            DLISMAux[k] = Op<T>::l(ISMi[k]);
            DUISMAux[k] = Op<T>::u(ISMi[k]);
          }
          img->add_cut( pOp, PolCut<T>::LE, 0., ISMEnv->ndiv(), POLISMAux[i].data(), DLISMAux.data(), vVar[i], -1. );
          img->add_cut( pOp, PolCut<T>::GE, 0., ISMEnv->ndiv(), POLISMAux[i].data(), DUISMAux.data(), vVar[i], -1. );
        }

        // polyhedral cut generation for MC  
        double rhs1 = -MCISMRes[j].cv(),  
               rhs2 = -MCISMRes[j].cc();
        for( unsigned i=0; i<nin; ++i ){ 
          rhs1 += MCISMRes[j].cvsub(i) * MCISMVar[i].cv();
          rhs2 += MCISMRes[j].ccsub(i) * MCISMVar[i].cc();
        }
        img->add_cut( pOp, PolCut<T>::LE, rhs1, nin, vVar, MCISMRes[j].cvsub(), vRes[j], -1. );
        img->add_cut( pOp, PolCut<T>::GE, rhs2, nin, vVar, MCISMRes[j].ccsub(), vRes[j], -1. );
      }
      return;
      
    // Affine superposition models
    case MLP<T>::Options::RELAXTYPE::ASM:
#ifdef MC__MCMLP_CHECK
      assert( ASMEnv->ndiv() == options.ISMDIV );
#endif
      // polyhedral cut generation
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__MCMLP_DEBUG
        std::cerr << "ASMRes[" << j << "] in " << ASMRes[j] << std::endl;
#endif
        switch(ASMRes[j].get_ASVar()){
          case 1:
            append_ASMcuts( img, pOp, vRes[j], vVar, ASMRes[j].get_cst() );
            break;    
                    
          case 2:
            append_ASMcuts( img, pOp, vRes[j], vVar, ASMRes[j].get_cst(), ASMRes[j].get_lnr() );
            break;

          case 3:
            append_ASMcuts( img, pOp, vRes[j], vVar, ASMRes[j].get_lst() );          
            break;

          case 4:
            append_ASMcuts( img, pOp, vRes[j], vVar, ASMRes[j].get_lst() );
            if( options.CUTSHADOW && options.ISMSHADOW )
              append_ASMcuts( img, pOp, vRes[j], vVar, ASMRes[j].get_shadow() );
            break;

          default:
           throw std::runtime_error("MLP::relax: **ERROR** invalid flag from get_ASVar()\n");
        }
      }
      return;
  }

#ifdef MC__MCMLP_DEBUG
  std::cerr << *img;
  {int dum; std::cout << "PAUSED, ENTER 1"; std::cin >> dum;}
#endif
}

template<typename T>
inline void
MLP<T>::append_ASMcuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& vRes, PolVar<T>* vVar,
  double const& rhsEst, std::vector<double> const& lnrEst )
{
  double rhs = rhsEst / (double)nin;
  for( unsigned i=0; i<nin; ++i ){
    double DX = mc::Op<T>::diam(vVar[i].range()),
           DY = !lnrEst.empty()? lnrEst[i] * DX: 0.,
           XL = mc::Op<T>::l(vVar[i].range()),
           YL = !lnrEst.empty()? lnrEst[i] * DX + rhs: rhs;
    T IY = !lnrEst.empty()? lnrEst[i] * vVar[i].range() + rhs: T(rhs);
    POLLASMAux[i].set( img, IY, true );
    img->add_cut( pOp, mc::PolCut<T>::EQ, DX*YL-DY*XL, POLLASMAux[i], DX, vVar[i], -DY );
  }
  img->add_cut( pOp, PolCut<T>::EQ, 0., nin, POLLASMAux.data(), 1., vRes, -1. );
}

template<typename T>
inline void
MLP<T>::append_ASMcuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& vRes, PolVar<T>* vVar,
  std::vector<UnivarPWL<T>> const& pwlEst )
{
  for( unsigned i=0; i<nin; ++i ){
    UnivarPWLE<double> const& uest = pwlEst[i].undEst;
    if( uest.empty() )
      POLLASMAux[i].set( img, T(0.), true );
    else{ 
      POLLASMAux[i].set( img, T(uest.get_lb(),uest.get_ub()), true );
      auto const [ucst,isuCst] = uest.get_cst();
      if( isuCst )
        img->add_cut( pOp, PolCut<T>::EQ, ucst, POLLASMAux[i], 1. );
      else{
        unsigned NK = uest.first.size()-1;
#ifdef MC__MCMLP_CHECK
        assert( uest.second.size() == uest.first.size() );
#endif
        if(NK==1){
          NK += 1;
          DXASMAux.assign( NK, 0. );
          DYASMAux.assign( NK, 0. );
          for( unsigned j=0; j<NK; ++j ){
            DXASMAux[j] = uest.first[j];
            DYASMAux[j] = uest.second[j];
          }            
        }
        else{
          DXASMAux.assign( NK, uest.first[0] );
          DYASMAux.assign( NK, uest.second[0] );
          for( unsigned j=0; j<NK; ++j ){
            DXASMAux[j] += uest.first[j+1];
            DYASMAux[j] += uest.second[j+1];
          }
        }
        img->add_semilinear_cuts( pOp, NK, vVar[i], DXASMAux.data(), POLLASMAux[i], DYASMAux.data(), mc::PolCut<T>::EQ );
      }
    }

    UnivarPWLE<double> const& oest = pwlEst[i].oveEst;
    if( oest.empty() )
      POLUASMAux[i].set( img, T(0.), true );
    else{
      POLUASMAux[i].set( img, T(oest.get_lb(),oest.get_ub()), true );
      auto const [ocst,isoCst] = oest.get_cst();
      if( isoCst )
        img->add_cut( pOp, PolCut<T>::EQ, ocst, POLUASMAux[i], 1. );
      else{
        unsigned NK = oest.first.size()-1;
#ifdef MC__MCMLP_CHECK
        assert( oest.second.size() == oest.first.size() );
#endif
        if(NK == 1){
          NK += 1;
          DXASMAux.assign( NK, 0. );
          DYASMAux.assign( NK, 0. );
          for( unsigned j=0; j<NK; ++j ){
            DXASMAux[j] = oest.first[j];
            DYASMAux[j] = oest.second[j];
          }
        }  
        else{  
          DXASMAux.assign( NK, oest.first[0] );
          DYASMAux.assign( NK, oest.second[0] );
          for( unsigned j=0; j<NK; ++j ){
            DXASMAux[j] += oest.first[j+1];
            DYASMAux[j] += oest.second[j+1];
          }
        }	      
        img->add_semilinear_cuts( pOp, NK, vVar[i], DXASMAux.data(), POLUASMAux[i], DYASMAux.data(), mc::PolCut<T>::EQ );
      }
    }
  }
  img->add_cut( pOp, PolCut<T>::LE, 0., nin, POLLASMAux.data(), 1., vRes, -1. );
  img->add_cut( pOp, PolCut<T>::GE, 0., nin, POLUASMAux.data(), 1., vRes, -1. );
} 

//! @brief C++ class defining neural networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLP is a C++ class for defining neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template<typename T, unsigned int ID>
class FFMLP
////////////////////////////////////////////////////////////////////////
: public FFOp
{
public:

  //! @brief Constructors
  FFMLP
    ()
    : FFOp( (int)EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      info = ID;
#ifdef MC__MCMLP_CHECK
      assert( nVar == pMLP->nin && idep < pMLP->nout );
#endif
      return *(insert_external_operation( *this, pMLP->nout, nVar, pVar )[idep]);
    }

  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      info = ID;
#ifdef MC__MCMLP_CHECK
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
      else if( idU == typeid( ISVar<T> ) )
        return eval( nRes, static_cast<ISVar<T>*>(vRes), nVar, static_cast<ISVar<T> const*>(vVar), mVar );
      else if( idU == typeid( ASVar<T> ) )
        return eval( nRes, static_cast<ASVar<T>*>(vRes), nVar, static_cast<ASVar<T> const*>(vVar), mVar );
      else if( idU == typeid( PolVar<T> ) )
        return eval( nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar );

      throw std::runtime_error( "FFMLP::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
    }

  template<typename U>
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
//! mc::FFMLP is a C++ class for defining gradient of neural networks
//! as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template<typename T, unsigned int ID>
class FFGRADMLP
////////////////////////////////////////////////////////////////////////
: public FFOp
{
public:

  //! @brief Constructors
  FFGRADMLP
    ()
    : FFOp( (int)EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      info = ID+1;
#ifdef MC__MCMLP_CHECK
      assert( nVar == pMLP->nin && idep < pMLP->nout );
#endif
      return *(insert_external_operation( *this, nVar*pMLP->nout, nVar, pVar )[idep]);
    }

  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar, MLP<T>* pMLP )
    const
    {
      data = pMLP;
      info = ID+1;
#ifdef MC__MCMLP_CHECK
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
    { std::ostringstream oss; oss << data; return "GRADMLP[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};


template<typename T, unsigned int ID>
template<typename U>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, U* vRes, unsigned const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__MLPOP_TRACE
  std::cout << "FFMLP::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  pMLP->evaluate( vRes, vVar );
}

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFDep\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  FFVar** pRes = operator()( nVar, vVar, pMLP );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: fadbad::F<FFVar>\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = operator()( nVar, vVarVal.data(), pMLP );
  FFGRADMLP<T,ID> ResDer;
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

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::deriv: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  FFGRADMLP<T,ID> ResDer;
  FFVar const*const* vResDer = ResDer( nVar, vVar, pMLP );
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
      vDer[k][i] = *vResDer[k+nRes*i];
}

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: SLiftVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template<typename T, unsigned int ID>
inline void
FFMLP<T,ID>::eval
( unsigned const nRes, PolVar<T>* vRes, unsigned const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: PolVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__MCMLP_CHECK
  assert( img && dag );
#endif
  FFVar** pRes = dag->curOp()->varout.data(); // ACCOUNT FOR MULTIPLE OUTPUTS
#ifdef MC__MCMLP_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  pMLP->resize_relax();
  pMLP->propagate_relax( img, pRes, vRes, vVar );
}

template<typename T, unsigned int ID>
inline bool
FFMLP<T,ID>::reval
( unsigned const nRes, PolVar<T> const* vRes, unsigned const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: PolVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nout );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__MCMLP_CHECK
  assert( img && pop );
#endif

  pMLP->backpropagate_relax( img, pop, vRes, vVar );
  return true;
}

template<typename T, unsigned int ID>
inline void
FFGRADMLP<T,ID>::eval
( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADMLP::eval: double\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
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

template<typename T, unsigned int ID>
inline void
FFGRADMLP<T,ID>::eval
( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGRADMLP::eval: FFDep\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template<typename T, unsigned int ID>
inline void
FFGRADMLP<T,ID>::eval
( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADMLP::eval: FFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  FFVar** pRes = operator()( nVar, vVar, pMLP );
  for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(pRes[j]);
}

template<typename T, unsigned int ID>
inline void
FFGRADMLP<T,ID>::eval
( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFODE_TRACE
  std::cout << "FFGRADMLP::eval: SLiftFVar\n";
#endif
  MLP<T>* pMLP = static_cast<MLP<T>*>( data );
#ifdef MC__MCMLP_CHECK
  assert( pMLP && nVar == pMLP->nin && nRes == pMLP->nin * pMLP->nout );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

} // end namespace mc

#endif
