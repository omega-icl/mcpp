#ifndef MC__FFEXTERN_HPP
#define MC__FFEXTERN_HPP

//#define MC__FFEXTERN_DEBUG
#define MC__FFEXTERN_CHECK

#include "mccormick.hpp"
#include "scmodel.hpp"
#include "supmodel.hpp"
#include "pwcu.hpp"
#include "pwlu.hpp"
#include "mcfunc.hpp"
#include "ffunc.hpp"
#include "ffexpr.hpp"
#include "ffdep.hpp"
#include "polimage.hpp"
#include "slift.hpp"

namespace mc
{

//! @brief C++ class defining external DAG operations in MC++ for template-type objects
////////////////////////////////////////////////////////////////////////
//! mc::FFEXTERN is a C++ class for defining an external DAG operations
//! in MC++ for a generic template-type object O. The other template
//! parameter T specifies the type for interval arithmetic.
////////////////////////////////////////////////////////////////////////
template< typename T, typename O >
class FFEXTERN
////////////////////////////////////////////////////////////////////////
: public FFOp
{

protected:

  // pointer to the object
  O*             _ptrObj;
  // whether this class owns the object
  bool           _ownObj;

  // set tracking any exception thrown during evaluation
  mutable std::set<int>  _excpEval;

  //! @brief Storage for Polyhedral relaxation
  mutable PolImg<T>*                             _POLEnv;
  mutable std::vector<PolVar<T>>                 _POLVar;
  mutable std::vector<PolVar<T>>                 _POLRes;
  mutable std::map<PolVar<T> const*,PolVar<T>,lt_PolVar<T>> _POLMap;

  //! @brief Storage for Interval bounds
  mutable std::vector<T>                         _IVar;
  mutable std::vector<T>                         _IRes;

  //! @brief Storage for McCormick relaxation
  mutable std::vector<McCormick<T>>              _MCVar;
  mutable std::vector<McCormick<T>>              _MCRes;

  //! @brief Storage for sparse Chebyshev models
  mutable SCModel<T>*                            _SCEnv;
  mutable std::vector<SCVar<T>>                  _SCVar;
  mutable std::vector<SCVar<T>>                  _SCRes;
  //std::map< t_mon, FFVar, lt_mon > _SCMon;

  //! @brief Storage for interval superposition models
  mutable SupModel<PWCU>*                        _PWCSEnv;
  mutable std::vector<SupVar<PWCU>>              _PWCSVar;
  mutable std::vector<SupVar<PWCU>>              _PWCSRes;
  mutable std::vector<std::vector<PolVar<T>>>    _POLPWCSAux;
  mutable std::vector<double>                    _DLPWCSAux;
  mutable std::vector<double>                    _DUPWCSAux;
  mutable std::vector<McCormick<SupVar<PWCU>>>   _MCPWCSVar;
  mutable std::vector<McCormick<SupVar<PWCU>>>   _MCPWCSRes;

  //! @brief Storage for affine superposition models
  mutable SupModel<PWLU>*                        _PWLSEnv;
  mutable std::vector<SupVar<PWLU>>              _PWLSVar;
  mutable std::vector<SupVar<PWLU>>              _PWLSRes;
  mutable std::vector<PolVar<T>>                 _POLPWLSAux;
  mutable std::vector<double>                    _DXPWLSAux;
  mutable std::vector<double>                    _DYPWLSAux;
  mutable std::vector<McCormick<SupVar<PWLU>>>   _MCPWLSVar;
  mutable std::vector<McCormick<SupVar<PWLU>>>   _MCPWLSRes;

  //! @brief Resize relaxation containers for the object
  void _resize_relax
    ( PolImg<T>* img )
    const;

  //! @brief Propagate polyhedral image for the object
  void _propagate_relax
    ( PolImg<T>* img, FFVar** pRes, PolVar<T>* vRes, PolVar<T> const* vVar )
    const;

  //! @brief Append polyhedral image cuts for the object
  void _backpropagate_relax
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const* vRes, PolVar<T>* vVar )
    const;

private:

  //! @brief Append polyhedral image cuts for McCormick relaxation
  template< typename U >
  void _append_MCcuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, McCormick<U> const& MCRes,
      PolVar<T> const* vVar, std::vector<McCormick<U>> const& vMCVar )
    const;

  //! @brief Append polyhedral image cuts for variables in piecewise-constant superposition relaxation
  void _append_PWCScuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T>& Var, std::vector<PolVar<T>>& vAux,
      PWCU const& uest, PWCU const& oest )
    const;

  //! @brief Append polyhedral image cuts for dependent in piecewise-constant superposition relaxation
  void _append_PWCScuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, std::vector<std::vector<PolVar<T>>> const& vAux,
      std::vector<PWCU> const& est, std::set<unsigned int> const& dep, bool const under )
    const;

  //! @brief Append polyhedral image cuts for variables and dependent in piecewise-linear superposition relaxation
  void _append_PWLScuts
    ( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, PolVar<T>* vVar,
      std::vector<PolVar<T>>& vAux, std::vector<double>& Xk, std::vector<double>& Yk, 
      std::vector<PWLU> const& est, std::set<unsigned int> const& dep, bool const under )
    const;

public:

  //! @brief FFEXTERN options
  struct Options
  {
    //! @brief Relaxation type
    enum RELAXTYPE{
      INT=0,   //!< Interval bounds
      AUX,     //!< Auxiliary variable polyhedral relaxation
      MC,      //!< McCormick relaxation with interval bounds
      PWCS,    //!< Piecewise-constant superposition model
      MCPWCS,  //!< McCormick relaxation with piecewise-constant superposition bounds
      PWLS,    //!< Piecewise-linear superposition model
      MCPWLS   //!< McCormick relaxation with piecewise-linear superposition bounds
    };

    //! @brief Default constructor
    Options()
      {
        reset();
      }

    //! @brief Assignment operator
    Options& operator=
      ( Options const& opt ){
        RELAX     = opt.RELAX;

        PWCDIV    = opt.PWCDIV;
        PWCREL    = opt.PWCREL;
        PWCSHADOW = opt.PWCSHADOW;
        PWCSUP    = opt.PWCSUP;

        PWLINI    = opt.PWLINI;
        PWLMAX    = opt.PWLMAX;
        PWLREL    = opt.PWLREL;
        PWLSHADOW = opt.PWLSHADOW;
        PWLSUP    = opt.PWLSUP;

        return *this;
      }

    //! @brief Reset options
    void reset
      ()
      {
        RELAX     = { MC };

        POLDEF    = 1;
        POLIMG.reset();

        PWCDIV    = 16;
        PWCREL    = 0;
        PWCSHADOW = 0;
        PWCSUP.reset();

        PWLINI    = 1;
        PWLMAX    = 0;
        PWLREL    = 0;
        PWLSHADOW = 0;
        PWLSUP.reset();
      }

    //! @brief Type of relaxation
    std::set<RELAXTYPE>               RELAX;
    //! @brief Whether for options to default to the parent polyhedral relaxation environment
    bool                              POLDEF;
    //! @brief Options for polyhedral relaxation environment
    typename PolImg<T>::Options       POLIMG;

    //! @brief Equipartition size in piecewise-constant superposition model
    size_t                            PWCDIV;
    //! @brief Representation of piecewise-constant univariates - 0: continuous relaxation; 1: binary encoding
    int                               PWCREL;
    //! @brief Whether to append cuts from shadow estimators in piecewise-constant superposition model relaxation
    bool                              PWCSHADOW;
    //! @brief Options for superposition model with piecewise-constant univariate estimators
    typename SupModel<PWCU>::Options  PWCSUP;

    //! @brief Initial partition size in piecewise-linear superposition model
    size_t                            PWLINI;
    //! @brief Maximal partition size in piecewise-linear superposition model
    size_t                            PWLMAX;
    //! @brief Representation of piecewise-linear univariates - 0: continuous relaxation; 1: binary encoding; 2: SOS2 encoding
    int                               PWLREL;
    //! @brief Whether to append cuts from shadow estimators in piecewise-linear superposition model relaxation
    bool                              PWLSHADOW;
    //! @brief Options for superposition model with adaptive piecewise-linear univariate estimators
    typename SupModel<PWLU>::Options  PWLSUP;
  }    options;

  //! @brief Default Constructor
  FFEXTERN
    ()
    : FFOp     ( EXTERN ),
      _ptrObj  ( nullptr ),
      _ownObj  ( false ),
      _POLEnv  ( nullptr ),
      _PWCSEnv ( nullptr ),
      _PWLSEnv ( nullptr )
    {}

  // Destructor
  virtual ~FFEXTERN
    ()
    {
      if( _ownObj && _ptrObj ) delete _ptrObj;
      if( _POLEnv ) delete _POLEnv;
      delete _PWCSEnv;
      delete _PWLSEnv;
    }

  // Copy constructor
  FFEXTERN
    ( FFEXTERN<T,O> const& other )
    : FFOp( other ),
      _POLEnv  ( nullptr ),
      _PWCSEnv ( nullptr ),
      _PWLSEnv ( nullptr ),
      options  ( other.options )
    {
#ifdef MC__FFEXTERN_TRACE
      std::cout << "FFEXTERN::copy constructor\n";
#endif
      if( !other._ptrObj )
        throw std::runtime_error( "FFEXTERN::copy constructor ** Null pointer to external object\n" );

      _ownObj = other._ownObj;      
      if( _ownObj )
        _ptrObj = new O( *other._ptrObj );
      else
        _ptrObj = other._ptrObj;
    }
};

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_resize_relax
( PolImg<T>* img )
const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();

  // Used as intermediate to intersect bounds computed with different arithmetics
  //std::cout << "\nnout = " << nout << std::endl;
  _IRes.resize( nout );

  // Set relaxation environment and containers
  for( auto relax : options.RELAX ){
    switch( relax ){

      case Options::INT:
      default:
        _IVar.resize( nin );
        //_IRes.resize( nout );
        break;

      case Options::AUX:
        if( !_POLEnv ) _POLEnv = new PolImg<T>;
        else           _POLEnv->reset();
        _POLEnv->options = (options.POLDEF? img->options: options.POLIMG);
        _POLVar.resize( nin );
        _POLRes.resize( nout );
        _POLMap.clear();
        break;

      case Options::MC:
        _MCVar.resize( nin );
        _MCRes.resize( nout );
        break;

      case Options::MCPWCS:
        _MCPWCSVar.resize( nin );
        _MCPWCSRes.resize( nout );
        // no break

      case Options::PWCS:
        if( !_PWCSEnv || _PWCSEnv->nvar() != nin ){
          delete _PWCSEnv;
          _PWCSEnv = new SupModel<PWCU>( nin );
        }
        _PWCSEnv->options = options.PWCSUP;
        _PWCSVar.resize( nin );
        _PWCSRes.resize( nout );
        _POLPWCSAux.resize( nin );
        _DLPWCSAux.resize( options.PWCDIV );
        _DUPWCSAux.resize( options.PWCDIV );
        break;

      case Options::MCPWLS:
        _MCPWLSVar.resize( nin );
        _MCPWLSRes.resize( nout );
        // no break

      case Options::PWLS:
        if( !_PWLSEnv || _PWLSEnv->nvar() != nin ){
          delete _PWLSEnv;
          _PWLSEnv = new SupModel<PWLU>( nin );
        }
        _PWLSEnv->options = options.PWLSUP;
        _PWLSVar.resize( nin );
        _PWLSRes.resize( nout );
        _POLPWLSAux.resize( nin );
        break;
    }
  }
}

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_propagate_relax
( PolImg<T>* img, FFVar** pRes, PolVar<T>* vRes, PolVar<T> const* vVar )
const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();

  // track any exceptions thrown
  _excpEval.clear();

  // Forward-propagate relaxations through object
  bool firstEval = true;
  for( auto relax : options.RELAX ){
    try{
      switch( relax ){

      // Interval bounds
      case Options::INT:
      default:
        for( unsigned i=0; i<nin; ++i )
          _IVar[i] = vVar[i].range();
        _ptrObj->eval( _IVar.data(), _IRes.data() );
#ifdef MC__FFEXTERN_DEBUG
        for( unsigned j=0; j<nout; ++j )
          std::cerr << "_IRes[" << j << "]: " << _IRes[j] << std::endl;
#endif
        break;

      // Polyhedral relaxation with auxiliary variables
      case Options::AUX:
        for( unsigned i=0; i<nin; ++i ){
          _POLVar[i].set( _POLEnv, _ptrObj->varin()[i], vVar[i].range(), true );
          _POLMap[&_POLVar[i]] = vVar[i];
        }
        _ptrObj->eval( _POLVar.data(), _POLRes.data() );
        // COULD DO CONSTRAINT PROPAGATION BASED ON vRes.range() - NEEDS A PRIORI BoUNDS IN vRes?
        for( unsigned j=0; j<nout; ++j ){
          if( firstEval )
            _IRes[j] = _POLRes[j].range();
          else if( !Op<T>::inter( _IRes[j], _IRes[j], _POLRes[j].range() ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << _POLRes[j].range() << std::endl;
//#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_POLRes[" << j << "]: " << _POLRes[j] << std::endl;
//#endif
        }
        break;

      // McCormick relaxations at mid-point with subgradient in each direction
      case Options::MC:
        for( unsigned i=0; i<nin; ++i )
          _MCVar[i] = McCormick<T>( vVar[i].range(), Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
        _ptrObj->eval( _MCVar.data(), _MCRes.data() );
        for( unsigned j=0; j<nout; ++j ){
          if( firstEval )
            _IRes[j] = _MCRes[j].I();
          else if( !Op<T>::inter( _IRes[j], _IRes[j], _MCRes[j].I() ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << _MCRes[j].I() << std::endl;
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_MCRes[" << j << "]: " << _MCRes[j] << std::endl;
#endif
        }
        break;

      // Piecewise-constant superposition models
      case Options::PWCS:
        for( unsigned i=0; i<nin; ++i )
          _PWCSVar[i].set( *_PWCSEnv, i, vVar[i].range(), options.PWCDIV );
        _ptrObj->eval( _PWCSVar.data(), _PWCSRes.data() );
        for( unsigned j=0; j<nout; ++j ){
          if( firstEval )
            _IRes[j] = T(_PWCSRes[j].l(),_PWCSRes[j].u());
          else if( !Op<T>::inter( _IRes[j], _IRes[j], T(_PWCSRes[j].l(),_PWCSRes[j].u()) ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << T(_PWCSRes[j].l(),_PWCSRes[j].u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_PWCSRes[" << j << "]: " << _PWCSRes[j] << std::endl;
#endif
        }
        break;
 
      // McCormick relaxation with piecewise-constant superposition bounds and subgradient calculated at mid-point
      case Options::MCPWCS:
        for( unsigned i=0; i<nin; ++i )
          _MCPWCSVar[i] = McCormick<SupVar<PWCU>>( SupVar<PWCU>( *_PWCSEnv, i, vVar[i].range(), options.PWCDIV ),
                                                   Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
        _ptrObj->eval( _MCPWCSVar.data(), _MCPWCSRes.data() );
        for( unsigned j=0; j<nout; ++j ){
          auto const& PWCSRes_j = _MCPWCSRes[j].I();
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_PWCSRes[" << j << "]: " << PWCSRes_j << std::endl;
#endif
          if( firstEval )
            _IRes[j] = T(PWCSRes_j.l(),PWCSRes_j.u());
          else if( !Op<T>::inter( _IRes[j], _IRes[j], T(PWCSRes_j.l(),PWCSRes_j.u()) ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << T(PWCSRes_j.l(),PWCSRes_j.u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_MCPWCSRes[" << j << "]: " << _MCPWCSRes[j] << std::endl;
#endif
        }
        break;

      // Piecewise-linear superposition models
      case Options::PWLS:
        //UnivarPWLE<double>::nbpsMax = options.ASMBPS;       
        for( unsigned i=0; i<nin; ++i )
          _PWLSVar[i].set( *_PWLSEnv, i, vVar[i].range(), options.PWLINI );
        _ptrObj->eval( _PWLSVar.data(), _PWLSRes.data() );
        for( unsigned j=0; j<nout; ++j ){
          if( firstEval )
            _IRes[j] = T(_PWLSRes[j].l(),_PWLSRes[j].u());
          else if( !Op<T>::inter( _IRes[j], _IRes[j], T(_PWLSRes[j].l(),_PWLSRes[j].u()) ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << T(_PWLSRes[j].l(),_PWLSRes[j].u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_PWLSRes[" << j << "]: " << _PWLSRes[j] << std::endl;
#endif
        }
        break;

      // McCormick relaxation with piecewise-linear superposition bounds and subgradient calculated at mid-point
      case Options::MCPWLS:
        for( unsigned i=0; i<nin; ++i )
          _MCPWLSVar[i] = McCormick<SupVar<PWLU>>( SupVar<PWLU>( *_PWLSEnv, i, vVar[i].range(), options.PWLINI ),
                                                   Op<T>::mid( vVar[i].range() ) ).sub( nin, i );
        _ptrObj->eval( _MCPWLSVar.data(), _MCPWLSRes.data() );
        for( unsigned j=0; j<nout; ++j ){
          auto const& PWLSRes_j = _MCPWLSRes[j].I();
          if( firstEval )
            _IRes[j] = T(PWLSRes_j.l(),PWLSRes_j.u());
          else if( !Op<T>::inter( _IRes[j], _IRes[j], T(PWLSRes_j.l(),PWLSRes_j.u()) ) )
            std::cerr << "Empty intersection of relaxation bounds: "
                      << _IRes[j] << "  " << T(PWLSRes_j.l(),PWLSRes_j.u()) << std::endl;
#ifdef MC__FFEXTERN_DEBUG
            std::cerr << "MCPWLSRes[" << j << "]: " << _MCPWLSRes[j] << std::endl;
#endif
        }
        break;
     }
     
     firstEval = false;
   }
   
   catch(...){
     _excpEval.insert( relax );
     continue;
   }
  }
  
  // Set intersected bounds
  for( unsigned j=0; j<nout; ++j )
    vRes[j].set( img, *pRes[j], _IRes[j] );

  // Track polyhedral relaxation dependents
  if( options.RELAX.count( Options::AUX ) )
    for( unsigned j=0; j<nout; ++j )
      _POLMap[&_POLRes[j]] = vRes[j];
}

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_backpropagate_relax
( PolImg<T>* img, FFOp* pOp, PolVar<T> const* vRes, PolVar<T>* vVar )
const
{
  size_t const nin  = varin.size();
  size_t const nout = varout.size();

  // Back-propagate relaxations through object
  for( auto relax : options.RELAX ){
    if( _excpEval.count( relax ) ) break; 
    switch( relax ){

    // Interval bounds
    case Options::INT:
    default:
      // No cuts apart from interval bounds on vRes
      break;

    // Polyhedral relaxation with auxiliary variables
    case Options::AUX:
      _POLEnv->generate_cuts( nout, _POLRes.data() );
#ifdef MC__FFEXTERN_DEBUG
      std::cerr << "_POLEnv:" << *_POLEnv << std::endl;
#endif
      img->insert_cuts( _POLEnv, _POLMap );
      break;
      
    // McCormick relaxations at mid-point with subgradient in each direction
    case Options::MC:
      for( unsigned j=0; j<nout; ++j ){
        // polyhedral cut generation for MC
        _append_MCcuts( img, pOp, vRes[j], _MCRes[j], vVar, _MCVar );
      }
      break;

    // Piecewise-constant superposition models
    case Options::PWCS:
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        _append_PWCScuts( img, pOp, vVar[i], _POLPWCSAux[i], _PWCSVar[i].uest()[i], _PWCSVar[i].oest()[i] );          
      }     
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__FFEXTERN_DEBUG
        std::cerr << "_PWCSRes["   << j << "] in " << _PWCSRes[j] << std::endl;
#endif
        // constant superposition relaxation
        if( _PWCSRes[j].sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, _PWCSRes[j].cst(), vRes[j], 1. );
          continue;
        }
        // polyhedral cuts for superposition relaxation
        _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, _PWCSRes[j].uest(0), _PWCSRes[j].sdep(), 1 );
        _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, _PWCSRes[j].oest(0), _PWCSRes[j].sdep(), 0 );
        if( options.PWCSHADOW ){
          _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, _PWCSRes[j].uest(1), _PWCSRes[j].sdep(), 1 );
          _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, _PWCSRes[j].oest(1), _PWCSRes[j].sdep(), 0 );
        }
      }
      break;

    // McCormick relaxation with piecewise-constant superposition bounds and subgradient calculated at mid-point
    case Options::MCPWCS:
      // define auxiliary variables 
      for( unsigned i=0; i<nin; ++i ){
        auto& PWCSVar_i = _MCPWCSVar[i].I();
        _append_PWCScuts( img, pOp, vVar[i], _POLPWCSAux[i], PWCSVar_i.uest()[i], PWCSVar_i.oest()[i] );          
      }     
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
        auto& PWCSRes_j = _MCPWCSRes[j].I();
#ifdef MC__FFEXTERN_DEBUG
        std::cerr << "_PWCSRes[" << j << "] in " << PWCSRes_j << std::endl;
        std::cerr << "_MCPWCSRes[" << j << "] in " << _MCPWCSRes[j] << std::endl;
#endif
        // constant superposition relaxation
        if( PWCSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWCSRes_j.cst(), vRes[j], 1. );
          continue;
        }
        // polyhedral cuts for superposition relaxation
        _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, PWCSRes_j.uest(0), PWCSRes_j.sdep(), 1 );          
        _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, PWCSRes_j.oest(0), PWCSRes_j.sdep(), 0 );          
        if( options.PWCSHADOW ){
          _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, PWCSRes_j.uest(1), PWCSRes_j.sdep(), 1 );
          _append_PWCScuts( img, pOp, vRes[j], _POLPWCSAux, PWCSRes_j.oest(1), PWCSRes_j.sdep(), 0 );
        }
        // polyhedral cut generation for MC  
        _append_MCcuts( img, pOp, vRes[j], _MCPWCSRes[j], vVar, _MCPWCSVar );
      }
      break;
      
    // Piecewise-linear superposition models
    case Options::PWLS:
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__FFEXTERN_DEBUG
        std::cerr << "PWLSRes[" << j << "]: " << _PWLSRes[j] << std::endl;
#endif
        // constant superposition relaxation
        if( _PWLSRes[j].sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, _PWLSRes[j].cst(), vRes[j], 1. );
          continue;
        }
        // polyhedral cuts for superposition relaxation
        for( auto& summand : _PWLSRes[j].uest(0) ) summand.reduce( true, options.PWLMAX );
        _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, _PWLSRes[j].uest(0), _PWLSRes[j].sdep(), true  );
        for( auto& summand : _PWLSRes[j].oest(0) ) summand.reduce( false, options.PWLMAX );
        _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, _PWLSRes[j].oest(0), _PWLSRes[j].sdep(), false );
        if( options.PWLSHADOW ){
          for( auto& summand : _PWLSRes[j].uest(1) ) summand.reduce( true, options.PWLMAX );
          _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, _PWLSRes[j].uest(1), _PWLSRes[j].sdep(), true  );
          for( auto& summand : _PWLSRes[j].oest(1) ) summand.reduce( false, options.PWLMAX );
          _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, _PWLSRes[j].oest(1), _PWLSRes[j].sdep(), false );
        }
      }
      break;

    // McCormick relaxation with piecewise-linear superposition bounds and subgradient calculated at mid-point
    case Options::MCPWLS:
      // add polyhedral cuts
      for( unsigned j=0; j<nout; ++j ){
#ifdef MC__FFEXTERN_DEBUG
        std::cerr << "MCPWLSRes[" << j << "]: " << _MCPWLSRes[j] << std::endl;
#endif
        auto& PWLSRes_j = _MCPWLSRes[j].I();
        // constant superposition relaxation
        if( PWLSRes_j.sdep().empty() ){
          *img->add_cut( pOp, PolCut<T>::EQ, PWLSRes_j.cst(), vRes[j], 1. );
          continue;
        }
        // polyhedral cuts for superposition relaxation
        for( auto& summand : PWLSRes_j.uest(0) ) summand.reduce( true, options.PWLMAX );
        _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, PWLSRes_j.uest(0), PWLSRes_j.sdep(), true  );
        for( auto& summand : PWLSRes_j.oest(0) ) summand.reduce( false, options.PWLMAX );
        _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, PWLSRes_j.oest(0), PWLSRes_j.sdep(), false );
        if( options.PWLSHADOW ){
          for( auto& summand : PWLSRes_j.uest(1) ) summand.reduce( true, options.PWLMAX );
          _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, PWLSRes_j.uest(1), PWLSRes_j.sdep(), true  );
          for( auto& summand : PWLSRes_j.oest(1) ) summand.reduce( false, options.PWLMAX );
          _append_PWLScuts( img, pOp, vRes[j], vVar, _POLPWLSAux, _DXPWLSAux, _DYPWLSAux, PWLSRes_j.oest(1), PWLSRes_j.sdep(), false );
        }
        // polyhedral cut generation for MC  
        _append_MCcuts( img, pOp, vRes[j], _MCPWLSRes[j], vVar, _MCPWLSVar );
      }
      break;
    } // end switch
  }

#ifdef MC__FFEXTERN_DEBUG
  std::cerr << *img;
  {int dum; std::cout << "PAUSED, ENTER 1"; std::cin >> dum;}
#endif
}

template< typename T, typename O >
template< typename U >
inline void
FFEXTERN<T,O>::_append_MCcuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, McCormick<U> const& MCRes,
  PolVar<T> const* vVar, std::vector<McCormick<U>> const& vMCVar )
const
{
  double rhs1 = -MCRes.cv(),
         rhs2 = -MCRes.cc();
  for( unsigned i=0; i<vMCVar.size(); ++i ){ 
    rhs1 += MCRes.cvsub(i) * vMCVar[i].cv();
    rhs2 += MCRes.ccsub(i) * vMCVar[i].cc();
  }
  img->add_cut( pOp, PolCut<T>::LE, rhs1, vMCVar.size(), vVar, MCRes.cvsub(), Res, -1. );
  img->add_cut( pOp, PolCut<T>::GE, rhs2, vMCVar.size(), vVar, MCRes.ccsub(), Res, -1. );
}

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_append_PWCScuts
( PolImg<T>* img, FFOp* pOp, PolVar<T>& Var, std::vector<PolVar<T>>& vAux,
  PWCU const& uest, PWCU const& oest )
const
{
#ifdef MC__FFEXTERN_CHECK
  assert( uest.yL().size() == options.PWCDIV
       && oest.yU().size() == options.PWCDIV );
#endif
  vAux.resize( options.PWCDIV );
  for( unsigned k=0; k<options.PWCDIV; ++k )
    vAux[k].set( img, Op<T>::zeroone(), options.PWCREL? true: false );
  // sum of auxiliaries equal to 1
  img->add_cut( pOp, PolCut<T>::EQ, 1., options.PWCDIV, vAux.data(), 1. );

  // link auxiliaries to independent variables
  img->add_cut( pOp, PolCut<T>::LE, 0., options.PWCDIV, vAux.data(), uest.yL().data(), Var, -1. );
  img->add_cut( pOp, PolCut<T>::GE, 0., options.PWCDIV, vAux.data(), oest.yU().data(), Var, -1. );
} 

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_append_PWCScuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, std::vector<std::vector<PolVar<T>>> const& vAux,
  std::vector<PWCU> const& est, std::set<unsigned int> const& dep, bool const under )
const
{
  if( est.empty() ) return;

  // polyhedral cuts for superposition relaxation
  if( under ){
    auto cut = *img->add_cut( pOp, PolCut<T>::LE, 0., Res, -1. );
    for( auto const& i : dep ){
#ifdef MC__FFEXTERN_CHECK
      assert( est[i].yL().size() == options.PWCDIV );
#endif
      cut->append( options.PWCDIV, vAux[i].data(), est[i].yL().data() );
    }
  }
  else{
    auto cut = *img->add_cut( pOp, PolCut<T>::GE, 0., Res, -1. );
    for( auto const& i : dep ){
#ifdef MC__FFEXTERN_CHECK
      assert( est[i].yU().size() == options.PWCDIV );
#endif
      cut->append( options.PWCDIV, vAux[i].data(), est[i].yU().data() );
    }
  }
}

template< typename T, typename O >
inline void
FFEXTERN<T,O>::_append_PWLScuts
( PolImg<T>* img, FFOp* pOp, PolVar<T> const& Res, PolVar<T>* vVar,
  std::vector<PolVar<T>>& vAux, std::vector<double>& Xk, std::vector<double>& Yk, 
  std::vector<PWLU> const& est, std::set<unsigned int> const& dep, bool const under )
const
{
  if( est.empty() ) return;
  
  for( auto const& i : dep ){
    vAux[i].set( img, T(est.at(i).l(),est.at(i).u()), true ); // continuous auxiliary

    // Generate breakpoints from segment-based representation in mc::PWLU
    auto const& est_i = est.at(i);
    size_t const NK = est_i.dx().size()+1;
    Xk.resize( NK ); Xk[0] = est_i.xL();
    Yk.resize( NK ); Yk[0] = est_i.yL();
    auto idx = est_i.dx().cbegin(), idy = est_i.dy().cbegin();
    for( unsigned k=1; k<NK; ++k, ++idx, ++idy ){
      Xk[k] = Xk[k-1] + *idx;
      Yk[k] = Yk[k-1] + *idx * *idy;
    }

    img->add_semilinear_cuts( pOp, NK, vVar[i], Xk.data(), vAux[i], Yk.data(), mc::PolCut<T>::EQ, options.PWLREL?0:1 );
  }
     
  img->add_cut( pOp, (under? PolCut<T>::LE: PolCut<T>::GE), 0., dep, vAux.data(), 1., Res, -1. );
}

} // end namespace mc

#endif
