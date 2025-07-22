#ifndef MC__FFEXTERN_HPP
#define MC__FFEXTERN_HPP

#define MC__FFEXTERN_DEBUG
#define MC__FFEXTERN_CHECK
#define MC__FFDAG_DEBUG
#define MC__FFDAG_CHECK

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

//! @brief C++ class for evaluation and relaxation of an operation defined by an expression tree
////////////////////////////////////////////////////////////////////////
//! mc::DAG is a C++ class for evaluation and relaxation of an
//! operation defined by an expression tree using MC++
////////////////////////////////////////////////////////////////////////
template <typename T> 
class DAG
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
#ifdef MC__DAG_DEBUG
      std::cout << "DAG:: Original DAG: " <<  dag << std::endl;
      std::cout << "DAG:: Copied DAG:   " << _dag << std::endl;
#endif
      _dag->insert( dag, varin,  _varin  );
      _dag->insert( dag, varout, _varout );
      _codelist.clear();
#ifdef MC__DAG_DEBUG
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
#ifdef MC__DAG_DEBUG
      std::cout << "DAG:: Work array: " << flag << " passes\n";
      for( unsigned i=0; i<_codelist.len_tap-_codelist.len_wrk; i++ )
        std::cout << "wk[" << i << "] = " << wk[i] << std::endl;
      for( unsigned i=0; i<_varin.size(); i++ )
        std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif
      return( flag<0? false: true );
    }

public:

  //! @brief DAG options
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
  DAG
    ():
    _dag( nullptr )
    {}

  //! @brief Data constructor
  DAG
    ( FFGraph* dag, std::vector<FFVar> const& varin, std::vector<FFVar> const& varout ):
    _dag( nullptr )
    {
      _set( dag, varin, varout );
    }

  //! @brief Copy constructor
  DAG
    ( DAG const& other ):
    _dag( nullptr )
    {
      _set( other._dag, other._varin, other._varout );
      options = other.options;
    }

  virtual ~DAG() 
    {
      delete _dag;
    }

  //! @brief Set DAG data
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

  //! @brief Differentiate DAG
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
          _POLVar[i].set( _POLEnv, _ptrObj->varin()[i], vVar[i].range() );
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
#ifdef MC__FFEXTERN_DEBUG
          std::cerr << "_POLRes[" << j << "]: " << _POLRes[j] << std::endl;
#endif
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
    vAux[k].set( img, Op<T>::zeroone(), !options.PWCREL );
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

    img->add_semilinear_cuts( pOp, NK, vVar[i], Xk.data(), vAux[i], Yk.data(), mc::PolCut<T>::EQ, options.PWLREL );
  }
     
  img->add_cut( pOp, (under? PolCut<T>::LE: PolCut<T>::GE), 0., dep, vAux.data(), 1., Res, -1. );
}

//! @brief C++ class defining expression tree as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFDAG is a C++ class for defining expression trees as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template< typename T >
class FFDAG
////////////////////////////////////////////////////////////////////////
: public FFEXTERN<T,DAG<T>>
{

protected:

  using FFEXTERN<T,DAG<T>>::_ptrObj;
  using FFEXTERN<T,DAG<T>>::_ownObj;
    
  // set the object and related operation in DAG
  FFVar** _set
    ( size_t const nVar, FFVar const* pVar, DAG<T>* pDAG, int policy )
    {
#ifdef MC__FFDAG_CHECK
      assert( nVar == pDAG->nin() );
#endif
      if( _ownObj && _ptrObj ) delete _ptrObj;
      _ownObj = ( policy>0? true: false ); //copy;
      //_ownObj = true;
      this->owndata = false;
      this->data = _ptrObj = pDAG;

      FFVar** ppRes = this->insert_external_operation( *this, pDAG->nout(), nVar, pVar );

      _ownObj = false;
      FFOp* pOp = (*ppRes)->opdef().first;
      if( policy > 0 )
        _ptrObj = static_cast<FFDAG*>(pOp)->_ptrObj; // set pointer to DAG copy
      else if( policy < 0 )
        static_cast<FFDAG*>(pOp)->_ownObj = true; // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFDAG_TRACE
      std::cerr << "DAG operation address: " << this << std::endl;
      std::cerr << "DAG address in DAG: " << _ptrObj << std::endl;
#endif
      return ppRes;
    }

public:

  //! @brief Enumeration type for copy policy of DAG object
  enum POLICY_TYPE{
    SHALLOW=0,  //!< Shallow copy of DAG object in FFGraph (without ownership)
    COPY=1,     //!< Deep copy of DAG object in FFGraph (with ownership)
    TRANSFER=-1 //!< Shallow copy of DAG object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFDAG
    ()
    : FFEXTERN<T,DAG<T>>()
    {}

  // Destructor
  virtual ~FFDAG
    ()
    {}

  // Copy constructor
  FFDAG
    ( FFDAG<T> const& Op )
    : FFEXTERN<T,DAG<T>>( Op )
    {}

  // Define operation
  FFVar** operator()
    ( std::vector<FFVar> const& vVar, DAG<T>* pDAG, int policy=COPY )
    {
      return _set( vVar.size(), vVar.data(), pDAG, policy );
    }

  FFVar& operator()
    ( size_t const idep, std::vector<FFVar> const& vVar, DAG<T>* pDAG, int policy=COPY )
    {
#ifdef MC__FFDAG_CHECK
      assert( idep < pDAG->nout() );
#endif
      return *(_set( vVar.size(), vVar.data(), pDAG, policy )[idep]);
    }

  // Define operation
  FFVar** operator()
    ( size_t const nVar, FFVar const* pVar, DAG<T>* pDAG, int policy=COPY )
    {
      return _set( nVar, pVar, pDAG, policy );
    }

  FFVar& operator()
    ( size_t const idep, size_t const nVar, FFVar const* pVar, DAG<T>* pDAG, int policy=COPY )
    {
#ifdef MC__FFDAG_CHECK
      assert( idep < pDAG->nout() );
#endif
      return *(_set( nVar, pVar, pDAG, policy )[idep]);
    }

  DAG<T>* pDAG
    ()
    const
    {
      //std::cerr << "DAG address retreived: " << _ptrObj << std::endl;
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

      throw std::runtime_error( "FFDAG::feval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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

      throw std::runtime_error( "FFDAG::reval: **ERROR** No evaluation method with type"+std::string(idU.name())+"\n" );
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
    { std::ostringstream oss; oss << this->data; return "DAG[" + oss.str() + "]"; }

  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

template< typename T >
inline void
FFDAG<T>::eval
( size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: FFDep\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  vRes[0] = 0;
  for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
  vRes[0].update( FFDep::TYPE::N );
  for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
}

template< typename T >
inline void
FFDAG<T>::eval
( size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: FFVar\n";
  std::cerr << "DAG operation address: " << this    << std::endl;
  std::cerr << "DAG address in DAG: "    << _ptrObj << std::endl;
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  FFVar** ppRes = this->insert_external_operation( *this, nRes, nVar, vVar );
  for( unsigned j=0; j<nRes; ++j )
    vRes[j] = *(ppRes[j]);
}

template< typename T >
template< typename U >
inline void
FFDAG<T>::eval
( size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: " << typeid( vRes[0] ).name() << " (generic)\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  _ptrObj->eval( vVar, vRes );
}

template< typename T >
inline void
FFDAG<T>::eval
( size_t const nRes, fadbad::F<FFVar>* vRes, size_t const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: fadbad::F<FFVar>\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  FFVar const*const* vResVal = this->insert_external_operation( *this, nRes, nVar, vVarVal.data() );

  auto&& sdervarout = _ptrObj->diff();
  DAG<T>* pDAGDer = new DAG<T>( _ptrObj->dag(), _ptrObj->varin(), std::get<2>(sdervarout) );
  FFDAG<T> ResDer;
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
      vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
    }
  }
}

template< typename T >
inline void
FFDAG<T>::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::deriv: FFVar\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  auto&& sdervarout = _ptrObj->diff();
  DAG<T>* pDAGDer = new DAG<T>( _ptrObj->dag(), _ptrObj->varin(), std::get<2>(sdervarout) );
  FFDAG<T> ResDer;
  FFVar const*const* vResDer = ResDer._set( nVar, vVar, pDAGDer, TRANSFER );

  // set all derivatives to zero
  for( unsigned k=0; k<nRes; ++k )
    for( unsigned i=0; i<nVar; ++i )
        vDer[k][i] = 0;

  // copy nonzero derivatives from sparse Jacobian
  for( unsigned int ie=0; ie<std::get<0>(sdervarout).size(); ++ie ){
    auto const& k = std::get<0>(sdervarout)[ie];
    auto const& i = std::get<1>(sdervarout)[ie];
    vDer[k][i] = *vResDer[k+nRes*i];
  }
}

template< typename T >
inline void
FFDAG<T>::eval
( size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: SLiftVar\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  vVar->env()->lift( nRes, vRes, nVar, vVar );
}

template< typename T >
inline void
FFDAG<T>::eval
( size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
  unsigned const* mVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::eval: PolVar\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag = vVar[0].var().dag();
#ifdef MC__FFDAG_CHECK
  assert( img && dag );
#endif
  FFVar** ppRes = dag->curOp()->varout.data(); // ACCOUNT FOR MULTIPLE OUTPUTS
#ifdef MC__FFDAG_CHECK
  assert( nRes == dag->curOp()->varout.size() );
#endif

  this->_resize_relax( img );
  this->_propagate_relax( img, ppRes, vRes, vVar );
}

template< typename T >
inline bool
FFDAG<T>::reval
( size_t const nRes, PolVar<T> const* vRes, size_t const nVar, PolVar<T>* vVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::reval: PolVar\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop = vVar[0].var().opdef().first;
#ifdef MC__FFDAG_CHECK
  assert( img && pop );
#endif

  this->_backpropagate_relax( img, pop, vRes, vVar );
  return true;
}

template< typename T >
inline bool
FFDAG<T>::reval
( size_t const nRes, T const* vRes, size_t const nVar, T* vVar )
const
{
#ifdef MC__FFDAG_TRACE
  std::cout << "FFDAG::reval: T\n";
#endif
#ifdef MC__FFDAG_CHECK
  assert( _ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin() );
#endif

  return _ptrObj->reval( vVar, vRes );
}

} // end namespace mc

#endif
