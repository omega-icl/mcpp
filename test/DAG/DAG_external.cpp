#undef  MC__FFUNC_CPU_EVAL
#undef  MC__FFUNC_EXTERN_DEBUG
#undef  MC__FFUNC_SBAD_DEBUG
#define MC__SELIM_DEBUG_PROCESS
#define MC__SELIM_DEBUG_MIP
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "mclapack.hpp"
#include "ffunc.hpp"
#include "slift.hpp"
#include "selim.hpp"
#include "fflin.hpp"
#include "ffspol.hpp"
#include "ffmlp.hpp"
#include "ffvect.hpp"
#include "ffdagext.hpp"
#include "ffcustom.hpp"

#if defined( MC__USE_PROFIL )
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#elif defined( MC__USE_FILIB )
 #include "mcfilib.hpp"
 typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
#elif defined( MC__USE_BOOST )
 #include "mcboost.hpp"
 typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
 typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
 typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
 typedef boost::numeric::interval<double,T_boost_policy> I;
#else
 #include "interval.hpp"
 typedef mc::Interval I;
#endif

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

#include "supmodel.hpp"
#include "pwcu.hpp"
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;
typedef mc::McCormick<PWCSV> MCPWCSV;

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

#include "polimage.hpp"
typedef mc::PolVar<I> POLV;


////////////////////////////////////////////////////////////////////////
// EXTERNAL OPERATIONS
////////////////////////////////////////////////////////////////////////
namespace mc
{

class FFnorm2
: public FFOp
{
public:
  // Constructors
  FFnorm2
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      return **insert_external_operation( *this, 1, nVar, pVar );
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( I ) )
        return eval( nRes, static_cast<I*>(vRes), nVar, static_cast<I const*>(vVar), mVar );
      else if( idU == typeid( McCormick<I> ) )
        return eval( nRes, static_cast<McCormick<I>*>(vRes), nVar, static_cast<McCormick<I> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFnorm2 with type"+std::string(idU.name())+"\n" );
    }

  template< typename T >
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      std::cout << "NORM2 generic instantiation\n";
      switch( nVar ){
        case 0: vRes[0] = T( 0. ); break;
        case 1: vRes[0] = vVar[0]; break;
        default: vRes[0] = Op<T>::sqr( vVar[0] );
                 for( unsigned i=1; i<nVar; ++i ) vRes[0] += Op<T>::sqr( vVar[i] );
                 vRes[0] = Op<T>::sqrt( vRes[0] ); break;
      }
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      std::cout << "NORM2 FFVar instantiation\n";
      vRes[0] = operator()( nVar, vVar );
    }
    
  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar && nRes == 1 );
      std::cout << "NORM2 SLiftVar instantiation\n";
      vVar->env()->lift( nRes, vRes, nVar, vVar );
    }

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const
    {
      assert( nRes == 1 );
      std::cout << "NORM2 FFVar differentiation\n";
      fadbad::F<FFVar> vFVar[nVar], vFRes[nRes];
      for( unsigned i=0; i<nVar; ++i ){
        vFVar[i] = vVar[i];
        vFVar[i].diff( i, nVar );
      }
      eval( nRes, vFRes, nVar, vFVar, nullptr );
      for( unsigned j=0; j<nRes; ++j )
        for( unsigned i=0; i<nVar; ++i )
          vDer[j][i] = vFRes[j].d(i);
    }

  // Properties
  std::string name
    ()
    const
    { return "NORM2"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return true; }
};

///////////////////////////////////////////////////////////////////////////////

class FFnorm12
: public FFOp
{
public:
  // Constructors
  FFnorm12
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar )
    const
    {
      return *(insert_external_operation( *this, 2, nVar, pVar )[idep]);
    }
  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      return insert_external_operation( *this, 2, nVar, pVar );
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( I ) )
        return eval( nRes, static_cast<I*>(vRes), nVar, static_cast<I const*>(vVar), mVar );
      else if( idU == typeid( McCormick<I> ) )
        return eval( nRes, static_cast<McCormick<I>*>(vRes), nVar, static_cast<McCormick<I> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFnorm12 with type"+std::string(idU.name())+"\n" );
    }

  template< typename T >
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 2 );
      std::cout << "NORM12 generic instantiation\n";
      switch( nVar ){
        case 0: vRes[0] = T( 0. );
                vRes[1] = T( 0. );
                break;
        case 1: vRes[0] = Op<T>::fabs( vVar[0] );
                vRes[1] = Op<T>::fabs( vVar[0] );
                break;
        default: vRes[0] = Op<T>::sqr( vVar[0] );
                 for( unsigned i=1; i<nVar; ++i ) vRes[0] += Op<T>::sqr( vVar[i] );
                 vRes[0] = Op<T>::sqrt( vRes[0] );
                 vRes[1] = Op<T>::fabs( vVar[0] );
                 for( unsigned i=1; i<nVar; ++i ) vRes[1] += Op<T>::fabs( vVar[i] );
                 break;
      }
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      std::cout << "NORM12 FFVar instantiation\n";
      FFVar** ppRes = operator()( nVar, vVar );
      for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
    }
    
  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar && nRes == 2 );
      std::cout << "NORM12 SLiftVar instantiation\n";
      vVar->env()->lift( nRes, vRes, nVar, vVar );
    }

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const
    {
      assert( nRes == 2 );
      std::cout << "NORM12 FFVar differentiation\n";
      fadbad::F<FFVar> vFVar[nVar], vFRes[nRes];
      for( unsigned i=0; i<nVar; ++i ){
        vFVar[i] = vVar[i];
        vFVar[i].diff( i, nVar );
      }
      eval( nRes, vFRes, nVar, vFVar, nullptr );
      for( unsigned j=0; j<nRes; ++j )
        for( unsigned i=0; i<nVar; ++i )
          vDer[j][i] = vFRes[j].d(i);
    }

  // Properties
  std::string name
    ()
    const
    { return "NORM12"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return true; }
};

///////////////////////////////////////////////////////////////////////////////

class FFxlog
: public FFOp
{
public:
  // Constructors
  FFxlog
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( FFVar const& Var )
    const
    {
      return *(insert_external_operation( *this, 1, Var )[0]);
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( I ) )
        return eval( nRes, static_cast<I*>(vRes), nVar, static_cast<I const*>(vVar), mVar );
      else if( idU == typeid( McCormick<I> ) )
        return eval( nRes, static_cast<McCormick<I>*>(vRes), nVar, static_cast<McCormick<I> const*>(vVar), mVar );
      else if( idU == typeid( fadbad::F<FFVar> ) )
        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( PolVar<I> ) )
        return eval( nRes, static_cast<PolVar<I>*>(vRes), nVar, static_cast<PolVar<I> const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFxlog with type"+std::string(idU.name())+"\n" );
    }

  template <typename T>
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog generic instantiation\n"; 
      vRes[0] = vVar[0] * Op<T>::log( vVar[0] );
    }

  template <typename T>
  void eval
    ( unsigned const nRes, McCormick<T>* vRes, unsigned const nVar, McCormick<T> const* vVar,
      unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog McCormick instantiation\n"; 
      vRes[0] = xlog( vVar[0] );
    }

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog fadbad::F<FFVar> instantiation\n"; 
      vRes[0] = operator()( vVar[0].val() );
      if( !vVar[0].depend() ) return;
      FFVar dxlog( log( vVar[0].val()) + 1 );
      vRes[0].setDepend( vVar[0] );
      for( unsigned int i=0; i<vRes[0].size(); ++i )
        vRes[0][i] = dxlog * vVar[0][i];
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      vRes[0] = operator()( vVar[0] );
    }

  template <typename T>
  void eval
    ( unsigned const nRes, PolVar<T>* vRes, unsigned const nVar, PolVar<T> const* vVar,
      unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog Polyhedral image instantiation\n"; 
      PolImg<T>* img = vVar[0].image();
      FFBase* dag = vVar[0].var().dag();
      assert( img && dag );
      FFVar* pRes = dag->curOp()->varout[0];
      T TRes = Op<I>::xlog( vVar[0].range() );
      vRes[0].set( img, *pRes, TRes );
      // vRes[0] = xlog( vVar[0] );
    }

  virtual bool reval
    ( std::type_info const& idU, unsigned const nRes, void const* vRes, unsigned const nVar, void* vVar )
    const
    {
      if( idU == typeid( PolVar<I> ) )
        return reval( nRes, static_cast<PolVar<I> const*>(vRes), nVar, static_cast<PolVar<I>*>(vVar) );

      throw std::runtime_error( "Error: No evaluation method for FFXlog with type"+std::string(idU.name())+"\n" );
    }

  template <typename T>
  bool reval
    ( unsigned const nRes, PolVar<T> const* vRes, unsigned const nVar, PolVar<T>* vVar )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog Polyhedral image generation\n"; 
      PolImg<T>* img = vVar[0].image();
      FFBase* dag = vVar[0].var().dag();
      FFOp* pop = vVar[0].var().opdef().first;
      assert( img && dag && pop );
      struct loc{ static std::pair<double,double> xlog
        ( const double x, const double*rusr, const int*iusr )
        { return std::make_pair( mc::xlog(x), std::log(x)+1. ); }
      };
      img->add_semilinear_cuts( pop, vVar[0], Op<T>::l(vVar[0].range()), Op<T>::u(vVar[0].range()),
        vRes[0], PolCut<T>::LE, loc::xlog );
      img->add_sandwich_cuts( pop, vVar[0], Op<T>::l(vVar[0].range()), Op<T>::u(vVar[0].range()),
        vRes[0], Op<T>::l(vRes[0].range()), Op<T>::u(vRes[0].range()), PolCut<T>::GE, loc::xlog );
      return true;
    }

  virtual void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const
    {
      assert( nVar == 1 && nRes == 1 );
      std::cout << "xlog FFVar differentiation\n";
      vDer[0][0] = log( vVar[0] ) + 1;
    }

  // Properties
  std::string name
    ()
    const
    { return "XLOG EXT"; }
};

///////////////////////////////////////////////////////////////////////////////

class FFdet
: public FFOp
{
public:
  // Constructors
  FFdet
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      return **insert_external_operation( *this, 1, nVar, pVar );
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( I ) )
        return eval( nRes, static_cast<I*>(vRes), nVar, static_cast<I const*>(vVar), mVar );
      else if( idU == typeid( McCormick<I> ) )
        return eval( nRes, static_cast<McCormick<I>*>(vRes), nVar, static_cast<McCormick<I> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFArrh with type"+std::string(idU.name())+"\n" );
    }

  template< typename T >
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      std::cout << "Det generic instantiation\n"; 
      const unsigned nDim = std::sqrt(nVar);
      switch( nDim ){
        case 0:  vRes[0] = T( 0. ); break;
	case 1:  vRes[0] = vVar[0]; break;
        default: vRes[0] = FFBase::det( nDim, vVar ); break;
      }
    }
  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      std::cout << "Det double instantiation\n"; 
      const unsigned nDim = std::sqrt(nVar);
      CPPL::dgematrix Amat( nDim, nDim );
      for( unsigned i=0; i<nDim; ++i )
        for( unsigned j=0; j<nDim; ++j )
          Amat(i,j) = vVar[i+j*nDim];
      if( dgeqrf( Amat, vRes[0] ) )
        throw FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    }
  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      vRes[0] = operator()( nVar, vVar );
    }
  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      vRes[0] = 0;
      for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
      switch( nVar ){
        case 0:
        case 1:  vRes[0].update( FFDep::TYPE::L ); break;
        case 2:  vRes[0].update( FFDep::TYPE::Q ); break;
        default: vRes[0].update( FFDep::TYPE::P ); break;
      }
    }

  // Properties
  std::string name
    ()
    const
    { return "DET"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

///////////////////////////////////////////////////////////////////////////////

struct FFDOptBase
{
  // Vector of atom matrices
  static std::vector< CPPL::dsymatrix > _A;

  // Read atom matrices from file
  static unsigned read
    ( unsigned const dim, std::string filename, bool const disp=false )
    {
      CPPL::dsymatrix Ai( dim );
      _A.clear();
      std::ifstream file( filename );
      if( !file ) throw std::runtime_error("Error: Could not open input file\n");
      std::string line;
      unsigned i = 0;
      bool empty = false;
      while( std::getline( file, line ) ){
        std::istringstream iss( line );
        for( unsigned j=0; j<dim; j++ ){
          if( !(iss >> Ai(i,j) ) ){
            if( j ) throw std::runtime_error("Error: Could not read input file\n");
            empty = true;
            break;
          }
          if( disp ) std::cout << "reading (" << i << "," << j << "): " << Ai(i,j) << std::endl;
        }
        i++;
        if( empty ){
          if( disp ) std::cout << "Atomic matrix #" << _A.size() << ":" << std::endl << Ai;
          _A.push_back( Ai );
          i = 0;
          empty = false;
        }
        if( i > dim ) throw std::runtime_error("Error: Could not read input file\n");
      }
      if( i ) _A.push_back( Ai );
      return _A.size();
    }
};

inline std::vector< CPPL::dsymatrix > FFDOptBase::_A;

class FFDOpt
: public FFOp,
  public FFDOptBase
{
public:
  // Constructors
  FFDOpt
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      return **insert_external_operation( *this, 1, nVar, pVar );
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
        return eval( nRes, static_cast<fadbad::F<FFVar>*>(vRes), nVar, static_cast<fadbad::F<FFVar> const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFDOpt with type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const
    {
      std::cout << "FFDOpt::eval: double\n"; 
      assert( nRes == 1 && nVar == _A.size() && _A.begin() != _A.end() );
      CPPL::dsymatrix Amat( _A[0].n );
      Amat.zero();
      for( unsigned i=0; i<nVar; ++i )
        if( !i ) Amat  = vVar[0] * _A[0];
        else     Amat += vVar[i] * _A[i];
      //std::cout << Amat;
      if( dgeqrf( Amat.to_dgematrix(), vRes[0] ) )
        throw FFBase::Exceptions( FFBase::Exceptions::EXTERN );
      vRes[0] = std::log( vRes[0] );
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      std::cout << "FFDOpt::eval: FFVar\n"; 
      vRes[0] = operator()( nVar, vVar );
    }

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == 1 );
      std::cout << "FFDOpt::eval: FFDep\n"; 
      vRes[0] = 0;
      for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
      vRes[0].update( FFDep::TYPE::N );
    }

  void eval
    ( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
      unsigned const* mVar )
    const;

  void deriv
    ( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
    const;

  // Properties
  std::string name
    ()
    const
    { return "DOPT"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

class FFDOptGrad
: public FFOp,
  public FFDOptBase
{
public:
  // Constructors
  FFDOptGrad
    ()
    : FFOp( EXTERN )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const idep, unsigned const nVar, FFVar const* pVar )
    const
    {
      return *(insert_external_operation( *this, nVar, nVar, pVar )[idep]);
    }
  FFVar** operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      return insert_external_operation( *this, nVar, nVar, pVar );
    }

  // Evaluation overloads
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

      throw std::runtime_error( "Error: No evaluation method for FFDOptGrad with type"+std::string(idU.name())+"\n" );
    }

  void eval
    ( unsigned const nRes, double* vRes, unsigned const nVar, double const* vVar, unsigned const* mVar )
    const
    {
      std::cout << "FFDOptGrad::eval: double\n"; 
      assert( nRes == nVar && nVar == _A.size() && _A.begin() != _A.end() );
      CPPL::dsymatrix Amat( _A[0].n );
      Amat.zero();
      for( unsigned i=0; i<nVar; ++i )
        if( !i ) Amat  = vVar[0] * _A[0];
        else     Amat += vVar[i] * _A[i];
      //std::cout << Amat;
      // Perform LDL' decomposition
      CPPL::dgematrix Lmat;
      std::vector<int> IPIV;
      if( dsytrf( Amat, Lmat, IPIV ) )
        throw FFBase::Exceptions( FFBase::Exceptions::EXTERN );
      CPPL::dgematrix Xmat;
      for( unsigned i=0; i<nVar; ++i ){
        if( dsytrs( Lmat, IPIV, _A[i].to_dgematrix(), Xmat ) )
          throw FFBase::Exceptions( FFBase::Exceptions::EXTERN );
        //std::cout << Xmat;
        for( int j=0; j<Xmat.n; ++j )
          if( !j ) vRes[i]  = Xmat(0,0);
          else     vRes[i] += Xmat(j,j);
      }
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == nVar );
      std::cout << "FFDOptGrad::eval: FFVar\n"; 
      FFVar** ppRes = operator()( nVar, vVar );
      for( unsigned j=0; j<nRes; ++j ) vRes[j] = *(ppRes[j]);
    }

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const
    {
      assert( nRes == nVar );
      std::cout << "FFDOpt::eval: FFDep\n"; 
      vRes[0] = 0;
      for( unsigned i=0; i<nVar; ++i ) vRes[0] += vVar[i];
      vRes[0].update( FFDep::TYPE::N );
      for( unsigned j=1; j<nRes; ++j ) vRes[j] = vRes[0];
    }

  // Properties
  std::string name
    ()
    const
    { return "DOPTGRAD"; }
  //! @brief Return whether or not operation is commutative
  bool commutative
    ()
    const
    { return false; }
};

inline void
FFDOpt::eval
( unsigned const nRes, fadbad::F<FFVar>* vRes, unsigned const nVar, fadbad::F<FFVar> const* vVar,
  unsigned const* mVar )
const
{
  assert( nRes == 1 && nVar == _A.size() && _A.begin() != _A.end() );
  std::cout << "FFDOpt::eval: fadbad::F<FFVar>\n";
  std::vector<FFVar> vVarVal( nVar );
  for( unsigned i=0; i<nVar; ++i )
    vVarVal[i] = vVar[i].val();
  vRes[0] = operator()( nVar, vVarVal.data() );
  FFDOptGrad DOptGrad;
  for( unsigned i=0; i<nVar; ++i )
    vRes[0].setDepend( vVar[i] );
  for( unsigned j=0; j<vRes[0].size(); ++j )
    for( unsigned i=0; i<nVar; ++i )
      if( !i ) vRes[0][j]  = DOptGrad( 0, nVar, vVarVal.data() ) * vVar[0][j];
      else     vRes[0][j] += DOptGrad( i, nVar, vVarVal.data() ) * vVar[i][j];
}

inline void
FFDOpt::deriv
( unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar, FFVar** vDer )
const
{
  assert( nRes == 1 && nVar == _A.size() && _A.begin() != _A.end() );
  std::cout << "FFDOpt::deriv: FFVar\n";
  FFDOptGrad DOptGrad;
  for( unsigned i=0; i<nVar; ++i )
    vDer[0][i] = DOptGrad( i, nVar, vVar );
}

///////////////////////////////////////////////////////////////////////////////

class FFArrh
: public FFOp
{
public:
  // Constructors
  FFArrh
    ()
    : FFOp( EXTERN )
    {}

  FFArrh
    ( FFArrh const& Op )
    : FFOp( Op )
    {}

  // Destructor
  virtual ~FFArrh
    ()
    {
//      std::cout << "FFArrh destructor invoked\n";
//      if( owndata ) delete static_cast<double*>( data );
    }

  // Define operation
  FFVar& operator()
    ( FFVar const& Var, double& r, bool const copy=false )
    const
    {
      data = &r; // this is assuming r isn't going out of scope
//      info = ID;
//      std::cout << "data: " << data << std::endl;
      return *(insert_external_operation( *this, 1, Var )[0]);
//      FFVar* pRes = insert_external_operation( *this, 1, Var )[0];
//      FFOp* pOp = pRes->opdef().first;
//      if( copy && !pOp->owndata ){
//        double* pr = new double( r );
//	auto opins = update_data( pOp, pr, true );
//	assert( opins.second );
//	pOp = opins.first;
////        data = pODE = static_cast<ODESLVS_CVODES<ExtOps...>*>( pOp->data );
//        data = pOp->data;
//      }
//      return *pRes;
    }

  // Evaluation overloads
  virtual void feval
    ( std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
      void const* vVar, unsigned const* mVar )
    const
    {
      if( idU == typeid( FFVar ) )
        return eval( nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar );
      else if( idU == typeid( FFDep ) )
        return eval( nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar );
      else if( idU == typeid( FFExpr ) )
        return eval( nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar );
      else if( idU == typeid( double ) )
        return eval( nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar );
      else if( idU == typeid( I ) )
        return eval( nRes, static_cast<I*>(vRes), nVar, static_cast<I const*>(vVar), mVar );
      else if( idU == typeid( McCormick<I> ) )
        return eval( nRes, static_cast<McCormick<I>*>(vRes), nVar, static_cast<McCormick<I> const*>(vVar), mVar );
      else if( idU == typeid( SLiftVar ) )
        return eval( nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar );

      throw std::runtime_error( "Error: No evaluation method for FFArrh with type"+std::string(idU.name())+"\n" );
    }

  template <typename T>
  void eval
    ( unsigned const nRes, T* vRes, unsigned const nVar, T const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 && data );
      std::cout << "FFArrh: generic instantiation\n"; 
      vRes[0] = Op<T>::exp( - *static_cast<double*>( data ) / vVar[0] );
    }

  template <typename T>
  void eval
    ( unsigned const nRes, McCormick<T>* vRes, unsigned const nVar, McCormick<T> const* vVar,
      unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 && data );
      std::cout << "FFArrh: McCormick instantiation\n"; 
      vRes[0] = arrh( vVar[0], *static_cast<double*>( data ) );
    }

  void eval
    ( unsigned const nRes, FFVar* vRes, unsigned const nVar, FFVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 && data );
      std::cout << "FFArrh: FFVar instantiation\n"; 
      vRes[0] = operator()( vVar[0], *static_cast<double*>( data ) );
    }

  void eval
    ( unsigned const nRes, FFDep* vRes, unsigned const nVar, FFDep const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 && data );
      std::cout << "FFArrh: FFDep instantiation\n";
      vRes[0] = vVar[0];
      vRes[0].update( FFDep::TYPE::N );
    }

  void eval
    ( unsigned const nRes, SLiftVar* vRes, unsigned const nVar, SLiftVar const* vVar, unsigned const* mVar )
    const
    {
      assert( nVar == 1 && nRes == 1 && data );
      std::cout << "FFArrh: SLiftVar instantiation\n";
      vVar[0].env()->lift( nRes, vRes, nVar, vVar );
    }

  // Properties
  std::string name
    ()
    const
    { //std::cout << "data: " << data << std::endl;
      return "ARRH[" + std::to_string(*static_cast<double*>( data )) + "]"; }
};

}

///////////////////////////////////////////////////////////////////////////////

int test_external0()
{
  std::cout << "\n==============================================\ntest_external0:\n";

  // Create DAG
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 3;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFnorm2  norm2;
  mc::FFnorm12 norm12;
  mc::FFVar F[NF] = { norm2( NX, X ), norm12( 0, NX, X ), norm12( 1, NX, X ) };
  std::cout << DAG;
  
  std::ofstream o_F( "external0_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  auto F_op  = DAG.subgraph( NF, F );
  DAG.output( F_op );

  // Evaluation in real arithmetic
  double dX[NX] = { 2., 3. }, dF[NF];
  std::vector<double> dwk;
  DAG.eval( F_op, dwk, NF, F, dF, NX, X, dX );
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  // Evaluation in McCormick arithmetic
  MC mcX[NX] = { MC(I(1.5,2.5),2.), MC(I(2.5,3.5),3.) }, mcF[NF];
  std::vector<MC> mcwk;
  DAG.eval( F_op, mcwk, NF, F, mcF, NX, X, mcX );
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << mcF[i] << std::endl;

  // Forward AD
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NX, X, true );
  std::ofstream o_dFdX( "external0_dFdX.dot", std::ios_base::out );
  DAG.dot_script( NX*NF, dFdX, o_dFdX );
  o_dFdX.close();

  auto dFdX_op  = DAG.subgraph( NX*NF, dFdX );
  DAG.output( dFdX_op );
  delete[] dFdX;

  auto F2_op  = DAG.subgraph( 1, F+1 );
  DAG.output( F2_op );

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external1()
{
  std::cout << "\n==============================================\ntest_external1:\n";

  // Create DAG
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFnorm2 norm2;
  mc::FFxlog  myxlog;
  mc::FFVar F[NF] = { xlog( norm2( NX, X) ), myxlog( norm2( NX, X ) ) };
  std::cout << DAG;
  
  std::ofstream o_F( "external1_F.dot", std::ios_base::out );
  DAG.dot_script( NF, F, o_F );
  o_F.close();

  auto F_op  = DAG.subgraph( NF, F );
  DAG.output( F_op );

  // Evaluation in real arithmetic
  double dX[NX] = { 2., 3. }, dF[NF];
  std::vector<double> dwk;
  DAG.eval( F_op, dwk, NF, F, dF, NX, X, dX );
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << dF[i] << std::endl;

  // Evaluation in McCormick arithmetic
  MC mcX[NX] = { MC(I(1.5,2.5),2.), MC(I(2.5,3.5),3.) }, mcF[NF];
  std::vector<MC> mcwk;
  DAG.eval( F_op, mcwk, NF, F, mcF, NX, X, mcX );
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << mcF[i] << std::endl;

  // Forward AD
  const mc::FFVar* dFdX = DAG.FAD( NF, F, NX, X, true );
  std::ofstream o_dFdX( "external1_dFdX.dot", std::ios_base::out );
  DAG.dot_script( NF*NX, dFdX, o_dFdX );
  o_dFdX.close();

  auto dFdX_op  = DAG.subgraph( NF*NX, dFdX );
  DAG.output( dFdX_op );
  delete[] dFdX;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external2()
{
  std::cout << "\n==============================================\ntest_external2:\n";

  mc::FFGraph DAG;
  const unsigned NX = 4, NF = NX*NX;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFVar F[NF];
  for( unsigned i=0; i<NX; ++i )
    for( unsigned j=0; j<NX; ++j )
      F[i+j*NX] = pow(X[i],(int)j);
  mc::FFdet det;
  mc::FFVar G = det( NF, F );
  std::cout << DAG;
  
  std::ofstream o_G( "external2_G.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_G );
  o_G.close();

  auto G_op  = DAG.subgraph( 1, &G );
  DAG.output( G_op );

  // Evaluation in real arithmetic
  double dX[NX], dG;
  for( unsigned i=0; i<NX; ++i ) dX[i] = i+1.;
  std::vector<double> dwk;
  DAG.eval( G_op, dwk, 1, &G, &dG, NX, X, dX );
  std::cout << "G = " << dG << std::endl;

  // Evaluation in interval arithmetic
  I IX[NX], IG;
  for( unsigned i=0; i<NX; ++i ) IX[i] = 1e-10*I(-1,1)+(double)i+1.;
  std::vector<I> iwk;
  DAG.eval( G_op, iwk, 1, &G, &IG, NX, X, IX );
  std::cout << "G = " << IG << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external3()
{
  std::cout << "\n==============================================\ntest_external3:\n";

  mc::FFGraph DAG;
  mc::FFVar X, Y, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  mc::FFxlog myxlog;
  //F = myxlog(X);
  F = exp(X);
  //G = sqr(Y)+F;
  G = myxlog(Y)+F;

  auto F_op  = DAG.subgraph( 1, &F );
  DAG.output( F_op, " F" );
  auto G_op  = DAG.subgraph( 1, &G );
  DAG.output( G_op, " G" );

  std::ofstream o_comp0( "external3_0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F );

  auto GoF_op  = DAG.subgraph( 1, GoF );
  DAG.output( GoF_op, " GoF" );

  std::ofstream o_comp1( "external3_1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external4()
{
  std::cout << "\n==============================================\ntest_external4:\n";

  mc::FFGraph DAG;
  mc::FFVar X;
  X.set( &DAG );
  mc::FFxlog myxlog;

  mc::FFVar F = myxlog(X);
  std::cout << DAG;
  auto F_op  = DAG.subgraph( 1, &F );
  DAG.output( F_op, " F" );

  // Evaluation in real arithmetic
  double dX = 2., dF;
  std::vector<double> dwk;
  DAG.eval( F_op, dwk, 1, &F, &dF, 1, &X, &dX );
  std::cout << "F = " << dF << std::endl;

  // Evaluation in interval arithmetic
  I IX = I(1.5,2.5), IF;
  std::vector<I> Iwk;
  DAG.eval( F_op, Iwk, 1, &F, &IF, 1, &X, &IX );
  std::cout << "F = " << IF << std::endl;

  // Evaluation in McCormick arithmetic
  MC mcX = MC(IX,2.), mcF;
  std::vector<MC> mcwk;
  DAG.eval( F_op, mcwk, 1, &F, &mcF, 1, &X, &mcX );
  std::cout << "F = " << mcF << std::endl;

  // Forward AD
  const mc::FFVar* dFdX = DAG.FAD( 1, &F, 1, &X );
  std::cout << DAG;
  DAG.output( DAG.subgraph( 1, dFdX ), " dFdX" );
  delete[] dFdX;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external5()
{
  std::cout << "\n==============================================\ntest_external5:\n";

  mc::FFGraph DAG;
  mc::FFVar X;
  X.set( &DAG );
  mc::FFxlog myxlog;

  mc::FFVar F = myxlog(X);
  std::cout << DAG;
  auto F_op  = DAG.subgraph( 1, &F );
  DAG.output( F_op, " F" );

  // Polyhedral relaxation
  mc::PolImg<I> IMG;
  I IX = { I(1,5) };
  POLV PX( &IMG, X, IX ), PF;
  std::vector<POLV> polwk;
  DAG.eval( F_op, polwk, 1, &F, &PF, 1, &X, &PX );
  IMG.generate_cuts( 1, &PF );
  std::cout << "F =" << IMG << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external6()
{
  std::cout << "\n==============================================\ntest_external6:\n";

  mc::FFGraph DAG;
  const unsigned NP = 4;
  const unsigned NS = mc::FFDOptBase::read( NP, "fims.txt", false ); 
  mc::FFVar S[NS];
  for( unsigned int i=0; i<NS; i++ ) S[i].set( &DAG );
  mc::FFDOpt DOpt;
  mc::FFVar F = DOpt( NS, S );
  //std::cout << DAG;

  // Forward AD
  const mc::FFVar* dFdS = DAG.FAD( 1, &F, NS, S );
  std::cout << DAG;
  auto dFdS_op = DAG.subgraph( NS, dFdS );
  DAG.output( dFdS_op, " dFdS" );

  // Evaluation in real arithmetic
  std::vector<double> dwk;
  double dS[NS];
  for( unsigned i=0; i<NS; ++i ) dS[i] = 1./NS;

  auto F_op = DAG.subgraph( 1, &F );
  //DAG.output( F_op, " dFdS" );
  double dF;
  DAG.eval( F_op, dwk, 1, &F, &dF, NS, S, dS );
  std::cout << "F = " << dF << std::endl;

  double ddFdS[NS];
  DAG.eval( dFdS_op, dwk, NS, dFdS, ddFdS, NS, S, dS );
  for( unsigned i=0; i<NS; ++i ) std::cout << "dFdS[" << i << "] = " << ddFdS[i] << std::endl;

  delete[] dFdS;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external7()
{
  std::cout << "\n==============================================\ntest_external7:\n";

  // Create DAG
  mc::FFGraph DAG;
  mc::FFVar X( &DAG );
  mc::FFArrh Arrh;
  double C1(2.), C2(3.);
  //std::cout << "C1: " << &C1 << "  C2: " << &C2 << std::endl;
  //mc::FFVar F[2] = { Arrh( X, C1 ) + Arrh( X, C2 ), Arrh( X, C1 ) - Arrh( X, C2 ) };
  mc::FFVar F[2] = { Arrh( X, C1, true ) + Arrh( X, C2, true ), Arrh( X, C1, true ) - Arrh( X, C2, true ) };
  std::cout << DAG;

  std::ofstream o_F( "external7_F.dot", std::ios_base::out );
  DAG.dot_script( 2, F, o_F );
  o_F.close();

  auto F_op  = DAG.subgraph( 2, F );
  DAG.output( F_op );
  std::cout << DAG;

  // Evaluation in real arithmetic
  double dX = 2., dF[2];
  DAG.eval( F_op, 2, F, dF, 1, &X, &dX );
  std::cout << "F[0] = " << dF[0] << std::endl;
  std::cout << "F[1] = " << dF[1] << std::endl;

  // Evaluation in McCormick arithmetic
  MC mcX = MC(I(1.5,2.5),2.), mcF[2];
  DAG.eval( F_op, 2, F, mcF, 1, &X, &mcX );
  std::cout << "F[0] = " << mcF[0] << std::endl;
  std::cout << "F[1] = " << mcF[1] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external8()
{
  std::cout << "\n==============================================\ntest_external8:\n";
/*
  // Create MLP
  mc::MLP<I> f;
  f.options.RELAX     = mc::MLP<I>::Options::AUX;//MC;//AUX;//MCISM;
  f.options.ISMDIV    = 16;
  f.options.ASMBPS    = 8;
  f.options.ISMCONT   = true;
  f.options.ISMSLOPE  = true;
  f.options.ISMSHADOW = false;//true;
  f.options.CUTSHADOW = false;

  //#include "ReLUANN_30L1.hpp"
  #include "ReLUANN_40L4.hpp"
  unsigned l=0;
  for( auto const& layer : MLPCOEF )
    f.append_data( layer, (++l)<MLPCOEF.size()? mc::MLP<I>::Options::RELU:
                                                mc::MLP<I>::Options::LINEAR );

  // Create DAG
  mc::FFGraph DAG;
  size_t NX = 2;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFMLP<I> MLP;
  mc::FFVar F = MLP( 0, NX, X, &f );
  std::cout << DAG;

  std::ofstream o_F( "external8_F.dot", std::ios_base::out );
  DAG.dot_script( 1, &F, o_F );
  o_F.close();

  auto F_op  = DAG.subgraph( 1, &F );
  DAG.output( F_op );
  std::cout << DAG;

  // Evaluation in real arithmetic
  double dX[NX] = { 2., 2. }, dF;
  DAG.eval( F_op, 1, &F, &dF, NX, X, dX );
  std::cout << "F = " << dF << std::endl;

  // Evaluation in McCormick arithmetic
  MC mcX[NX] = { MC(I(-3.,3.),2.), MC(I(-3.,3.),2.) }, mcF;
  DAG.eval( F_op, 1, &F, &mcF, NX, X, mcX );
  std::cout << "F = " << mcF << std::endl;

  // Polyhedral relaxation
  mc::PolImg<I> IMG;
  IMG.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::CONT;//BIN;//SOS2;
  IMG.options.AGGREG_LQ       = true;
  IMG.options.BREAKPOINT_RTOL =
  IMG.options.BREAKPOINT_ATOL = 0e0;
  IMG.options.ALLOW_DISJ      = { mc::FFOp::FABS, mc::FFOp::MAXF };
  IMG.options.ALLOW_NLIN      = { mc::FFOp::TANH, mc::FFOp::EXP  };
  I IX[NX] = { I(-3.,3.), I(-3.,3.) };
  POLV polX[NX] = { POLV( &IMG, X[0], IX[0] ), POLV( &IMG, X[1], IX[1] ) }, polF;
  DAG.eval( F_op, 1, &F, &polF, NX, X, polX );
  IMG.generate_cuts( 1, &polF );
  std::cout << "F =" << IMG << std::endl;

  // Evaluation of forward symbolic derivatives in real arithmetic
  const mc::FFVar* dFdX_F = DAG.FAD( 1, &F, NX, X );
  std::ofstream o_dFdX_F( "external8_dFdX_F.dot", std::ios_base::out );
  DAG.dot_script( NX, dFdX_F, o_dFdX_F );
  o_dFdX_F.close();

  auto op_dFdX_F = DAG.subgraph( NX, dFdX_F );
  DAG.output( op_dFdX_F );
  //std::cout << DAG;

  double ddFdX_F[NX];
  DAG.eval( op_dFdX_F, NX, dFdX_F, ddFdX_F, NX, X, dX );
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dFdX_F[" << i << "] = " << ddFdX_F[i] << std::endl;
  delete[] dFdX_F;

  // Evaluation of forward automatic derivatives in real arithmetic
  fadbad::F<double> fdX[NX], fdF;
  for( unsigned i=0; i<NX; ++i ){
    fdX[i] = dX[i];
    fdX[i].diff(i,NX);
  }
  DAG.eval( F_op, 1, &F, &fdF, NX, X, fdX );
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dFdX[" << i << "] = " << fdF.d(i) << std::endl;

  // Evaluation of forward symbolic derivatives in real arithmetic
  const mc::FFVar* dFdX_B = DAG.BAD( 1, &F, NX, X );
  std::ofstream o_dFdX_B( "external8_dFdX_B.dot", std::ios_base::out );
  DAG.dot_script( NX, dFdX_B, o_dFdX_B );
  o_dFdX_B.close();

  auto op_dFdX_B = DAG.subgraph( NX, dFdX_B );
  DAG.output( op_dFdX_B );
  //std::cout << DAG;

  double ddFdX_B[NX];
  DAG.eval( op_dFdX_B, NX, dFdX_B, ddFdX_B, NX, X, dX );
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dFdX_B[" << i << "] = " << ddFdX_B[i] << std::endl;
  delete[] dFdX_B;

  // Evaluation of backward automatic derivatives in real arithmetic
  fadbad::B<double> bdX[NX], bdF;
  for( unsigned i=0; i<NX; ++i )
    bdX[i] = dX[i];
  DAG.eval( F_op, 1, &F, &bdF, NX, X, bdX );
  bdF.diff(0,1);
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dFdX[" << i << "] = " << bdX[i].d(0) << std::endl;
*/
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external9()
{
  std::cout << "\n==============================================\ntest_external9:\n";

  // Create sparse polynomial
  typedef mc::SPoly<mc::FFVar const*,mc::lt_FFVar> t_SPoly;
  t_SPoly::options.BASIS = t_SPoly::Options::CHEB;//MONOM; //
  t_SPoly::options.DISPLINE = true;

  size_t const NX = 3;
  t_SPoly X[NX], P;

  mc::FFGraph DAG;
  mc::FFVar DAGX[NX];
  for( size_t i=0; i<NX; i++ ) X[i].var( &DAGX[i].set( &DAG ) );
  P = pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 );
  //P = pow( X[0], 3 );//2 * sqr( X[0] ) - 1;

  std::cout << "\nSparse multivariate polynomial:\n";
  std::cout << "P = " << P << std::endl;

  // Evaluate polynomial
  std::map<mc::FFVar const*,double,mc::lt_FFVar> Xval{ { &DAGX[0], 0.5 }, { &DAGX[1], -0.5 }, { &DAGX[2], 0.75 } };
  std::cout << "\nSparse multivariate polynomial value:\n";
  std::cout << "P = " << P.eval( Xval ) << std::endl;

  // Insert sparse polynomial in DAG
  mc::FFVar DAGP = mc::FFSPoly<>::insert( P, &DAG );
  auto P_op = DAG.subgraph( 1, &DAGP );
  DAG.output( P_op, " OF P" );
  //std::cout << DAG;

  // Incorporate sparse polynomial as DAG operation
  mc::FFSPoly<I> SPoly;
  DAGP = SPoly( P );
  //std::cout << DAG;

  std::ofstream o_P( "external9_P.dot", std::ios_base::out );
  DAG.dot_script( 1, &DAGP, o_P );
  o_P.close();

  P_op = DAG.subgraph( 1, &DAGP );
  DAG.output( P_op, " OF P" );

  // Evaluation in real arithmetic
  double dX[NX] = { 0.5, -0.5, 0.75 }, dP;
  DAG.eval( P_op, 1, &DAGP, &dP, NX, DAGX, dX );
  std::cout << "\nSparse multivariate polynomial value from DAG:\n";
  std::cout << "P = " << dP << std::endl;

  // Evaluation in McCormick arithmetic
  I IX[NX] = { I(0.4,0.6), I(-0.6,-0.4), I(0.7,0.8) };
  MC mcX[NX] = { MC(IX[0],0.5), MC(IX[1],-0.5), MC(IX[2],0.75) }, mcP;
  DAG.eval( P_op, 1, &DAGP, &mcP, NX, DAGX, mcX );
  std::cout << "\nSparse multivariate polynomial McCormick relaxation from DAG:\n";
  std::cout << "P = " << mcP << std::endl;

  // Differentiate polynomial and evaluate derivative
  for( unsigned i=0; i<NX; ++i ){
    t_SPoly dPdXi = P.diff( &DAGX[i] );
    std::cout << "\nSparse multivariate polynomial derivative:\n";
    std::cout << "dPdX[" << i << "] = " << dPdXi << std::endl;

    // Evaluate polynomial derivative
    std::cout << "\nSparse multivariate polynomial derivative value:\n";
    std::cout << "dPdX[" << i << "] = " << dPdXi.eval( Xval ) << std::endl;
  }
  
  // Evaluation of forward symbolic derivatives in real arithmetic
  const mc::FFVar* DAGdPdX_FAD = DAG.FAD( 1, &DAGP, NX, DAGX );
  std::ofstream o_DAGdPdX_FAD( "external9_dPdX_FAD.dot", std::ios_base::out );
  DAG.dot_script( NX, DAGdPdX_FAD, o_DAGdPdX_FAD );
  o_DAGdPdX_FAD.close();

  auto op_DAGdPdX_FAD = DAG.subgraph( NX, DAGdPdX_FAD );
  DAG.output( op_DAGdPdX_FAD, " OF dPdX" );

  double ddPdX_FAD[NX];
  DAG.eval( op_DAGdPdX_FAD, NX, DAGdPdX_FAD, ddPdX_FAD, NX, DAGX, dX );
  std::cout << "\nSparse multivariate polynomial derivative value from DAG forward symbolic differentiation:\n";
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dPdX[" << i << "] = " << ddPdX_FAD[i] << std::endl;
  delete[] DAGdPdX_FAD;

  // Evaluation of forward automatic derivatives in real arithmetic
  fadbad::F<double> fdX[NX], fdP;
  for( unsigned i=0; i<NX; ++i ){
    fdX[i] = dX[i];
    fdX[i].diff(i,NX);
  }
  DAG.eval( P_op, 1, &DAGP, &fdP, NX, DAGX, fdX );
  std::cout << "\nSparse multivariate polynomial derivative value from DAG forward automatic differentiation:\n";
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dPdX[" << i << "] = " << fdP.d(i) << std::endl;

  // Evaluation of backward symbolic derivatives in real arithmetic
  const mc::FFVar* DAGdPdX_BAD = DAG.BAD( 1, &DAGP, NX, DAGX );
  std::ofstream o_DAGdPdX_BAD( "external9_dPdX_BAD.dot", std::ios_base::out );
  DAG.dot_script( NX, DAGdPdX_BAD, o_DAGdPdX_BAD );
  o_DAGdPdX_BAD.close();

  auto op_DAGdPdX_BAD = DAG.subgraph( NX, DAGdPdX_BAD );
  DAG.output( op_DAGdPdX_BAD, " OF dPdX" );

  double ddPdX_BAD[NX];
  DAG.eval( op_DAGdPdX_BAD, NX, DAGdPdX_BAD, ddPdX_BAD, NX, DAGX, dX );
  std::cout << "\nSparse multivariate polynomial derivative value from DAG backward symbolic differentiation:\n";
  for( unsigned i=0; i<NX; ++i )
    std::cout << "dPdX[" << i << "] = " << ddPdX_BAD[i] << std::endl;
  delete[] DAGdPdX_BAD;

  // Evaluation of backward automatic derivatives in real arithmetic
  if( t_SPoly::options.BASIS == t_SPoly::Options::MONOM ){
    fadbad::B<double> bdX[NX], bdF;
    for( unsigned i=0; i<NX; ++i )
      bdX[i] = dX[i];
    DAG.eval( P_op, 1, &DAGP, &bdF, NX, DAGX, bdX );
    bdF.diff(0,1);
    std::cout << "\nSparse multivariate polynomial derivative value from DAG backward automatic differentiation:\n";
    for( unsigned i=0; i<NX; ++i )
      std::cout << "dFdX[" << i << "] = " << bdX[i].d(0) << std::endl;
  
    // Lifting of polynomial subexpression
    mc::SLiftEnv SPE( &DAG );
    SPE.options.KEEPFACT = false;
    SPE.options.LIFTIPOW = false;
    SPE.process( 1, &DAGP, true );
    std::cout << "\nLifting of sparse multivariate polynomial in DAG:\n";
    std::cout << SPE;
  }
/*
  // Polyhedral relaxation
  mc::PolImg<I> IMG;
  IMG.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::CONT;//BIN;//SOS2;
  IMG.options.AGGREG_LQ       = true;
  IMG.options.BREAKPOINT_RTOL =
  IMG.options.BREAKPOINT_ATOL = 0e0;
  IMG.options.ALLOW_DISJ      = { mc::FFOp::FABS, mc::FFOp::MAXF };
  IMG.options.ALLOW_NLIN      = { mc::FFOp::TANH, mc::FFOp::EXP  };
  I IX[NX] = { I(-3.,3.), I(-3.,3.) };
  POLV polX[NX] = { POLV( &IMG, X[0], IX[0] ), POLV( &IMG, X[1], IX[1] ) }, polF;
  DAG.eval( F_op, 1, &F, &polF, NX, X, polX );
  IMG.generate_cuts( 1, &polF );
  std::cout << "F =" << IMG << std::endl;
*/

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external9_1()
{
  std::cout << "\n==============================================\ntest_external9_1:\n";

  // Create graph
  size_t const NX = 3;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(NX);
  for( size_t i=0; i<NX; i++ ) X[i].set( &DAG );

  // Insert linear expression in DAG
  mc::FFLin<I> Sum;
  std::vector<mc::FFVar> F{ sqrt( Sum( X, 1. ) ) - 0., 0.5*X[0]+2.*X[1]-1. };
  //std::vector<mc::FFVar> F{ Sum( X, 1. ), Sum( X, {0.5,2.,3.}, 1. ) - 0.5*X[0] };
  //std::vector<mc::FFVar> F{ exp( Sum( X, 1. ) ), Sum( X, {0.5,2.,3.}, 1. ) + 0.5*X[0] };
  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );
  auto strout = mc::FFExpr::subgraph( &DAG, opF );
  for( unsigned i=0; i<strout.size(); ++i )
    std::cout << "F[" << i << "]: " << strout[i] << std::endl;

  std::ofstream oF( "external9_1.dot", std::ios_base::out );
  DAG.dot_script( F, oF );
  oF.close();

  // Evaluation in real arithmetic
  std::vector<double> DX{ 0.5, -0.5, 0.75 }, DF;
  DAG.eval( opF, F, DF, X, DX );
  std::cout << "\nLinear subexpression values from DAG:\n";
  for( unsigned i=0; i<DF.size(); ++i )
    std::cout << "F[" << i << "]: " << DF[i] << std::endl;

  // Evaluation in interval arithmetic
  std::vector<I> IX{ I(0.4,0.6), I(-0.6,-0.4), I(0.7,0.8) }, IF;
  DAG.eval( opF, F, IF, X, IX );
  std::cout << "\nLinear subexpression interval bounds from DAG:\n";
  for( unsigned i=0; i<IF.size(); ++i )
    std::cout << "F[" << i << "]: " << IF[i] << std::endl;

  // Constraint propagation
  IF = { I(2.3), I(2.75) };
  auto flag = DAG.reval( opF, F, IF, X, IX, I(-1e10,1e10), 10 );
  std::cout << "\nLinear subexpression constraint propagation from DAG: " << flag << " passes\n";
  for( unsigned i=0; i<IX.size(); ++i )
    std::cout << "X[" << i << "]: " << IX[i] << std::endl;
  for( unsigned i=0; i<IF.size(); ++i )
    std::cout << "F[" << i << "]: " << IF[i] << std::endl;

  // Evaluation in McCormick arithmetic
  std::vector<MC> MCX{ MC(IX[0],0.5), MC(IX[1],-0.5), MC(IX[2],0.75) }, MCF;
  DAG.eval( opF, F, MCF, X, MCX );
  std::cout << "\nLinear subexpression McCormick relaxations from DAG:\n";
  for( unsigned i=0; i<MCF.size(); ++i )
    std::cout << "F[" << i << "]: " << MCF[i] << std::endl;

  // Polyhedral relaxation
  mc::PolImg<I> IMG;
  IMG.options.AGGREG_LQ = true;
  std::vector<POLV> PX{ POLV( &IMG, X[0], IX[0] ), POLV( &IMG, X[1], IX[1] ), POLV( &IMG, X[2], IX[2] ) }, PF;
  DAG.eval( opF, F, PF, X, PX );
  IMG.generate_cuts( PF );
  std::cout << "\nLinear subexpression polyhedral relaxation from DAG:\n";
  std::cout << IMG << std::endl;
  
  // Forward differentiation
  std::vector<mc::FFVar> dFdX_FAD = DAG.FAD( F, X );
  auto opdFdX_FAD = DAG.subgraph( dFdX_FAD );
  DAG.output( opdFdX_FAD, " OF dFdX" );
  
  // Lifting of polynomial subexpression
  mc::SLiftEnv SLE( &DAG );
  SLE.options.KEEPFACT = false;
  SLE.process( F, true );
  std::cout << "\nLifting of linear subexpression:\n";
  std::cout << SLE;

  // Variable elimination
  mc::SElimEnv SEE( &DAG );
  SEE.options.SLIFT.KEEPFACT = false;
  SEE.options.MIPDISPLEVEL   = 0;
  SEE.options.MIPOUTPUTFILE  = "test_external9_1.lp";
  SEE.options.DISPFULL       = false;
  SEE.process( F );
  std::cout << "\nVariable elimination in linear subexpression:\n";
  std::cout << SEE;


  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external9_2()
{
  std::cout << "\n==============================================\ntest_external9_2:\n";

  // Create graph
  size_t const NX = 3;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X = DAG.add_vars( NX );

  // Insert linear expression in DAG
  mc::FFLin<I> Sum;
  mc::FFVar Sum1;
  mc::FFSubgraph opSum1;
  std::vector<mc::FFExpr> strSum1;

  Sum1 = Sum( X, 1. );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  Sum1 = Sum( X, 1., 2. );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  Sum1 = Sum( X, 2. );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  Sum1 = Sum( X, 2., 1. );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  std::cout << DAG;

  Sum1 = Sum( X, std::vector<double>(NX,-1.) );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  std::cout << DAG;

  Sum1 = Sum( X, -1. );
  opSum1  = DAG.subgraph( {Sum1} );
  strSum1 = mc::FFExpr::subgraph( &DAG, opSum1 );
  std::cout << "Sum1: " << strSum1[0] << std::endl << std::endl;

  std::cout << DAG;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external10()
{
  std::cout << "\n==============================================\ntest_external10:\n";

  // Create DAG
  mc::FFGraph DAG;
  size_t const NX = 3;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> Y{ pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 ), X[0] };
  std::vector<std::vector<mc::FFVar>> vY( 3, Y );
  
  mc::Vect vF( &DAG, X, vY );
  mc::FFVect<I> OpVect;
  
  //mc::FFVar** ppF = OpVect( &vF );
  std::vector<mc::FFVar> F;
  for( size_t i=0; i<vF.nFun(); ++i )
    F.push_back( OpVect( i, &vF ) );

  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );

  // Evaluation in real arithmetic
  std::vector<double> dX{ 0.5, -0.5, 0.75 }, dF;
  DAG.eval( opF, F, dF, X, dX );
  std::cout << "\nVectorized function value:\n";
  for( auto const& dFi : dF )
    std::cout << dFi << std::endl;

  // Evaluation of symbolic derivatives in real arithmetic
  std::vector<mc::FFVar>&& dFdX = DAG.FAD( F, X );
  auto opdFdX = DAG.subgraph( dFdX );
  DAG.output( opdFdX, " OF dFdX" );

  std::vector<double> ddFdX;
  DAG.eval( opdFdX, dFdX, ddFdX, X, dX );
  std::cout << "\nVectorized function derivative:\n";
  for( auto const& dFdXi : ddFdX )
    std::cout << dFdXi << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external11()
{
  std::cout << "\n==============================================\ntest_external11:\n";

  // Create DAG
  mc::FFGraph DAG;
  size_t const NX = 3;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> C{ mc::FFVar(&DAG) };
  std::vector<mc::FFVar> Y{ pow( X[0] + sqr( X[1] ) - C[0] * X[2], 3 ), X[0] };
  std::vector<std::vector<mc::FFVar>> vY( 3, Y );
  
  mc::Vect vF( &DAG, X, C, vY );
  mc::FFVect<I> OpVect;
  
  //mc::FFVar** ppF = OpVect( &vF );
  std::vector<mc::FFVar> F;
  for( size_t i=0; i<vF.nFun(); ++i )
    F.push_back( OpVect( i, &vF ) );

  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );

  // Evaluation in real arithmetic
  std::vector<double> dX{ 0.5, -0.5, 0.75 }, dC{ 2. }, dF;
  DAG.eval( opF, F, dF, X, dX, C, dC );
  std::cout << "\nVectorized function value:\n";
  for( auto const& dFi : dF )
    std::cout << dFi << std::endl;

  // Evaluation of symbolic derivatives in real arithmetic
  std::vector<mc::FFVar>&& dFdX = DAG.FAD( F, X );
  auto opdFdX = DAG.subgraph( dFdX );
  DAG.output( opdFdX, " OF dFdX" );

  std::vector<double> ddFdX;
  DAG.eval( opdFdX, dFdX, ddFdX, X, dX, C, dC );
  std::cout << "\nVectorized function derivative:\n";
  for( auto const& dFdXi : ddFdX )
    std::cout << dFdXi << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external12()
{
  std::cout << "\n==============================================\ntest_external12:\n";

  // Create DAG
  mc::FFGraph DAG;
  size_t const NX = 4;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> Y{ pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 ), sqrt( X[3] ) };

  mc::DAGEXT<I> DAGY0( &DAG, X, Y );
  DAGY0.options.AUTODIFF = DAGY0.options.F;
  mc::FFDAGEXT<I> Expr( false );//true );
  std::vector<mc::FFVar> F{ Expr( 0, X, &DAGY0, 1 ), Expr( 1, X, &DAGY0, 1 ) }; // with DAG copy

  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );
  auto strout = mc::FFExpr::subgraph( &DAG, opF );
  for( unsigned i=0; i<strout.size(); ++i )
    std::cout << "F: " << strout[i] << std::endl;

  // Evaluation in real arithmetic
  std::vector<double> dX{ 0.5, -0.5, 0.75, 4. }, dF;
  DAG.eval( opF, F, dF, X, dX );
  std::cout << "\nFunction value: [ ";
  for( auto const& Fi : dF )
    std::cout << Fi << " ";
  std::cout << "]" << std::endl;

  // Evaluation of symbolic derivatives in real arithmetic
  std::vector<mc::FFVar>&& dFdX = DAG.FAD( F, X );
  auto opdFdX = DAG.subgraph( dFdX );
  DAG.output( opdFdX, " OF dFdX" );
  strout = mc::FFExpr::subgraph( &DAG, opdFdX );
  for( unsigned i=0; i<strout.size(); ++i )
    std::cout << "dFdX[" << i << "]: " << strout[i] << std::endl;

  // Evaluation in real arithmetic
  std::vector<double> ddFdX;
  DAG.eval( opdFdX, dFdX, ddFdX, X, dX );
  std::cout << "\nFunction derivative: [ ";
  for( auto const& dFdXi : ddFdX )
    std::cout << dFdXi << " ";
  std::cout << "]" << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external13()
{
  std::cout << "\n==============================================\ntest_external13:\n";

  // Create DAG
  mc::FFGraph DAG;
  size_t const NX = 3;
  std::vector<mc::FFVar> X(NX);
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> Y{ pow( X[0] + sqr( X[1] ) - 2 * X[2], 3 ) };

  mc::DAGEXT<I> DAGY0( &DAG, X, Y );
  mc::FFDAGEXT<I> Expr;
  Expr.options.RELAX  = { Expr.options.MCPWCS };//AUX };//MC };//INT };
  Expr.options.PWLINI = 4;
  std::vector<mc::FFVar> F{ Expr( 0, X, &DAGY0, 1 ) }; // with DAG copy

  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );
  auto strout = mc::FFExpr::subgraph( &DAG, opF );
  std::cout << "F: " << strout[0] << std::endl;

  // Evaluation in interval arithmetic
  std::vector<I> IX{ I(-1.,1.), I(0.,1.), I(-1.,0.) }, IF;
  DAG.eval( opF, F, IF, X, IX );
  std::cout << "\nFunction enclosure: " << IF[0] << std::endl;

  // Evaluation in McCormick arithmetic
  std::vector<MC> MCX{ MC( IX[0], mc::Op<I>::mid(IX[0]) ).sub(NX,0),
                       MC( IX[1], mc::Op<I>::mid(IX[1]) ).sub(NX,1),
                       MC( IX[2], mc::Op<I>::mid(IX[2]) ).sub(NX,2) }, MCF;
  DAG.eval( opF, F, MCF, X, MCX );
  std::cout << "\nFunction enclosure: " << MCF[0] << std::endl;

  // Evaluation in McCormick arithmetic
  PWCSM PWCSenv( NX );
  std::vector<PWCSV> PWCSX{ PWCSV( PWCSenv, 0, IX[0], 16 ),
                            PWCSV( PWCSenv, 1, IX[1], 16 ),
                            PWCSV( PWCSenv, 2, IX[2], 16 ) }, PWCSF;
  DAG.eval( opF, F, PWCSF, X, PWCSX );
  std::cout << "\nFunction enclosure: " << PWCSF[0] << std::endl;

  // Evaluation in McCormick arithmetic with PWC support bounds
  std::vector<MCPWCSV> MCPWCSX{ MCPWCSV( PWCSX[0], mc::Op<I>::mid(IX[0]) ).sub(NX,0),
                                MCPWCSV( PWCSX[1], mc::Op<I>::mid(IX[1]) ).sub(NX,1),
                                MCPWCSV( PWCSX[2], mc::Op<I>::mid(IX[2]) ).sub(NX,2) }, MCPWCSF;
  DAG.eval( opF, F, MCPWCSF, X, MCPWCSX );
  std::cout << "\nFunction enclosure: " << MCPWCSF[0] << std::endl;
  //std::cout << "\nFunction enclosure: " << MCPWCSF[0].I() << std::endl;
  
  // Polyhedral relaxation
  mc::PolImg<I> IMG;
  IMG.options.BREAKPOINT_TYPE = mc::PolImg<I>::Options::CONT;//BIN;//SOS2;
  IMG.options.AGGREG_LQ       = true;
  IMG.options.BREAKPOINT_RTOL =
  IMG.options.BREAKPOINT_ATOL = 0e0;
  IMG.options.ALLOW_DISJ      = { mc::FFOp::FABS, mc::FFOp::MAXF };
  IMG.options.ALLOW_NLIN      = { mc::FFOp::TANH, mc::FFOp::EXP  };
  std::vector<POLV> polX{ POLV( &IMG, X[0], IX[0] ), POLV( &IMG, X[1], IX[1] ), POLV( &IMG, X[2], IX[2] ) }, polF;
  DAG.eval( opF, F, polF, X, polX );
  IMG.generate_cuts( polF );
  std::cout << "\n Polyhedral relaxation:" << IMG << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external14()
{
  std::cout << "\n==============================================\ntest_external14:\n";

  // Create DAG
  mc::FFGraph DAG;
  size_t const NX = 2, NF = 1;
  std::vector<mc::FFVar> X( NX );
  for( auto& Xi : X ) Xi.set( &DAG );
  std::vector<mc::FFVar> Y{ sqr(X[0])+X[0]*X[1]+4 };
  auto opY  = DAG.subgraph( Y );

  mc::DAGEXT<I> DAGY0( &DAG, X, Y );
  DAGY0.options.CPMAX = 1;
  
  mc::FFDAGEXT<I> Expr;
  Expr.options.RELAX  = { Expr.options.AUX };
  std::vector<mc::FFVar> F{ Expr( 0, X, &DAGY0, 1 ) }; // with DAG copy

  std::cout << DAG;

  auto opF  = DAG.subgraph( F );
  DAG.output( opF, " OF F" );
  auto strout = mc::FFExpr::subgraph( &DAG, opF );
  std::cout << "F: " << strout[0] << std::endl;

  // Evaluate in interval arithmetic, both forward and backward
  std::vector<I> IWKF;
  std::vector<I> IX{ I(-0.8,-0.3), I(6.,9.) }, IF;
  std::cout << "\nInterval hull:\n";
  std::cout << "X[0] = "
            << I(-mc::Op<I>::l(IX[1])/2.+sqrt(mc::sqr(mc::Op<I>::l(IX[1])/2.)-4),
                 -mc::Op<I>::u(IX[1])/2.+sqrt(mc::sqr(mc::Op<I>::u(IX[1])/2.)-4))
            << std::endl;
  std::cout << "X[1] = " << IX[1] << std::endl;  

  DAG.eval( opY, IWKF, Y, IF, X, IX );
  std::cout << "\nInterval forward evaluation:\n";
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "Y[" << i << "] = " << IF[i] << std::endl;

  DAG.eval( opF, IWKF, F, IF, X, IX );
  std::cout << "\nInterval forward evaluation:\n";
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;

  IF[0] = 0.;
  IX = { I(-0.8,-0.3), I(6.,9.) };
  int flag = DAG.reval( opY, IWKF, Y, IF, X, IX, DAGY0.options.CPINF*I(-1,1), 10 );
  std::cout << "\nInterval forward/backward evaluation: " << flag << std::endl;
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "Y[" << i << "] = " << IF[i] << std::endl;

  IX = { I(-0.8,-0.3), I(6.,9.) };
  flag = DAG.reval( opF, IWKF, F, IF, X, IX, DAGY0.options.CPINF*I(-1,1), 10 );
  std::cout << "\nInterval forward/backward evaluation: " << flag << std::endl;
  for( unsigned i=0; i<NX; i++ ) std::cout << "X[" << i << "] = " << IX[i] << std::endl;
  for( unsigned i=0; i<NF; i++ ) std::cout << "F[" << i << "] = " << IF[i] << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

std::vector<double> addition( std::vector<double> const& x ) { return {x[0]+x[1]}; }

int test_external15()
{
  std::cout << "\n==============================================\ntest_external15:\n";

  // Create graph
  size_t const NX = 2;
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X = DAG.add_vars( NX );

  // Insert custom expressions in DAG
  mc::FFCustom<I> Custom;
  std::function<std::vector<double>( std::vector<double> const& )> DEval
    ( [](std::vector<double> x)
      -> std::vector<double>
      { return { x[0]*x[1] }; } );
  std::function<std::vector<I>( std::vector<I> const& )> IEval
    ( [](std::vector<I> x)
      -> std::vector<I>
      { return { x[0]*x[1] }; } );
  //std::function<std::vector<double>( std::vector<double> const& )> DEval( addition );
  Custom.set_eval( DEval );
  Custom.set_eval( IEval );

  mc::FFCustom<I> DCustom;
  std::function<std::vector<double>( std::vector<double> const& )> DDerEval
    ( [](std::vector<double> x)
      -> std::vector<double>
      { return { x[1], x[0] }; } );
  DCustom.set_eval( DDerEval );
  Custom.set_deriv( DCustom, 1 ); // uid=1

  mc::FFVar Y = Custom( X, 0 ); // uid=0
  mc::FFSubgraph opY = DAG.subgraph( {Y} );
  std::vector<mc::FFExpr> strY = mc::FFExpr::subgraph( &DAG, opY );
  std::cout << "\nY: " << strY[0] << std::endl << std::endl;

  std::cout << DAG;

  // Evaluation in real arithmetic
  std::vector<double> DX{2,3}, DY;
  DAG.eval( {Y}, DY, X, DX );
  std::cout << "Y = " << DY[0] << std::endl; 

  // Evaluation in interval arithmetic
  std::vector<I> IX{I(1,3), I(2,4)}, IY;
  DAG.eval( {Y}, IY, X, IX );
  std::cout << "Y = " << IY[0] << std::endl; 

  // Differentiation and evaluation in real arithmetic
  std::vector<mc::FFVar>&& dYdX = DAG.FAD( {Y}, X );
  mc::FFSubgraph opdYdX = DAG.subgraph( dYdX );
  std::vector<mc::FFExpr> strdYdX = mc::FFExpr::subgraph( &DAG, opdYdX );
  std::cout << std::endl;
  for( unsigned i=0; i<strdYdX.size(); ++i )
    std::cout << "dYdX[" << i << "]: " << strdYdX[i] << std::endl;

  std::cout << DAG;

  // Derivative evaluation in real arithmetic
  std::vector<double> DdYdX;
  DAG.eval( dYdX, DdYdX, X, DX );
  for( unsigned i=0; i<DdYdX.size(); ++i )
    std::cout << "dYdX[" << i << "] = " << DdYdX[i] << std::endl; 

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_slift_external0()
{
  std::cout << "\n==============================================\ntest_slift_external0:\n";

  // Create DAG
  mc::FFGraph DAG;
  const unsigned NX = 2, NF = 2;
  mc::FFVar X[NX];
  for( unsigned int i=0; i<NX; i++ ) X[i].set( &DAG );
  mc::FFnorm2 norm2;
  mc::FFnorm12 norm12;
  //mc::FFVar F[NF] = { max( X[0], X[1] ), norm12( 0, NX, X ), norm12( 1, NX, X ) };
  //mc::FFVar F[NF] = { norm2( NX, X ), norm12( 0, NX, X ), norm12( 1, NX, X ) };
  mc::FFVar F[NF] = { norm12( 0, NX, X ), norm12( 1, NX, X ) };
  std::cout << DAG;

  mc::SLiftEnv SPE( &DAG );
  SPE.process( 2, F, true );
  std::cout << SPE;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_slift_external1()
{
  std::cout << "\n==============================================\ntest_slift_external1:\n";

  // Create DAG
  mc::FFGraph DAG;
  mc::FFVar X( &DAG );
  mc::FFArrh Arrh;
  double C1(2.), C2(3.);
  //std::cout << "C1: " << &C1 << "  C2: " << &C2 << std::endl;
  mc::FFVar F = 1 - Arrh( X, C1 ) / Arrh( X, C2 );
  std::cout << DAG;

  mc::SLiftEnv SPE( &DAG );
  SPE.process( 1, &F, true );
  std::cout << SPE;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
//    test_external0();
//    test_external1();
//    test_external2();
//    test_external3();
//    test_external4();
//    test_external5();
//    test_external6();
//    test_external7();
//    test_external8();
//    test_external9_1();
//    test_external9_2();
//    test_external10();
//    test_external11();
    test_external12();
//    test_external13();
//    test_external14();
//    test_external15();
//    test_slift_external0();
//    test_slift_external1();

  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
	          << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
  catch( MC::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in McCormick relaxation:" << std::endl
	          << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( SCM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in sparse Chebyshev model arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( CM::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in dense Chebyshev model arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
  catch( mc::PolImg<I>::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in polyhedral image arithmetic:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

