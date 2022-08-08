#undef  MC__FFUNC_CPU_EVAL
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "mctime.hpp"
#include "mclapack.hpp"
#include "ffunc.hpp"

#ifdef MC__USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
 #else
  #ifdef MC__USE_BOOST
   #include "mcboost.hpp"
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

#include "cmodel.hpp"
typedef mc::CModel<I> CM;
typedef mc::CVar<I> CV;

#include "scmodel.hpp"
typedef mc::SCModel<I> SCM;
typedef mc::SCVar<I> SCV;

///////////////////////////////////////////////////////////////////////////////
namespace mc
{
class FFnorm2
: public FFOp
{
public:
  // Constructors
  FFnorm2
    ()
    : FFOp( (int)EXTERN )
    {}
  FFnorm2
    ( int const type, unsigned const nVar, FFVar** pVar, FFVar* pRes )
    : FFOp( type, nVar, pVar, pRes )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      auto dep = FFDep();
      for( unsigned i=0; i<nVar; ++i ) dep += pVar[i].dep();
      dep.update( FFDep::TYPE::N );
      return insert_external_operation( *this, dep, nVar, pVar );
    }

  // Evaluation overloads
  template< typename T > void eval
    ( T& vRes, unsigned const nVar, T const* vVar )
    {
      switch( nVar ){
        case 0: vRes = T( 0. ); break;
        case 1: vRes = vVar[0]; break;
        default: vRes = Op<T>::sqr( vVar[0] );
                 for( unsigned i=1; i<nVar; ++i ) vRes += Op<T>::sqr( vVar[i] );
                 vRes = Op<T>::sqrt( vRes ); break;
      }
    }
  void eval
    ( FFVar& vRes, unsigned const nVar, FFVar const* pVar )
    const
    {
      vRes = operator()( nVar, pVar );
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
  //! @brief Test type-info
  bool sameid
    ( std::type_info const& id )
    const
    { return( typeid(*this) == id ); };
};

class FFxlog
: public FFOp
{
public:
  // Constructors
  FFxlog
    ()
    : FFOp( (int)EXTERN+1 )
    {}
  FFxlog
    ( int const type, FFVar* pVar, FFVar* pRes )
    : FFOp( type, pVar, pRes )
    {}

  // Functor
  FFVar& operator()
    ( FFVar const& Var )
    const
    {
      auto dep = Var.dep();
      dep.update( FFDep::TYPE::N );
      return insert_external_operation( *this, dep, Var );
    }

  // Evaluation overloads
  template< typename T > void eval
    ( T& vRes, unsigned const nVar, T const* vVar )
    const
    { 
      assert( nVar == 1 );
      std::cout << "Generic instantiation!\n"; 
      vRes = vVar[0] * Op<T>::log( vVar[0] );
    }
  template< typename T > void eval
    ( McCormick<T>& vRes, unsigned const nVar, McCormick<T> const* vVar )
    const
    {
      assert( nVar == 1 );
      std::cout << "McCormick instantiation!\n"; 
      vRes = xlog( vVar[0] );
    }
  void eval
    ( FFVar& vRes, unsigned const nVar, FFVar const* vVar )
    const
    {
      assert( nVar == 1 );
      vRes = operator()( vVar[0] );
    }

  // Properties
  std::string name
    ()
    const
    { return "XLOG"; }

  //! @brief Test type-info
  bool sameid
    ( std::type_info const& id )
    const
    { return( typeid(*this) == id ); }
};

class FFdet
: public FFOp
{
public:
  // Constructors
  FFdet
    ()
    : FFOp( (int)EXTERN+2 )
    {}
  FFdet
    ( int const type, unsigned const nVar, FFVar** pVar, FFVar* pRes )
    : FFOp( type, nVar, pVar, pRes )
    {}

  // Functor
  FFVar& operator()
    ( unsigned const nVar, FFVar const* pVar )
    const
    {
      auto dep = FFDep();
      for( unsigned i=0; i<nVar; ++i ) dep += pVar[i].dep();
      switch( nVar ){
        case 0:
        case 1:  dep.update( FFDep::TYPE::L ); break;
        case 2:  dep.update( FFDep::TYPE::Q ); break;
        default: dep.update( FFDep::TYPE::P ); break;
      }
      return insert_external_operation( *this, dep, nVar, pVar );
    }

  // Evaluation overloads
  template< typename T > void eval
    ( T& vRes, unsigned const nVar, T const* vVar )
    const
    {
      const unsigned nDim = std::sqrt(nVar);
      switch( nDim ){
        case 0:  vRes = T( 0. ); break;
        case 1:  vRes = vVar[0]; break;
        default: vRes = FFBase::det( nDim, vVar ); break;
      }
    }
  void eval
    ( double& vRes, unsigned const nVar, double const* pVar )
    const
    {
      const unsigned nDim = std::sqrt(nVar);
      CPPL::dgematrix Amat( nDim, nDim );
      for( unsigned i=0; i<nDim; ++i )
        for( unsigned j=0; j<nDim; ++j )
          Amat(i,j) = pVar[i+j*nDim];
      if( dgeqrf( Amat, vRes ) )
        throw FFBase::Exceptions( FFBase::Exceptions::EXTERN );
    }
  void eval
    ( FFVar& vRes, unsigned const nVar, FFVar const* pVar )
    const
    {
      vRes = operator()( nVar, pVar );
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
    { return true; }
  //! @brief Test type-info
  bool sameid
    ( std::type_info const& id )
    const
    { return( typeid(*this) == id ); };
};
}

///////////////////////////////////////////////////////////////////////////////

int test_external1()
{
  std::cout << "\n==============================================\ntest_external1:\n";

  // Create DAG
  //mc::FFGraphExt DAG;
  mc::FFGraph<mc::FFnorm2,mc::FFxlog> DAG;
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

  //mc::FFGraphExt DAG;
  mc::FFGraph<mc::FFdet> DAG;
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
  for( unsigned i=0; i<NX; ++i ) IX[i] = i+1+1e-10*I(-1,1);
  std::vector<I> iwk;
  DAG.eval( G_op, iwk, 1, &G, &IG, NX, X, IX );
  std::cout << "G = " << IG << std::endl;

  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int test_external3()
{
  std::cout << "\n==============================================\ntest_external2:\n";

  //mc::FFGraphExt DAG;
  mc::FFGraph<mc::FFnorm2,mc::FFxlog> DAG;
  mc::FFVar X, Y, F, G;
  X.set( &DAG );
  Y.set( &DAG );
  mc::FFxlog myxlog;
  //F = myxlog(X);
  F = exp(X);
  //G = sqr(Y)+F;
  G = myxlog(Y)+F;
  std::cout << DAG;

  std::ofstream o_comp0( "external3_0.dot", std::ios_base::out );
  DAG.dot_script( 1, &G, o_comp0 );
  o_comp0.close();

  const mc::FFVar* GoF = DAG.compose( 1, &G, 1, &Y, &F );
  std::cout << DAG;

  std::ofstream o_comp1( "external3_1.dot", std::ios_base::out );
  DAG.dot_script( 1, GoF, o_comp1 );
  o_comp1.close();

  delete[] GoF;
  return 0;
}

///////////////////////////////////////////////////////////////////////////////

int main()
{
  try{
    test_external1();
    test_external2();
    test_external3();
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
}

