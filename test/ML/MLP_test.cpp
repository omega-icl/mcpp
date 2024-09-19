#define MC__FFMLPREG_CHECK
#define MC__FFMLPREG_DEBUG
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

#include "ffmlpreg.hpp"

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

///////////////////////////////////////////////////////////////////////////////

int test_MLP0()
{
  std::cout << "\n==============================================\ntest_MLP0:\n";

  // Create MLP
  mc::MLPREG<I> NN;
  std::pair<size_t,std::vector<std::pair<size_t,int>>> NNstruct
    { 3,
      { {5,mc::MLPREG<I>::Options::TANH},
        {5,mc::MLPREG<I>::Options::TANH}, 
        {1,mc::MLPREG<I>::Options::LINEAR} }
    };
  NN.set_structure( NNstruct );

  // Evaluate MLP
  std::vector<double> vin( NN.nin() );
  for( unsigned i=0; i<NN.nin(); ++i ) vin[i] = 0.1+0.1*i;
  std::vector<double> vwei( NN.nwei() );
  for( unsigned i=0; i<NN.nwei(); ++i ) vwei[i] = 0.01+0.01*i;
  std::vector<double> vout( NN.nout() );
  NN.eval( vout.data(), vin.data(), vwei.data() );

  std::cout << "y =";
  for( unsigned i=0; i<NN.nout(); ++i )
    std::cout << " " << vout[i];
  std::cout << std::endl;

  // Evaluate MLP derivatives
  std::vector<fadbad::F<double>> Fvin( NN.nin() );
  for( unsigned i=0; i<NN.nin(); ++i ){
    Fvin[i] = vin[i];
    Fvin[i].diff(i,NN.nin()+NN.nwei());    
  }
  std::vector<fadbad::F<double>> Fvwei( NN.nwei() );
  for( unsigned i=0; i<NN.nwei(); ++i ){
    Fvwei[i] = vwei[i];
    Fvwei[i].diff(NN.nin()+i,NN.nin()+NN.nwei());
  }
  std::vector<fadbad::F<double>> Fvout( NN.nout() );
  NN.eval( Fvout.data(), Fvin.data(), Fvwei.data() );

  std::vector<fadbad::B<double>> Bvin( NN.nin() );
  for( unsigned i=0; i<NN.nin(); ++i ) Bvin[i] = vin[i];
  std::vector<fadbad::B<double>> Bvwei( NN.nwei() );
  for( unsigned i=0; i<NN.nwei(); ++i ) Bvwei[i] = vwei[i];
  std::vector<fadbad::B<double>> Bvout( NN.nout() );
  NN.eval( Bvout.data(), Bvin.data(), Bvwei.data() );
  for( unsigned j=0; j<NN.nout(); ++j ) Bvout[j].diff(j,NN.nout());

  for( unsigned i=0; i<NN.nin(); ++i ){
    std::cout << "dydin[" << i << "] =";
    for( unsigned j=0; j<NN.nout(); ++j )
      std::cout << " F: " << Fvout[j].d(i) << "  B: " << Bvin[i].d(j) << std::endl;
  }
  for( unsigned i=0; i<NN.nwei(); ++i ){
    std::cout << "dydwei[" << i << "] =";
    for( unsigned j=0; j<NN.nout(); ++j )
      std::cout << " F: " << Fvout[j].d(NN.nin()+i) << "  B: " << Bvwei[i].d(j) << std::endl;
  }

  // Create DAG
  mc::FFGraph DAG;
  std::vector<mc::FFVar> Xin(NN.nin());
  for( unsigned int i=0; i<NN.nin(); i++ ) Xin[i].set( &DAG );
  std::vector<mc::FFVar> Xwei(NN.nwei());
  for( unsigned int i=0; i<NN.nwei(); i++ ) Xwei[i].set( &DAG );
  mc::FFMLPREG<I> OpNN;
  std::vector<mc::FFVar> Xout(NN.nout());
  for( unsigned int j=0; j<NN.nout(); j++ )
    Xout[j] = OpNN( j, NN.nin(), Xin.data(), NN.nwei(), Xwei.data(), &NN );

  std::ofstream o_F( "MLP0_Y.dot", std::ios_base::out );
  DAG.dot_script( NN.nout(), Xout.data(), o_F );
  o_F.close();

  auto Xout_op  = DAG.subgraph( NN.nout(), Xout.data() );
  DAG.output( Xout_op );

  // Evaluation in real arithmetic
  DAG.eval( Xout_op, NN.nout(), Xout.data(), vout.data(),
            NN.nin(), Xin.data(), vin.data(),
            NN.nwei(), Xwei.data(), vwei.data() );

  std::cout << "y =";
  for( unsigned i=0; i<NN.nout(); ++i )
    std::cout << " " << vout[i];
  std::cout << std::endl;

  // Evaluation of forward automatic derivatives in real arithmetic
  DAG.eval( Xout_op, NN.nout(), Xout.data(), Fvout.data(),
            NN.nin(), Xin.data(), Fvin.data(),
            NN.nwei(), Xwei.data(), Fvwei.data() );

  // Evaluation of forward symbolic derivatives in real arithmetic
  const mc::FFVar* dXoutdXin = DAG.FAD( NN.nout(), Xout.data(), NN.nin(), Xin.data() );

  std::ofstream o_dXoutdXin( "MLP0_dYdXin.dot", std::ios_base::out );
  DAG.dot_script( NN.nout()*NN.nin(), dXoutdXin, o_dXoutdXin );
  o_dXoutdXin.close();

  auto op_dXoutdXin = DAG.subgraph( NN.nout()*NN.nin(), dXoutdXin );
  DAG.output( op_dXoutdXin );
  //std::cout << DAG;

  std::vector<double> vdXoutdXin(NN.nout()*NN.nin());
  DAG.eval( op_dXoutdXin, NN.nout()*NN.nin(), dXoutdXin, vdXoutdXin.data(),
            NN.nin(), Xin.data(), vin.data(),
            NN.nwei(), Xwei.data(), vwei.data() );

  for( unsigned i=0; i<NN.nin(); ++i ){
    std::cout << "dydin[" << i << "] =";
    for( unsigned j=0; j<NN.nout(); ++j )
      std::cout << " F: " << Fvout[j].d(i) << "  S: " << vdXoutdXin[j*NN.nin()+i] << std::endl;
  }
  for( unsigned i=0; i<NN.nwei(); ++i ){
    std::cout << "dydwei[" << i << "] =";
    for( unsigned j=0; j<NN.nout(); ++j )
      std::cout << " F: " << Fvout[j].d(NN.nin()+i) << std::endl;//"  B: " << Bvwei[i].d(j) << std::endl;
  }

  delete[] dXoutdXin;


/*
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

int main()
{
  try{
    test_MLP0();
  }
  catch( mc::FFBase::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
}

