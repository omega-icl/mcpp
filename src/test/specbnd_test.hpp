// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.


#ifndef MC__SPECBND_TEST_HPP
#define MC__SPECBND_TEST_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "interval.hpp"
#include "specbnd.hpp"
#include "mcfadbad.hpp"

typedef mc::Interval I;
typedef mc::Specbnd<I> SB;
typedef fadbad::F<I> FI;
typedef fadbad::F<FI> FFI;

namespace mc
{
//! @brief C++ class for test of mc::TModel / mc::TVar class using CppUnit
////////////////////////////////////////////////////////////////////////
//! TModelTest is a C++ class for testing the multivariate Taylor
//! arithmetic class mc::Taylor using CppUnit.
////////////////////////////////////////////////////////////////////////
class SpecbndTest : public CppUnit::TestFixture
////////////////////////////////////////////////////////////////////////
{
// Implementation of: static CppUnit::TestSuite *suite();
CPPUNIT_TEST_SUITE( SpecbndTest );
CPPUNIT_TEST( testExpression1 );
CPPUNIT_TEST( testExpression2 );
CPPUNIT_TEST( testExpression3 );
CPPUNIT_TEST( testExpression4 );
CPPUNIT_TEST( testExpression5 );
CPPUNIT_TEST_EXCEPTION( testSize, SB::Exceptions );
CPPUNIT_TEST_SUITE_END();

private:

  bool Eq( const double&D1, const double&D2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return isequal( D1, D2, atol, rtol );
  }

  bool Eq( const I&I1, const I&I2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return( isequal( I1.l(), I2.l(), atol, rtol )
         && isequal( I1.u(), I2.u(), atol, rtol ) );
  }

  bool Eq( const SB&SB1, const SB&SB2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    if( SB1.dep() != SB2.dep() ) return false;
    for( unsigned int i=0; i<SB1.dep(); i++ )
      if( !Eq( SB1.FI().deriv(i), SB2.FI().deriv(i), atol, rtol ) )
        return false;
    return( Eq( SB1.I(),  SB2.I(),  atol, rtol )
         && Eq( SB1.SI(), SB2.SI(), atol, rtol ) );
  }
 
  SB SBX1;
  SB SBX2;
  SB SBX3;
  SB SBF;
  
public:

  void setUp(){}
  
  void tearDown(){}

  void testExpression1(){
    I::options.DISPLAY_DIGITS = 15;
    SBX1.set( I(-0.3,0.2), 0, 3 );
    SBX2.set( I(-0.1,0.6), 1, 3 );
    SBX3.set( I(-0.4,0.5), 2, 3 );
    SBF = exp( SBX1 - 2*sqr(SBX2) + 3*pow(SBX3,3) );
    //std::cout << SBF << std::endl;
    CPPUNIT_ASSERT( Eq( SBF.I(), I( 2.976014808681888e-01, 1.777130526914038e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(0), I( 2.976014808681888e-01, 1.777130526914038e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(1), I(-4.265113264593691e+00, 7.108522107656153e-01) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(2), I( 0.000000000000000e+00, 3.998543685556586e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.SI(), I( -1.990386190143723e+01, 3.700430039666756e+01) ) );

    FI FX1 = SBX1.I(); FX1.diff(0,3);
    FI FX2 = SBX2.I(); FX2.diff(1,3);
    FI FX3 = SBX3.I(); FX3.diff(2,3);
    FFI FFX1 = FX1; FFX1.diff(0,3);
    FFI FFX2 = FX2; FFX2.diff(1,3);
    FFI FFX3 = FX3; FFX3.diff(2,3);
    FFI FFF = exp( FFX1 - 2*pow(FFX2,2) + 3*pow(FFX3,3) );

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    //std::cout << I(spbndG.first,spbndG.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndG.first,  -2.639038832467347e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndG.second,  3.858594656562106e+01 ) );

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    //std::cout << I(spbndHR.first,spbndHR.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndHR.first,  -2.119728596072178e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndHR.second,  3.052719845857329e+01 ) );
  }

  void testExpression2(){
    I::options.DISPLAY_DIGITS = 15;
    SBX1.set( I(-2., 0.), 0, 2 );
    SBX2.set( I( 1., 3.), 1, 2 );
    SBF = -1./(pow(SBX1-4.,2)+pow(SBX2-4.,2)+0.1)
          -1./(pow(SBX1-1.,2)+pow(SBX2-1.,2)+0.2)
          -1./(pow(SBX1-8.,2)+pow(SBX2-8.,2)+0.2);
    //std::cout << SBF << std::endl;
    CPPUNIT_ASSERT( Eq( SBF.I(), I( -9.030236278289147e-01,-1.046329376284258e-01) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(0), I(-4.210218558692969e+00,-1.613029532115490e-02) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(1), I(-2.227867078139488e-02, 2.776345275060252e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.SI(), I( -6.025886165342616e+01, 2.791959925274383e+00) ) );

    FI FX1 = SBX1.I(); FX1.diff(0,3);
    FI FX2 = SBX2.I(); FX2.diff(1,3);
    FI FX3 = SBX3.I(); FX3.diff(2,3);
    FFI FFX1 = FX1; FFX1.diff(0,3);
    FFI FFX2 = FX2; FFX2.diff(1,3);
    FFI FFX3 = FX3; FFX3.diff(2,3);
    FFI FFF = -1./(pow(FFX1-4.,2)+pow(FFX2-4.,2)+0.1)
              -1./(pow(FFX1-1.,2)+pow(FFX2-1.,2)+0.2)
              -1./(pow(FFX1-8.,2)+pow(FFX2-8.,2)+0.2);

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    //std::cout << I(spbndG.first,spbndG.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndG.first,  -6.949017247140594e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndG.second,  2.917316513057383e+01 ) );

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    //std::cout << I(spbndHR.first,spbndHR.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndHR.first,  -6.021748572413666e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndHR.second,  2.917072504731187e+01) );
  }

  void testExpression3(){
    I::options.DISPLAY_DIGITS = 15;
    SBX1.set( I(-.5, .5), 0, 2 );
    SBX2.set( I(-.5, .5), 1, 2 );
    SBF = 1.+SBX1*sin(2.*SBX1+3.*SBX2)-cos(3.*SBX1-5.*SBX2);
    //std::cout << SBF << std::endl;
    CPPUNIT_ASSERT( Eq( SBF.I(), I( -5e-01, 2.5e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(0), I(-5e+00, 5e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(1), I(-6.5e+00, 6.5e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.SI(), I(-4.570783850655786e+01, 4.610555127546399e+01) ) );

    FI FX1 = SBX1.I(); FX1.diff(0,3);
    FI FX2 = SBX2.I(); FX2.diff(1,3);
    FI FX3 = SBX3.I(); FX3.diff(2,3);
    FFI FFX1 = FX1; FFX1.diff(0,3);
    FFI FFX2 = FX2; FFX2.diff(1,3);
    FFI FFX3 = FX3; FFX3.diff(2,3);
    FFI FFF = 1.+FFX1*sin(2.*FFX1+3.*FFX2)-cos(3.*FFX1-5.*FFX2);

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    //std::cout << I(spbndG.first,spbndG.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndG.first,  -5.05e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndG.second,  5.05e+01 ) );

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    //std::cout << I(spbndHR.first,spbndHR.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndHR.first,  -4.420150445666332e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndHR.second,  4.446626656303890e+01 ) );
  }

  void testExpression4(){
    I::options.DISPLAY_DIGITS = 15;
    SBX1.set( I(-.5, .5), 0, 2 );
    SBX2.set( I(-.5, .5), 1, 2 );
    SBF = tan(SBX1*SBX2) + asin(SBX1)*acos(SBX2);
    //std::cout << SBF << std::endl;
    CPPUNIT_ASSERT( Eq( SBF.I(), I( -1.351964632453187e+00, 1.351964632453187e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(0), I( 5.145978028301729e-01, 2.950998900678716e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(1), I(-1.137199536444498e+00, 1.137199536444498e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.SI(), I(-4.680165331176305e+00, 4.680165331176305e+00) ) );

    FI FX1 = SBX1.I(); FX1.diff(0,3);
    FI FX2 = SBX2.I(); FX2.diff(1,3);
    FI FX3 = SBX3.I(); FX3.diff(2,3);
    FFI FFX1 = FX1; FFX1.diff(0,3);
    FFI FFX2 = FX2; FFX2.diff(1,3);
    FFI FFX3 = FX3; FFX3.diff(2,3);
    FFI FFF = tan(FFX1*FFX2) + asin(FFX1)*acos(FFX2);

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    //std::cout << I(spbndG.first,spbndG.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndG.first,  -2.217589520854308e+00 ) );
    CPPUNIT_ASSERT( Eq( spbndG.second,  2.217589520854308e+00 ) );

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    //std::cout << I(spbndHR.first,spbndHR.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndHR.first,  -1.909043632551918e+00 ) );
    CPPUNIT_ASSERT( Eq( spbndHR.second,  1.909043632551918e+00 ) );
  }

  void testExpression5(){
    I::options.DISPLAY_DIGITS = 15;
    SBX1.set( I( .5, 2.), 0, 2 );
    SBX2.set( I( .5, 2.), 1, 2 );
    SBF = sqrt(pow(SBX1,2)+pow(SBX2,2));
    //std::cout << SBF << std::endl;
    CPPUNIT_ASSERT( Eq( SBF.I(), I( 7.071067811865476e-01, 2.828427124746190e+00) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(0), I( 1.767766952966369e-01, 2.828427124746190e+00 ) ) );
    CPPUNIT_ASSERT( Eq( SBF.FI().deriv(1), I( 1.767766952966369e-01, 2.828427124746190e+00 ) ) );
    CPPUNIT_ASSERT( Eq( SBF.SI(), I(-3.2e+01, 4e+00) ) );

    FI FX1 = SBX1.I(); FX1.diff(0,3);
    FI FX2 = SBX2.I(); FX2.diff(1,3);
    FI FX3 = SBX3.I(); FX3.diff(2,3);
    FFI FFX1 = FX1; FFX1.diff(0,3);
    FFI FFX2 = FX2; FFX2.diff(1,3);
    FFI FFX3 = FX3; FFX3.diff(2,3);
    FFI FFF = sqrt(pow(FFX1,2)+pow(FFX2,2));

    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> spbndG = SB::spectral_bound( FFF );
    //std::cout << I(spbndG.first,spbndG.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndG.first,  -2.121320343559642e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndG.second,  1.268372788753369e+01 ) );

    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> spbndHR = SB::spectral_bound( FFF );
    //std::cout << I(spbndHR.first,spbndHR.second) << std::endl;
    CPPUNIT_ASSERT( Eq( spbndHR.first,  -2.121320343559642e+01 ) );
    CPPUNIT_ASSERT( Eq( spbndHR.second,  1.268372788753369e+01 ) );
  }
 
  void testSize(){
    // The following line should throw an instance of SB::Exceptions
    SBX1.set( I(0., 1.), 0, 1 );
    SBX2.set( I(0., 1.), 1, 3 );
    SBX1 + SBX2;
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( mc::SpecbndTest );

} // end namespace mc

#endif
