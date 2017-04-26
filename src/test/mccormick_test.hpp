// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.


#ifndef MC__MCCORMICK_TEST_HPP
#define MC__MCCORMICK_TEST_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "interval.hpp"
#include "mccormick.hpp"

namespace mc
{
//! @brief C++ class for test of mc::McCormick class using CppUnit
////////////////////////////////////////////////////////////////////////
//! McCormickTest is a C++ class for testing the convex/concave
//! relaxation arithmetic class mc::McCormick using CppUnit.
////////////////////////////////////////////////////////////////////////
class McCormickTest : public CppUnit::TestFixture
////////////////////////////////////////////////////////////////////////
{
// Implementation of: static CppUnit::TestSuite *suite();
CPPUNIT_TEST_SUITE( McCormickTest );
CPPUNIT_TEST( testAddition );
CPPUNIT_TEST( testSubtraction );
CPPUNIT_TEST( testProduct );
CPPUNIT_TEST( testDivision );
CPPUNIT_TEST( testExpression1 );
CPPUNIT_TEST( testExpression2 );
CPPUNIT_TEST( testExpression3 );
CPPUNIT_TEST( testExpression4 );
CPPUNIT_TEST( testExpression5 );
CPPUNIT_TEST_EXCEPTION( testDivisionByZero, McCormick<Interval>::Exceptions );
CPPUNIT_TEST_EXCEPTION( testNegativeLog,    McCormick<Interval>::Exceptions );
CPPUNIT_TEST_EXCEPTION( testNegativeSqrt,   McCormick<Interval>::Exceptions );
CPPUNIT_TEST_EXCEPTION( testAsinOutrange,   McCormick<Interval>::Exceptions );
CPPUNIT_TEST_EXCEPTION( testTanOutrange,    McCormick<Interval>::Exceptions );
CPPUNIT_TEST_SUITE_END();

private:

  bool Eq( const McCormick<Interval>&MC1, const McCormick<Interval>&MC2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    bool subissame = true;
    if( MC1.nsub() || MC2.nsub() )
      for( unsigned int isub=0; subissame && isub<max(MC1.nsub(),MC2.nsub()); isub++ )
        subissame = isequal( MC1.cvsub(isub), MC2.cvsub(isub), atol, rtol )
                && isequal( MC1.cvsub(isub), MC2.cvsub(isub), atol, rtol );
    return( isequal( MC1.l(), MC2.l(), atol, rtol )
         && isequal( MC1.u(), MC2.u(), atol, rtol )
         && isequal( MC1.cv(), MC2.cv(), atol, rtol )
         && isequal( MC1.cc(), MC2.cc(), atol, rtol )
         && subissame );
  }

  bool Eq( const double&D1, const double&D2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return isequal( D1, D2, atol, rtol );
  }
 
  McCormick<Interval> X;
  McCormick<Interval> Y;
  McCormick<Interval> Z;
  
public:

  void setUp(){
    X = McCormick<Interval>( Interval(1., 2.), 1.5 );
    X.sub( 3, 0 );
    Y = McCormick<Interval>( Interval(0.5, 0.7), 0.6, 0.65 );
    Y.sub( 3, 1 );
    Z = McCormick<Interval>( Interval(-1., 2.), 0.2, 1. );
    Z.sub( 3, 2 );
  }
  
  void tearDown(){}

  void testAddition(){
    CPPUNIT_ASSERT( Eq( X+Y, McCormick<Interval>( X.I()+Y.I(), X.cv()+Y.cv(), X.cc()+Y.cc() ) ) );
    CPPUNIT_ASSERT( Eq( X+Z, McCormick<Interval>( X.I()+Z.I(), X.cv()+Z.cv(), X.cc()+Z.cc() ) ) );
    CPPUNIT_ASSERT( Eq( Y+Z, McCormick<Interval>( Y.I()+Z.I(), Y.cv()+Z.cv(), Y.cc()+Z.cc() ) ) );
    CPPUNIT_ASSERT( Eq( McCormick<Interval>(X.I())+Y, McCormick<Interval>( X.I()+Y.I(), X.l()+Y.cv(), X.u()+Y.cc() ) ) );
    CPPUNIT_ASSERT( Eq( 1.+Y, McCormick<Interval>( 1.+Y.I(), 1.+Y.cv(), 1.+Y.cc() ) ) );
  }

  void testSubtraction(){
    CPPUNIT_ASSERT( Eq( X-Y, McCormick<Interval>( X.I()-Y.I(), X.cv()-Y.cc(), X.cc()-Y.cv() ) ) );
    CPPUNIT_ASSERT( Eq( X-Z, McCormick<Interval>( X.I()-Z.I(), X.cv()-Z.cc(), X.cc()-Z.cv() ) ) );
    CPPUNIT_ASSERT( Eq( Y-Z, McCormick<Interval>( Y.I()-Z.I(), Y.cv()-Z.cc(), Y.cc()-Z.cv() ) ) );
    CPPUNIT_ASSERT( Eq( X-McCormick<Interval>(Y.I()), McCormick<Interval>( X.I()-Y.I(), X.cv()-Y.u(), X.cc()-Y.l() ) ) );
    CPPUNIT_ASSERT( Eq( 1.-Y, McCormick<Interval>( 1.-Y.I(), 1.-Y.cc(), 1.-Y.cv() ) ) );
  }

  void testProduct(){
    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << X*Y << std::endl;
    McCormick<Interval> XY( X.I()*Y.I(), .85, 1. );
    double XY_cvsub[3] = { .5, 1., 0. };
    double XY_ccsub[3] = { .7, 1., 0. };
    XY.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      XY.cvsub( isub ) = XY_cvsub[isub];
      XY.ccsub( isub ) = XY_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*Y, XY ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << X*Y << std::endl;
    McCormick<Interval> XY2( X.I()*Y.I(), .85, 1. );
    double XY2_cvsub[3] = { .6, 1.5, 0. };
    double XY2_ccsub[3] = { .7, 1., 0. };
    XY2.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      XY2.cvsub( isub ) = XY2_cvsub[isub];
      XY2.ccsub( isub ) = XY2_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*Y, XY2 ) );

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << X*Z << std::endl;
    McCormick<Interval> XZ( X.I()*Z.I(), -.3, 2. );
    double XZ_cvsub[3] = { -1., 0., 1. };
    double XZ_ccsub[3] = { 2., 0., 1. };
    XZ.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      XZ.cvsub( isub ) = XZ_cvsub[isub];
      XZ.ccsub( isub ) = XZ_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*Z, XZ ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << X*Z << std::endl;
    McCormick<Interval> XZ2( X.I()*Z.I(), -.3, 2. );
    double XZ2_cvsub[3] = { -1., 0., 1. };
    double XZ2_ccsub[3] = { 2., 0., 1. };
    XZ2.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      XZ2.cvsub( isub ) = XZ2_cvsub[isub];
      XZ2.ccsub( isub ) = XZ2_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*Z, XZ2 ) );

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << Y*Z << std::endl;
    McCormick<Interval> YZ( Y.I()*Z.I(), -.05, .8 );
    double YZ_cvsub[3] = { 0., -1., .5 };
    double YZ_ccsub[3] = { 0., -1., .7 };
    YZ.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      YZ.cvsub( isub ) = YZ_cvsub[isub];
      YZ.ccsub( isub ) = YZ_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( Y*Z, YZ ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << Y*Z << std::endl;
    McCormick<Interval> YZ2( Y.I()*Z.I(), -.02, 2.3/3. );
    double YZ2_cvsub[3] = { 0., 0., 1.7/3. };
    double YZ2_ccsub[3] = { 0., 0., 1.9/3. };
    YZ2.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      YZ2.cvsub( isub ) = YZ2_cvsub[isub];
      YZ2.ccsub( isub ) = YZ2_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( Y*Z, YZ2 ) );
  }

  void testDivision(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << Y/X << std::endl;
    McCormick<Interval> YdX( Y.I()/X.I(), 2.3/6., 0.5 );
    double YdX_cvsub[3] = { -2./9., .5, 0. };
    double YdX_ccsub[3] = { -.35, .5, 0. };
    YdX.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      YdX.cvsub( isub ) = YdX_cvsub[isub];
      YdX.ccsub( isub ) = YdX_ccsub[isub];
    }
    //std::cout << YdX << std::endl;
    CPPUNIT_ASSERT( Eq( Y/X, YdX ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << Y/X << std::endl;
    McCormick<Interval> YdX2( Y.I()/X.I(), .3972026594366538, 0.475 );
    double YdX2_cvsub[3] = { -.2648017729577692, 2./3., 0. };
    double YdX2_ccsub[3] = { -.25, 1.0, 0. };
    YdX2.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      YdX2.cvsub( isub ) = YdX2_cvsub[isub];
      YdX2.ccsub( isub ) = YdX2_ccsub[isub];
    }
    //std::cout << YdX2 << std::endl;
    CPPUNIT_ASSERT( Eq( Y/X, YdX2 ) );
  }

  void testExpression1(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << X*exp(-pow(X,2)) << std::endl;
    McCormick<Interval> F11( X.I()*exp(-pow(X.I(),2)), 9.124281806826588e-02, 4.061675774727018e-01 );
    double F11_cvsub[3] = { -2.279393569829622e-01, 0., 0. };
    double F11_ccsub[3] = {  1.831563888873416e-02, 0., 0. };
    F11.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F11.cvsub( isub ) = F11_cvsub[isub];
      F11.ccsub( isub ) = F11_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*exp(-pow(X,2)), F11 ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << X*exp(-pow(X,2)) << std::endl;
    McCormick<Interval> F12( X.I()*exp(-pow(X.I(),2)), 9.124281806826588e-02, 4.061675774727018e-01 );
    double F12_cvsub[3] = { -2.279393569829622e-01, 0., 0. };
    double F12_ccsub[3] = {  1.831563888873416e-02, 0., 0. };
    F12.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F12.cvsub( isub ) = F12_cvsub[isub];
      F12.ccsub( isub ) = F12_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( X*exp(-pow(X,2)), F12 ) );
  }

  void testExpression2(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << sin(pow(Y,-3))*cos(pow(Y,2)) << std::endl;
    McCormick<Interval> F21( sin(pow(Y.I(),-3))*cos(pow(Y.I(),2)), -9.358968236779348e-01, 6.095699354841704e-01 );
    double F21_cvsub[3] = { 0.,  4.227290799301080e-01, 0. };
    double F21_ccsub[3] = { 0., -4.004437741413581e+00, 0. };
    F21.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F21.cvsub( isub ) = F21_cvsub[isub];
      F21.ccsub( isub ) = F21_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( sin(pow(Y,-3))*cos(pow(Y,2)), F21 ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << sin(pow(Y,-3))*cos(pow(Y,2)) << std::endl;
    McCormick<Interval> F22( sin(pow(Y.I(),-3))*cos(pow(Y.I(),2)), -9.358968236779348e-01, 6.095699354841704e-01 );
    double F22_cvsub[3] = { 0.,  4.227290799301080e-01, 0. };
    double F22_ccsub[3] = { 0., -4.004437741413581e+00, 0. };
    F22.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F22.cvsub( isub ) = F22_cvsub[isub];
      F22.ccsub( isub ) = F22_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( sin(pow(Y,-3))*cos(pow(Y,2)), F22 ) );
  }

  void testExpression3(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << pow(-fabs(Z-0.5)-sqrt(X),3) << std::endl;
    McCormick<Interval> F31( pow(-fabs(Z.I()-0.5)-sqrt(X.I()),3), -2.239865823691492e+01, -1.758883476483184e+00 );
    double F31_cvsub[3] = { -5.065077037389579e+00, 0., 0. };
    double F31_ccsub[3] = { -1.810660171779821e+00, 0., 0. };
    F31.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F31.cvsub( isub ) = F31_cvsub[isub];
      F31.ccsub( isub ) = F31_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( pow(-fabs(Z-0.5)-sqrt(X),3), F31 ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << pow(-fabs(Z-0.5)-sqrt(X),3) << std::endl;
    McCormick<Interval> F32( pow(-fabs(Z.I()-0.5)-sqrt(X.I()),3), -2.239865823691492e+01, -1.758883476483184e+00 );
    double F32_cvsub[3] = { -5.065077037389579e+00, 0., 0. };
    double F32_ccsub[3] = { -1.810660171779821e+00, 0., 0. };
    F32.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F32.cvsub( isub ) = F32_cvsub[isub];
      F32.ccsub( isub ) = F32_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( pow(-fabs(Z-0.5)-sqrt(X),3), F32 ) );
  }

  void testExpression4(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << (asin(Y)+acos(Y))*atan(Z) << std::endl;
    McCormick<Interval> F41( (asin(Y.I())+acos(Y.I()))*atan(Z.I()), -3.030073040676173e-01, 1.391051187094292e+00 );
    double F41_cvsub[3] = { 0., -7.063560898229274e-03, 8.396783008670946e-01 };
    double F41_ccsub[3] = { 0.,  9.957258313039266e-03, 6.594988028912212e-01 };
    F41.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F41.cvsub( isub ) = F41_cvsub[isub];
      F41.ccsub( isub ) = F41_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( (asin(Y)+acos(Y))*atan(Z), F41 ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << (asin(Y)+acos(Y))*atan(Z) << std::endl;
    McCormick<Interval> F42( (asin(Y.I())+acos(Y.I()))*atan(Z.I()), -2.905738429932485e-01, 1.391051187094292e+00 );
    double F42_cvsub[3] = { 0., 0., 9.727224172270860e-01 };
    double F42_ccsub[3] = { 0.,  9.957258313039266e-03, 6.594988028912212e-01 };
    F42.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F42.cvsub( isub ) = F42_cvsub[isub];
      F42.ccsub( isub ) = F42_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( (asin(Y)+acos(Y))*atan(Z), F42 ) );
  }

  void testExpression5(){
    McCormick<Interval>::options.DISPLAY_DIGITS=15;

    McCormick<Interval>::options.MVCOMP_USE = false;
    //std::cout << tan(Y)*erf(X*Y*Z) << std::endl;
    McCormick<Interval> F51( tan(Y.I())*erf(X.I()*Y.I()*Z.I()), -5.141693048301641e-01, 7.497584127994054e-01 );
    double F51_cvsub[3] = { -1.857260989733198e-01, -1.674637795015533e+00, 1.326614992666570e-01 };
    double F51_ccsub[3] = {  9.096067836899530e-02,  1.609762265101710e+00, 3.248595656035547e-02 };
    F51.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F51.cvsub( isub ) = F51_cvsub[isub];
      F51.ccsub( isub ) = F51_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( tan(Y)*erf(X*Y*Z), F51 ) );

    McCormick<Interval>::options.MVCOMP_USE = true;
    //std::cout << tan(Y)*erf(X*Y*Z) << std::endl;
    McCormick<Interval> F52( tan(Y.I())*erf(X.I()*Y.I()*Z.I()), -5.141693048301641e-01, 7.497584127994054e-01 );
    double F52_cvsub[3] = { -1.857260989733198e-01, -1.674637795015533e+00, 1.326614992666570e-01 };
    double F52_ccsub[3] = {  9.096067836899530e-02,  1.609762265101710e+00, 3.248595656035547e-02 };
    F52.sub( 3 );
    for( unsigned int isub=0; isub<3; isub++ ){
      F52.cvsub( isub ) = F52_cvsub[isub];
      F52.ccsub( isub ) = F52_ccsub[isub];
    }
    CPPUNIT_ASSERT( Eq( tan(Y)*erf(X*Y*Z), F52 ) );
  }

  void testDivisionByZero(){
    // The following line should throw an instance of McCormick<Interval>::Exceptions
    inv(Z);
  }

  void testNegativeLog(){
    // The following line should throw an instance of McCormick<Interval>::Exceptions
    log( Z );
  }

  void testNegativeSqrt(){
    // The following line should throw an instance of McCormick<Interval>::Exceptions
    sqrt( Z );
  }

  void testAsinOutrange(){
    // The following line should throw an instance of McCormick<Interval>::Exceptions
    asin( Z );
  }

  void testTanOutrange(){
    // The following line should throw an instance of McCormick<Interval>::Exceptions
    tan( Z );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( mc::McCormickTest );

} // end namespace mc

#endif
