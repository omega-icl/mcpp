// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.


#ifndef MC__INTERVAL_TEST_HPP
#define MC__INTERVAL_TEST_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "interval.hpp"


namespace mc
{
//! @brief C++ class for test of mc::Interval class using CppUnit
////////////////////////////////////////////////////////////////////////
//! IntervalTest is a C++ class for testing the interval arithmetic
//! class mc::Interval using CppUnit.
////////////////////////////////////////////////////////////////////////
class IntervalTest : public CppUnit::TestFixture
////////////////////////////////////////////////////////////////////////
{
// Implementation of: static CppUnit::TestSuite *suite();
CPPUNIT_TEST_SUITE( IntervalTest );
CPPUNIT_TEST( testAddition );
CPPUNIT_TEST( testSubtraction );
CPPUNIT_TEST( testProduct );
CPPUNIT_TEST( testDivision );
CPPUNIT_TEST( testFunction );
CPPUNIT_TEST( testExpression );
CPPUNIT_TEST_EXCEPTION( testDivisionByZero, mc::Interval::Exceptions );
CPPUNIT_TEST_EXCEPTION( testNegativeLog,    mc::Interval::Exceptions );
CPPUNIT_TEST_EXCEPTION( testNegativeSqrt,   mc::Interval::Exceptions );
CPPUNIT_TEST_EXCEPTION( testAsinOutrange,   mc::Interval::Exceptions );
CPPUNIT_TEST_EXCEPTION( testAcosOutrange,   mc::Interval::Exceptions );
CPPUNIT_TEST_EXCEPTION( testTanOutrange,    mc::Interval::Exceptions );
CPPUNIT_TEST_SUITE_END();

private:

  bool Eq( const Interval&I1, const Interval&I2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return( isequal( I1.l(), I2.l(), atol, rtol )
         && isequal( I1.u(), I2.u(), atol, rtol ) );
  }

  bool Eq( const double&D1, const double&D2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return isequal( D1, D2, atol, rtol );
  }
 
  Interval X;
  Interval Y;
  Interval Z;
  
public:

  void setUp(){
    X = Interval( 1., 2. );
    Y = Interval( 0.5, 0.7 );
    Z = Interval( -1., 2. );
  }
  
  void tearDown(){}

  void testAddition(){
    CPPUNIT_ASSERT( Eq( Interval( 0, 1 ) + Interval( -1, 2 ), Interval( -1, 3 ) ) );
    CPPUNIT_ASSERT( Eq( Interval( 0, 1 ) + 1, Interval( 1, 2 ) ) );
    CPPUNIT_ASSERT( Eq( 1 + Interval( 0, 1 ), Interval( 1, 2 ) ) );
  }

  void testSubtraction(){
    CPPUNIT_ASSERT( Eq( Interval( 0, 1 ) - Interval( -1, 2 ), Interval( -2, 2 ) ) );
    CPPUNIT_ASSERT( Eq( Interval( 1, 2 ) - 1, Interval( 0, 1 ) ) );
    CPPUNIT_ASSERT( Eq( 1 - Interval( -1, 3 ), Interval( -2, 2 ) ) );
  }

  void testProduct(){
    //std::cout << *I_1_2 * *I_n1_2 << std::endl;
    CPPUNIT_ASSERT( Eq( Interval( 1, 2 ) * Interval( -1, 2 ), Interval( -2, 4 ) ) );
    CPPUNIT_ASSERT( Eq( Interval( -1, 2 ) * 2., Interval( -2, 4 ) ) );
    CPPUNIT_ASSERT( Eq( 2. * Interval( -1, 2 ), Interval( -2, 4 ) ) );
  }

  void testDivision(){
    //std::cout << 2. / *I_1_2 << std::endl;
    CPPUNIT_ASSERT( Eq( Interval( -2, 4 ) / Interval( 2, 4 ), Interval( -1, 2 ) ) );
    CPPUNIT_ASSERT( Eq( Interval( -2, 4 ) / 2., Interval( -1, 2 ) ) );
    CPPUNIT_ASSERT( Eq( 2. / Interval( 1, 2 ), Interval( 1, 2 ) ) );
  }

  void testFunction(){
    //std::cout << inv( Interval( 1, 2 ) ) << std::endl;
    CPPUNIT_ASSERT( Eq( exp( Interval( -1, 2 ) ), Interval( std::exp(-1), std::exp(2) ) ) );
    CPPUNIT_ASSERT( Eq( log( Interval( std::exp(-1), std::exp(2) ) ), Interval( -1, 2 ) ) );
    CPPUNIT_ASSERT( Eq( sqrt( Interval( 1, 2 ) ), Interval( 1, std::sqrt(2) ) ) );
    CPPUNIT_ASSERT( Eq( sqr( Interval( 1, 2 ) ), Interval( 1, mc::sqr(2) ) ) );
    CPPUNIT_ASSERT( Eq( sqr( Interval( -1, 2 ) ), Interval( 0, mc::sqr(2) ) ) );
    CPPUNIT_ASSERT( Eq( mid( Interval( -4, 2 ) ), -1. ) );
    CPPUNIT_ASSERT( Eq( diam( Interval( -4, 2 ) ), 6. ) );
    CPPUNIT_ASSERT( Eq( abs( Interval( -4, 2 ) ), 4. ) );
    CPPUNIT_ASSERT( Eq( inv( Interval( 1, 2 ) ), Interval( 0.5, 1 ) ) );
  }

  void testExpression(){
    mc::Interval::options.DISPLAY_DIGITS=15;
    //std::cout << X*exp(-pow(X,2)) << std::endl;
    CPPUNIT_ASSERT( Eq( X*exp(-pow(X,2)),
                    Interval( 1.831563888873418e-02, 7.357588823428847e-01 ) ));
    //std::cout << sin(pow(Y,-3))*cos(pow(Y,2)) << std::endl;
    CPPUNIT_ASSERT( Eq( sin(pow(Y,-3))*cos(pow(Y,2)),
                    Interval( -9.689124217106447e-01, 9.689124217106447e-01 ) ));
    //std::cout << pow(-fabs(Z-0.5)-sqrt(X),3) << std::endl;
    CPPUNIT_ASSERT( Eq( pow(-fabs(Z-0.5)-sqrt(X),3),
                    Interval( -2.474936867076458e+01, -1.000000000000000e+00 ) ));
    //std::cout << (asin(Y)+acos(Y))*atan(Z) << std::endl;
    CPPUNIT_ASSERT( Eq( (asin(Y)+acos(Y))*atan(Z),
                    Interval( -1.431462803165178e+00, 2.017883770237767e+00 ) ));
  }

  void testDivisionByZero(){
    // The following line should throw an instance of mc::Interval::Exceptions
    1. / Interval( -1, 2 );
  }

  void testNegativeLog(){
    // The following line should throw an instance of mc::Interval::Exceptions
    log( Interval( -1, 2 ) );
  }

  void testNegativeSqrt(){
    // The following line should throw an instance of mc::Interval::Exceptions
    sqrt( Interval( -1, 2 ) );
  }

  void testAsinOutrange(){
    // The following line should throw an instance of mc::Interval::Exceptions
    asin( Interval( -1, 2 ) );
  }

  void testAcosOutrange(){
    // The following line should throw an instance of mc::Interval::Exceptions
    acos( Interval( -1, 2 ) );
  }

  void testTanOutrange(){
    // The following line should throw an instance of mc::Interval::Exceptions
    acos( Interval( -1, 2 ) );
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( mc::IntervalTest );

} // end namespace mc

#endif
