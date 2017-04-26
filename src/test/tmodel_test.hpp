// Copyright (C) 2009-2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.


#ifndef MC__TMODEL_TEST_HPP
#define MC__TMODEL_TEST_HPP

#include <iostream>
#include <iomanip>
#include <string>
#include <cppunit/extensions/HelperMacros.h>
#include "interval.hpp"
#include "tmodel.hpp"

namespace mc
{
//! @brief C++ class for test of mc::TModel / mc::TVar class using CppUnit
////////////////////////////////////////////////////////////////////////
//! TModelTest is a C++ class for testing the multivariate Taylor
//! arithmetic class mc::Taylor using CppUnit.
////////////////////////////////////////////////////////////////////////
class TModelTest : public CppUnit::TestFixture
////////////////////////////////////////////////////////////////////////
{
// Implementation of: static CppUnit::TestSuite *suite();
CPPUNIT_TEST_SUITE( TModelTest );
CPPUNIT_TEST( testExpression1 );
CPPUNIT_TEST( testExpression2 );
CPPUNIT_TEST( testExpression3 );
CPPUNIT_TEST( testExpression4 );
CPPUNIT_TEST( testExpression5 );
CPPUNIT_TEST_EXCEPTION( testDivisionByZero, TModel<Interval>::Exceptions );
CPPUNIT_TEST_EXCEPTION( testTMEnv, TModel<Interval>::Exceptions );
CPPUNIT_TEST_SUITE_END();

private:

  bool Eq( const Interval&I1, const Interval&I2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return( isequal( I1.l(), I2.l(), atol, rtol )
         && isequal( I1.u(), I2.u(), atol, rtol ) );
  }

  bool Eq( const TVar<Interval>&TV1, const TVar<Interval>&TV2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    if( TV1.env() != TV2.env() ) return false;

    std::pair< unsigned int, const double* > TV1_coef = TV1.coefmon();
    std::pair< unsigned int, const double* > TV2_coef = TV2.coefmon();
    if( TV1_coef.first != TV2_coef.first ) return false;

    bool polissame = true;
    for( unsigned int imon=0; polissame && imon<TV1_coef.first; imon++ )
      if( !isequal( TV1_coef.second[imon], TV2_coef.second[imon], atol, rtol ) )
        polissame = false;
    return( isequal( TV1.R().l(), TV2.R().l(), atol, rtol )
         && isequal( TV1.R().u(), TV2.R().u(), atol, rtol )
         && polissame );
  }

  bool Eq( const double&D1, const double&D2,
           const double atol=1e2*machprec(), const double rtol=1e2*machprec() ) const
  {
    return isequal( D1, D2, atol, rtol );
  }
 
  TModel<Interval>* TM_1d;
  TModel<Interval>* TM_2d;
  TVar<Interval> TVX;
  TVar<Interval> TVX1;
  TVar<Interval> TVX2;
  double* TM_1d_coef;
  double* TM_2d_coef;
  
public:

  void setUp(){
    TM_1d  = new TModel<Interval>( 1, 4 );
    TM_2d  = new TModel<Interval>( 2, 4 );
    TM_1d_coef = new double[TM_1d->nmon()];
    TM_2d_coef = new double[TM_2d->nmon()];
  }
  
  void tearDown(){
    delete TM_1d;
    delete TM_2d;
    delete[] TM_1d_coef;
    delete[] TM_2d_coef;
  }

  void testExpression1(){
    Interval::options.DISPLAY_DIGITS=15;
    TM_1d->options.DISPLAY_DIGITS=15;
    TM_1d->options.BOUNDER_ORDER = 0;
    TM_1d->options.BERNSTEIN_USE = false;
    TM_1d->options.REF_MIDPOINT  = true;
    TVX  = TVar<Interval>( TM_1d,  0, Interval(1., 2.) );

    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TVar<Interval> TVF;
    TVF.set( TM_1d );
    TM_1d_coef[0] = -6.5;
    TM_1d_coef[1] = -4.25;
    TM_1d_coef[2] =  1.;
    TM_1d_coef[3] =  1./3.;
    TM_1d_coef[4] =  0.;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(0.) );
    //std::cout << 1.-5.*TVX-pow(TVX,2)/2.+pow(TVX,3)/3. << std::endl;
    CPPUNIT_ASSERT( Eq( 1.-5.*TVX-pow(TVX,2)/2.+pow(TVX,3)/3., TVF ) );
    CPPUNIT_ASSERT( Eq( (1.-5.*TVX-pow(TVX,2)/2.+pow(TVX,3)/3.).B(),
                        Interval(-8.416666666666666,-4.083333333333333) ) );

    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::BERNSTEIN;
    TM_1d->options.BOUNDER_ORDER = 0;
    //std::cout << 1.-5.*TVX-pow(TVX,2)/2.+pow(TVX,3)/3. << std::endl;
    CPPUNIT_ASSERT( Eq( (1.-5.*TVX-pow(TVX,2)/2.+pow(TVX,3)/3.).B(),
                        Interval(-8.333333333333334,-4.166666666666667) ) );
  }

  void testExpression2(){
    Interval::options.DISPLAY_DIGITS=15;
    TM_1d->options.DISPLAY_DIGITS=15;
    TM_1d->options.BOUNDER_ORDER = 0;
    TM_1d->options.BERNSTEIN_USE = false;
    TM_1d->options.REF_MIDPOINT  = true;
    TVX  = TVar<Interval>( TM_1d,  0, Interval(1., 2.) );

    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TVar<Interval> TVF;
    TVF.set( TM_1d );
    TM_1d_coef[0] =  1.580977915390473e-01;
    TM_1d_coef[1] = -3.688347260758422e-01;
    TM_1d_coef[2] =  2.356737265178344e-01;
    TM_1d_coef[3] =  1.485652969885668e-01;
    TM_1d_coef[4] = -3.443722207893254e-01;
    TVF.set( TM_1d_coef );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TVF.set( Interval(-4.570756442941450e-02,5.140720829574581e-02) );
#else
    TVF.set( Interval(-8.400835035769887e-02,8.565165355280621e-02) );
#endif
//    std::cout << TVX*exp(-pow(TVX,2)) << std::endl;
//    std::cout << TVF << std::endl;
    CPPUNIT_ASSERT( Eq( TVX*exp(-pow(TVX,2)), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (TVX*exp(-pow(TVX,2))).B(),
                        Interval(-5.320263022173342e-02,4.714114566257437e-01) ) );
#else
    CPPUNIT_ASSERT( Eq( (TVX*exp(-pow(TVX,2))).B(),
                        Interval(-9.150341615001778e-02,5.056559018828041e-01) ) );
#endif

    TM_1d->options.BERNSTEIN_USE = true;
    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::BERNSTEIN;
    TM_1d->options.BOUNDER_ORDER = 0;
    TM_1d_coef[0] =  2.044161987792845e-01;
    TM_1d_coef[1] = -4.167261676526294e-01;
    TM_1d_coef[2] =  8.008080273494700e-03;
    TM_1d_coef[3] =  3.721766955641853e-01;
    TM_1d_coef[4] = -6.714016978805343e-02;
    TVF.set( TM_1d_coef );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TVF.set( Interval(-7.811311073209398e-02,7.624573321767887e-03) );
#else
    TVF.set( Interval(-7.811311073209398e-02,7.624573321767887e-03) );
#endif
    //std::cout << TVX*exp(-pow(TVX,2)) << std::endl;
    CPPUNIT_ASSERT( Eq( TVX*exp(-pow(TVX,2)), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (TVX*exp(-pow(TVX,2))).B(),
                        Interval(-3.773214937698066e-02,3.716875284384643e-01) ) );
#else
    CPPUNIT_ASSERT( Eq( (TVX*exp(-pow(TVX,2))).B(),
                        Interval(-3.773214937698066e-02,3.716875284384643e-01) ) );
#endif
  }

  void testExpression3(){
    Interval::options.DISPLAY_DIGITS=15;
    TM_1d->options.DISPLAY_DIGITS=15;
    TM_1d->options.BOUNDER_ORDER = 0;
    TM_1d->options.BERNSTEIN_USE = false;
    TM_1d->options.REF_MIDPOINT  = true;
    TVX  = TVar<Interval>( TM_1d, 0, Interval(PI/6., PI/3.) );

    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TVar<Interval> TVF;
    TVF.set( TM_1d );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TM_1d_coef[0] =  5.108695620951995e-01;
    TM_1d_coef[1] =  8.396384309211719e-01;
    TM_1d_coef[2] = -3.292783145553981e+01;
    TM_1d_coef[3] =  9.390926766970267e+01;
    TM_1d_coef[4] =  2.014706944034799e+01;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.737024530017051e+01,1.817503069112297e+01) );
#else
    TM_1d_coef[0] =  5.124433481497594e-01;
    TM_1d_coef[1] =  8.705227224518149e-01;
    TM_1d_coef[2] = -3.278767773072964e+01;
    TM_1d_coef[3] =  9.355611183090004e+01;
    TM_1d_coef[4] =  1.895354210210460e+01;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-3.250078574465327e+01,3.314610609542391e+01) );
#endif
    //std::cout << sin(pow(TVX,-3))*cos(sqrt(TVX)) << std::endl;
    CPPUNIT_ASSERT( Eq( sin(pow(TVX,-3))*cos(sqrt(TVX)), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (sin(pow(TVX,-3))*cos(sqrt(TVX))).B(),
                        Interval(-2.102108653472737e+01,2.047095133478485e+01) ) );
#else
    CPPUNIT_ASSERT( Eq( (sin(pow(TVX,-3))*cos(sqrt(TVX))).B(),
                        Interval(-3.614219586378493e+01,3.543208259844760e+01) ) );
#endif

    TM_1d->options.BERNSTEIN_USE = true;
    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::BERNSTEIN;
    TM_1d->options.BOUNDER_ORDER = 0;
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TM_1d_coef[0] =  4.707121568629750e-01;
    TM_1d_coef[1] =  1.136642710335418e-01;
    TM_1d_coef[2] = -3.580366547598283e+01;
    TM_1d_coef[3] =  1.025227823675693e+02;
    TM_1d_coef[4] =   3.921204569420465e+01;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.963527146726179e+01,2.059725243295820e+01) );
#else
    TM_1d_coef[0] =  4.558260020807168e-01;
    TM_1d_coef[1] =  2.231950017408094e+00;
    TM_1d_coef[2] = -4.044683947297301e+01;
    TM_1d_coef[3] =  6.222119101576310e+01;
    TM_1d_coef[4] =  3.197963635244135e+02;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-2.096755493263365e+01,2.088659163331111e+01) );
#endif
    //std::cout << sin(pow(TVX,-3))*cos(sqrt(TVX)) << std::endl;
    CPPUNIT_ASSERT( Eq( sin(pow(TVX,-3))*cos(sqrt(TVX)), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (sin(pow(TVX,-3))*cos(sqrt(TVX))).B(),
                        Interval(-2.330367136044478e+01,2.207014796259254e+01) ) );
#else
    CPPUNIT_ASSERT( Eq( (sin(pow(TVX,-3))*cos(sqrt(TVX))).B(),
                        Interval(-2.348242729007603e+01,2.376874864044477e+01) ) );
#endif
  }

  void testExpression4(){
    Interval::options.DISPLAY_DIGITS=15;
    TM_1d->options.DISPLAY_DIGITS=15;
    TM_1d->options.BOUNDER_ORDER = 0;
    TM_1d->options.BERNSTEIN_USE = false;
    TM_1d->options.REF_MIDPOINT = true;
    TVX  = TVar<Interval>( TM_1d, 0, Interval(0., PI/3.) );

    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TVar<Interval> TVF;
    TVF.set( TM_1d );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TM_1d_coef[0] =  1.453761383551061e+00;
    TM_1d_coef[1] = -6.943048491745324e-01;
    TM_1d_coef[2] = -1.458374062951880e+00;
    TM_1d_coef[3] = -3.012545882573837e-01;
    TM_1d_coef[4] =  2.168153310674373e+00;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.088244798268075e+00,1.033322426312697e+00) );
#else
    TM_1d_coef[0] =  1.453761384163116e+00;
    TM_1d_coef[1] = -6.943048542770530e-01;
    TM_1d_coef[2] = -1.458374060439291e+00;
    TM_1d_coef[3] = -3.012545325252489e-01;
    TM_1d_coef[4] =  2.168153306781894e+00;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.320509119518708e+00,1.254534489319857e+00) );
#endif
    //std::cout << tan(cos(TVX*atan(TVX))) << std::endl;
    CPPUNIT_ASSERT( Eq( tan(cos(TVX*atan(TVX))), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (tan(cos(TVX*atan(TVX)))).B(),
                        Interval(-4.410864796348388e-01,2.775925891013603e+00) ) );
#else
    CPPUNIT_ASSERT( Eq( (tan(cos(TVX*atan(TVX)))).B(),
                        Interval(-6.733507942560369e-01,2.997137947697021e+00) ) );
#endif

    TM_1d->options.BERNSTEIN_USE = true;
    TM_1d->options.BOUNDER_TYPE = TModel<Interval>::Options::BERNSTEIN;
    TM_1d->options.BOUNDER_ORDER = 0;
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TM_1d_coef[0] =  1.453765264415975e+00;
    TM_1d_coef[1] = -6.943420069793143e-01;
    TM_1d_coef[2] = -1.458303602210573e+00;
    TM_1d_coef[3] = -3.010649766275364e-01;
    TM_1d_coef[4] =  2.167972166686784e+00;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.085903017585574e+00,1.026253890852634e+00) );
#else
    TM_1d_coef[0] =  1.453765264989323e+00;
    TM_1d_coef[1] = -6.943420118511072e-01;
    TM_1d_coef[2] = -1.458303599468113e+00;
    TM_1d_coef[3] = -3.010649230187787e-01;
    TM_1d_coef[4] =  2.167972159667138e+00;
    TVF.set( TM_1d_coef );
    TVF.set( Interval(-1.203887915715804e+00,1.137788607381716e+00) );
#endif
    //std::cout << tan(cos(TVX*atan(TVX))) << std::endl;
    CPPUNIT_ASSERT( Eq( tan(cos(TVX*atan(TVX))), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (tan(cos(TVX*atan(TVX)))).B(),
                        Interval(-2.757660546130731e-01,2.776234243192992e+00) ) );
#else
    CPPUNIT_ASSERT( Eq( (tan(cos(TVX*atan(TVX)))).B(),
                        Interval(-3.937509468011575e-01,2.887768959517196e+00) ) );
#endif
  }

  void testExpression5(){
    Interval::options.DISPLAY_DIGITS=15;
    TM_2d->options.DISPLAY_DIGITS=15;
    TM_2d->options.BOUNDER_ORDER = 0;
    TM_2d->options.BERNSTEIN_USE = false;
    TM_2d->options.REF_MIDPOINT  = true;
    TVX1 = TVar<Interval>( TM_2d, 0, Interval(-2., 0.) );
    TVX2 = TVar<Interval>( TM_2d, 1, Interval( 1., 3.) );

    TM_2d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TVar<Interval> TVF;
    TVF.set( TM_2d );
    TM_2d_coef[0]  = -2.348862971128921e-01;
    TM_2d_coef[1]  =  6.665511124293279e-02;
    TM_2d_coef[2]  = -1.576224354004002e-01;
    TM_2d_coef[3]  =  1.203361699643438e-02;
    TM_2d_coef[4]  =  9.492821616831086e-02;
    TM_2d_coef[5]  = -6.530516531879975e-02;
    TM_2d_coef[6]  = -1.808066707272379e-02;
    TM_2d_coef[7]  =  1.165640446481486e-02;
    TM_2d_coef[8]  =  4.946967344108813e-02;
    TM_2d_coef[9]  = -1.601171183587026e-03;
    TM_2d_coef[10] =  2.463676238986581e-03;
    TM_2d_coef[11] = -3.085553162279163e-02;
    TM_2d_coef[12] =  1.503488328577997e-02;
    TM_2d_coef[13] = -1.147814066490426e-02;
    TM_2d_coef[14] =  1.834342816689918e-02;
    TVF.set( TM_2d_coef );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TVF.set( Interval(-4.837665578992553e-01,1.650102051910515e-01) );
#else
    TVF.set( Interval(-2.604336002272144e+03,2.604322374334143e+03) );
#endif
    //std::cout << -1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
    //             -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
    //         -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2) << std::endl;
    CPPUNIT_ASSERT( Eq( -1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                        -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (-1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                         -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                 -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2)).B(),
                        Interval(-1.214271754596066e+00,3.550416987090134e-01) ) );
#else
    CPPUNIT_ASSERT( Eq( (-1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                         -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                 -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2)).B(),
                        Interval(-2.605066507468841e+03,2.604512405827661e+03) ) );
#endif

    TM_2d->options.BERNSTEIN_USE = true;
    TM_2d->options.BOUNDER_TYPE = TModel<Interval>::Options::LSB;
    TM_2d->options.BOUNDER_ORDER = 0;
    TM_2d_coef[0]  = -3.555231260976054e-01;
    TM_2d_coef[1]  =  1.474029473937707e-01;
    TM_2d_coef[2]  = -3.217035195213951e-01;
    TM_2d_coef[3]  =  3.613251092666555e-02;
    TM_2d_coef[4]  =  1.626063439176184e-01;
    TM_2d_coef[5]  = -9.111044391153904e-02;
    TM_2d_coef[6]  = -3.551455782018680e-02;
    TM_2d_coef[7]  =  4.932700853093925e-02;
    TM_2d_coef[8]  =  2.593170712881997e-02;
    TM_2d_coef[9]  =  3.774619342673860e-02;
    TM_2d_coef[10] = -2.193704358801821e-03;
    TM_2d_coef[11] = -3.131188928431736e-02;
    TM_2d_coef[12] =  1.437515771981686e-02;
    TM_2d_coef[13] = -2.400729777855814e-02;
    TM_2d_coef[14] =  1.872155617061304e-02;
    TVF.set( TM_2d_coef ); 
#ifdef MC__TMODEL_TIGHT_REMAINDER
    TVF.set( Interval(-1.515370247137803e-01,3.590323492160009e-01) );
#else
    TVF.set( Interval(-1.515739454501995e-01,3.590738223750071e-01) );
#endif
    //std::cout << -1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
    //             -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
    //             -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2) << std::endl;
    CPPUNIT_ASSERT( Eq( -1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                        -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2), TVF ) );
#ifdef MC__TMODEL_TIGHT_REMAINDER
    CPPUNIT_ASSERT( Eq( (-1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                         -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                 -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2)).B(),
                        Interval(-1.399783252957405e+00,8.171794688262966e-01) ) );
#else
    CPPUNIT_ASSERT( Eq( (-1./(pow(TVX1-4.,2)+pow(TVX2-4.,2)+0.1)
                         -1./(pow(TVX1-1.,2)+pow(TVX2-1.,2)+0.2)
	                 -1./(pow(TVX1-8.,2)+pow(TVX2-8.,2)+0.2)).B(),
                        Interval(-1.399820173693825e+00,8.172209419853024e-01) ) );
#endif
  }

  void testDivisionByZero(){
    // The following line should throw an instance of TModel<Interval>::Exceptions
    TVX1/0.;
  }

  void testTMEnv(){
    // The following line should throw an instance of TModel<Interval>::Exceptions
    TVX  = TVar<Interval>( TM_1d, 0, Interval(0., PI/3.) );
    TVX2 = TVar<Interval>( TM_2d, 1, Interval(0., 1.) );
    TVX+TVX2;
  }

};

CPPUNIT_TEST_SUITE_REGISTRATION( mc::TModelTest );

} // end namespace mc

#endif
