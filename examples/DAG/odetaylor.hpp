// Copyright (C) 2013 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// ***************************************************************
// ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS OF THE
//                         COPYRIGHT NOTICE
// ***************************************************************

#ifndef MC__ODETAYLOR_HPP
#define MC__ODETAYLOR_HPP

#include <cassert>

#include "odestruct.hpp"
#include "mcop.hpp"
#include "mcfadbad.hpp"

//#undef MC__ODETAYLOR_DEBUG

namespace mc
{
//! @brief C++ class doing Taylor expansion of the solutions of parametric ODEs using FADBAD++
////////////////////////////////////////////////////////////////////////
//! mc::ODETAYLOR is a C++ class for Taylor expansion of the solutions
//! of parametric initial value problems in ODEs.
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP>
class ODETAYLOR : public virtual ODESTRUCT
////////////////////////////////////////////////////////////////////////
{
public:

  typedef fadbad::T< U > TU;
  typedef fadbad::F< U > FU;
  typedef fadbad::T< FU > TFU;
  typedef fadbad::F< FU > FFU;
  typedef fadbad::T< FFU > TFFU;

  //! @brief Constructor
  ODETAYLOR();
  //! @brief Destructor
  virtual ~ODETAYLOR();

  //! @brief Compute the Taylor coefficient values of the ODE solution
  template< typename V > void T_expand
    ( const unsigned int is, const double t, const U*x, const V*p,
      const unsigned int nts, U**f );
  //! @brief Compute the Taylor coefficient values as well as first derivatives of the ODE solution
  template< typename V > void TF_expand
    ( const unsigned int is, const double t, const U*x, const V*p,
      const unsigned int nts, U**f, U**dfdx, U**dfdp=0 );
  //! @brief Compute the Taylor coefficient values as well as directional second-order derivatives of the ODE solution
  template< typename V > void TFFd_expand
    ( const unsigned int is, const double t, const U*x, const V*p,
      const unsigned int nts, U**f, const U*d, U**d2fdx2dd );
  //! @brief Compute the invariant values of the ODE system
  void invariant
    ( const unsigned int is, const double t, const U*x, const U*p, U*inv );
  //! @brief Compute the invariant values as well as first derivatives of the ODE system
  void F_invariant
    ( const unsigned int is, const double t, const U*x, const U*p,
      U*inv, U*dinvdx, U*dinvdp=0 );
  //! @brief Compute the invariant values as well as directional second-order derivatives of the ODE system
  void FFd_invariant
    ( const unsigned int is, const double t, const U*x, const U*p,
      U*inv, const U*d, U*d2invdx2dd );

  //! @brief Retrieve/set T counter
  unsigned int& T_calls()
    { return _numT; }
  //! @brief Retrieve T counter
  unsigned int T_calls() const
    { return _numT; }
  //! @brief Retrieve/set TF counter
  unsigned int& TF_calls()
    { return _numTF; }
  //! @brief Retrieve TF counter
  unsigned int TF_calls() const
    { return _numTF; }
  //! @brief Retrieve/set TFFd counter
  unsigned int& TFFd_calls()
    { return _numTFFd; }
  //! @brief Retrieve TFFd counter
  unsigned int TFFd_calls() const
    { return _numTFFd; }
  //! @brief Reset counters
  void reset_calls()
    { _numT = _numTF = _numTFFd = 0; }

private:

  //! @brief Number of TC evaluations
  unsigned int _numT;
  //! @brief Number of TC Jacobian evaluations
  unsigned int _numTF;
  //! @brief Number of TC (directional) Hessian evaluations
  unsigned int _numTFFd;

  //! @brief Parameter values
  TU *Tp;
  //! @brief State values
  TU *Tx;
  //! @brief State time derivatives
  TU *Tdx;

  //! @brief Parameter values
  TFU *TFp;
  //! @brief State values
  TFU *TFx;
  //! @brief State time derivatives
  TFU *TFdx;
  //! @brief Invariant
  TFU *TFinv;

  //! @brief Parameter values
  TFFU *TFFp;
  //! @brief State values
  TFFU *TFFx;
  //! @brief State time derivatives
  TFFU *TFFdx;

  //! @brief Parameter values
  FU *Fp;
  //! @brief State values
  FU *Fx;
  //! @brief Invariant
  FU *Finv;

  //! @brief Parameter values
  FFU *FFp;
  //! @brief State values
  FFU *FFx;
  //! @brief Invariant
  FFU *FFinv;
};

////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP>
inline
ODETAYLOR<U,IVP>::ODETAYLOR():
  ODESTRUCT( IVP() )
////////////////////////////////////////////////////////////////////////   
{
  // Allocate Tx, Tp, Tdx
  Tp = Tx = Tdx = 0;

  // Allocate TFx, TFp, TFdx
  TFp = TFx = TFdx = 0;

  // Allocate TFFx, TFFp, TFFdx
  TFFp = TFFx = TFFdx = 0;

  // Allocate Fx, Fp, Finv
  Fp = Fx = Finv = 0;

  // Allocate FFx, FFp, FFinv
  FFp = FFx = FFinv = 0;
}

////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP> inline
ODETAYLOR<U,IVP>::~ODETAYLOR()
////////////////////////////////////////////////////////////////////////   
{
  delete[] Tp;
  delete[] Tx;
  delete[] Tdx;

  delete[] TFp;
  delete[] TFx;
  delete[] TFdx;

  delete[] TFFp;
  delete[] TFFx;
  delete[] TFFdx;

  delete[] Fp;
  delete[] Fx;
  delete[] Finv;

  delete[] FFp;
  delete[] FFx;
  delete[] FFinv;
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::T_expand(
//! const unsigned int is, const double t, const U*x, const U*p,
//! const unsigned int nts, U**f )
//!
//! This function computes the first \f$n_{ts}\f$ Taylor coefficients
//! \f$f^{[1]},\ldots,f^{[n_{ts}]}\f$ of the ODE solutions \f$x\f$ at
//! time \f$t\f$ and in stage \f$is\f$, for the parameter values \f$p\f$
//!
//! - Input Arguments: is; t; x; p; nts
//! - Output Arguments: f[q], for q=0,...,nts
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP> template <typename V>
inline void
ODETAYLOR<U,IVP>::T_expand
( const unsigned int is, const double t, const U*x, const V*p,
  const unsigned int nts, U**f )
////////////////////////////////////////////////////////////////////////
{
  assert( x && f );

  // Initialize Tt, Tp, Tx and Tdx
  if( !Tp )  Tp  = new TU[_np];
  for( unsigned int i=0; i<_np; i++ ) Tp[i] = p[i];
  if( !Tx )  Tx  = new TU[_nx];
  for( unsigned int i=0; i<_nx; i++ ) Tx[i] = f[0][i] = x[i];
  TU Tt = U(t); Tt[1] = 1.;
  if( !Tdx ) Tdx = new TU[_nx];
  for( unsigned int i=0; i<_nx; i++ ) Tdx[i] = IVP().RHS( i, Tp, Tx, Tt, is );
  _numT++;

  // Taylor coefficients of the solution x (up to order nts)
#ifdef MC__ODETAYLOR_DEBUG
  for( unsigned int i=0; i<_nx; i++ )
    std::cout << "ODETAYLOR::T_expand *** f(" << i << ")[0] = "
              << f[0][i] << std::endl;
#endif
  for( unsigned int q=0; q<nts; q++ ){
    assert( f[q+1] );
    for( unsigned int i=0; i<_nx; i++ ){
      // Evaluate q-th Taylor coefficient for f[i]
      Tdx[i].eval(q);
      // Set result as (q+1)-th Taylor coefficient for x[i]
      f[q+1][i] = Tx[i][q+1] = Tdx[i][q]/double(q+1);
#ifdef MC__ODETAYLOR_DEBUG
      std::cout << "ODETAYLOR::T_expand *** f(" << i << ")[" << q+1 << "] = "
                << f[q+1][i] << std::endl;
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::TF_expand(
//! const unsigned int is, const double t, const U*x, const U*p,
//! const unsigned int nts, U**f, U**dfdx, U**dfdp=0 )
//!
//! This function computes the first \f$k\f$ Taylor coefficients
//! \f$f^{[1]},\ldots,f^{[k]}\f$ of the ODE solutions \f$x\f$ at
//! \f$t\f$ for parameter values \f$p\f$, as well as their Jacobian
//! w.r.t. \f$x\f$ and \f$p\f$
//!
//! - Input Arguments: is; t; x; p; nts
//! - Output Arguments: f[q]; dfdx[q]; dfdp[q], for q=0,...,nts
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP> template <typename V>
inline void
ODETAYLOR<U,IVP>::TF_expand
( const unsigned int is, const double t, const U*x, const V*p,
  const unsigned int nts, U**f, U**dfdx, U**dfdp )
////////////////////////////////////////////////////////////////////////
{
  assert( x && f && dfdx );

  // Initialize TFt, TFp, TFx and TFdx
  const unsigned int nf = (dfdp? _nx+_np: _nx);
  if( !TFp )  TFp  = new TFU[_np];
  for( unsigned int i=0; i<_np; i++ ){
    TFp[i] = p[i];
    if( dfdp ){
      assert( dfdp[0] );
      TFp[i][0].diff(_nx+i,nf);
      for( unsigned int j=0; j<_nx; j++ ) dfdp[0][i*_nx+j] = 0.;
    }
  }
  if( !TFx )  TFx  = new TFU[_nx];
  for( unsigned int i=0; i<_nx; i++ ){
    assert( f[0] && dfdx[0] );
    TFx[i] = f[0][i] = x[i];
    TFx[i][0].diff(i,nf);
    for( unsigned int j=0; j<_nx; j++ ) dfdx[0][i*_nx+j] = TFx[i][0].d(j);
  }
  TFU TFt = U(t); TFt[1] = 1.;
  if( !TFdx ) TFdx = new TFU[_nx];
  for( unsigned int i=0; i<_nx; i++ ) TFdx[i] = IVP().RHS( i, TFp, TFx, TFt, is );
  _numTF++;

  // Taylor coefficients of the solution x (up to order nts)
#ifdef MC__ODETAYLOR_DEBUG
  for( unsigned int i=0; i<_nx; i++ ){
    std::cout << "ODETAYLOR::TF_expand *** f(" << i << ")[0] = "
              << f[0][i] << std::endl;
    std::cout << "ODETAYLOR::TF_expand *** dfdx(" << i << ",*)[0] = ";
    for( unsigned int j=0; j<_nx; j++ ) std::cout << dfdx[0][j*_nx+i] << "  ";
    std::cout << std::endl;
    std::cout << "ODETAYLOR::TF_expand *** dfdp(" << i << ",*)[0] = ";
    for( unsigned int j=0; dfdp && j<_np; j++ ) std::cout << dfdp[0][j*_nx+i] << "  ";
    std::cout << std::endl;
  }
#endif
  for( unsigned int q=0; q<nts; q++ ){
    assert( f[q+1] && dfdx[q+1] && (!dfdp || dfdp[q+1]) );
    for( unsigned int i=0; i<_nx; i++ ){
      // Evaluate q-th Taylor coefficient for f[i]
      TFdx[i].eval(q);
      TFx[i][q+1] = TFdx[i][q]/double(q+1);
      // Set result as (q+1)-th Taylor coefficient for x[i]
      f[q+1][i] = TFx[i][q+1].x();
      for( unsigned int j=0; j<_nx; j++ )
        dfdx[q+1][j*_nx+i] = TFx[i][q+1].d(j);
      for( unsigned int j=0; dfdp && j<_np; j++ )
        dfdp[q+1][j*_nx+i] = TFx[i][q+1].d(_nx+j);
#ifdef MC__ODETAYLOR_DEBUG
      std::cout << "ODETAYLOR::TF_expand *** f(" << i << ")[" << q+1 << "] = "
                << f[q+1][i] << std::endl;
      std::cout << "ODETAYLOR::TF_expand *** dfdx(" << i << ",*)[" << q+1 << "] = ";
      for( unsigned int j=0; j<_nx; j++ ) std::cout << dfdx[q+1][j*_nx+i] << "  ";
      std::cout << std::endl;
      std::cout << "ODETAYLOR::TF_expand *** dfdp(" << i << ",*)[" << q+1 << "] = ";
      for( unsigned int j=0; dfdp && j<_np; j++ ) std::cout << dfdp[q+1][j*_nx+i] << "  ";
      std::cout << std::endl;
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::TFFd_expand(
//! const unsigned int is, const double t, const U*x, const U*p,
//! const unsigned int nts, U**f, U*d, U**d2fdx2dd )
//!
//! This function computes the first \f$k\f$ Taylor coefficients
//! \f$f^{[1]},\ldots,f^{[k]}\f$ of the ODE solutions \f$x\f$ at
//! \f$t\f$ for parameter values \f$p\f$, as well as their directional
//! second-order derivatives \f$d^T\frac{\partial^2 f_i^{[k]}}{\partial x^2}d\f$
//! w.r.t. \f$x\f$ in the direction \f$d\f$.
//!
//! - Input Arguments: is; t; x; p; d; nts
//! - Output Arguments: f[q]; d2fdx2dd[q], for q=0,...,nts
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP> template <typename V>
inline void
ODETAYLOR<U,IVP>::TFFd_expand
( const unsigned int is, const double t, const U*x, const V*p,
  const unsigned int nts, U**f, const U*d, U**d2fdx2dd )
////////////////////////////////////////////////////////////////////////
{
  assert( x && f && d && d2fdx2dd );

  // Initialize TFt, TFp, TFx and TFdx
  if( !TFFp ) TFFp = new TFFU[_np];
  for( unsigned int i=0; i<_np; i++ )
    TFFp[i] = p[i];   
  if( !TFFx ) TFFx = new TFFU[_nx];
  assert( f[0] && d2fdx2dd[0] );
  for( unsigned int i=0; i<_nx; i++ ){
    TFFx[i] = f[0][i] = x[i];
    TFFx[i][0].diff(0,1) = d[i];
    TFFx[i][0].x().diff(0,1) = d[i];
    d2fdx2dd[0][i] = TFFx[i][0].d(0).d(0);
  }
  TFFU TFFt = U(t); TFFt[1] = 1.;
  if( !TFFdx ) TFFdx = new TFFU[_nx];
  for( unsigned int i=0; i<_nx; i++ )
    TFFdx[i] = IVP().RHS( i, TFFp, TFFx, TFFt, is );
  _numTFFd++;

  // Taylor coefficients of the solution x (up to order nts)
#ifdef MC__ODETAYLOR_DEBUG
  for( unsigned int i=0; i<_nx; i++ ){
    std::cout << "ODETAYLOR::TFFd_expand *** f(" << i << ")[0] = "
              << f[0][i] << std::endl;
    std::cout << "ODETAYLOR::TFFd_expand *** d2fdx2dd(" << i << ")[0] = "
              << d2fdx2dd[0][i] << std::endl;
  }
#endif
  for( unsigned int q=0; q<nts; q++ ){
    assert( f[q+1] && d2fdx2dd[q+1] );
    for( unsigned int i=0; i<_nx; i++ ){
      // Evaluate q-th Taylor coefficient for f[i]
      TFFdx[i].eval(q);
      TFFx[i][q+1] = TFFdx[i][q]/double(q+1);
      // Set result as (q+1)-th Taylor coefficient for x[i]
      f[q+1][i] = TFFx[i][q+1].x().x();
      d2fdx2dd[q+1][i] = TFFx[i][q+1].d(0).d(0);
#ifdef MC__ODETAYLOR_DEBUG
      std::cout << "ODETAYLOR::TFFd_expand *** f(" << i << ")[" << q+1 << "] = "
                << f[q+1][i] << std::endl;
      std::cout << "ODETAYLOR::TFFd_expand *** d2fdx2dd(" << i << ")[" << q+1 << "] = "
                << d2fdx2dd[q+1][i] << std::endl;
#endif
    }
  }
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::invariant(
//! const unsigned int is, const double t, const U*x, const U*p, U*inv )
//!
//! This function computes the invariants \f$\Omega\f$ of the ODE system
//! at \f$t\f$ and in stage \f$is\f$, for parameter values \f$p\f$ and
//! state values \f$x\f$
//!
//! - Input Arguments: is; t; x; p
//! - Output Arguments: inv
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP>
inline void
ODETAYLOR<U,IVP>::invariant
( const unsigned int is, const double t, const U*x, const U*p, U*inv )
////////////////////////////////////////////////////////////////////////
{
  assert( x && inv );

  // Invariant value and directional 2nd-order derivative
  for( unsigned int i=0; i<_ni; i++ ){
    inv[i] = IVP().INV( i, p, x, t, is );
#ifdef MC__ODETAYLOR_DEBUG
    std::cout << "ODETAYLOR::invariant *** inv(" << i << ") = "
              << inv[i] << std::endl;
#endif
  }
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::F_invariant(
//! const unsigned int is, const double t, const U*x, const U*p,
//! U*inv, U*dinvdx, U*dinvdp=0 )
//!
//! This function computes the invariants \f$\Omega\f$ of the ODE system
//! at \f$t\f$ for parameter values \f$p\f$ and state values \f$x\f$,
//! as well as their Jacobian w.r.t. \f$x\f$ and \f$p\f$.
//! 
//! - Input Arguments: is; t; x; p
//! - Output Arguments: inv; dinvdx; dinvdp
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP>
inline void
ODETAYLOR<U,IVP>::F_invariant
( const unsigned int is, const double t, const U*x, const U*p,
  U*inv, U*dinvdx, U*dinvdp )
////////////////////////////////////////////////////////////////////////
{
  assert( x && inv && dinvdx );

  // Initialize Ft, Fp, Fx
  const unsigned int nf = (dinvdp? _nx+_np: _nx);
  if( !Fp ) Fp = new FU[_np];
  for( unsigned int i=0; i<_np; i++ ){
    Fp[i] = p[i];
    if( dinvdp ) Fp[i].diff(_nx+i,nf);
  }
  if( !Fx ) Fx = new FU[_nx];
  for( unsigned int i=0; i<_nx; i++ ){
    Fx[i] = x[i];
    Fx[i].diff(i,nf);
  }
  FU Ft = U(t);
  if( !Finv ) Finv = new FU[_nx];

  // Invariant value and 1st-order derivative
  for( unsigned int i=0; i<_ni; i++ ){
    Finv[i] = IVP().INV( i, Fp, Fx, Ft, is );
    inv[i] = Finv[i].x();
    for( unsigned int j=0; j<_nx; j++ )
      dinvdx[j*_ni+i] = Finv[i].d(j);
    for( unsigned int j=0; dinvdp && j<_np; j++ )
      dinvdp[j*_ni+i] = Finv[i].d(_nx+j);

#ifdef MC__ODETAYLOR_DEBUG
    std::cout << "ODETAYLOR::F_invariant *** inv(" << i << ") = "
              << inv[i] << std::endl;
    std::cout << "ODETAYLOR::F_invariant *** dinvdx(" << i << ",*) = ";
    for( unsigned int j=0; j<_nx; j++ ) std::cout << dinvdx[j*_ni+i] << "  ";
    std::cout << std::endl;
    std::cout << "ODETAYLOR::F_invariant *** dinvdp(" << i << ",*) = ";
    for( unsigned int j=0; dinvdp && j<_np; j++ ) std::cout << dinvdp[j*_ni+i] << "  ";
#endif
  }
  //{ int dum; std::cin >> dum; }
}

////////////////////////////////////////////////////////////////////////
//! @fn template <typename U, typename IVP> void ODETAYLOR<U,IVP>::FFd_invariant(
//! const unsigned int is, const double t, const U*x, const U*p,
//! U*inv, U*d, U*d2invdx2dd )
//!
//! This function computes the invariants \f$\Omega\f$ of the ODE system
//! at \f$t\f$ for parameter values \f$p\f$ and state values \f$x\f$,
//! as well as their directional second-order derivatives
//! \f$d^T\frac{\partial^2 \Omega_i}{\partial x^2}d\f$
//! w.r.t. \f$x\f$ in the direction \f$d\f$.
//! 
//! - Input Arguments: is; t; x; p; d
//! - Output Arguments: inv; d2invdx2dd
////////////////////////////////////////////////////////////////////////
template <typename U, typename IVP>
inline void
ODETAYLOR<U,IVP>::FFd_invariant
( const unsigned int is, const double t, const U*x, const U*p,
  U*inv, const U*d, U*d2invdx2dd )
////////////////////////////////////////////////////////////////////////
{
  assert( x && inv && d && d2invdx2dd );

  // Initialize TFt, TFp, TFx
  if( !FFp ) FFp = new FFU[_np];
  for( unsigned int i=0; i<_np; i++ )
    FFp[i] = p[i];   
  if( !FFx ) FFx = new FFU[_nx];
  for( unsigned int i=0; i<_nx; i++ ){
    FFx[i] = x[i];
    FFx[i].diff(0,1) = d[i];
    FFx[i].x().diff(0,1) = d[i];
  }
  FFU FFt = U(t);
  if( !FFinv ) FFinv = new FFU[_ni];

  // Invariant value and directional 2nd-order derivative
  for( unsigned int i=0; i<_ni; i++ ){
    FFinv[i] = IVP().INV( i, FFp, FFx, FFt, is );
    inv[i] = FFinv[i].x().x();
    d2invdx2dd[i] = FFinv[i].d(0).d(0);
#ifdef MC__ODETAYLOR_DEBUG
    std::cout << "ODETAYLOR::FFd_invariant *** inv(" << i << ") = "
              << inv[i] << std::endl;
    std::cout << "ODETAYLOR::FFd_invariant *** d2invdx2dd(" << i << ") = "
              << d2invdx2dd[i] << std::endl;
#endif
  }
}

} // namespace mc

#endif
