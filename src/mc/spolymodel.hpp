// Copyright (C) 2021 Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__SPOLYMODEL_H
#define MC__SPOLYMODEL_H

#include <iostream>
#include <iomanip>
#include <typeinfo>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <stdarg.h>
#include <cassert>
#include <climits>
#include <limits>
#include <stdlib.h>
#include <complex>
#include <numeric>

#include "mcfunc.hpp"
#include "mclapack.hpp"
#include "mcop.hpp"
#include "spolymon.hpp"

#undef  MC__SPOLYMODEL_DEBUG
#define MC__SPOLYMODEL_DEBUG_POLYBOUND
#define MC__SPOLYMODEL_CHECK

namespace mc
{

//! @brief C++ base class for the computation of sparse polynomial models for factorable functions: Variable
////////////////////////////////////////////////////////////////////////
//! mc::SPolyVar is a C++ base class for definition of polynomial model
//! variables in sparse format.
////////////////////////////////////////////////////////////////////////
template <typename T>
class SPolyVar
////////////////////////////////////////////////////////////////////////
{
public:

  
public:
  /** @addtogroup POLYMOD Polynomial Model Arithmetic for Factorable Functions
   *  @{
   */

  //! @brief Set multivariate polynomial coefficients in variable as <tt>coefmon</tt>
  virtual SPolyVar<T>& set
    ( t_coefmon& coefmon )
    { _coefmon = coefmon; _unset_bndT(); _unset_bndpol();
      return *this; } // this is assuming the same order and number of variables

  //! @brief Set remainder term in variable as <tt>bndrem</tt>
  virtual SPolyVar<T>& set
    ( T const& bndrem )
    { _bndrem = bndrem; return *this; }

  //! @brief Set multivariate polynomial coefficients and remainder term equal to those in variable <tt>var</tt>, possibly defined in another polynomial model environment with fewer variables or with a different expansion order. Coefficients involving other variables or higher order are initialized to 0 if <tt>reset=true</tt> (default), otherwise they are left unmodified. Higher-order terms in TV are bounded and added to the remainder bound.
  virtual SPolyVar<T>& set
    ( SPolyVar<T> const& var, bool const reset=true );
/*
  //! @brief Copy multivariate polynomial coefficients from current variable into variable <tt>var</tt>, possibly defined in another polynomial model environment with less variables or with a lower expansion order. Copied coefficients are reset to 0 in current Taylor variable if <tt>reset=true</tt>, otherwise they are left unmodified (default).
  virtual SPolyVar<T>& get
    ( SPolyVar<T>&var, const bool reset=false );
*/




  /** @} */
};

///////////////////////////////// SPolyVar /////////////////////////////////////

template <typename T> inline SPolyVar<T>&
SPolyVar<T>::set
( const SPolyVar<T>& var, const bool reset )
{
  if( reset ){
    _ndxvar.clear();
    _coefmon.clear();
  }
  auto&& itord = var._coefmon.upper_bound( SPolyMon( maxord()+1, std::map<unsigned,unsigned>() ) );
  for( auto&& it=var._coefmon.begin(); it!=itord; ++it ){
    for( auto&& [var,order] : it->first.expr ) _ndxvar.insert( var );
    _coefmon.insert( *it );
  }
  _bndrem  = var._bndrem + var.bndord(maxord()+1);  // var may have a higher order than *this
  _unset_bndT();
  _unset_bndpol();

  return *this;
}
/*
template <typename T> inline SPolyVar<T>&
SPolyVar<T>::get
( SPolyVar<T>& var, const bool reset )
{

  if( !_PM ){
    var._ndxmon.clear();
    if( var._PM && var._PM->sparse() ) var._ndxmon.insert(0);
    var._coefmon[0] = _coefmon[0];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[0] = 0.;
    return *this;
  }
  if( !var._PM || nvar() < var.nvar() || nord() < var.nord() ) // looks fishy...
    return *this;

  // Copy monomial coefficients from *this into var
  unsigned*iexp = new unsigned[nvar()];
  // Case: sparse representation
  for( auto it=_ndxmon.begin(); it!=_ndxmon.end() && *it<_posord(nord()+1); ++it ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(*it)[ivar]: 0 );
    //unsigned imon = _loc_expmon(iexp);
    auto pmon = _ndxmon.insert( _loc_expmon(iexp) );
    var._coefmon[*it] = _coefmon[*(pmon.first)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[*(pmon.first)] = 0.;
  }
  // Case: dense representation
  for( unsigned jmon=0; _ndxmon.empty() && jmon<var.nmon(); jmon++ ){
    for( unsigned ivar=0; ivar<nvar(); ivar++ )
      iexp[ivar] = ( ivar<var.nvar()? var._expmon(jmon)[ivar]: 0 );
    var._coefmon[jmon] = _coefmon[_loc_expmon(iexp)];
    // Reset monomial coefficients to 0
    if( reset ) _coefmon[_loc_expmon(iexp)] = 0.;
  }
  delete[] iexp;
  *var._bndrem = 0.;
  var._unset_bndT();
  var._unset_bndpol();
  var._bndord_uptd = false;

  if( reset ){
    _unset_bndT();
    _unset_bndpol();
    _bndord_uptd = false;
  }

  return *this;
}
*/

} // namespace mc

#endif

