// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

/*!
\page page_FFPWLU Piecewise-linear univariates as external DAG operations
\author Benoit Chachuat

The class mc::FFPWLU defines continuous piecewise-linear functions of a single variable as an
external operation in a DAG: one or several functions of the same operand, each an mc::PWLU
extended linearly beyond its domain by its first and last segments, so that it is defined on the
whole real line. Such functions arise as the edges of Kolmogorov-Arnold networks, as tabulated data
or as cost curves.

\code
  #include "ffpwlu.hpp"
  typedef mc::Interval I;

  mc::FFGraph DAG;
  mc::FFVar X( &DAG );
  mc::PWLU F1( { -1., 0., 2. }, { 1., 0., 1. }, 0., 0.5 );    // breakpoints, values, end slopes
  mc::PWLU F2 = mc::PWLU::relu( 0.1, -0.5, { 1.2, -0.4 }, { 1., -2. }, { 0.3, 0.2 } );
  mc::FFPWLU<I> OpF;
  std::vector<mc::FFVar> Y = OpF( X, { F1, F2 } );            // Y[0] = F1(X), Y[1] = F2(X)
\endcode

mc::PWLU::relu converts the form \f$c + v x + \sum_r a_r \max(w_r x + b_r, 0)\f$ exactly.

The operation is evaluated:
- in double precision and in the forward and reverse modes of AD, with the right derivative at a
  breakpoint (mc::PWLU::val, mc::PWLU::slope);
- symbolically, its derivative being \f$s_0 + \sum_k \Delta s_k\, {\rm fstep}(x - x_k)\f$, with
  \f$s_0\f$ the slope of the first segment and \f$\Delta s_k\f$ the change of slope at breakpoint
  \f$x_k\f$;
- in interval arithmetic, exactly, as the range of the function restricted to the range of the
  operand (mc::PWLU::restrict); backward, the operand is contracted to the hull of the preimage
  of the range of each result;
- in McCormick relaxations, with the convex and concave envelopes of each function on the range of
  the operand (mc::PWLU::hull);
- in polyhedral relaxations (mc::PolImg), by an exact encoding of each function on the range of its
  operand, chosen by FFPWLU::Options::ENCODING: the convex hull of the graph by cuts, without
  auxiliary variables (HULL); or on a subdivision of the operand at the breakpoints, with
  continuous (CONT), binary (BIN) or SOS2 (SOS2) auxiliary variables (see
  mc::PolImg::add_subdivision). The functions of an operation share one subdivision, at the union
  of their breakpoints, unless FFPWLU::Options::SEPARATE is set;
- in superposition models (mc::SupVar), exactly where the operand is a single-dependency
  estimator whose under- and overestimators coincide, such as a variable of the model: each result
  is then the composite of the function with that estimator, as an mc::PWLU summand with the
  crossings of its breakpoints added (mc::PWLU::compose), or as an mc::PWCU summand with the exact
  range of the function on each cell; as after any composition, the summands are then reduced to
  mc::SupModel::Options::MAX_SUBDIV segments, which leaves them exact as long as the composite has
  no more segments; other operands go through the form below;
- in Chebyshev models, the other superposition models and the arithmetics of mc::OCVar, through
  the form \f$y_0 + s_0 (x - x_0) + \sum_k \Delta s_k \max(x - x_k, 0)\f$.

In a polyhedral relaxation, breakpoints closer than mc::PolImg::Options::BREAKPOINT_ATOL and
BREAKPOINT_RTOL to one another or to a bound of the operand are merged by mc::PolImg, and the
encoding is then that of the interpolant at the remaining breakpoints.
*/

#ifndef MC__FFPWLU_HPP
#define MC__FFPWLU_HPP

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

#include "ffunc.hpp"
#include "ffdep.hpp"
#include "ffexpr.hpp"
#include "ffinv.hpp"
#include "slift.hpp"
#include "interval.hpp"
#include "mccormick.hpp"
#include "polimage.hpp"
#include "pwcu.hpp"
#include "pwlu.hpp"
#include "scmodel.hpp"
#include "supmodel.hpp"
#include "ocbase.hpp"

namespace mc
{

//! @brief C++ class defining piecewise-linear univariates as external DAG
//! operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFPWLU is a C++ class defining one or several continuous
//! piecewise-linear functions of the same variable, each an mc::PWLU
//! extended linearly beyond its domain, as an external DAG operation in
//! MC++. The template parameter specifies the type for interval
//! arithmetic. See \ref page_FFPWLU.
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class FFPWLU : public FFOp
////////////////////////////////////////////////////////////////////////
{
 public:
  //! @brief Options of mc::FFPWLU
  struct Options
  {
    //! @brief Constructor
    Options() { reset(); }
    //! @brief Encodings in a polyhedral relaxation, numbered as
    //! mc::PolImg::Options::REFINE
    enum ENCODING_TYPE
    {
      IMAGE = -1,  //!< As set by mc::PolImg::Options::BREAKPOINT_TYPE, with HULL for NONE
      HULL  = 0,   //!< Convex hull of the graph by cuts, without auxiliary variables
      CONT,        //!< Subdivision with continuous auxiliary variables
      BIN,         //!< Subdivision with binary auxiliary variables
      SOS2         //!< Subdivision with SOS2 auxiliary variables
    };
    //! @brief Reset to default options
    void
    reset()
    {
      ENCODING = IMAGE;
      SEPARATE = false;
    }
    //! @brief Encoding in a polyhedral relaxation - Default: IMAGE
    int ENCODING;
    //! @brief Whether each function of an operation is encoded on its own
    //! subdivision, at its own breakpoints, rather than all on one, at the
    //! union of their breakpoints - Default: false
    bool SEPARATE;
  } options;

  //! @brief Default constructor
  FFPWLU() : FFOp(EXTERN) {}

  //! @brief Copy constructor
  FFPWLU(FFPWLU<T> const& other) : FFOp(other), options(other.options), _F(other._F) {}

  //! @brief Destructor
  virtual ~FFPWLU() {}

  //! @brief Define the functions <a>F</a> of <a>X</a>, one result each
  std::vector<FFVar>
  operator()(FFVar const& X, std::vector<PWLU> const& F)
  {
    FFVar** ppRes = _set(X, F);
    std::vector<FFVar> vRes(F.size());
    for (size_t j = 0; j < F.size(); ++j) vRes[j] = *ppRes[j];
    return vRes;
  }

  //! @brief Define the function <a>F</a> of <a>X</a>
  FFVar&
  operator()(FFVar const& X, PWLU const& F)
  {
    return *_set(X, std::vector<PWLU>(1, F))[0];
  }

  //! @brief Define the functions <a>F</a> of <a>X</a> and return result <a>idep</a>
  FFVar&
  operator()(size_t const idep, FFVar const& X, std::vector<PWLU> const& F)
  {
    if (idep >= F.size()) throw std::runtime_error("mc::FFPWLU ** Index of result out of range\n");
    return *_set(X, F)[idep];
  }

  //! @brief Functions of the operation
  std::vector<PWLU> const&
  functions() const
  {
    return _F;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type(idU, typeid(FFVar)))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar, static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FADType<FFVar>)))
      return eval(nRes, static_cast<FADType<FFVar>*>(vRes), nVar,
                  static_cast<FADType<FFVar> const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFDep)))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar, static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFInv)))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar, static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(double)))
      return eval(nRes, static_cast<double*>(vRes), nVar, static_cast<double const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FADType<double>)))
      return _evallin(nRes, static_cast<FADType<double>*>(vRes),
                      static_cast<FADType<double> const*>(vVar));
    else if (mc::same_type(idU, typeid(BADType<double>)))
      return _evallin(nRes, static_cast<BADType<double>*>(vRes),
                      static_cast<BADType<double> const*>(vVar));
    else if (mc::same_type(idU, typeid(T)))
      return eval(nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(SCVar<T>)))
      return _eval(nRes, static_cast<SCVar<T>*>(vRes), static_cast<SCVar<T> const*>(vVar));
    else if (mc::same_type(idU, typeid(McCormick<T>)))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(SupVar<PWCU>)))
      return _evalsup(nRes, static_cast<SupVar<PWCU>*>(vRes), static_cast<SupVar<PWCU> const*>(vVar));
    else if (mc::same_type(idU, typeid(SupVar<PWLU>)))
      return _evalsup(nRes, static_cast<SupVar<PWLU>*>(vRes), static_cast<SupVar<PWLU> const*>(vVar));
    else if (mc::same_type(idU, typeid(PolVar<T>)))
      return eval(nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T> const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(OCVar<double>)))
      return _eval(nRes, static_cast<OCVar<double>*>(vRes), static_cast<OCVar<double> const*>(vVar));
    else if (mc::same_type(idU, typeid(OCVar<FADType<double>>)))
      return _eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes),
                   static_cast<OCVar<FADType<double>> const*>(vVar));
    else if (mc::same_type(idU, typeid(OCVar<FFDep>)))
      return _eval(nRes, static_cast<OCVar<FFDep>*>(vRes), static_cast<OCVar<FFDep> const*>(vVar));
    else if (mc::same_type(idU, typeid(SLiftVar)))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar, static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFExpr)))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar, static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error("FFPWLU::feval: **ERROR** No evaluation method with type" +
                             std::string(idU.name()) + "\n");
  }

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
            FADType<FFVar> const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar, FFInv const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, double* vRes, size_t const nVar, double const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, T* vRes, size_t const nVar, T const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, McCormick<T>* vRes, size_t const nVar,
            McCormick<T> const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
            unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes, unsigned const nVar,
        void* vVar) const
  {
    if (mc::same_type(idU, typeid(T)))
      return reval(nRes, static_cast<T*>(vRes), nVar, static_cast<T*>(vVar));
    else if (mc::same_type(idU, typeid(PolVar<T>)))
      return reval(nRes, static_cast<PolVar<T>*>(vRes), nVar, static_cast<PolVar<T>*>(vVar));

    throw std::runtime_error("FFPWLU::reval: **ERROR** No evaluation method with type" +
                             std::string(idU.name()) + "\n");
  }

  bool reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const;

  bool reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T>* vVar) const;

  // Derivatives
  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
             FFVar** vDer) const;

  // Ordering
  bool lt(FFOp const* op) const;

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << _F.data();
    return "PWLU[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }

 private:
  //! @brief Functions of the operation
  std::vector<PWLU> _F;

  // Define the operation in the DAG
  FFVar**
  _set(FFVar const& X, std::vector<PWLU> const& F)
  {
    if (F.empty()) throw std::runtime_error("mc::FFPWLU ** No function\n");
    for (auto const& f : F)
      if (!f.size() || f.size() != f.dy().size())
        throw std::runtime_error("mc::FFPWLU ** Function without segment\n");
    _F            = F;
    this->owndata = false;
    this->data    = nullptr;
    return this->insert_external_operation(*this, F.size(), X);
  }

  // Evaluation in arithmetic U through the form y0 + s0 (x-x0) + sum_k ds_k max(x-xk,0), with
  // max(z,0) = (z+|z|)/2: Chebyshev models and mc::OCVar have fabs but no max with a constant, and
  // in superposition models fabs (a convex non-monotonic composition) gives tighter estimators
  // than max with a constant where the operand depends on several variables
  template <typename U>
  void
  _eval(size_t const nRes, U* vRes, U const* vVar) const
  {
    for (size_t j = 0; j < nRes; ++j)
    {
      PWLU const& f = _F[j];
      vRes[j]       = f.yL() + f.dy().front() * (vVar[0] - f.xL());
      double xk     = f.xL();
      for (size_t k = 0; k + 1 < f.size(); ++k)
      {
        xk += f.dx()[k];
        double const ds = f.dy()[k + 1] - f.dy()[k];
        if (ds == 0.) continue;
        U const z = vVar[0] - xk;
        vRes[j] += (0.5 * ds) * (z + fabs(z));
      }
    }
  }

  // Summand of a function of an exact univariate estimator: the PWLU composite, exact
  static PWLU
  _summand(PWLU const& g, PWLU const& u)
  {
    return PWLU(u).compose(g);
  }

  // Summand of a function of an exact univariate estimator: in each cell, the exact range of g on
  // the range of the estimator, and slope bounds from the range of the slopes of g there
  static PWCU
  _summand(PWLU const& g, PWCU const& u)
  {
    PWCU r(u);
    size_t const N = u.yL().size();
    bool const slope = !u.y0().empty();
    for (size_t k = 0; k < N; ++k)
    {
      PWLU const R = g.restrict(u.yL()[k], u.yU()[k]);
      r.yL()[k] = R.l();
      r.yU()[k] = R.u();
      if (!slope) continue;
      auto const [itL, itU] = std::minmax_element(R.dy().cbegin(), R.dy().cend());
      double const sL = *itL, sU = *itU;
      double const sLL_sL = u.sLL()[k] * sL, sLL_sU = u.sLL()[k] * sU, sUL_sL = u.sUL()[k] * sL,
                   sUL_sU = u.sUL()[k] * sU;
      r.sLL()[k] = std::min(std::min(sLL_sL, sLL_sU), std::min(sUL_sL, sUL_sU));
      r.sUL()[k] = std::max(std::max(sLL_sL, sLL_sU), std::max(sUL_sL, sUL_sU));
      double const sLU_sL = u.sLU()[k] * sL, sLU_sU = u.sLU()[k] * sU, sUU_sL = u.sUU()[k] * sL,
                   sUU_sU = u.sUU()[k] * sU;
      r.sLU()[k] = std::max(std::max(sLU_sL, sLU_sU), std::max(sUU_sL, sUU_sU));
      r.sUU()[k] = std::min(std::min(sLU_sL, sLU_sU), std::min(sUU_sL, sUU_sU));
    }
    if (slope)
      for (size_t k = 0; k <= N; ++k) r.y0()[k] = g.val(u.y0()[k], true);
    return r;
  }

  // Whether two summands are the same estimator
  static bool
  _same(PWLU const& a, PWLU const& b)
  {
    return a.xL() == b.xL() && a.yL() == b.yL() && a.dx() == b.dx() && a.dy() == b.dy();
  }
  static bool
  _same(PWCU const& a, PWCU const& b)
  {
    return a.xL() == b.xL() && a.xU() == b.xU() && a.yL() == b.yL() && a.yU() == b.yU() &&
           a.y0() == b.y0() && a.sLL() == b.sLL() && a.sLU() == b.sLU() && a.sUL() == b.sUL() &&
           a.sUU() == b.sUU();
  }

  // Evaluation in superposition arithmetic: where the operand is a single-dependency estimator
  // whose under- and overestimators coincide (a variable of the superposition model, say), each
  // result is the summand of g composed with it; otherwise through the form of _eval
  template <typename S>
  void
  _evalsup(size_t const nRes, SupVar<S>* vRes, SupVar<S> const* vVar) const
  {
    SupVar<S> const& X = vVar[0];
    if (X.ndep() == 1 && X.uest(1).empty() && X.oest(1).empty() && X.l(2) == X.l(0) &&
        X.u(2) == X.u(0))
    {
      unsigned const i = *X.sdep().cbegin();
      if (_same(X.uest()[i], X.oest()[i]))
      {
        size_t const nseg = X.mod()->options.MAX_SUBDIV;
        for (size_t j = 0; j < nRes; ++j)
        {
          // reduced to SupModel::Options::MAX_SUBDIV segments, as after any composition: the
          // under- and overestimators then differ where the composite has more segments
          S const est = _summand(_F[j], X.uest()[i]);
          vRes[j]     = X;
          vRes[j].uest()[i] = est;
          vRes[j].uest()[i].reduce(true, nseg);
          vRes[j].oest()[i] = est;
          vRes[j].oest()[i].reduce(false, nseg);
          vRes[j].l(3);  // bounds recomputed from the new summands
          vRes[j].u(3);
        }
        return;
      }
    }
    _eval(nRes, vRes, vVar);
  }

  // Evaluation in AD arithmetic U: tangent of each function at the value of the operand
  template <typename U>
  void
  _evallin(size_t const nRes, U* vRes, U const* vVar) const
  {
    double const x = vVar[0].val();
    for (size_t j = 0; j < nRes; ++j)
    {
      double const s = _F[j].slope(x);
      vRes[j]        = s * vVar[0] + (_F[j].val(x, true) - s * x);
    }
  }

  // Value and slope of the envelope H at t; at a breakpoint, the slope of the piece towards ref
  static std::pair<double, double>
  _hull_eval(PWLU const& H, double const t, double const ref)
  {
    double const s = (ref < t ? H.slope(std::nextafter(t, -INFINITY)) : H.slope(t));
    return {H.val(t, true), s};
  }

  // Encoding type in a polyhedral image
  int
  _encoding(PolImg<T> const* img) const
  {
    if (options.ENCODING != Options::IMAGE) return options.ENCODING;
    switch (img->options.BREAKPOINT_TYPE)
    {
      case PolImg<T>::Options::CONT: return Options::CONT;
      case PolImg<T>::Options::BIN:  return Options::BIN;
      case PolImg<T>::Options::SOS2: return Options::SOS2;
      default:                        return Options::HULL;
    }
  }

  // Cuts of the subdivision encoding of the functions R[j], j in J, restricted to [l,u]
  void _subdiv_cuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& X, double const l,
                    double const u, std::vector<PWLU> const& R, std::vector<size_t> const& J,
                    PolVar<T> const* vRes, int const enc) const;
};

////////////////////////////////////////////////////////////////////////

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
                unsigned const* mVar) const
{
  FFVar** ppRes = this->insert_external_operation(*this, nRes, vVar[0]);
  for (size_t j = 0; j < nRes; ++j) vRes[j] = *ppRes[j];
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
                FADType<FFVar> const* vVar, unsigned const* mVar) const
{
  FFVar const X = vVar[0].val();
  FFVar** ppRes = this->insert_external_operation(*this, nRes, X);
  std::vector<FFVar> vDer(nRes);
  std::vector<FFVar*> pDer(nRes);
  for (size_t j = 0; j < nRes; ++j) pDer[j] = &vDer[j];
  deriv(nRes, nullptr, 1, &X, pDer.data());
  for (size_t j = 0; j < nRes; ++j)
  {
    vRes[j] = *ppRes[j];
    vRes[j].setDepend(vVar[0]);
    for (unsigned k = 0; k < vRes[j].size(); ++k)
    {
      if (vVar[0][k].cst() && vVar[0][k].num().val() == 0.)
        vRes[j][k] = 0.;
      else
        vRes[j][k] = vDer[j] * vVar[0][k];
    }
  }
}

template <typename T>
inline void
FFPWLU<T>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar, FFVar const* vVar,
                 FFVar** vDer) const
{
  for (size_t j = 0; j < nRes; ++j)
  {
    PWLU const& f = _F[j];
    FFVar D       = f.dy().front();
    double xk     = f.xL();
    for (size_t k = 0; k + 1 < f.size(); ++k)
    {
      xk += f.dx()[k];
      double const ds = f.dy()[k + 1] - f.dy()[k];
      if (ds == 0.) continue;
      D += ds * fstep(vVar[0] - xk);
    }
    vDer[j][0] = D;
  }
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
                unsigned const* mVar) const
{
  for (size_t j = 0; j < nRes; ++j)
  {
    vRes[j] = vVar[0];
    if (!_F[j].bkpts().empty()) vRes[j].update(FFDep::TYPE::N);
  }
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar, FFInv const* vVar,
                unsigned const* mVar) const
{
  for (size_t j = 0; j < nRes; ++j)
  {
    vRes[j] = vVar[0];
    vRes[j].update(FFInv::TYPE::U);  // Not a candidate for inversion
  }
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, double* vRes, size_t const nVar, double const* vVar,
                unsigned const* mVar) const
{
  for (size_t j = 0; j < nRes; ++j) vRes[j] = _F[j].val(vVar[0], true);
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, T* vRes, size_t const nVar, T const* vVar,
                unsigned const* mVar) const
{
  double const l = Op<T>::l(vVar[0]), u = Op<T>::u(vVar[0]);
  for (size_t j = 0; j < nRes; ++j)
  {
    PWLU const R = _F[j].restrict(l, u);
    vRes[j]      = T(R.l(), R.u());
  }
}

template <typename T>
inline bool
FFPWLU<T>::reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const
{
  // Operand contracted to the hull of the preimage of the range of each result
  for (size_t j = 0; j < nRes; ++j)
  {
    double const yl = Op<T>::l(vRes[j]), yu = Op<T>::u(vRes[j]);
    PWLU const R = _F[j].restrict(Op<T>::l(vVar[0]), Op<T>::u(vVar[0]));
    double pl = INFINITY, pu = -INFINITY, xa = R.xL(), ya = R.yL();
    for (size_t k = 0; k < R.size(); ++k)
    {
      double const xb = xa + R.dx()[k], yb = ya + R.dx()[k] * R.dy()[k];
      if (std::max(ya, yb) >= yl && std::min(ya, yb) <= yu)
      {
        double a = xa, b = xb;
        if (R.dy()[k] != 0.)
        {
          // Points of the segment where the function is yl and yu
          double const tl = xa + (yl - ya) / R.dy()[k], tu = xa + (yu - ya) / R.dy()[k];
          a = std::max(a, std::min(tl, tu));
          b = std::min(b, std::max(tl, tu));
        }
        pl = std::min(pl, a);
        pu = std::max(pu, b);
      }
      xa = xb;
      ya = yb;
    }
    if (pl > pu) return false;
    if (!Op<T>::inter(vVar[0], vVar[0], T(pl, pu))) return false;
    // Range of the result on the contracted operand
    PWLU const Rc = _F[j].restrict(Op<T>::l(vVar[0]), Op<T>::u(vVar[0]));
    if (!Op<T>::inter(vRes[j], vRes[j], T(Rc.l(), Rc.u()))) return false;
  }
  return true;
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, McCormick<T>* vRes, size_t const nVar,
                McCormick<T> const* vVar, unsigned const* mVar) const
{
  McCormick<T> const& MX = vVar[0];
  unsigned const nsub    = MX.nsub();
  std::vector<double> cvsub(nsub), ccsub(nsub);
  auto mid = [](double const a, double const b, double const c) {
    return std::max(std::min(a, b), std::min(std::max(a, b), c));
  };
  for (size_t j = 0; j < nRes; ++j)
  {
    PWLU const R = _F[j].restrict(MX.l(), MX.u());

    // Minimizer and maximizer among the vertices
    double xmin = R.xL(), xmax = R.xL(), fmin = R.yL(), fmax = R.yL(), xk = R.xL(), yk = R.yL();
    for (size_t k = 0; k < R.size(); ++k)
    {
      xk += R.dx()[k];
      yk += R.dx()[k] * R.dy()[k];
      if (yk < fmin) { fmin = yk; xmin = xk; }
      if (yk > fmax) { fmax = yk; xmax = xk; }
    }

    // Convex envelope at mid(cv,cc,xmin), concave envelope at mid(cv,cc,xmax)
    double const zcv      = mid(MX.cv(), MX.cc(), xmin);
    auto const [fcv, scv] = _hull_eval(R.hull(true), zcv, xmin);
    double const zcc      = mid(MX.cv(), MX.cc(), xmax);
    auto const [fcc, scc] = _hull_eval(R.hull(false), zcc, xmax);
    for (unsigned i = 0; i < nsub; ++i)
    {
      cvsub[i] = (zcv == xmin ? 0. : zcv == MX.cv() ? scv * MX.cvsub(i) : scv * MX.ccsub(i));
      ccsub[i] = (zcc == xmax ? 0. : zcc == MX.cv() ? scc * MX.cvsub(i) : scc * MX.ccsub(i));
    }
    vRes[j] = McCormick<T>(T(fmin, fmax), fcv, fcc);
    vRes[j].sub(nsub, cvsub.data(), ccsub.data()).cut();
  }
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T> const* vVar,
                unsigned const* mVar) const
{
  PolImg<T>* img = vVar[0].image();
  double const l = Op<T>::l(vVar[0].range()), u = Op<T>::u(vVar[0].range());
  for (size_t j = 0; j < nRes; ++j)
  {
    PWLU const R = _F[j].restrict(l, u);
    vRes[j].set(img, *this->varout[j], T(R.l(), R.u()));
  }
}

template <typename T>
inline void
FFPWLU<T>::_subdiv_cuts(PolImg<T>* img, FFOp* pOp, PolVar<T> const& X, double const l,
                        double const u, std::vector<PWLU> const& R, std::vector<size_t> const& J,
                        PolVar<T> const* vRes, int const enc) const
{
  // Subdivision of the operand at the union of the breakpoints of the functions in J
  std::vector<double> Xk(1, l);
  for (size_t j : J)
  {
    auto const bk = R[j].bkpts();
    Xk.insert(Xk.end(), bk.begin(), bk.end());
  }
  Xk.push_back(u);
  std::sort(Xk.begin(), Xk.end());
  Xk.erase(std::unique(Xk.begin(), Xk.end()), Xk.end());
  int const type = (enc == Options::SOS2 ? 2 : enc == Options::BIN ? 1 : 0);
  std::vector<double> const XKNOT = img->add_subdivision(pOp, Xk.size(), X, Xk.data(), type);

  // Each function as the interpolant of its values at the knots of the subdivision
  std::vector<double> Yk(XKNOT.size());
  for (size_t j : J)
  {
    for (size_t k = 0; k < XKNOT.size(); ++k) Yk[k] = R[j].val(XKNOT[k], true);
    img->add_subdivision_cut(pOp, X, vRes[j], Yk.data(), PolCut<T>::EQ, type);
  }
}

template <typename T>
inline bool
FFPWLU<T>::reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar, PolVar<T>* vVar) const
{
  PolImg<T>* img = vVar[0].image();
  FFOp* pOp      = vRes[0].var().opdef().first;
  double const l = Op<T>::l(vVar[0].range()), u = Op<T>::u(vVar[0].range());

  // Operand fixed: each result is constant
  if (!(l < u))
  {
    for (size_t j = 0; j < nRes; ++j)
      img->add_cut(pOp, PolCut<T>::EQ, _F[j].val(l, true), vRes[j], 1.);
    return true;
  }

  // Functions affine on [l,u]: equality; the others are encoded below
  std::vector<PWLU> R(nRes);
  std::vector<size_t> J;
  for (size_t j = 0; j < nRes; ++j)
  {
    R[j] = _F[j].restrict(l, u);
    if (!R[j].bkpts().empty())
    {
      J.push_back(j);
      continue;
    }
    double const s = R[j].dy().front();
    img->add_cut(pOp, PolCut<T>::EQ, R[j].yL() - s * l, vRes[j], 1., vVar[0], -s);
  }
  if (J.empty()) return true;

  int const enc = _encoding(img);
  if (enc == Options::HULL)
  {
    // Convex hull of the graph of each function: the segments of its convex and concave envelopes
    for (size_t j : J)
      for (bool under : {true, false})
      {
        PWLU const H = R[j].hull(under);
        double xk = H.xL(), yk = H.yL();
        for (size_t k = 0; k < H.size(); ++k)
        {
          double const s = H.dy()[k];
          img->add_cut(pOp, under ? PolCut<T>::GE : PolCut<T>::LE, yk - s * xk, vRes[j], 1.,
                       vVar[0], -s);
          xk += H.dx()[k];
          yk += H.dx()[k] * s;
        }
      }
    return true;
  }

  // Subdivision of the operand: one for all functions, or one for each
  if (!options.SEPARATE)
    _subdiv_cuts(img, pOp, vVar[0], l, u, R, J, vRes, enc);
  else
    for (size_t j : J) _subdiv_cuts(img, pOp, vVar[0], l, u, R, {j}, vRes, enc);
  return true;
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
                unsigned const* mVar) const
{
  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFPWLU<T>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
                unsigned const* mVar) const
{
  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (size_t j = 0; j < nRes; ++j)
      {
        std::ostringstream os;
        os << name() << "[" << j << "]";
        vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
      }
      break;
    case FFExpr::Options::GAMS:
    default:
      throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
  }
}

template <typename T>
inline bool
FFPWLU<T>::lt(FFOp const* op) const
{
  auto const* other = mc::type_cast<FFPWLU<T> const>(op);
  if (options.ENCODING != other->options.ENCODING)
    return options.ENCODING < other->options.ENCODING;
  if (options.SEPARATE != other->options.SEPARATE) return options.SEPARATE < other->options.SEPARATE;
  if (_F.size() != other->_F.size()) return _F.size() < other->_F.size();
  for (size_t j = 0; j < _F.size(); ++j)
  {
    PWLU const &f = _F[j], &g = other->_F[j];
    if (f.xL() != g.xL()) return f.xL() < g.xL();
    if (f.yL() != g.yL()) return f.yL() < g.yL();
    if (f.dx() != g.dx()) return f.dx() < g.dx();
    if (f.dy() != g.dy()) return f.dy() < g.dy();
  }
  return false;
}

}  // namespace mc

#endif
