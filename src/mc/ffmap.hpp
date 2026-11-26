// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__FFMAP_HPP
#define MC__FFMAP_HPP

#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

#include "ffextern.hpp"

// #define MC__FFMAP_DEBUG
#define MC__FFMAP_CHECK

namespace mc
{

// For the accessors MLP() and pMLP() of the operations on an mc::MLP
template <typename T>
class MLP;

//! @brief C++ class defining networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMAP is a C++ class for defining a network -- a map from its
//! inputs to its outputs held by an object of class MAP -- as an external
//! DAG operation in MC++: mc::FFMLP<T> is mc::FFMAP<T, mc::MLP<T>>, for a
//! multilayer perceptron, and mc::FFKAN<T> is mc::FFMAP<T, mc::KAN<T>>, for
//! a Kolmogorov-Arnold network. As any mc::FFEXTERN, the operation has the
//! network as its home object in the DAG, and is evaluated by the network
//! in the arithmetics of MC++. The class MAP provides, as mc::MLP and
//! mc::KAN do:
//! - nin(), nout(), nwei() and uid(), and the static tag(), the name of the
//!   operation;
//! - eval(U const* x, U* y) in the arithmetics U of mc::FFEXTERN, and
//!   reval(T* x, T* y), the backward propagation of interval bounds;
//! - grad(x, J, y), tang(x, dx, dy, y) and adj(x, ay, ax, y) in double
//!   precision, for the operations mc::FFGradMAP, mc::FFTanMAP and
//!   mc::FFAdjMAP that the differentiation of the DAG inserts -- the last
//!   two where its options TANGENT and ADJOINT are set;
//! - copy construction, for the policy COPY.
////////////////////////////////////////////////////////////////////////
template <typename T, typename MAP>
class FFMAP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MAP>
{
 protected:
  using FFEXTERN<T, MAP>::_ptrObj;
  using FFEXTERN<T, MAP>::_ownObj;

  // set the object and related operation in DAG
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MAP* pNet, int policy)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin());
#endif
    // home of the network in the DAG of the operands: made on its first
    // insertion there, with a copy of the network (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pNet, policy);
    FFVar** ppRes =
        this->insert_external_operation(*this, _ptrObj->nout(), nVar, pVar);
    return ppRes;
  }

  // set the operation in DAG on the network of another operation, whose home
  // it shares
  FFVar**
  _set(size_t const nVar, FFVar const* pVar,
       FFEXTERN<T, MAP> const& parent)
  {
    this->owndata = false;
    this->_join(parent);
    return this->insert_external_operation(*this, _ptrObj->nout(), nVar, pVar);
  }

 public:
  //! @brief Enumeration type for copy policy of the network object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< The DAG refers to the network object, which it does not own
    COPY    = 1,  //!< The DAG holds a copy of the network object: one per DAG
    TRANSFER =
        -1  //!< The DAG refers to the network object, and deletes it in the end
  };

  //! @brief Default constructor
  FFMAP() : FFEXTERN<T, MAP>() {}

  // Destructor
  virtual ~FFMAP() {}

  // Copy constructor
  FFMAP(FFMAP<T, MAP> const& Other) : FFEXTERN<T, MAP>(Other) {}

  // Define operation
  // FFVar** operator()
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin());
#endif
    // return _set( vVar.size(), vVar.data(), pNet, policy );
    FFVar** ppDer = _set(vVar.size(), vVar.data(), pNet, policy);
    std::vector<FFVar> vDer(pNet->nout());
    for (size_t i = 0; i < vDer.size(); ++i) vDer[i] = *ppDer[i];
    return vDer;  // std::move( vDer );
  }

  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vVar, MAP* pNet,
             int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin() && idep < pNet->nout());
#endif
    return *(_set(vVar.size(), vVar.data(), pNet, policy)[idep]);
  }

  FFVar**
  operator()(size_t const nVar, FFVar const* pVar, MAP* pNet,
             int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin());
#endif
    return _set(nVar, pVar, pNet, policy);
  }

  FFVar&
  operator()(size_t const idep, size_t const nVar, FFVar const* pVar,
             MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin() && idep < pNet->nout());
#endif
    return *(_set(nVar, pVar, pNet, policy)[idep]);
  }

  //! @brief Network (map) object of the DAG: after an insertion, the copy that the
  //! DAG holds (policy COPY) or the object given, shared by every operation
  //! on that network in that DAG -- see mc::FFEXTERN::object. Throws before
  //! any insertion
  MAP&
  map() const
  {
    return this->object();
  }

  //! @brief The same under its name for an mc::MLP (mc::FFMLP): MLP()
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  N&
  MLP() const
  {
    return this->object();
  }

  //! @brief Pointer to the MLP object of the DAG, null before any insertion
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  [[deprecated("use map(), which returns a reference")]] N*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FADType<FFVar>) ))
      return eval(nRes, static_cast<FADType<FFVar>*>(vRes), nVar,
                  static_cast<FADType<FFVar> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FADType<double>) ))
      return eval(nRes, static_cast<FADType<double>*>(vRes), nVar,
                  static_cast<FADType<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(BADType<double>) ))
      return evalbad(nRes, static_cast<BADType<double>*>(vRes), nVar,
                     static_cast<BADType<double> const*>(vVar));
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(PolVar<T>) ))
      return eval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                  static_cast<PolVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFMAP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  //! @brief Evaluation in the reverse mode of AD: each output is a single
  //! node whose partial derivatives are the Jacobian of the network, computed
  //! in double precision by passes through the layers (see the member grad of the network), in
  //! place of a tape of every neuron
  void evalbad(size_t const nRes, BADType<double>* vRes, size_t const nVar,
               BADType<double> const* vVar) const
  {
#ifdef MC__FFMAP_CHECK
    assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif
    std::vector<double> x(nVar), y(nRes), J(nRes * nVar);
    for (size_t i = 0; i < nVar; ++i) x[i] = vVar[i].val();
    _ptrObj->grad(x.data(), J.data(), y.data());
    for (size_t k = 0; k < nRes; ++k)
      vRes[k] = bad_node(y[k], nVar, vVar, J.data() + k * nVar);
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
            FADType<FFVar> const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
            PolVar<T> const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    if (mc::same_type( idU, typeid(T) ))
      return reval(nRes, static_cast<T*>(vRes), nVar, static_cast<T*>(vVar));
    else if (mc::same_type( idU, typeid(PolVar<T>) ))
      return reval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                   static_cast<PolVar<T>*>(vVar));

    throw std::runtime_error(
        "FFMAP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  bool reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const;

  bool reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
             PolVar<T>* vVar) const;

  // Derivatives
  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
             FFVar const* vVar, FFVar** vDer) const;

  // Whether nSweep passes through the network, each costing npass
  // evaluations, are estimated to be cheaper in a DAG than the Jacobian
  // operation and the products that the differentiation forms with its
  // entries in each sweep. In units of one weight in an evaluation of the
  // network, with W weights, nin inputs and nout outputs, timings of the
  // derivative DAG in double arithmetic give about
  //  - a pass:       npass W + 100, the latter for the operation itself
  //  - the Jacobian: W (1/2 + min(nin,nout)/4)
  //  - the products: 35 nin nout, when the Jacobian is evaluated and again
  //                  in each sweep
  // The passes are always the cheaper where the Jacobian has many entries
  // for the size of the network, and never, even for a single sweep, where
  // it has few -- a network with a single output above all. Over networks
  // of 2 to 200 inputs and 1 to 50 outputs and 1 to 32 sweeps, the choice
  // came within 4% of the better of the two in total; where it was the
  // wrong one it lost tens of microseconds, and 0.2 to 0.7 ms (16 to 34%)
  // at worst, with 200 inputs for 5 outputs and 12 to 16 sweeps
  bool
  passes_pay(unsigned const nSweep, double const npass) const
  {
    double const W = _ptrObj->nwei(), E = _ptrObj->nin() * _ptrObj->nout(),
                 R = std::min(_ptrObj->nin(), _ptrObj->nout());
    return nSweep * (npass * W + 100. - 35. * E) <=
           W * (0.5 + 0.25 * R) + 35. * E;
  }

  // Tangents: with MAP::Options::TANGENT, the forward symbolic
  // differentiation inserts the directional derivative operation FFTanMAP
  // along the tangents of the operands, in place of the Jacobian FFGradMAP
  // -- where that is estimated to be the cheaper for nSweep sweeps, a
  // tangent pass costing about 1.1 evaluations of the network
  bool
  has_tangent(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.TANGENT && passes_pay(nSweep, 1.1);
  }
  void tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vTan,
               FFVar* vResTan) const;

  // Adjoints: with MAP::Options::ADJOINT, the backward symbolic
  // differentiation inserts the adjoint operation FFAdjMAP for the adjoints
  // of the results, in place of the Jacobian FFGradMAP -- where that is
  // estimated to be the cheaper for nSweep sweeps, an adjoint pass costing
  // about 1.3 evaluations of the network
  bool
  has_adjoint(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.ADJOINT && passes_pay(nSweep, 1.3);
  }
  void adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vResAdj,
               FFVar* vVarAdj) const;

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return MAP::tag() + "[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

//! @brief C++ class defining gradient of neural networks as external DAG
//! operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradMAP is a C++ class for defining gradient of neural
//! networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template <typename T, typename MAP>
class FFGradMAP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MAP>
{
  friend class FFMAP<T, MAP>;

 protected:
  using FFEXTERN<T, MAP>::_ptrObj;
  using FFEXTERN<T, MAP>::_ownObj;

  // set the object and related operation in DAG
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MAP* pNet, int policy)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin());
#endif
    // home of the network in the DAG of the operands: made on its first
    // insertion there, with a copy of the network (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pNet, policy);
    FFVar** ppRes =
        this->insert_external_operation(*this, nVar * _ptrObj->nout(), nVar, pVar);
    return ppRes;
  }

  // set the operation in DAG on the network of another operation, whose home
  // it shares
  FFVar**
  _set(size_t const nVar, FFVar const* pVar,
       FFEXTERN<T, MAP> const& parent)
  {
    this->owndata = false;
    this->_join(parent);
    return this->insert_external_operation(*this, nVar * _ptrObj->nout(), nVar, pVar);
  }

 public:
  //! @brief Enumeration type for copy policy of the network object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< The DAG refers to the network object, which it does not own
    COPY    = 1,  //!< The DAG holds a copy of the network object: one per DAG
    TRANSFER =
        -1  //!< The DAG refers to the network object, and deletes it in the end
  };

  //! @brief Default constructor
  FFGradMAP() : FFEXTERN<T, MAP>() {}

  // Destructor
  virtual ~FFGradMAP() {}

  // Copy constructor
  FFGradMAP(FFGradMAP<T, MAP> const& Other) : FFEXTERN<T, MAP>(Other) {}

  // Define operation
  FFVar**
  operator()(std::vector<FFVar> const& vVar, MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin());
#endif
    return _set(vVar.size(), vVar.data(), pNet, policy);
  }

  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vVar, MAP* pNet,
             int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin() && idep < pNet->nout() * pNet->nin());
#endif
    return *(_set(vVar.size(), vVar.data(), pNet, policy)[idep]);
  }

  FFVar**
  operator()(size_t const nVar, FFVar const* pVar, MAP* pNet,
             int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin());
#endif
    return _set(nVar, pVar, pNet, policy);
  }

  FFVar&
  operator()(size_t const idep, size_t const nVar, FFVar const* pVar,
             MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin() && idep < pNet->nout() * pNet->nin());
#endif
    return *(_set(nVar, pVar, pNet, policy)[idep]);
  }

  //! @brief Network (map) object of the DAG: after an insertion, the copy that the
  //! DAG holds (policy COPY) or the object given, shared by every operation
  //! on that network in that DAG -- see mc::FFEXTERN::object. Throws before
  //! any insertion
  MAP&
  map() const
  {
    return this->object();
  }

  //! @brief The same under its name for an mc::MLP (mc::FFMLP): MLP()
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  N&
  MLP() const
  {
    return this->object();
  }

  //! @brief Pointer to the MLP object of the DAG, null before any insertion
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  [[deprecated("use map(), which returns a reference")]] N*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFGradMAP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFGradMAP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return std::string("Grad") + MAP::tag() + "[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

//! @brief C++ class defining directional derivatives of neural networks as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFTanMAP is a C++ class for defining the directional derivative
//! of a neural network as an external DAG operation in MC++: its
//! operands are the inputs of the network followed by a direction of
//! the same size, and its results are the derivatives of the outputs
//! of the network along that direction, evaluated in a single forward
//! pass. It is inserted by the forward symbolic differentiation of
//! mc::FFMAP when the option TANGENT of the network is set (FFOp::tangent),
//! and can also be defined directly. It has no derivative of its own.
////////////////////////////////////////////////////////////////////////
template <typename T, typename MAP>
class FFTanMAP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MAP>
{
  friend class FFMAP<T, MAP>;

 protected:
  using FFEXTERN<T, MAP>::_ptrObj;
  using FFEXTERN<T, MAP>::_ownObj;

  // set the object and related operation in DAG; the operands are the inputs
  // followed by the direction
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MAP* pNet, int policy)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == 2 * pNet->nin());
#endif
    // home of the network in the DAG of the operands: made on its first
    // insertion there, with a copy of the network (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pNet, policy);
    FFVar** ppRes =
        this->insert_external_operation(*this, _ptrObj->nout(), nVar, pVar);
    return ppRes;
  }

  // set the operation in DAG on the network of another operation, whose home
  // it shares
  FFVar**
  _set(size_t const nVar, FFVar const* pVar,
       FFEXTERN<T, MAP> const& parent)
  {
    this->owndata = false;
    this->_join(parent);
    return this->insert_external_operation(*this, _ptrObj->nout(), nVar, pVar);
  }

 public:
  //! @brief Enumeration type for copy policy of the network object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< The DAG refers to the network object, which it does not own
    COPY    = 1,  //!< The DAG holds a copy of the network object: one per DAG
    TRANSFER =
        -1  //!< The DAG refers to the network object, and deletes it in the end
  };

  //! @brief Default constructor
  FFTanMAP() : FFEXTERN<T, MAP>() {}

  // Destructor
  virtual ~FFTanMAP() {}

  // Copy constructor
  FFTanMAP(FFTanMAP<T, MAP> const& Other) : FFEXTERN<T, MAP>(Other) {}

  // Define operation: derivatives of the outputs at vVar along vDir
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vDir,
             MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin() && vDir.size() == pNet->nin());
#endif
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vDir.begin(), vDir.end());
    FFVar** ppTan = _set(vOpd.size(), vOpd.data(), pNet, policy);
    std::vector<FFVar> vTan(pNet->nout());
    for (size_t k = 0; k < vTan.size(); ++k) vTan[k] = *ppTan[k];
    return vTan;
  }

  //! @brief Network (map) object of the DAG: after an insertion, the copy that the
  //! DAG holds (policy COPY) or the object given, shared by every operation
  //! on that network in that DAG -- see mc::FFEXTERN::object. Throws before
  //! any insertion
  MAP&
  map() const
  {
    return this->object();
  }

  //! @brief The same under its name for an mc::MLP (mc::FFMLP): MLP()
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  N&
  MLP() const
  {
    return this->object();
  }

  //! @brief Pointer to the MLP object of the DAG, null before any insertion
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  [[deprecated("use map(), which returns a reference")]] N*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFTanMAP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFTanMAP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return std::string("Tan") + MAP::tag() + "[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T, typename MAP>
template <typename U>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                  unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFTanMAP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  _ptrObj->tang(vVar, vVar + _ptrObj->nin(), vRes);
}

template <typename T, typename MAP>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                  FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                  FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                  FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  FFVar** ppRes = this->_reinsert(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T, typename MAP>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                  SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T, typename MAP>
inline void
FFTanMAP<T, MAP>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                  FFExpr const* vVar, unsigned const* mVar) const
{
  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
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

//! @brief C++ class defining adjoints of neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFAdjMAP is a C++ class for defining the adjoint of a neural
//! network as an external DAG operation in MC++: its operands are the
//! inputs of the network followed by one weight per output, and its
//! results are the derivatives of the weighted sum of the outputs with
//! respect to the inputs -- the product of the transposed Jacobian with
//! the weights -- evaluated in a single backward pass. It is inserted
//! by the backward symbolic differentiation of mc::FFMAP when
//! the option ADJOINT of the network is set (FFOp::adjoint), and can also be
//! defined directly. It has no derivative of its own.
////////////////////////////////////////////////////////////////////////
template <typename T, typename MAP>
class FFAdjMAP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MAP>
{
  friend class FFMAP<T, MAP>;

 protected:
  using FFEXTERN<T, MAP>::_ptrObj;
  using FFEXTERN<T, MAP>::_ownObj;

  // set the object and related operation in DAG; the operands are the inputs
  // followed by the weights of the outputs
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MAP* pNet, int policy)
  {
#ifdef MC__FFMAP_CHECK
    assert(nVar == pNet->nin() + pNet->nout());
#endif
    // home of the network in the DAG of the operands: made on its first
    // insertion there, with a copy of the network (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pNet, policy);
    FFVar** ppRes =
        this->insert_external_operation(*this, _ptrObj->nin(), nVar, pVar);
    return ppRes;
  }

  // set the operation in DAG on the network of another operation, whose home
  // it shares
  FFVar**
  _set(size_t const nVar, FFVar const* pVar,
       FFEXTERN<T, MAP> const& parent)
  {
    this->owndata = false;
    this->_join(parent);
    return this->insert_external_operation(*this, _ptrObj->nin(), nVar, pVar);
  }

 public:
  //! @brief Enumeration type for copy policy of the network object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< The DAG refers to the network object, which it does not own
    COPY    = 1,  //!< The DAG holds a copy of the network object: one per DAG
    TRANSFER =
        -1  //!< The DAG refers to the network object, and deletes it in the end
  };

  //! @brief Default constructor
  FFAdjMAP() : FFEXTERN<T, MAP>() {}

  // Destructor
  virtual ~FFAdjMAP() {}

  // Copy constructor
  FFAdjMAP(FFAdjMAP<T, MAP> const& Other) : FFEXTERN<T, MAP>(Other) {}

  // Define operation: derivatives in vVar of the sum of the outputs weighted
  // by vWei
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vWei,
             MAP* pNet, int policy = COPY)
  {
#ifdef MC__FFMAP_CHECK
    assert(vVar.size() == pNet->nin() && vWei.size() == pNet->nout());
#endif
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vWei.begin(), vWei.end());
    FFVar** ppAdj = _set(vOpd.size(), vOpd.data(), pNet, policy);
    std::vector<FFVar> vAdj(pNet->nin());
    for (size_t i = 0; i < vAdj.size(); ++i) vAdj[i] = *ppAdj[i];
    return vAdj;
  }

  //! @brief Network (map) object of the DAG: after an insertion, the copy that the
  //! DAG holds (policy COPY) or the object given, shared by every operation
  //! on that network in that DAG -- see mc::FFEXTERN::object. Throws before
  //! any insertion
  MAP&
  map() const
  {
    return this->object();
  }

  //! @brief The same under its name for an mc::MLP (mc::FFMLP): MLP()
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  N&
  MLP() const
  {
    return this->object();
  }

  //! @brief Pointer to the MLP object of the DAG, null before any insertion
  template <typename N = MAP, typename = typename std::enable_if<std::is_same<N, mc::MLP<T>>::value>::type>
  [[deprecated("use map(), which returns a reference")]] N*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFAdjMAP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFAdjMAP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return std::string("Adj") + MAP::tag() + "[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T, typename MAP>
template <typename U>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                  unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFAdjMAP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  _ptrObj->adj(vVar, vVar + _ptrObj->nin(), vRes);
}

template <typename T, typename MAP>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                  FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                  FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                  FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  FFVar** ppRes = this->_reinsert(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T, typename MAP>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                  SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T, typename MAP>
inline void
FFAdjMAP<T, MAP>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                  FFExpr const* vVar, unsigned const* mVar) const
{
  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
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

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                  FFVar const* vVar, FFVar const* const* vTan,
                  FFVar* vResTan) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::tangent\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  // Operands of the directional derivative: the inputs, then their tangents
  // -- zero where the tangent of an input is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(2 * nVar);
  for (unsigned i = 0; i < nVar; ++i)
    vOpd.push_back(vTan[i] ? *vTan[i] : FFVar(0.));

  FFTanMAP<T, MAP> ResTan;
  FFVar const* const* vResDer = ResTan._set(vOpd.size(), vOpd.data(), *this);
  for (unsigned k = 0; k < nRes; ++k) vResTan[k] = *vResDer[k];
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                  FFVar const* vVar, FFVar const* const* vResAdj,
                  FFVar* vVarAdj) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::adjoint\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  // Operands of the adjoint operation: the inputs, then the adjoints of the
  // outputs -- zero where the adjoint of an output is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(nVar + nRes);
  for (unsigned k = 0; k < nRes; ++k)
    vOpd.push_back(vResAdj[k] ? *vResAdj[k] : FFVar(0.));

  FFAdjMAP<T, MAP> ResAdj;
  FFVar const* const* vResDer = ResAdj._set(vOpd.size(), vOpd.data(), *this);
  for (unsigned i = 0; i < nVar; ++i) vVarAdj[i] = *vResDer[i];
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
               FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: FFDep\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
               FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: FFInv\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                   FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: FFDep\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                   FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: FFInv\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
               FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: FFVar\n";
  std::cerr << "FFMAP operation address: " << this << std::endl;
  std::cerr << "Network address in DAG: " << _ptrObj << std::endl;
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  FFVar** ppRes = this->_reinsert(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T, typename MAP>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                   FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: FFVar\n";
  std::cerr << "FFGradMAP operation address: " << this << std::endl;
  std::cerr << "Network address in DAG: " << _ptrObj << std::endl;
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  FFVar** ppRes = this->_reinsert(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T, typename MAP>
template <typename U>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
               unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  _ptrObj->eval(vVar, vRes);
}

template <typename T, typename MAP>
template <typename U>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                   unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  _ptrObj->grad(vVar, vRes);
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
               FADType<FFVar> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: FADType<FFVar>\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  std::vector<FFVar> vVarVal(nVar);
  for (unsigned i = 0; i < nVar; ++i) vVarVal[i] = vVar[i].val();
  FFVar const* const* vResVal =
      this->_reinsert(*this, nRes, nVar, vVarVal.data());

  FFGradMAP<T, MAP> ResDer;
  FFVar const* const* vResDer =
      ResDer._set(nVar, vVarVal.data(), *this);

  for (unsigned k = 0; k < nRes; ++k)
  {
    vRes[k] = *vResVal[k];
    for (unsigned i = 0; i < nVar; ++i) vRes[k].setDepend(vVar[i]);
    for (unsigned j = 0; j < vRes[k].size(); ++j)
    {
      vRes[k][j] = 0.;
      for (unsigned i = 0; i < nVar; ++i)
      {
        if (vVar[i][j].cst() && vVar[i][j].num().val() == 0.) continue;
        // vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
        vRes[k][j] += *vResDer[k * nVar + i] * vVar[i][j];
      }
    }
  }
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                FFVar const* vVar, FFVar** vDer) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::deriv\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  FFGradMAP<T, MAP> ResDer;
  FFVar const* const* vResDer = ResDer._set(nVar, vVar, *this);
  for (unsigned k = 0; k < nRes; ++k)
    for (unsigned i = 0; i < nVar; ++i)
      // vDer[k][i] = *vResDer[k+nRes*i];
      vDer[k][i] = *vResDer[k * nVar + i];
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
               SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T, typename MAP>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                   SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
               FFExpr const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: FFExpr\n";
#endif

  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
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

template <typename T, typename MAP>
inline void
FFGradMAP<T, MAP>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                   FFExpr const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFGradMAP::eval: FFExpr\n";
#endif

  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
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

template <typename T, typename MAP>
inline void
FFMAP<T, MAP>::eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
               PolVar<T> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::eval: PolVar\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag    = vVar[0].var().dag();
#ifdef MC__FFMAP_CHECK
  assert(img && dag);
#endif
  FFVar** ppRes = dag->curOp()->varout.data();
#ifdef MC__FFMAP_CHECK
  assert(nRes == dag->curOp()->varout.size());
#endif

  this->_resize_relax(img);
  this->_propagate_relax(img, ppRes, vRes, vVar);
}

template <typename T, typename MAP>
inline bool
FFMAP<T, MAP>::reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                PolVar<T>* vVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::reval: PolVar\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop      = vVar[0].var().opdef().first;
#ifdef MC__FFMAP_CHECK
  assert(img && pop);
#endif

  this->_backpropagate_relax(img, pop, vRes, vVar);
  return true;
}

template <typename T, typename MAP>
inline bool
FFMAP<T, MAP>::reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const
{
#ifdef MC__FFMAP_TRACE
  std::cout << "FFMAP::reval: T\n";
#endif
#ifdef MC__FFMAP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  return _ptrObj->reval(vVar, vRes);
}

}  // end namespace mc

#endif
