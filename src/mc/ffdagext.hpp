// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__FFDAGEXT_HPP
#define MC__FFDAGEXT_HPP

// #define MC__FFDAGEXT_DEBUG
#define MC__FFDAGEXT_CHECK

#include <cstddef>  // std::ptrdiff_t
#include <limits>   // std::numeric_limits
#include <new>      // std::bad_alloc

#include "interval.hpp"
#include "ffextern.hpp"

namespace mc
{

//! @brief C++ class for evaluation and relaxation of an operation defined by an
//! expression tree
////////////////////////////////////////////////////////////////////////
//! mc::DAGEXT is a C++ class for evaluation and relaxation of an
//! operation defined by an expression tree using MC++
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class DAGEXT
////////////////////////////////////////////////////////////////////////
{
 private:
  //! @brief Expression tree
  FFGraph* _dag;
  //! @brief Variables
  std::vector<FFVar> _varin;
  //! @brief Dependents
  std::vector<FFVar> _varout;
  //! @brief Codelist
  FFSubgraph _codelist;
  //! @brief Stamp of the expression tree: see uid()
  size_t _uid;

  //! @brief Intermediate storage for DAG evaluation
  std::vector<double>                    _wkD;
  std::vector<FADType<double>>           _wkFD;
  std::vector<BADType<double>>           _wkBD;
  std::vector<T>                         _wkI;
  std::vector<T>                         _wkCPI;
  std::vector<McCormick<T>>              _wkMC;
  std::vector<Specbnd<T>>                _wkSB;
  std::vector<SCVar<T>>                  _wkSC;
  std::vector<SupVar<PWCU>>              _wkPWCS;
  std::vector<SupVar<PWLU>>              _wkPWLS;
  std::vector<PolVar<T>>                 _wkPOL;
  std::vector<OCVar<double>>             _wkOCD;
  std::vector<OCVar<FADType<double>>>    _wkOCFD;

  //! @brief Set expression tree
  void
  _set(FFGraph* dag, std::vector<FFVar> const& varin,
       std::vector<FFVar> const& varout)
  {
    delete _dag;
    _dag = new FFGraph;
    _uid = ffextern_uid();
#ifdef MC__DAGEXT_DEBUG
    std::cout << "DAGEXT:: Original DAG: " << dag << std::endl;
    std::cout << "DAGEXT:: Copied DAG:   " << _dag << std::endl;
#endif
    _dag->insert(dag, varin, _varin);
    _dag->insert(dag, varout, _varout);
#ifdef MC__DAGEXT_DEBUG
    _codelist                  = _dag->subgraph(_varout);
    std::vector<FFExpr> strout = FFExpr::subgraph(_dag, _codelist);
    for (size_t i = 0; i < strout.size(); ++i)
      std::cout << "F" << i << ": " << strout[i] << std::endl;
#else
    _codelist.clear();
#endif
  }

  //! @brief Set expression tree
  void
  _set(FFGraph* dag, std::set<unsigned> const& ndxin,
       std::vector<FFVar> const& varin, std::vector<FFVar> const& varout)
  {
    delete _dag;
    _dag = new FFGraph;
    _uid = ffextern_uid();
#ifdef MC__DAGEXT_DEBUG
    std::cout << "DAGEXT:: Original DAG: " << dag << std::endl;
    std::cout << "DAGEXT:: Copied DAG:   " << _dag << std::endl;
#endif
    // Set inputs
    _dag->insert(dag, ndxin, varin, _varin);
    size_t i = 0;
    for (auto const& j : ndxin) _varin[i++] = _varin[j];
    _varin.resize(i);

    // Set outputs
    _dag->insert(dag, varout, _varout);
#ifdef MC__DAGEXT_DEBUG
    _codelist                  = _dag->subgraph(_varout);
    std::vector<FFExpr> strout = FFExpr::subgraph(_dag, _codelist);
    for (size_t i = 0; i < strout.size(); ++i)
      std::cout << "F" << i << ": " << strout[i] << std::endl;
#else
    _codelist.clear();
#endif
  }

  //! @brief Set expression tree
  void
  _set(FFGraph* dag, std::vector<FFVar> const& varin,
       std::set<unsigned> const& ndxout, std::vector<FFVar> const& varout)
  {
    delete _dag;
    _dag = new FFGraph;
    _uid = ffextern_uid();
#ifdef MC__DAGEXT_DEBUG
    std::cout << "DAGEXT:: Original DAG: " << dag << std::endl;
    std::cout << "DAGEXT:: Copied DAG:   " << _dag << std::endl;
#endif
    // Set inputs
    _dag->insert(dag, varin, _varin);

    // Set outputs
    _dag->insert(dag, ndxout, varout, _varout);
    size_t i = 0;
    for (auto const& j : ndxout) _varout[i++] = _varout[j];
    _varout.resize(i);
#ifdef MC__DAGEXT_DEBUG
    _codelist                  = _dag->subgraph(_varout);
    std::vector<FFExpr> strout = FFExpr::subgraph(_dag, _codelist);
    for (size_t i = 0; i < strout.size(); ++i)
      std::cout << "F" << i << ": " << strout[i] << std::endl;
#else
    _codelist.clear();
#endif
  }

  //! @brief Set expression tree
  void
  _set(FFGraph* dag, std::set<unsigned> const& ndxin,
       std::vector<FFVar> const& varin, std::set<unsigned> const& ndxout,
       std::vector<FFVar> const& varout)
  {
    delete _dag;
    _dag = new FFGraph;
    _uid = ffextern_uid();
#ifdef MC__DAGEXT_DEBUG
    std::cout << "DAGEXT:: Original DAG: " << dag << std::endl;
    std::cout << "DAGEXT:: Copied DAG:   " << _dag << std::endl;
#endif
    // Set inputs
    _dag->insert(dag, ndxin, varin, _varin);
    size_t i = 0;
    for (auto const& j : ndxin) _varin[i++] = _varin[j];
    _varin.resize(i);

    // Set outputs
    _dag->insert(dag, ndxout, varout, _varout);
    i = 0;
    for (auto const& j : ndxout) _varout[i++] = _varout[j];
    _varout.resize(i);
#ifdef MC__DAGEXT_DEBUG
    _codelist                  = _dag->subgraph(_varout);
    std::vector<FFExpr> strout = FFExpr::subgraph(_dag, _codelist);
    for (size_t i = 0; i < strout.size(); ++i)
      std::cout << "F" << i << ": " << strout[i] << std::endl;
#else
    _codelist.clear();
#endif
  }

  //! @brief Evaluate expression tree
  template <typename U>
  void
  _eval(U const* valin, U* valout, std::vector<U>& wk)
  //    const
  {
#ifdef MC__DAGEXT_DEBUG
    std::cout << std::endl;
    auto exout = mc::FFExpr::subgraph(_dag, _codelist);
    for (unsigned i = 0; i < exout.size(); ++i)
      std::cout << "_varout[" << i << "] = " << exout[i] << std::endl;
#endif
    _dag->eval(_codelist, wk, _varout.size(), _varout.data(), valout,
               _varin.size(), _varin.data(), valin);
  }

  //! @brief Evaluate expression tree
  template <typename U>
  bool
  _reval(U* valin, U* valout, std::vector<U>& wk, U const& inf)
  //    const
  {
    int flag = _dag->reval(_codelist, wk, _varout.size(), _varout.data(),
                           valout, _varin.size(), _varin.data(), valin, inf,
                           options.CPMAX, options.CPTHRES);
#ifdef MC__DAGEXT_DEBUG
    std::cout << "DAGEXT:: Work array: " << flag << " passes\n";
    for (unsigned i = 0; i < _codelist.len_tap - _codelist.len_wrk; i++)
      std::cout << "wk[" << i << "] = " << wk[i] << std::endl;
    for (unsigned i = 0; i < _varin.size(); i++)
      std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif
    return (flag < 0 ? false : true);
  }

 public:
  //! @brief DAGEXT options
  struct Options
  {
    //! @brief Default constructor
    Options() { reset(); }

    //! @brief Assignment operator
    Options&
    operator=(Options const& other)
    {
      AUTODIFF = other.AUTODIFF;
      CPMAX    = other.CPMAX;
      CPTHRES  = other.CPTHRES;
      CPINF    = other.CPINF;

      return *this;
    }

    //! @brief Reset options
    void
    reset()
    {
      AUTODIFF = F;
      CPMAX =
          1;  // Leave it to the outer DAG to iterate as necessary by default
      CPTHRES = 1e4 * DBL_EPSILON;
      CPINF   = 1e30;
    }

    //! @brief Enumeration type for AD strategy
    enum AD_TYPE
    {
      F = 0,  //!< Forward differentiation
      B       //!< Backward differentiation
    };

    //! @brief Whether to apply forward or reverse automatic differentiation
    int AUTODIFF;
    //! @brief Maximum rounds of constraint propagation
    size_t CPMAX;
    //! @brief Threshold for repeating constraint propagation (minimum relative
    //! reduction in any variable)
    double CPTHRES;
    //! @brief Infinite value for unbounded variables in constraint propagation
    double CPINF;
  } options;

  //! @brief Default constructor
  DAGEXT() : _dag(nullptr), _uid(ffextern_uid()) {}

  //! @brief Data constructor
  DAGEXT(FFGraph* dag, std::vector<FFVar> const& varin,
         std::vector<FFVar> const& varout)
      : _dag(nullptr), _uid(0)
  {
    _set(dag, varin, varout);
  }

  //! @brief Data constructor
  DAGEXT(FFGraph* dag, std::set<unsigned> const& ndxin,
         std::vector<FFVar> const& varin, std::set<unsigned> const& ndxout,
         std::vector<FFVar> const& varout)
      : _dag(nullptr), _uid(0)
  {
    _set(dag, ndxin, varin, ndxout, varout);
  }

  //! @brief Copy constructor
  DAGEXT(DAGEXT const& other) : _dag(nullptr), _uid(0)
  {
    _set(other._dag, other._varin, other._varout);
    options = other.options;
    _uid    = other._uid;  // the same expression tree
  }

  virtual ~DAGEXT() { delete _dag; }

  //! @brief Set DAGEXT data
  void
  set(FFGraph* dag, std::vector<FFVar> const& varin,
      std::vector<FFVar> const& varout)
  {
    _set(dag, varin, varout);
  }

  //! @brief DAG automatic differentiation
  void
  grad(double const* valin, double* gradout, double* const valout = nullptr)
  {
    size_t const nin  = _varin.size();
    size_t const nout = _varout.size();

    switch (options.AUTODIFF)
    {
      default:
      case Options::AD_TYPE::F:
      {
        std::vector<FADType<double>> vFvalin(nin);
        for (size_t i = 0; i < nin; ++i)
        {
          vFvalin[i] = valin[i];
          vFvalin[i].diff(i, nin);
        }
        std::vector<FADType<double>> vFvalout(nout);
        _eval(vFvalin.data(), vFvalout.data(), _wkFD);
        for (size_t k = 0; k < nout; ++k)
        {
          if (valout) valout[k] = vFvalout[k].x();
          for (size_t i = 0; i < nin; ++i)
            gradout[k * nin + i] = vFvalout[k].d(i);
        }
        break;
      }

      case Options::AD_TYPE::B:
      {
        std::vector<BADType<double>> vBvalin(nin);
        for (size_t i = 0; i < nin; ++i) vBvalin[i] = valin[i];
        std::vector<BADType<double>> vBvalout(nout);
        _eval(vBvalin.data(), vBvalout.data(), _wkBD);
        for (size_t k = 0; k < nout; ++k) vBvalout[k].diff(k, nout);
        // AD reverse mode propagates adjoints when expression nodes
        // are released. The persistent work array holds
        // references to intermediatess, so release it before reading input
        // derivatives.
        _wkBD.clear();
        for (size_t k = 0; k < nout; ++k)
        {
          if (valout) valout[k] = vBvalout[k].x();
          for (size_t i = 0; i < nin; ++i)
            gradout[k * nin + i] = vBvalin[i].d(k);
        }
        break;
      }
    }
  }

  //! @brief DAG evaluation
  void
  eval(double const* valin, double* valout)
  {
    _eval(valin, valout, _wkD);
  }
  void
  eval(FADType<double> const* valin, FADType<double>* valout)
  {
    _eval(valin, valout, _wkFD);
  }
  void
  eval(T const* valin, T* valout)
  {
    _eval(valin, valout, _wkI);
  }
  void
  eval(McCormick<T> const* valin, McCormick<T>* valout)
  {
#ifdef MC__DAGEXT_DEBUG
    for (unsigned i = 0; i < _varin.size(); ++i)
      std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif
    _eval(valin, valout, _wkMC);
  }
  void
  eval(Specbnd<T> const* valin, Specbnd<T>* valout)
  {
    _eval(valin, valout, _wkSB);
  }
  void
  eval(SCVar<T> const* valin, SCVar<T>* valout)
  {
    _eval(valin, valout, _wkSC);
  }
  void
  eval(SupVar<PWCU> const* valin, SupVar<PWCU>* valout)
  {
    _eval(valin, valout, _wkPWCS);
  }
  void
  eval(SupVar<PWLU> const* valin, SupVar<PWLU>* valout)
  {
    _eval(valin, valout, _wkPWLS);
  }
  void
  eval(PolVar<T> const* valin, PolVar<T>* valout)
  {
    _eval(valin, valout, _wkPOL);
  }
  void
  eval(OCVar<double> const* valin, OCVar<double>* valout)
  {
    _eval(valin, valout, _wkOCD);
  }
  void
  eval(OCVar<FADType<double>> const* valin, OCVar<FADType<double>>* valout)
  {
    _eval(valin, valout, _wkOCFD);
  }
  template <typename U>
  void
  eval(U const* valin, U* valout)
  {
    static thread_local std::vector<U> wkU;
    _eval(valin, valout, wkU);
  }

  //! @brief DAG reverse evaluation
  bool
  reval(T* valin, T* valout)
  {
    return _reval(valin, valout, _wkCPI, options.CPINF * T(-1, 1));
  }

  //! @brief DAG differentiation
  std::tuple<std::vector<unsigned>, std::vector<unsigned>, std::vector<FFVar>>
  diff() const
  {
    switch (options.AUTODIFF)
    {
      case Options::F:
      default:
        return _dag->SFAD(_varout, _varin);
      case Options::B:
        return _dag->SBAD(_varout, _varin);
    }
  }

  //! @brief Stamp of the expression tree: a number that no other expression
  //! tree has had, drawn again whenever the tree is set and kept by the
  //! copies. The options are not part of it. A DAGEXT operation that holds a
  //! copy of the object is identified in a DAG by this stamp: see
  //! mc::FFEXTERN::lt
  size_t
  uid() const
  {
    return _uid;
  }

  //! @brief Query members
  size_t
  nin() const
  {
    return _varin.size();
  }
  size_t
  nout() const
  {
    return _varout.size();
  }
  std::vector<FFVar> const&
  varin() const
  {
    return _varin;
  }
  std::vector<FFVar> const&
  varout() const
  {
    return _varout;
  }
  FFGraph*
  dag() const
  {
    return _dag;
  }
};

//! @brief C++ structure for ordering of operation defined by an expression tree
template <typename T>
struct lt_DAGEXT
{
  bool
  operator()(DAGEXT<T> const* pDAGEXT1, DAGEXT<T> const* pDAGEXT2) const
  {
    // Order based on their number of input variables first
    if (pDAGEXT1->nin() < pDAGEXT2->nin()) return true;
    if (pDAGEXT1->nin() > pDAGEXT2->nin()) return false;
    // Order based on their number of ougput variables second
    if (pDAGEXT1->nout() < pDAGEXT2->nout()) return true;
    if (pDAGEXT1->nout() > pDAGEXT2->nout()) return false;
    // Order based on their input variables next
    for (unsigned i = 0; i < pDAGEXT1->nin(); i++)
    {
      if (lt_FFVar()(pDAGEXT1->varin()[i], pDAGEXT2->varin()[i])) return true;
      if (lt_FFVar()(pDAGEXT2->varin()[i], pDAGEXT1->varin()[i])) return false;
    }
    // Order based on their output variables last
    for (unsigned i = 0; i < pDAGEXT1->nout(); i++)
    {
      if (lt_FFVar()(pDAGEXT1->varout()[i], pDAGEXT2->varout()[i])) return true;
      if (lt_FFVar()(pDAGEXT2->varout()[i], pDAGEXT1->varout()[i]))
        return false;
    }
    return false;
  }
};

//! @brief C++ class defining expression tree as external DAG operations in
//! MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFDAGEXT is a C++ class for defining expression trees as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class FFDAGEXT
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, DAGEXT<T>>
{
 protected:
  using FFEXTERN<T, DAGEXT<T>>::_ptrObj;
  using FFEXTERN<T, DAGEXT<T>>::_ownObj;

  using typename FFEXTERN<T, DAGEXT<T>>::Home;

  // Jacobian of the expression tree of a home, kept there for all the
  // operations on that tree
  struct Jacobian
  {
    // stamp of the expression tree that it was made from, and its AD strategy
    size_t stamp;
    int mode;
    // row and column of the entries that are not structural zeros
    std::vector<unsigned> row, col;
    // home of the expression tree of these entries
    std::weak_ptr<Home> home;
  };

  // attach this handle to the home of the Jacobian of the expression tree of
  // parent. The Jacobian is differentiated on the first call for that tree --
  // and again if the tree or its AD strategy has changed, or once no
  // operation on the previous one remains -- in a home that has the options
  // of its parent; its expression tree takes the options of the parent tree
  // as they are at each call
  std::shared_ptr<Jacobian>
  _attach_jacobian(FFDAGEXT<T> const& parent)
  {
    if (!parent._home)
      throw std::runtime_error(
          "FFDAGEXT::_attach_jacobian ** Operation without a home\n");
    std::lock_guard<std::recursive_mutex> lock(this->_mutex());
    Home& home = *parent._home;
    std::shared_ptr<Jacobian> jac = std::static_pointer_cast<Jacobian>(home.aux);
    std::shared_ptr<Home> jachome;
    if (jac) jachome = jac->home.lock();
    if (!jachome || jac->stamp != home.obj->uid() ||
        jac->mode != home.obj->options.AUTODIFF)
    {
      auto&& sdervarout = home.obj->diff();
      jac               = std::make_shared<Jacobian>();
      jac->stamp        = home.obj->uid();
      jac->mode         = home.obj->options.AUTODIFF;
      jac->row          = std::get<0>(sdervarout);
      jac->col          = std::get<1>(sdervarout);
      this->_attach_derived(
          parent, new DAGEXT<T>(home.obj->dag(), home.obj->varin(),
                                std::get<2>(sdervarout)));
      jac->home = this->_home;
      home.aux  = jac;
    }
    else
      this->_set_home(jachome);
    _ptrObj->options = home.obj->options;
    return jac;
  }

  // operation on the Jacobian of the expression tree of this operation, on
  // the operands vVar: its results are the entries listed in jac
  FFVar const* const*
  _jacobian(unsigned const nVar, FFVar const* vVar,
            std::shared_ptr<Jacobian>& jac) const
  {
    FFDAGEXT<T> ResDer;
    jac = ResDer._attach_jacobian(*this);
    return ResDer.insert_external_operation(ResDer, jac->row.size(), nVar,
                                            vVar);
  }

  // set the object and related operation in DAG
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, DAGEXT<T>* pDAG, int policy)
  {
#ifdef MC__FFDAGEXT_CHECK
    assert(nVar == pDAG->nin());
#endif
    // home of the object in the DAG of the operands: made on its first
    // insertion there, with a copy of the object (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pDAG, policy);
    FFVar** ppRes = this->insert_external_operation(*this, _ptrObj->nout(),
                                                    nVar, pVar);
    return ppRes;
  }

  // set the object and related operation in DAG
  FFVar**
  _set(std::set<unsigned> const& ndxVar, FFVar const* pVar, DAGEXT<T>* pDAG,
       int policy)
  {
#ifdef MC__FFDAGEXT_CHECK
    assert(ndxVar.size() == pDAG->nin());
#endif
    FFBase* dag = nullptr;
    for (auto const& i : ndxVar)
      if (!dag) dag = pVar[i].dag();
    // home of the object in the DAG of the operands: made on its first
    // insertion there, with a copy of the object (policy > 0), found again
    // afterwards -- see mc::FFEXTERN::_attach
    this->owndata = false;
    this->_attach(dag, pDAG, policy);
    FFVar** ppRes = this->insert_external_operation(*this, _ptrObj->nout(),
                                                    ndxVar, pVar);
    return ppRes;
  }

 public:
  //! @brief Enumeration type for copy policy of DAGEXT object
  enum POLICY_TYPE
  {
    SHALLOW  = 0,  //!< The DAG refers to the DAGEXT object, which it does not own
    COPY     = 1,  //!< The DAG holds a copy of the DAGEXT object: one per DAG
    TRANSFER = -1  //!< The DAG refers to the DAGEXT object, and deletes it in the
                   //!< end
  };

  //! @brief Default constructor
  FFDAGEXT(bool const sparse = true) : FFEXTERN<T, DAGEXT<T>>()
  {
    this->sparse = sparse;
  }

  // Destructor
  virtual ~FFDAGEXT() {}

  // Copy constructor
  FFDAGEXT(FFDAGEXT<T> const& Other) : FFEXTERN<T, DAGEXT<T>>(Other) {}

  // Define operation
  std::vector<FFVar>
  operator()(std::set<unsigned> const& ndxVar, std::vector<FFVar> const& vVar,
             DAGEXT<T>* pDAG, int policy = COPY)
  {
    // return _set( vVar.size(), vVar.data(), pDAG, policy );
    FFVar** ppDep = _set(ndxVar, vVar.data(), pDAG, policy);
    std::vector<FFVar> vDep(pDAG->nout());
    for (size_t i = 0; i < vDep.size(); ++i) vDep[i] = *ppDep[i];
    return vDep;  // std::move( vDer );
  }

  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, DAGEXT<T>* pDAG, int policy = COPY)
  {
    // return _set( vVar.size(), vVar.data(), pDAG, policy );
    FFVar** ppDep = _set(vVar.size(), vVar.data(), pDAG, policy);
    std::vector<FFVar> vDep(pDAG->nout());
    for (size_t i = 0; i < vDep.size(); ++i) vDep[i] = *ppDep[i];
    return vDep;  // std::move( vDer );
  }

  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vVar, DAGEXT<T>* pDAG,
             int policy = COPY)
  {
#ifdef MC__FFDAGEXT_CHECK
    assert(idep < pDAG->nout());
#endif
    return *(_set(vVar.size(), vVar.data(), pDAG, policy)[idep]);
  }

  FFVar**
  operator()(size_t const nVar, FFVar const* pVar, DAGEXT<T>* pDAG,
             int policy = COPY)
  {
    return _set(nVar, pVar, pDAG, policy);
  }

  FFVar**
  operator()(std::set<unsigned> const& ndxVar, FFVar const* pVar,
             DAGEXT<T>* pDAG, int policy = COPY)
  {
    return _set(ndxVar, pVar, pDAG, policy);
  }

  FFVar&
  operator()(size_t const idep, size_t const nVar, FFVar const* pVar,
             DAGEXT<T>* pDAG, int policy = COPY)
  {
#ifdef MC__FFDAGEXT_CHECK
    assert(idep < pDAG->nout());
#endif
    return *(_set(nVar, pVar, pDAG, policy)[idep]);
  }

  // DAGEXT pointer: after an insertion, the DAGEXT object of the DAG -- the
  // copy that it holds (COPY), or the object given
  DAGEXT<T>*
  pDAGEXT() const
  {
    // std::cerr << "DAGEXT address retreived: " << _ptrObj << std::endl;
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
        "FFDAGEXT::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

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
        "FFDAGEXT::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  bool reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const;

  bool reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
             PolVar<T>* vVar) const;

  // Derivatives
  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
             FFVar const* vVar, FFVar** vDer) const;

  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
             FFVar const* vVar, FFVar** vDer, size_t* nnz,
             size_t** colnz) const;

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return "DAGEXT[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T>
inline void
FFDAGEXT<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                  FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: FFVar\n";
  std::cerr << "DAGEXT operation address: " << this << std::endl;
  std::cerr << "DAGEXT address in DAG: " << _ptrObj << std::endl;
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  FFVar** ppRes = this->_reinsert(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T>
inline void
FFDAGEXT<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                  FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: FFInv\n";
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
template <typename U>
inline void
FFDAGEXT<T>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                  unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  _ptrObj->eval(vVar, vRes);
}

template <typename T>
inline void
FFDAGEXT<T>::eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
                  FADType<FFVar> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: FADType<FFVar>\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  std::vector<FFVar> vVarVal(nVar);
  for (unsigned i = 0; i < nVar; ++i) vVarVal[i] = vVar[i].val();
  FFVar const* const* vResVal =
      this->_reinsert(*this, nRes, nVar, vVarVal.data());

  std::shared_ptr<Jacobian> jac;
  FFVar const* const* vResDer = _jacobian(nVar, vVarVal.data(), jac);

  // initialize forward AD and set all derivatives to zero
  for (unsigned k = 0; k < nRes; ++k)
  {
    vRes[k] = *vResVal[k];
    for (unsigned i = 0; i < nVar; ++i) vRes[k].setDepend(vVar[i]);
    for (unsigned j = 0; j < vRes[k].size(); ++j) vRes[k][j] = 0.;
  }

  // accumulate nonzero derivatives from sparse Jacobian
  for (unsigned int ie = 0; ie < jac->row.size(); ++ie)
  {
    for (unsigned j = 0; j < vRes[jac->row[ie]].size(); ++j)
    {
      auto const& k = jac->row[ie];
      auto const& i = jac->col[ie];
      if (vVar[i][j].cst() && vVar[i][j].num().val() == 0.) continue;
      vRes[k][j] += *vResDer[ie] * vVar[i][j];
      // vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
    }
  }
}

template <typename T>
inline void
FFDAGEXT<T>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                   FFVar const* vVar, FFVar** vDer, size_t* nnz,
                   size_t** colnz) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::deriv (sparse)\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
  assert(this->sparse && nnz && colnz);
#endif

  // std::tuple< std::vector<unsigned>, std::vector<unsigned>,
  // std::vector<FFVar> > diff
  std::shared_ptr<Jacobian> jac;
  FFVar const* const* vResDer = _jacobian(nVar, vVar, jac);

  // count non-zero dependencies for each output
  for (size_t k = 0; k < nRes; ++k) nnz[k] = 0;
  for (size_t ie = 0; ie < jac->row.size(); ++ie)
    ++nnz[jac->row[ie]];

  // resize each derivative and component arrays
  for (size_t k = 0; k < nRes; ++k)
  {
    assert(!vDer[k] && !colnz[k]);
    size_t const n = nnz[k];
    if (!n) continue;
    // n is bounded by the number of nonzero Jacobian entries; this assert
    // documents the real semantic bound (active only in debug builds).
    assert(n <= jac->row.size());
    // Explicit cap against a compile-time limit: hardens the allocation and
    // gives the compiler a provable upper bound so that n*sizeof(T) cannot
    // exceed the maximum object size (silences -Walloc-size-larger-than=).
    if (n > std::numeric_limits<std::ptrdiff_t>::max() / sizeof(FFVar))
      throw std::bad_alloc();
    vDer[k]  = new FFVar[n];
    colnz[k] = new size_t[n];
#ifdef MC__FFDAGEXT_DEBUG
    std::cout << "NNZ[" << k << "] = " << n << std::endl;
#endif
  }

  // copy nonzero derivatives from sparse Jacobian
  for (size_t k = 0; k < nRes; ++k) nnz[k] = 0;
  for (size_t ie = 0; ie < jac->row.size(); ++ie)
  {
    auto const& k    = jac->row[ie];
    auto const& i    = jac->col[ie];
    colnz[k][nnz[k]] = i;
    vDer[k][nnz[k]]  = *vResDer[ie];
    ++nnz[k];
  }
}

template <typename T>
inline void
FFDAGEXT<T>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                   FFVar const* vVar, FFVar** vDer) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::deriv (dense)\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(!this->sparse && _ptrObj && nRes == _ptrObj->nout() &&
         nVar == _ptrObj->nin());
#endif

  // std::tuple< std::vector<unsigned>, std::vector<unsigned>,
  // std::vector<FFVar> > diff
  std::shared_ptr<Jacobian> jac;
  FFVar const* const* vResDer = _jacobian(nVar, vVar, jac);

  // set all derivatives to zero
  for (size_t k = 0; k < nRes; ++k)
    for (size_t i = 0; i < nVar; ++i) vDer[k][i] = 0;

  // copy nonzero derivatives from sparse Jacobian
  for (size_t ie = 0; ie < jac->row.size(); ++ie)
  {
    auto const& k = jac->row[ie];
    auto const& i = jac->col[ie];
    vDer[k][i]    = *vResDer[ie];
  }
}

template <typename T>
inline void
FFDAGEXT<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                  SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: SLiftVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFDAGEXT<T>::eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                  PolVar<T> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::eval: PolVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag    = vVar[0].var().dag();
#ifdef MC__FFDAGEXT_CHECK
  assert(img && dag);
#endif
  FFVar** ppRes = dag->curOp()->varout.data();  // ACCOUNT FOR MULTIPLE OUTPUTS
#ifdef MC__FFDAGEXT_CHECK
  assert(nRes == dag->curOp()->varout.size());
#endif

  // at each relaxation, as FFMLP does: the options may have changed
  this->_resize_relax(img);
  this->_propagate_relax(img, ppRes, vRes, vVar);
}

template <typename T>
inline bool
FFDAGEXT<T>::reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                   PolVar<T>* vVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::reval: PolVar\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop      = vVar[0].var().opdef().first;
#ifdef MC__FFDAGEXT_CHECK
  assert(img && pop);
#endif

  this->_backpropagate_relax(img, pop, vRes, vVar);
  return true;
}

template <typename T>
inline bool
FFDAGEXT<T>::reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const
{
#ifdef MC__FFDAGEXT_TRACE
  std::cout << "FFDAGEXT::reval: T\n";
#endif
#ifdef MC__FFDAGEXT_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  return _ptrObj->reval(vVar, vRes);
}

}  // end namespace mc

#endif
