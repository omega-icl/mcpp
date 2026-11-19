// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__FFHOME_HPP
#define MC__FFHOME_HPP

#include <atomic>
#include <cstdint>
#include <map>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <utility>

#include "ffunc.hpp"

namespace mc
{

//! @brief Stamp for the content of an object held by external operations: a
//! number that no other content has had. An object that defines
//! <tt>size_t uid() const</tt> -- a stamp drawn here on construction and
//! again whenever its content changes, and kept by its copies -- has the
//! operations holding a copy of it (see mc::FFHOME) identified by that
//! content, whatever becomes of the original
inline size_t
ffhome_uid()
{
  static std::atomic<size_t> count{0};
  return ++count;
}

//! @brief How mc::FFHOME copies an object for a DAG (policy COPY, and the
//! copies of a DAG). The default uses the copy constructor; a class that
//! cannot be copied that way -- one that needs a DAG of its own for its
//! copy, say -- specializes this template, the shared pointer then owning
//! whatever the copy needs besides the object
template <typename O>
struct ffhome_traits
{
  static std::shared_ptr<O>
  copy(O const& obj)
  {
    return std::make_shared<O>(obj);
  }
};

//! @brief Options of the operations of a class without options of their own
struct ffhome_nooptions
{
};

//! @brief C++ base class of the external DAG operations on an object, which
//! gives the object one home per DAG
////////////////////////////////////////////////////////////////////////
//! mc::FFHOME is the base class of the external operations of MC++ and of
//! the libraries built on it that hold an object of type O -- a network
//! (mc::FFMLP), an expression tree (mc::FFDAGEXT), a solver (CANON's
//! FFNLPSLV) -- and options of type OPT for the operations on it.
//!
//! <b>Home of an object in a DAG.</b> A class of operations that calls
//! _attach before an insertion gives its object a home in the DAG: the
//! object -- a copy made for that DAG, or the object itself, according to
//! the policy of the insertion -- and the options of the operations on it.
//! There is one home per object and per DAG, shared by every operation on
//! that object in that DAG: on whichever operands, from whichever handle,
//! and the operations that differentiation makes on the same object (_join)
//! or on an object derived from it (_attach_derived). Hence
//! - one copy of the object per DAG, however many operations use it, made
//!   at its first insertion there (with mc::ffhome_traits); a copy of the
//!   DAG takes a copy of its own, and of the options;
//! - the member options of a handle, given to the homes that the handle
//!   makes, <i>are</i> the options of those homes: changed after the
//!   insertion, they change every operation of the home from the next
//!   evaluation that reads them. The options are not part of what
//!   identifies an operation. A handle that inserts an object which another
//!   handle has inserted in the DAG before is attached to the existing home,
//!   whose options are then reached with home_options();
//! - the object of the home is returned by object() -- MLP() for mc::FFMLP,
//!   DAGEXT() for mc::FFDAGEXT -- for the options that it holds itself:
//!   they too apply to every operation of the home, from the next
//!   evaluation or differentiation that reads them. Should the content of a
//!   copy be changed that way, the operations of its home follow it, and
//!   the original object, inserted again, gets a home of its own: a home
//!   whose copy no longer has the stamp it was filed under is not found any
//!   more (_current);
//! - the home keeps, for the class of the operation, whatever the operations
//!   on the object share besides (Home::aux): the expression tree of a
//!   Jacobian, a cache of solutions, ...
//!
//! An object is found again in a DAG by its key: the stamp of its content
//! where the DAG holds a copy and the object has a stamp (see
//! mc::ffhome_uid), so that an object changed since, or another one at the
//! same address, gets a home and operations of its own; its address
//! otherwise. The register of the homes is protected against concurrent
//! use; the objects are not: the operations of a home, like any single
//! operation, are to be evaluated by one thread at a time unless the object
//! allows otherwise.
//!
//! A class that sets _ptrObj and _ownObj itself, without _attach, keeps an
//! object and a copy of the options per operation.
////////////////////////////////////////////////////////////////////////
template <typename O, typename OPT = ffhome_nooptions>
class FFHOME
    ////////////////////////////////////////////////////////////////////////
    : public FFOp
{
 public:
  //! @brief Options of the operations on the object
  typedef OPT Options;
  struct Home;

  //! @brief What identifies an object in a DAG: the stamp of its content
  //! where the DAG holds a copy of an object that has one (see
  //! mc::ffhome_uid), the address of the object otherwise
  typedef std::pair<size_t, std::uintptr_t> t_Key;

 protected:
  // pointer to the object: the one of the home where the operation has one
  O* _ptrObj;
  // whether this operation owns the object -- operations without a home only
  bool _ownObj;
  // ownership of a copy of the object -- operations without a home only
  std::shared_ptr<void> _keep;
  // home of the object in the DAG, shared with every operation on that
  // object in that DAG
  std::shared_ptr<Home> _home;
  // options that the member options refers to: those of the home for an
  // operation in a DAG, those given to the homes it makes for a handle
  std::shared_ptr<OPT> _popt;

  //! @brief Stamp of the content of an object: its uid() where it has one, 0
  //! otherwise
  template <typename OBJ>
  static auto
  _stamp_of(OBJ const* pObj, int) -> decltype(size_t(pObj->uid()))
  {
    return pObj->uid();
  }
  template <typename OBJ>
  static size_t
  _stamp_of(OBJ const*, long)
  {
    return 0;
  }

  //! @brief Key of an object, whether or not the DAG holds a copy of it
  static t_Key
  _key(O const* pObj, bool const copy)
  {
    size_t const stamp = (copy ? _stamp_of(pObj, 0) : 0);
    return stamp ? t_Key(stamp, 0)
                 : t_Key(0, reinterpret_cast<std::uintptr_t>(pObj));
  }

  //! @brief Homes by DAG and key. A home leaves the register when the last
  //! operation or handle that shares it is destroyed
  typedef std::map<std::pair<FFBase const*, t_Key>, std::weak_ptr<Home>>
      t_Register;
  // (neither is ever destroyed: a home may outlive them at exit)
  static t_Register&
  _register()
  {
    static t_Register* reg = new t_Register;
    return *reg;
  }
  static std::recursive_mutex&
  _mutex()
  {
    static std::recursive_mutex* mtx = new std::recursive_mutex;
    return *mtx;
  }

  //! @brief Whether a home still holds the content it is filed under: false
  //! for a copy whose content was changed in the DAG, through object()
  static bool
  _current(Home const& home)
  {
    return !home.key.first || _stamp_of(home.obj, 0) == home.key.first;
  }

  //! @brief Home with key <a>key</a> in <a>dag</a>, if any. A home that no
  //! longer holds the content of its key leaves the register: its operations
  //! keep it, and the next insertion of that content makes another one
  static std::shared_ptr<Home>
  _home_find(FFBase const* dag, t_Key const& key)
  {
    std::lock_guard<std::recursive_mutex> lock(_mutex());
    auto it = _register().find(std::make_pair(dag, key));
    if (it == _register().end()) return nullptr;
    std::shared_ptr<Home> home = it->second.lock();
    if (home && !_current(*home))
    {
      _register().erase(it);
      return nullptr;
    }
    return home;
  }

  //! @brief New home in the register, for the object <a>pObj</a>, owned
  //! through <a>hold</a> if not null
  static std::shared_ptr<Home>
  _home_make(FFBase const* dag, t_Key const& key, O* pObj,
             std::shared_ptr<void> const& hold, std::shared_ptr<OPT> const& popt,
             std::shared_ptr<Home> const& parent)
  {
    std::shared_ptr<Home> home(new Home(dag, key, pObj, hold, popt, parent));
    std::lock_guard<std::recursive_mutex> lock(_mutex());
    _register()[std::make_pair(dag, key)] = home;
    return home;
  }

  //! @brief Home in <a>dag</a> of the object of <a>home</a>, a home of
  //! another DAG: found by its key, or made on the first visit. An object
  //! that its home owns is copied, an object that is referred to is not; the
  //! options are copied, so that each DAG has its own
  static std::shared_ptr<Home>
  _home_image(std::shared_ptr<Home> const& home, FFBase const* dag)
  {
    bool const own     = bool(home->hold);
    size_t const stamp = (own ? _stamp_of(home->obj, 0) : 0);
    t_Key const key    = (stamp ? t_Key(stamp, 0) : home->key);
    std::lock_guard<std::recursive_mutex> lock(_mutex());
    std::shared_ptr<Home> image = _home_find(dag, key);
    if (image) return image;
    std::shared_ptr<Home> parent;
    if (home->parent) parent = _home_image(home->parent, dag);
    std::shared_ptr<O> copy;
    if (own) copy = ffhome_traits<O>::copy(*home->obj);
    return _home_make(dag, key, own ? copy.get() : home->obj, copy,
                      parent ? parent->popt : std::make_shared<OPT>(*home->popt),
                      parent);
  }

  //! @brief Take a home
  void
  _set_home(std::shared_ptr<Home> const& home)
  {
    if (_ownObj && _ptrObj && !_keep) delete _ptrObj;
    _ownObj    = false;
    _keep.reset();
    _home      = home;
    _ptrObj    = home->obj;
    this->data = home->obj;
  }

  //! @brief DAG of the operands, null if none of them is in a DAG
  static FFBase*
  _dag_of(size_t const nVar, FFVar const* pVar)
  {
    for (size_t i = 0; i < nVar; ++i)
      if (pVar[i].dag()) return pVar[i].dag();
    return nullptr;
  }

  //! @brief Attach a handle to the home of the object <a>pObj</a> in
  //! <a>dag</a>, before the operation is inserted. The home is made on the
  //! first insertion of the object in that DAG -- with a copy of the object
  //! (<a>policy</a> > 0), the object itself (0), or the object itself, to be
  //! deleted with the home (< 0) -- and with the options of this handle,
  //! which stay those of the home. A later insertion finds the home again:
  //! by the stamp of the object where a copy was made and the object has one
  //! (see mc::ffhome_uid), by its address otherwise
  void
  _attach(FFBase const* dag, O* pObj, int const policy)
  {
    if (!dag) throw FFBase::Exceptions(FFBase::Exceptions::DAG);
    if (!pObj)
      throw std::runtime_error(
          "FFHOME::_attach ** Null pointer to external object\n");
    t_Key const key = _key(pObj, policy > 0);
    std::shared_ptr<Home> home;
    if (_home && _home->dag == dag && _home->key == key && _current(*_home))
      home = _home;
    else
    {
      std::lock_guard<std::recursive_mutex> lock(_mutex());
      home = _home_find(dag, key);
      if (!home)
      {
        std::shared_ptr<O> copy;
        if (policy > 0) copy = ffhome_traits<O>::copy(*pObj);
        home = _home_make(dag, key, policy > 0 ? copy.get() : pObj, copy,
                          _popt, nullptr);
      }
    }
    // an object handed over is deleted with its home, or now if the home
    // holds another object
    if (policy < 0)
    {
      if (home->obj == pObj)
      {
        if (!home->hold) home->hold = std::shared_ptr<O>(pObj);
      }
      else
        delete pObj;
    }
    _set_home(home);
  }

  //! @brief Attach a handle to the home of another operation: for the
  //! operations that differentiation makes on the same object
  void
  _join(FFHOME<O, OPT> const& parent)
  {
    if (!parent._home)
      throw std::runtime_error("FFHOME::_join ** Operation without a home\n");
    _set_home(parent._home);
  }

  //! @brief Attach a handle to a new home for <a>pObj</a>, an object made
  //! from the object of <a>parent</a> -- its derivative, say -- and deleted
  //! with its home. The home is in the DAG of the parent, has its options,
  //! and keeps the home of the parent alive
  void
  _attach_derived(FFHOME<O, OPT> const& parent, O* pObj)
  {
    if (!parent._home)
      throw std::runtime_error(
          "FFHOME::_attach_derived ** Operation without a home\n");
    _set_home(_home_make(parent._home->dag, _key(pObj, true), pObj,
                         std::shared_ptr<O>(pObj), parent._home->popt,
                         parent._home));
  }

  //! @brief Insert the operation <a>self</a>, this one, on other operands
  //! -- as its evaluation in the arithmetic of DAG variables does. Where the
  //! operands belong to another DAG, the operation is inserted there with
  //! the home of its object in that DAG: see _home_image
  template <typename ExtOp>
  FFVar**
  _reinsert(ExtOp const& self, unsigned const nRes, unsigned const nVar,
            FFVar const* pVar) const
  {
    FFBase const* dag = _dag_of(nVar, pVar);
    if (!_home || !dag || _home->dag == dag)
      return this->insert_external_operation(self, nRes, nVar, pVar);
    ExtOp other(self);
    other._set_home(_home_image(_home, dag));
    return other.insert_external_operation(other, nRes, nVar, pVar);
  }

 public:
  //! @brief Home of an object in a DAG: the object and the options of the
  //! operations on it, shared by every operation on that object in that DAG
  //! -- on whichever operands, and the operations that differentiation makes
  //! on the same object -- and by the handles attached to it
  struct Home
  {
    //! @brief DAG of the home
    FFBase const* dag;
    //! @brief Key of the object in the DAG: how the home is found
    t_Key key;
    //! @brief Number of the home, unique: how the operations of the home are
    //! told from those of any other (see lt), even one with the same key
    size_t serial;
    //! @brief Object: a copy made for the DAG, or the object given
    O* obj;
    //! @brief Ownership of the object, null where it is referred to: what is
    //! released with the home
    std::shared_ptr<void> hold;
    //! @brief Options of the operations on the object
    std::shared_ptr<OPT> popt;
    //! @brief Home of the object that this one was made from, if any
    std::shared_ptr<Home> parent;
    //! @brief For the class of the operation: what it keeps about the
    //! object for all its operations -- the home of its derivative, say
    std::shared_ptr<void> aux;

    Home(FFBase const* dag_, t_Key const& key_, O* obj_,
         std::shared_ptr<void> const& hold_, std::shared_ptr<OPT> const& popt_,
         std::shared_ptr<Home> const& parent_)
        : dag(dag_), key(key_), serial(ffhome_uid()), obj(obj_), hold(hold_),
          popt(popt_), parent(parent_)
    {
    }

    ~Home()
    {
      {
        std::lock_guard<std::recursive_mutex> lock(_mutex());
        auto it = _register().find(std::make_pair(dag, key));
        if (it != _register().end() && it->second.expired())
          _register().erase(it);
      }
      // what the class keeps about the object goes before the object
      aux.reset();
      hold.reset();
    }

   private:
    Home(Home const&);
    Home& operator=(Home const&);
  };

  //! @brief Options of the operation. For a handle, they are given to the
  //! homes that it makes -- the same options, not a copy: they may be
  //! changed after the insertion, and the operations follow. For an
  //! operation in a DAG, they are those of its home
  OPT& options;

  //! @brief Options of the home that the handle is attached to -- the object
  //! it inserted last: the member options, unless another handle had inserted
  //! that object in the DAG first, the options of a home being those of the
  //! handle that made it. The member options before any insertion
  OPT&
  home_options()
  {
    return _home ? *_home->popt : *_popt;
  }

  //! @brief Options of home_options(), shared: they are kept alive for as
  //! long as the pointer is held -- for interfaces to other languages
  std::shared_ptr<OPT>
  share_options() const
  {
    return _home ? _home->popt : _popt;
  }

  //! @brief Object of the home, shared: the home is kept alive for as long
  //! as the pointer is held -- for interfaces to other languages. Null
  //! without a home
  std::shared_ptr<O>
  share_object() const
  {
    return _home ? std::shared_ptr<O>(_home, _home->obj) : nullptr;
  }

  //! @brief Object of the operation: that of its home where it has one --
  //! the copy made for the DAG (policy COPY), or the object given. For a
  //! handle, the object of the home of its last insertion. Throws before
  //! any insertion
  O&
  object() const
  {
    if (!_ptrObj)
      throw std::runtime_error(
          "FFHOME::object ** No object: the operation was not inserted in a "
          "DAG\n");
    return *_ptrObj;
  }

  //! @brief Whether the handle is attached to a home, or the operation has
  //! one
  bool
  attached() const
  {
    return bool(_home);
  }

  //! @brief Default Constructor
  FFHOME()
      : FFOp(EXTERN),
        _ptrObj(nullptr),
        _ownObj(false),
        _popt(std::make_shared<OPT>()),
        options(*_popt)
  {
  }

  // Destructor
  virtual ~FFHOME()
  {
    if (_ownObj && _ptrObj && !_keep) delete _ptrObj;
  }

  // Copy constructor: the copy of an operation with a home shares it
  FFHOME(FFHOME<O, OPT> const& other)
      : FFOp(other),
        _ptrObj(nullptr),
        _ownObj(false),
        _home(other._home),
        _popt(other._home ? other._home->popt
                          : std::make_shared<OPT>(*other._popt)),
        options(*_popt)
  {
    if (_home)
    {
      _ptrObj = _home->obj;
      return;
    }

    // Operation without a home: its own object, or a copy of it
    if (!other._ptrObj)
      throw std::runtime_error(
          "FFHOME::copy constructor ** Null pointer to external object\n");

    _ownObj = other._ownObj;
    if (_ownObj)
    {
      std::shared_ptr<O> copy = ffhome_traits<O>::copy(*other._ptrObj);
      _ptrObj                 = copy.get();
      _keep                   = copy;
    }
    else
      _ptrObj = other._ptrObj;
  }

  FFHOME& operator=(FFHOME<O, OPT> const&) = delete;

  //! @brief Ordering of operations of the same class on the same operands:
  //! by their home in the DAG, that is by the object they were attached to
  //! by its key (see _attach). With a home that holds a copy of an object
  //! with a stamp, the key is that stamp, that is the content copied: the
  //! original may since have changed, or been destroyed and another object
  //! created at its address, and an insertion of either is another
  //! operation; an unchanged original, or a copy of it, is the same one.
  //! Otherwise the key is the address of the object. Two homes are told
  //! apart by their number -- also two with the same key, which arises when
  //! the copy of a home was changed (see _current). Operations without a
  //! home come first, ordered by the address in the field data, as with
  //! mc::FFOp::lt
  virtual bool
  lt(FFOp const* op) const
  {
    // lt_FFOp compares operations of the same class only
    auto const* other = static_cast<FFHOME<O, OPT> const*>(op);
    if (!_home || !other->_home)
      return (_home || other->_home) ? !_home : this->data < op->data;
    return _home->serial < other->_home->serial;
  }
};

}  // namespace mc

#endif
