// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__ELLIPSOID_HPP
#define MC__ELLIPSOID_HPP

#include <stdarg.h>

#include <cassert>
#include <iomanip>
#include <iostream>
#ifdef MC__ELLIPSOID_USE_SDPA
#include <sdpa_call.h>
#endif

#include <armadillo>
#include <cfloat>

#include "mcfunc.hpp"
#include "mcop.hpp"

#undef MC__ELLIPSOID_DEBUG_SQRT
#undef MC__ELLIPSOID_DEBUG_HPINTERSECTION
#undef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
#undef MC__ELLIPSOID_DEBUG_ELLINTERSECTION_IA
#undef MC__ELLIPSOID_DEBUG

namespace mc
{

class Ellipsoid;

template <typename T>
Ellipsoid minksum_ea(const Ellipsoid& E, const T* I, const double& TOL,
                     const double EPS = machprec());
Ellipsoid minksum_ea(const Ellipsoid& E1, const Ellipsoid& E2,
                     const double EPS = machprec());
Ellipsoid ellintersection_ia(
    const std::vector<std::pair<arma::vec, double> >& HP,
    const double& TOL = 1e-4, const unsigned& MAXIT = 100);

//! @brief C++ class for ellipsoidal calculus
////////////////////////////////////////////////////////////////////////
//! mc::Ellipsoid is a C++ class for ellipsoidal calculus. Round-off
//! errors are not accounted for in the computations (non-verified
//! implementation).
////////////////////////////////////////////////////////////////////////
class Ellipsoid
////////////////////////////////////////////////////////////////////////
{
  friend std::ostream& operator<<(std::ostream&, const Ellipsoid&);
  friend Ellipsoid ell_unitball(const unsigned n);
  friend Ellipsoid mtimes(const Ellipsoid& E, const arma::mat& A,
                          const arma::vec& b);
  friend Ellipsoid minksum_ea(const std::vector<Ellipsoid>& E,
                              const arma::vec& D);
  friend std::vector<Ellipsoid> minksum_ea(const std::vector<Ellipsoid>& E,
                                           const std::vector<arma::vec>& D);
  friend Ellipsoid minksum_ea(const Ellipsoid& E,
                              const std::pair<arma::vec, arma::vec>& I,
                              const double& TOL, const double EPS);
  template <typename T>
  friend Ellipsoid minksum_ea(const Ellipsoid& E, const T* I, const double& TOL,
                              const double EPS);
  friend Ellipsoid minksum_ea(const Ellipsoid& E1, const Ellipsoid& E2,
                              const double EPS);
  friend Ellipsoid inv(const Ellipsoid& E);
  friend std::vector<Ellipsoid> inv(const std::vector<Ellipsoid>& E);
  friend double dist(const Ellipsoid& E,
                     const std::pair<arma::vec, double>& HP);
  friend Ellipsoid hpintersection(const Ellipsoid& E,
                                  const std::pair<arma::vec, double>& HP);
  friend Ellipsoid hpintersection(
      const Ellipsoid& E, const std::vector<std::pair<arma::vec, double> >& HP);
  friend Ellipsoid intersection_ea(const Ellipsoid& E,
                                   const std::pair<arma::vec, double>& HP,
                                   const double& TOL);
  friend double dist(const Ellipsoid& E1, const Ellipsoid& E2);
  friend Ellipsoid intersection_ea(const Ellipsoid& E1, const Ellipsoid& E2,
                                   const double& TOL);
  friend Ellipsoid ellintersection_ia(
      const std::vector<std::pair<arma::vec, double> >& HP, const double& TOL,
      const unsigned& MAXIT);

 private:
  //! @brief Ellipsoid shape matrix (dense symmetric format)
  arma::mat _Q;
  //! @brief Ellipsoid center
  arma::vec _c;
  //! @brief Pointer to eigenvalues
  std::pair<arma::vec, arma::mat> _eigQ;
  //! @brief Pointer to square root of shape matrix
  arma::mat _sqrtQ;
  //! @brief Pointer to singular values
  std::pair<arma::vec, std::pair<arma::mat, arma::mat> > _svdQ;
  //! @brief Pointer to inverse of shape matrix
  arma::mat _invQ;

  //! @brief Whether shape matrix was checked for postive semi-definiteness
  bool _PSDchecked;

  //! @brief Delete pointers in class
  void
  _reset_auxiliary()
  {
    _eigQ.first.reset();
    _eigQ.second.reset();
    _sqrtQ.reset();
    _svdQ.first.reset();
    _svdQ.second.first.reset();
    _svdQ.second.second.reset();
    _invQ.reset();
    _PSDchecked = false;
  }

 public:
  /** @defgroup ELLIPSOID Ellipsoidal Calculus and Ellipsoidal Arithmetic for
   * Factorable Functions
   *  @{
   */
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in
  //! \mathbb R^n\f$, returns \f$(\mathcal E + d)\f$
  Ellipsoid&
  operator+=(const arma::vec& d)
  {
    _c += d;
    return *this;
  }
  //! @brief Given vector \f$d\in\mathbb R^n\f$ and ellipsoid \f$\mathcal E\in
  //! \mathbb R^n\f$, returns \f$(\mathcal E - d)\f$
  Ellipsoid&
  operator-=(const arma::vec& d)
  {
    _c -= d;
    return *this;
  }
  //! @brief Assignment operator
  Ellipsoid&
  operator=(const Ellipsoid& E)
  {
    _Q = E._Q;
    _c = E._c;
    _reset_auxiliary();
    return *this;
  }

  //! @brief Ellipsoid options
  static struct Options
  {
    //! @brief Constructor
    Options()
        : PSDCHK(false),
          PSDTOL(machprec() * 1e2),
          RKTOLA(machprec()),
          RKTOLR(mc::machprec()),
          ROOTTOL(1e-10),
          ROOTSECANT(false),
          ROOTMAXIT(0)
    {
    }
    //! @brief Whether or not to check positive semi-definiteness of shape
    //! matrices (default=false)
    bool PSDCHK;
    //! @brief Absolute tolerance for positive semi-definiteness check of shape
    //! matrix (default=1e2*mc::machprec())
    double PSDTOL;
    //! @brief Absolute tolerance for rank and regularization of shape matrix
    //! (default=mc::machprec())
    double RKTOLA;
    //! @brief Relative tolerance for rank and regularization of shape matrix
    //! (default=1e6*mc::machprec())
    double RKTOLR;
    //! @brief Absolute stopping tolerance for root-finding method (objective
    //! function value less than ROOTTOL; default=1e-10)
    double ROOTTOL;
    //! @brief Whether to use the secant method for root finding (default=false)
    bool ROOTSECANT;
    //! @brief Maximum number of iteration for root-finding method (default=0 -
    //! no maximum)
    unsigned ROOTMAXIT;
  } options;

  //! @brief Ellipsoid exceptions
  class Exceptions
  {
   public:
    //! @brief Enumeration type for Interval exception handling
    enum TYPE
    {
      NONPSD = 1,  //!< Non positive-semi definite shape matrix
      LAPACK,      //!< Linear algebra routine in LAPACK failed
      ROOT,        //!< Root-finding routine failed
    };
    //! @brief Constructor for error <a>ierr</a>
    Exceptions(TYPE ierr) : _ierr(ierr) {}
    //! @brief Return error flag
    int
    ierr()
    {
      return _ierr;
    }
    //! @brief Return error description
    std::string
    what()
    {
      switch (_ierr)
      {
        case NONPSD:
          return "Not positive semi-definite shape matrix";
        case LAPACK:
          return "Failure in LAPACK components";
        case ROOT:
        default:
          return "Failure in root-finding method";
      }
    }

   private:
    TYPE _ierr;
  };

  //! @brief Default constructor (needed to declare arrays of Ellipsoid class)
  Ellipsoid() : _Q(), _c(), _PSDchecked(false) {}

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with center \f$c\f$
  //! and shape matrix \f$Q\f$
  Ellipsoid(const arma::mat& Q, const arma::vec& c = arma::vec())
      : _Q(0.5 * (Q + Q.t())), _c(c), _PSDchecked(false)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!c.n_elem || Q.n_cols == c.n_elem);
#endif
    if (options.PSDCHK)
    {
      bool PSD;
      _isPSD(Q, _eigQ.first, PSD);
      if (!PSD)
      {
        std::cout << std::scientific << std::setprecision(10);
        std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !"
                  << std::endl;
        throw Exceptions(Exceptions::NONPSD);
      }
      _PSDchecked = true;
    }
    if (_c.n_elem != _Q.n_cols)
    {
      _c.resize(_Q.n_cols);
      _c.zeros();
    }
  }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ with center \f$c\f$
  //! and shape matrix \f$Q\f$ (lower triangular part stored contiguously and
  //! columnwise)
  Ellipsoid(const unsigned n, const double* Q, const double* c = 0)
      : _PSDchecked(false)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(n && Q);
#endif
    _Q.set_size(n, n);
    _Q.zeros();
    for (unsigned j = 0, k = 0; j < n; j++)
      for (unsigned i = j; i < n; i++) _Q(i, j) = Q[k++];
    _Q = arma::symmatl(_Q);
    _c.resize(_Q.n_cols);
    for (unsigned i = 0; i < n; i++) _c(i) = (c ? c[i] : 0.);

    if (options.PSDCHK)
    {
      bool PSD;
      _isPSD(_Q, _eigQ.first, PSD);
      if (!PSD)
      {
        std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !"
                  << std::endl;
        throw Exceptions(Exceptions::NONPSD);
      }
      _PSDchecked = true;
    }
  }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ enclosing interval
  //! vector of radius \f$r\f$ centered at \f$c\f$
  Ellipsoid(const arma::vec& r, const arma::vec& c = arma::vec())
      : _Q(r.n_elem, r.n_elem, arma::fill::zeros), _c(c), _PSDchecked(false)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!c.n_elem || r.n_elem == c.n_elem);
#endif
    if (options.PSDCHK)
    {
      for (arma::uword i = 0; i < r.n_elem; i++)
        if (r(i) < 0)
        {
          std::cout << "Interval radius: " << r(i) << " < 0 !" << std::endl;
          throw Exceptions(Exceptions::NONPSD);
        }
      _PSDchecked = true;
    }
    _Q.zeros();
    double nrm1_r = arma::norm(r, 1);
    for (arma::uword i = 0; i < r.n_elem; i++) _Q(i, i) = r(i) * nrm1_r;
    if (_c.n_elem != _Q.n_cols)
    {
      _c.resize(_Q.n_cols);
      _c.zeros();
    }
  }

  //! @brief Constructor for ellipsoid of dimension \f$n\f$ enclosing interval
  //! vector \f$B\f$
  template <typename T>
  Ellipsoid(const unsigned n, const T* B)
      : _Q(n, n, arma::fill::zeros), _c(n), _PSDchecked(false)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!n || B);
#endif
    if (options.PSDCHK) _PSDchecked = true;
    _Q.zeros();
    arma::vec r(n);
    for (unsigned i = 0; i < n; i++)
    {
      r(i)  = 0.5 * Op<T>::diam(B[i]);
      _c(i) = Op<T>::mid(B[i]);
    }
    // Min-trace axis-aligned ellipsoid enclosing the box [c-r,c+r]:
    // Q_ii = r_i * ||r||_1, so a box corner satisfies
    // sum_i r_i^2/Q_ii = ||r||_1/||r||_1 = 1 (it lies ON the boundary).
    // The previous ||r||_2 scaling under-sized the ellipsoid (corners gave
    // ||r||_1/||r||_2 > 1 for n>=2) and so did NOT enclose the box.
    double nrm1_r = arma::norm(r, 1);
    for (unsigned i = 0; i < n; i++) _Q(i, i) = r(i) * nrm1_r;
  }

  //! @brief Copy constructor
  Ellipsoid(const Ellipsoid& E) : _Q(E._Q), _c(E._c), _PSDchecked(false) {}

  //! @brief Destructor
  virtual ~Ellipsoid() {}

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape
  //! matrix \f$Q\f$
  Ellipsoid&
  unitball(const unsigned n)
  {
    _reset_auxiliary();
    _Q.eye(n, n);
    if (_c.n_elem != n)
    {
      _c.resize(_Q.n_cols);
      _c.zeros();
    }
    return *this;
  }

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape
  //! matrix \f$Q\f$
  Ellipsoid&
  set(const arma::mat& Q = arma::mat(), const arma::vec& c = arma::vec())
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!c.n_elem || Q.n_cols == c.n_elem);
#endif
    _reset_auxiliary();
    if (Q.n_cols && options.PSDCHK)
    {
      bool PSD;
      _isPSD(Q, _eigQ.first, PSD);
      if (!PSD)
      {
        std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !"
                  << std::endl;
        throw Exceptions(Exceptions::NONPSD);
      }
      _PSDchecked = true;
    }
    _Q = 0.5 * (Q + Q.t());
    _c = c;
    if (_c.n_elem != _Q.n_cols)
    {
      _c.resize(_Q.n_cols);
      _c.zeros();
    }
    return *this;
  }

  //! @brief Define ellipsoid of dimension \f$n\f$ with center \f$c\f$ and shape
  //! matrix \f$Q\f$ (lower triangular part stored contiguously and columnwise)
  Ellipsoid&
  set(const unsigned n, const double* Q, const double* c = 0)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(n && Q);
#endif
    _Q.set_size(n, n);
    _Q.zeros();
    for (unsigned j = 0, k = 0; j < n; j++)
      for (unsigned i = j; i < n; i++) _Q(i, j) = Q[k++];
    _Q = arma::symmatl(_Q);
    _c.resize(_Q.n_cols);
    for (unsigned i = 0; i < n; i++) _c(i) = (c ? c[i] : 0.);

    _reset_auxiliary();
    if (options.PSDCHK)
    {
      bool PSD;
      _isPSD(_Q, _eigQ.first, PSD);
      if (!PSD)
      {
        std::cout << "Min. eigenvalue: " << _eigQ.first(0) << " < 0 !"
                  << std::endl;
        throw Exceptions(Exceptions::NONPSD);
      }
      _PSDchecked = true;
    }
    return *this;
  }

  //! @brief Define ellipsoid of dimension \f$n\f$ enclosing interval vector of
  //! radius \f$r\f$ centered at \f$c\f$
  Ellipsoid&
  set(const arma::vec& r, const arma::vec& c = arma::vec())
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!c.n_elem || r.n_elem == c.n_elem);
#endif
    _reset_auxiliary();
    if (options.PSDCHK)
    {
      for (arma::uword i = 0; i < r.n_elem; i++)
        if (r(i) < 0)
        {
          std::cout << "Interval radius: " << r(i) << " < 0 !" << std::endl;
          throw Exceptions(Exceptions::NONPSD);
        }
      _PSDchecked = true;
    }
    _Q.set_size(r.n_elem, r.n_elem);
    _Q.zeros();
    double nrm1_r = arma::norm(r, 1);
    for (arma::uword i = 0; i < r.n_elem; i++) _Q(i, i) = r(i) * nrm1_r;
    _c = c;
    if (_c.n_elem != _Q.n_cols)
    {
      _c.resize(_Q.n_cols);
      _c.zeros();
    }
    return *this;
  }

  //! @brief Define ellipsoid of dimension \f$n\f$ enclosing interval vector of
  //! radius \f$r\f$ centered at \f$c\f$
  template <typename T>
  Ellipsoid&
  set(const unsigned n, const T* B)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(!n || B);
#endif
    arma::vec r(n), c(n);
    for (unsigned i = 0; i < n; i++)
    {
      r(i) = 0.5 * Op<T>::diam(B[i]);
      c(i) = Op<T>::mid(B[i]);
    }
    return set(r, c);
  }

  //! @brief Reset ellipsoid
  Ellipsoid&
  reset()
  {
    _Q.reset();
    _c.reset();
    _reset_auxiliary();
    return *this;
  }

  //! @brief Reset ellipsoid
  Ellipsoid&
  reset_aux()
  {
    _reset_auxiliary();
    return *this;
  }

  //! @brief Extend dimension by one, by appending row <a>Qi</a> to shape matrix
  //! and entry <a>ci</a> to center
  Ellipsoid&
  extend(const arma::rowvec& Qi, const double& ci = 0.)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(Qi.n_elem == _Q.n_cols + 1);
#endif
    _reset_auxiliary();
    arma::mat Qext(_Q.n_cols + 1, _Q.n_cols + 1, arma::fill::zeros);
    arma::vec cext(_Q.n_cols + 1);
    for (arma::uword i = 0; i < _Q.n_cols; i++)
    {
      cext(i) = _c(i);
      for (arma::uword j = 0; j <= i; j++) Qext(i, j) = _Q(i, j);
    }
    cext(_Q.n_cols) = ci;
    for (arma::uword i = 0; i <= _Q.n_cols; i++) Qext(_Q.n_cols, i) = Qi(i);
    _c = cext;
    _Q = arma::symmatl(Qext);
    return *this;
  }

  //! @brief Recenter ellipsoid at the origin by canceling out the centre
  Ellipsoid
  O() const
  {
    return Ellipsoid(_Q);
  }
  //! @brief Return dimension of ellipsoid
  unsigned
  n() const
  {
    return _Q.n_cols;
  }
  //! @brief Return center of ellipsoid
  const arma::vec&
  c() const
  {
    return _c;
  }
  //! @brief Return center of ellipsoid
  arma::vec&
  c()
  {
    return _c;
  }
  //! @brief Return center coefficient
  double
  c(unsigned i) const
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i < _Q.n_cols);
#endif
    return _c(i);
  }
  //! @brief Return shape matrix of ellipsoid
  const arma::mat&
  Q() const
  {
    return _Q;
  }
  //! @brief Return shape matrix coefficient
  double
  Q(unsigned i, unsigned j) const
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i < _Q.n_cols && j < _Q.n_cols);
#endif
    return _Q(i, j);
  }
  //! @brief Return/set shape matrix coefficient
  double&
  Q(unsigned i, unsigned j)
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i < _Q.n_cols && j < _Q.n_cols);
#endif
    return _Q(i, j);
  }

  //! @brief Return trace of shape matrix
  double
  trQ() const
  {
    if (!_Q.n_cols) return 0.;
    double tr(_Q(0, 0));
    for (arma::uword i = 1; i < _Q.n_cols; i++) tr += _Q(i, i);
    return tr;
  }

  //! @brief Return eigenvalues and eigenvectors of shape matrix
  const std::pair<arma::vec, arma::mat>&
  eigQ()
  {
    if (!_eigQ.first.n_elem || !_eigQ.second.n_cols)
      _eigen(_Q, _eigQ.first, _eigQ.second);
    return _eigQ;
  }

  //! @brief Return whether or not shape matrix is positive-semi definite
  bool
  psdQ()
  {
    bool PSD;
    _isPSD(_Q, _eigQ.first, PSD);
    return PSD;
  }

  //! @brief Return square root of shape matrix
  const arma::mat&
  sqrtQ(const bool complete = false)
  {
    if (!_sqrtQ.n_cols)
    {
      if (!_eigQ.first.n_elem || !_eigQ.second.n_cols)
        _eigen(_Q, _eigQ.first, _eigQ.second);
      _sqrt(_Q, _eigQ.first, _eigQ.second, _sqrtQ, _PSDchecked);
    }
    if (complete) _sqrtQ = 0.5 * (_sqrtQ + _sqrtQ.t());
    return _sqrtQ;
  }

  //! @brief Return rank of shape matrix
  unsigned
  rankQ()
  {
    svdQ();
    unsigned rank = _Q.n_cols;
    for (arma::uword i = 0; i < _Q.n_cols; i++, rank--)
      if (_svdQ.first(_Q.n_cols - i - 1) > options.RKTOLA &&
          _svdQ.first(_Q.n_cols - i - 1) > options.RKTOLR * _svdQ.first(0))
        break;
    return rank;
  }

  //! @brief Return singular value decomposition of shape matrix
  const std::pair<arma::vec, std::pair<arma::mat, arma::mat> >&
  svdQ()
  {
    if (!_svdQ.first.n_elem || !_svdQ.second.first.n_cols ||
        !_svdQ.second.second.n_cols)
      _svd(_Q, _svdQ.first, _svdQ.second.first, _svdQ.second.second);
    return _svdQ;
  }

  //! @brief Return pointer to regularized shape matrix
  const arma::mat&
  regQ()
  {
    const arma::uword r = rankQ();
#ifdef MC__ELLIPSOID_DEBUG_REGQ
    std::cout << "***Ellipsoid::regQ -- r = " << r << std::endl;
#endif
    if (r == _Q.n_cols) return _Q;

    arma::mat E(_Q.n_cols, _Q.n_cols, arma::fill::zeros);
    const double eps0 =
        std::max(options.RKTOLA, options.RKTOLR * svdQ().first(0));
    const double eps = std::nextafter(eps0, DBL_MAX);
    for (arma::uword i = 0; i < _Q.n_cols - r; i++) E(r + i, r + i) = eps;
    arma::mat UEUT = svdQ().second.first * E * svdQ().second.first.t();
#ifdef MC__ELLIPSOID_DEBUG_REGQ
    std::cout << "***Ellipsoid::regQ -- UEUT = " << UEUT << std::endl;
#endif
    for (arma::uword i = 0; i < _Q.n_cols; i++)
    {
      for (arma::uword j = 0; j <= i; j++)
        _Q(i, j) += 0.5 * (UEUT(i, j) + UEUT(j, i));
    }
    _Q = arma::symmatl(_Q);
    // The shape matrix has just been mutated, so every cached factorization
    // (eigQ, sqrtQ, svdQ, invQ) is now stale and MUST be invalidated;
    // otherwise their lazy getters would return decompositions of the
    // pre-regularization matrix.  This reset previously sat behind
    // MC__ELLIPSOID_DEBUG_REGQ, which silently returned stale caches in
    // ordinary builds.
    _reset_auxiliary();
    return _Q;
  }

  //! @brief Return pointer to inverse shape matrix
  const arma::mat&
  invQ()
  {
    if (!_invQ.n_cols) _inv(_Q, _invQ);
    return _invQ;
  }

  //! @brief Return pointer to inverse shape matrix
  const double&
  invQ(unsigned i, unsigned j)
  {
    if (!_invQ.n_cols) _inv(_Q, _invQ);
    return _invQ(i, j);
  }

  //! @brief Computes an orthogonal matrix rotating the vector <a>x</a> so that
  //! it is parallel to the vector <a>v</a>
  arma::mat
  align(const arma::vec& v, const arma::vec& x) const
  {
    if (v.n_elem != x.n_elem) return arma::mat();

    arma::mat vmat(v.n_elem, 1), xmat(x.n_elem, 1);
    for (arma::uword i = 0; i < v.n_elem; i++)
    {
      vmat(i, 0) = v(i);
      xmat(i, 0) = x(i);
    }
    arma::vec Sv, Sx;
    arma::mat Uv, Ux, VTv, VTx;
    _svd(vmat, Sv, Uv, VTv);
    _svd(xmat, Sx, Ux, VTx);
    return (Uv * VTv(0, 0)) * (Ux * VTx(0, 0)).t();
  }

  //! @brief Return lower bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double
  l(const unsigned i) const
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i >= 0 && i < _Q.n_cols);
#endif
    return _c(i) - std::sqrt(std::max(0., _Q(i, i)));
  }
  //! @brief Return upper bound for \f$x_i\f$ for index \f$i\in\{0,...,n-1\}\f$
  double
  u(const unsigned i) const
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i >= 0 && i < _Q.n_cols);
#endif
    return _c(i) + std::sqrt(std::max(0., _Q(i, i)));
  }
  //! @brief Return maximum radius for \f$x_i\f$ for index
  //! \f$i\in\{0,...,n-1\}\f$
  double
  r(const unsigned i) const
  {
#ifdef MC__ELLIPSOID_DEBUG
    assert(i >= 0 && i < _Q.n_cols);
#endif
    return std::sqrt(std::max(0., _Q(i, i)));
  }
  /** @} */

 private:
  //! @brief Wrapper to LAPACK function <TT>dsyev_</TT> doing eigenvalue
  //! decomposition of a symmetric matrix
  static void _eigen(const arma::mat& Q, arma::vec& D, arma::mat& U);

  //! @brief Wrapper to LAPACK function <TT>dsyev_</TT> doing eigenvalue
  //! decomposition of a symmetric matrix (eigenvalues only)
  static void _eigen(const arma::mat& Q, arma::vec& D);

  //! @brief Check if a symmetric matrix is positive definite
  static void _isPSD(const arma::mat& Q, arma::vec& D, bool& PSD);

  //! @brief Compute the square-root of a symmetric matrix
  static void _sqrt(const arma::mat& Q, arma::vec& D, arma::mat& U,
                    arma::mat& R, bool& PSDchecked);

  //! @brief Wrapper to LAPACK function <TT>_dgesvd</TT> doing singular value
  //! decomposition of a general matrix
  static void _svd(const arma::mat& Q, arma::vec& S, arma::mat& U,
                   arma::mat& VT);

  //! @brief Wrapper to LAPACK functions <TT>_dsytrf</TT> and <TT>_dsytri</TT>
  //! doing factorization and inversion of a symmetric indefinite matrix
  static void _inv(const arma::mat& Q, arma::mat& Qinv);

  //! @brief Prototype function for finding junction points in convex/concave
  //! envelopes of univariate terms
  typedef double(puniv)(const double x, const arma::vec& c1,
                        const arma::mat& W1, const arma::vec& c2,
                        const arma::mat& W2, const unsigned n);
  //! @brief Secant method for root finding
  static double _secant(const double x0, const double x1, const double xL,
                        const double xU, puniv f, const arma::vec& c1,
                        const arma::mat& W1, const arma::vec& c2,
                        const arma::mat& W2, const unsigned n);
  //! @brief Golden section search method for root finding
  static double _goldsect(const double xL, const double xU, puniv f,
                          const arma::vec& c1, const arma::mat& W1,
                          const arma::vec& c2, const arma::mat& W2,
                          const unsigned n);
  //! @brief Golden section search iterations
  static double _goldsect_iter(const bool init, const double a, const double fa,
                               const double b, const double fb, const double c,
                               const double fc, puniv f, const arma::vec& c1,
                               const arma::mat& W1, const arma::vec& c2,
                               const arma::mat& W2, const unsigned n);
  //! @brief Function whose root in the interval (0, 1) determines the minimal
  //! volume ellipsoid overapproximating the intersection of two ellipsoids
  static double _ell_fusionlambda(const double a, const arma::vec& c1,
                                  const arma::mat& W1, const arma::vec& c2,
                                  const arma::mat& W2, const unsigned n);
  //! @brief Function to find parameter value for minimal volume ellipsoid in
  //! ellispoid intersection
  static double _ell_get_lambda(const arma::vec& c1, const arma::mat& W1,
                                const arma::vec& c2, const arma::mat& W2,
                                const bool isHP);

  //! @brief Function to calculate the trace of matrix <a>W</a>
  template <typename M>
  static double
  _trace(const M& W)
  {
    const arma::mat A(W);
    double trace = 0.;
    for (arma::uword i = 0; i < A.n_rows && i < A.n_cols; i++) trace += A(i, i);
    return trace;
  }

  //! @brief Function to calculate the determinant of symmetric matrix given its
  //! eignevalues <a>eigW</a>
  static double
  _det(const arma::vec& eigW)
  {
    double det = 1.;
    for (arma::uword i = 0; i < eigW.n_elem; i++) det *= eigW(i);
    return det;
  }

  //! @brief Pause the program execution and prompt the user
  static void
  _pause()
  {
    double tmp;
    std::cout << "ENTER <1> TO CONTINUE" << std::endl;
    std::cin >> tmp;
  }
};

////////////////////////////////////////////////////////////////////////

inline Ellipsoid::Options Ellipsoid::options;

inline void
Ellipsoid::_eigen(const arma::mat& Q, arma::vec& D)
{
  if (Q.is_empty())
  {
    D.reset();
    return;
  }
  if (!arma::eig_sym(D, 0.5 * (Q + Q.t())))
    throw Exceptions(Exceptions::LAPACK);
}

inline void
Ellipsoid::_eigen(const arma::mat& Q, arma::vec& D, arma::mat& U)
{
  if (Q.is_empty())
  {
    D.reset();
    U.reset();
    return;
  }
  if (!arma::eig_sym(D, U, 0.5 * (Q + Q.t())))
    throw Exceptions(Exceptions::LAPACK);
}

inline void
Ellipsoid::_isPSD(const arma::mat& Q, arma::vec& D, bool& PSD)
{
  if (Q.is_empty())
  {
    PSD = false;
    D.reset();
    return;
  }
  _eigen(Q, D);

  PSD = true;
  if (D(0) < -std::fabs(options.PSDTOL))
  {
#ifdef MC__ELLIPSOID_DEBUG
    std::cout << "Min. eigenvalue: " << D(0) << " < 0 !" << std::endl;
    _pause();
#endif
    PSD = false;
  }
  for (arma::uword i = 0; i < D.n_elem; i++)
  {
    if (!options.PSDCHK && D(i) < 0.)
      D(i) = 0.;
    else
      break;  // b/c eigenvalues are returned in increasing order
  }
}

inline void
Ellipsoid::_sqrt(const arma::mat& Q, arma::vec& D, arma::mat& U,
                 arma::mat& sqrtQ, bool& PSDchecked)
{
  bool PSD;
  _isPSD(Q, D, PSD);
  if (options.PSDCHK && !PSDchecked && !PSD)
  {
    PSDchecked = true;
    throw Exceptions(Exceptions::NONPSD);
  }
  if (D.is_empty() || U.is_empty())
  {
    sqrtQ.reset();
    return;
  }

  arma::vec Dpos     = arma::clamp(D, 0.0, arma::datum::inf);
  arma::mat sqrtQ_ge = U * arma::diagmat(arma::sqrt(Dpos)) * U.t();
  sqrtQ              = 0.5 * (sqrtQ_ge + sqrtQ_ge.t());
#ifdef MC__ELLIPSOID_DEBUG_SQRT
  std::cout << "sqrtQ*sqrtQ:\n" << sqrtQ * sqrtQ;
#endif
}

inline void
Ellipsoid::_svd(const arma::mat& Q, arma::vec& S, arma::mat& U, arma::mat& VT)
{
  if (Q.is_empty())
  {
    S.reset();
    U.reset();
    VT.reset();
    return;
  }
  arma::mat V;
  if (!arma::svd(U, S, V, Q)) throw Exceptions(Exceptions::LAPACK);
  VT = V.t();
}

inline void
Ellipsoid::_inv(const arma::mat& Q, arma::mat& Qinv)
{
  if (Q.is_empty())
  {
    Qinv.reset();
    return;
  }
  const arma::mat S = 0.5 * (Q + Q.t());
  arma::mat R;
  if (!arma::inv_sympd(R, S) && !arma::inv(R, S)) R = arma::pinv(S);
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "\n***Ellispoid: Q\n" << Q;
  std::cout << "***Ellispoid: inv(Q) (ill-conditioned)\n" << R;
#endif
  arma::mat invRQ_R;
  if (!arma::inv(invRQ_R, R * S)) invRQ_R = arma::pinv(R * S);
  invRQ_R *= R;
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "***Ellispoid: inv(R*Q)*R\n" << invRQ_R;
#endif
  Qinv = 0.5 * (invRQ_R + invRQ_R.t());
#ifdef MC__ELLIPSOID_DEBUG_INV
  std::cout << "***Ellispoid: inv(Q) (proper)\n" << invRQ_R;
#endif
}

inline Ellipsoid
mtimes(const Ellipsoid& E, const arma::mat& A, const arma::vec& b = arma::vec())
{
#ifdef MC__ELLIPSOID_DEBUG
  assert(E._Q.n_cols == A.n_cols);
#endif
  Ellipsoid AE;
  if (!A.n_rows) return AE;

  // Transformed center
  AE._c = A * E._c;
  if (b.n_elem == A.n_rows) AE._c += b;

  // Transformed shape
  arma::mat AQAT = A * E._Q * A.t();
  AE._Q.set_size(A.n_rows, A.n_rows);
  for (arma::uword i = 0; i < A.n_rows; i++)
    for (arma::uword j = 0; j <= i; j++)
      AE._Q(i, j) = 0.5 * (AQAT(i, j) + AQAT(j, i));
  AE._Q = arma::symmatl(AE._Q);
#ifdef MC__ELLIPSOID_DEBUG
  bool PSD;
  AE._isPSD(AE._Q, AE._eigQ.first, PSD);
  assert(PSD);
#endif
  return AE;
}

inline Ellipsoid
minksum_ea(const std::vector<Ellipsoid>& E, const arma::vec& D)
{
  Ellipsoid EA;
  if (!E.size()) return EA;

  // Check size
  const unsigned n                                   = (*E.begin())._Q.n_cols;
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
#ifdef MC__ELLIPSOID_DEBUG
  for (; Ei != E.end(); ++Ei) assert((*Ei)._Q.n_cols == n);
  assert(D.n_elem == n);
#endif

  // Construct center and shape matrix
  Ei = E.begin();
  arma::vec c_EA((*Ei)._c);
  for (++Ei; Ei != E.end(); ++Ei) c_EA += (*Ei)._c;
  arma::mat Q_EA(n, n, arma::fill::zeros);
  double sum_sqrt_lTQil = 0.;
  for (Ei = E.begin(); Ei != E.end(); ++Ei)
  {
    double sqrt_lTQil = std::sqrt(arma::dot(D, ((*Ei)._Q * D)));
    sum_sqrt_lTQil += sqrt_lTQil;
    if (Ei == E.begin())
      Q_EA = (*Ei)._Q / sqrt_lTQil;
    else
      Q_EA += (*Ei)._Q / sqrt_lTQil;
  }
  Q_EA *= sum_sqrt_lTQil;
  EA.set(Q_EA, c_EA);

  return EA;
}

inline std::vector<Ellipsoid>
minksum_ea(const std::vector<Ellipsoid>& E, const std::vector<arma::vec>& D)
{
  std::vector<Ellipsoid> EA;
  if (!E.size() || !D.size()) return EA;

  // Check size
  const unsigned n                                   = (*E.begin())._Q.n_cols;
  typename std::vector<arma::vec>::const_iterator Di = D.begin();
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
#ifdef MC__ELLIPSOID_DEBUG
  for (; Ei != E.end(); ++Ei) assert((*Ei)._Q.n_cols == n);
  for (; Di != D.end(); ++Di) assert((*Di).n_elem == n);
#endif

  // Construct centers and shape matrices
  Ei = E.begin();
  arma::vec c_EA((*Ei)._c);
  for (++Ei; Ei != E.end(); ++Ei) c_EA += (*Ei)._c;
  arma::mat Q_EA(n, n, arma::fill::zeros);
  for (Di = D.begin(); Di != D.end(); ++Di)
  {
    double sum_sqrt_lTQil = 0.;
    for (Ei = E.begin(); Ei != E.end(); ++Ei)
    {
      double sqrt_lTQil = std::sqrt(arma::dot((*Di), ((*Ei)._Q * (*Di))));
      sum_sqrt_lTQil += sqrt_lTQil;
      if (Ei == E.begin())
        Q_EA = (*Ei)._Q / sqrt_lTQil;
      else
        Q_EA += (*Ei)._Q / sqrt_lTQil;
    }
    Q_EA *= sum_sqrt_lTQil;
    EA.push_back(Ellipsoid(Q_EA, c_EA));
  }

  return EA;
}

template <typename T>
inline Ellipsoid
minksum_ea(const Ellipsoid& E, const T* I, const double& TOL, const double EPS)
{
  if (!E._Q.n_cols) return E;

  arma::vec c(E._Q.n_cols), r(E._Q.n_cols);
  for (arma::uword i = 0; i < E._Q.n_cols; i++)
  {
    c(i) = Op<T>::mid(I[i]);
    r(i) = 0.5 * Op<T>::diam(I[i]);
  }
  return minksum_ea(E, std::make_pair(r, c), TOL, EPS);
}

inline Ellipsoid
minksum_ea(const Ellipsoid& E, const std::pair<arma::vec, arma::vec>& I,
           const double& TOL, const double EPS = machprec())
{
#ifdef MC__ELLIPSOID_DEBUG
  assert(E._Q.n_cols == I.first.n_elem &&
         (!I.second.n_elem || E._Q.n_cols == I.second.n_elem));
#endif
  Ellipsoid EA;
  if (!E._Q.n_cols) return EA;

  // Construct center and shape matrix
  arma::vec c_EA(E._c);
  if (I.second.n_elem) c_EA += I.second;

  double trQ_scal = 0.0;
  arma::vec rad_scal(E._Q.n_cols);
  for (arma::uword i = 0; i < E._Q.n_cols; i++)
  {
    trQ_scal += E._Q(i, i) / (E._Q(i, i) + TOL);
    rad_scal(i) = I.first(i) / std::sqrt(E._Q(i, i) + TOL) + EPS;
  }
  double sqrt_trQ_scal = std::sqrt(trQ_scal) + EPS;
  double sum           = sqrt_trQ_scal;
  for (arma::uword i = 0; i < E._Q.n_cols; i++) sum += rad_scal(i);

  arma::mat Q_EA(E._Q / sqrt_trQ_scal);
  for (arma::uword i = 0; i < E._Q.n_cols; i++)
    Q_EA(i, i) += sqr(I.first(i)) / rad_scal(i);
  Q_EA *= sum;
  /*
    double trQ = 0.0;
    arma::vec sqrR(E._Q.n_cols);
    for( arma::uword i=0; i<E._Q.n_cols; i++ ){
      trQ += E._Q(i,i)/(E._Q(i,i)+TOL);
      sqrR(i) = I.first(i)/std::sqrt(E._Q(i,i)+TOL);
      // std::cout << "std::sqrt(E._Q(i,i)+TOL) =" << std::sqrt(E._Q(i,i)+TOL)
      //           << "  " << E._Q(i,i)+TOL << std::endl;
    }
    trQ = std::sqrt(trQ);
    double kappa = trQ;
    for( arma::uword i=0; i<E._Q.n_cols; i++ ) kappa += sqrR(i);

    arma::mat Q_EA( E._Q*kappa/(trQ+EPS) );
    for( arma::uword i=0; i<E._Q.n_cols; i++ )
      Q_EA(i,i) += I.first(i)*I.first(i)*kappa/(sqrR(i)+sqr(EPS))+sqr(EPS);
  */
  return EA.set(Q_EA, c_EA);
}

inline Ellipsoid
minksum_ea(const Ellipsoid& E1, const Ellipsoid& E2, const double EPS)
{
#ifdef MC__ELLIPSOID_DEBUG
  assert(E1._Q.n_cols == E2._Q.n_cols);
#endif
  Ellipsoid EA;
  const int n = E1._Q.n_cols;
  if (!n) return EA;

  // Construct center and shape matrix
  arma::vec c_EA = E1._c + E2._c;

  double trQ1 = 0., trQ2 = 0.;
  for (int i = 0; i < n; i++)
  {
    trQ1 += E1._Q(i, i);
    trQ2 += E2._Q(i, i);
  }
  double s1 = std::sqrt(trQ1) + EPS, s2 = std::sqrt(trQ2) + EPS;
  arma::mat Q_EA = (E1._Q / s1 + E2._Q / s2) * (s1 + s2);

  return EA.set(Q_EA, c_EA);
}

inline Ellipsoid
inv(const Ellipsoid& E)
{
  Ellipsoid Einv(E);
  if (!E._Q.n_cols) return Einv;
  Einv.regQ();
  Einv._Q = Einv.invQ();
  Einv._reset_auxiliary();
  return Einv;
}

inline std::vector<Ellipsoid>
inv(const std::vector<Ellipsoid>& E)
{
  std::vector<Ellipsoid> Einv;
  typename std::vector<Ellipsoid>::const_iterator Ei = E.begin();
  for (; Ei != E.end(); ++Ei) Einv.push_back(inv(*Ei));
  return Einv;
}

inline double
dist(const Ellipsoid& E, const std::pair<arma::vec, double>& HP)
{
  if (E._Q.n_cols != HP.first.n_elem)
  {
    std::cerr << "mc::dist: ellipsoid and hyperplane must be of the same "
                 "dimension.\n";
    return 0. / 0.;  // return NaN
  }

  return (std::fabs(HP.second - arma::dot(E._c, HP.first)) -
          std::sqrt(arma::dot(HP.first, (E._Q * HP.first)))) /
         (std::sqrt(arma::dot(HP.first, HP.first)));
}

inline double
dist(const Ellipsoid& E1, const Ellipsoid& E2)
{
  if (E1._Q.n_cols != E2._Q.n_cols)
  {
    std::cerr << "mc::dist: ellipsoids must be of the same dimension.\n";
    return 0. / 0.;  // return NaN
  }
  const arma::uword n = E1._Q.n_cols;
  if (!n) return 0.;

  // Separation test for two (possibly degenerate) ellipsoids E_i = { x :
  // (x-c_i)^T Q_i^{-1} (x-c_i) <= 1 }, following the 1-D parametrisation of
  // Gilitschenski & Hanebeck (2012).  With
  //     K(s) = 1 - s(1-s) (c2-c1)^T ( s Q2 + (1-s) Q1 )^{-1} (c2-c1),
  // the two ellipsoids overlap iff min_{s in [0,1]} K(s) >= 0.  We therefore
  // return -min_s K(s): the value is > 0 iff the ellipsoids are separated,
  // which is exactly the convention expected by intersection_ea() and the
  // hyperplane dist() overload above.  This replaces the previous stub that
  // unconditionally returned 0 (i.e. always assumed intersection).
  const arma::vec dc = E2._c - E1._c;
  if (arma::dot(dc, dc) <= 0.) return -1.;  // coincident centres -> overlap

  const arma::mat& Q1 = E1._Q;
  const arma::mat& Q2 = E2._Q;
  auto Kval           = [&](const double s) -> double
  {
    const arma::mat M = s * Q2 + (1. - s) * Q1;
    arma::vec x;
    if (!arma::solve(x, M, dc)) return 1.;  // singular blend -> skip
    return 1. - s * (1. - s) * arma::dot(dc, x);
  };

  // Coarse scan over the open interval, then a ternary refinement around the
  // best node (K is smooth on (0,1) with K(0)=K(1)=1, so the minimum is
  // interior).
  const unsigned NS = 256;
  double sbest = 0.5, Kmin = 1.;
  for (unsigned k = 1; k < NS; ++k)
  {
    const double s = double(k) / double(NS);
    const double K = Kval(s);
    if (K < Kmin)
    {
      Kmin  = K;
      sbest = s;
    }
  }
  double lo = std::max(0., sbest - 1. / double(NS));
  double hi = std::min(1., sbest + 1. / double(NS));
  for (unsigned it = 0; it < 60 && hi - lo > 1e-12; ++it)
  {
    const double m1 = lo + (hi - lo) / 3., m2 = hi - (hi - lo) / 3.;
    const double f1 = Kval(m1), f2 = Kval(m2);
    if (f1 < f2)
    {
      hi = m2;
      if (f1 < Kmin) Kmin = f1;
    }
    else
    {
      lo = m1;
      if (f2 < Kmin) Kmin = f2;
    }
  }
  return -Kmin;  // > 0  <=>  ellipsoids are separated
}

inline Ellipsoid
intersection_ea(const Ellipsoid& E1, const Ellipsoid& E2,
                const double& TOL = machprec())
{
  if (E1._Q.n_cols != E2._Q.n_cols)
  {
    std::cerr
        << "mc::intersection_ea: ellipsoids must be of the same dimension.\n";
    return Ellipsoid();
  }
  const unsigned n = E1._Q.n_cols;
  if (!n || dist(E1, E2) > 0.) return Ellipsoid();

  Ellipsoid EM1       = inv(E1);
  const arma::vec& c1 = EM1._c;
  const arma::mat& W1 = EM1.invQ();
  Ellipsoid EM2       = inv(E2);
  const arma::vec& c2 = EM2._c;
  const arma::mat& W2 = EM2.invQ();

  double a = Ellipsoid::_ell_get_lambda(c1, W1, c2, W2, false);
  // #ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "a (ell_get_lambda):" << a << std::endl;
  // #endif
  arma::mat W = a * W1 + (1 - a) * W2;
  arma::mat Winv(n, n, arma::fill::zeros);
  Ellipsoid::_inv(W, Winv);
  double k =
      1. - a * (1 - a) * arma::dot((c2 - c1), (W2 * Winv * W1 * (c2 - c1)));
  arma::vec q = Winv * (a * W1 * c1 + (1 - a) * W2 * c2);
  arma::mat Q = (1 + TOL) * k * Winv;
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "W:\n" << W << std::endl;
  std::cout << "Winv:\n" << Winv << std::endl;
  std::cout << "k:" << k << std::endl;
  std::cout << "q:" << q << std::endl;
  std::cout << "Q:\n" << Q << std::endl;
#endif
  return Ellipsoid(Q, q);
}

inline Ellipsoid
intersection_ea(const Ellipsoid& E, const std::pair<arma::vec, double>& HP,
                const double& TOL = machprec())
{
  if (E._Q.n_cols != HP.first.n_elem)
  {
    std::cerr << "mc::intersection_ea: ellipsoid and hyperplane must be of the "
                 "same dimension.\n";
    return Ellipsoid();
  }
  const unsigned n = E._Q.n_cols;
  if (!n) return Ellipsoid();

  arma::vec v = -HP.first / std::sqrt(arma::dot(HP.first, HP.first));
  double c    = -HP.second / std::sqrt(arma::dot(HP.first, HP.first));
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "dist(E,HP): " << dist(E, HP) << std::endl;
  std::cout << "HP.vT E.c > 0 ? " << (arma::dot(v, E._c) > c ? 'T' : 'F')
            << std::endl;
#endif
  if (dist(E, HP) > 0 && arma::dot(v, E._c) > c) return E;
  if (dist(E, HP) > 0 && arma::dot(v, E._c) < c) return Ellipsoid();
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "E:\n" << E << std::endl;
  std::cout << "HP:\n" << HP.first << "; " << HP.second << std::endl;
#endif

  Ellipsoid EM(E);
  arma::mat W2(n, n, arma::fill::zeros);
  arma::mat vmat(n, 1);
  for (unsigned i = 0; i < n; i++) vmat(i, 0) = v(i);
  arma::mat vmatvmatT = vmat * vmat.t();
  for (unsigned i = 0; i < n; i++)
    for (unsigned j = 0; j <= i; j++)
      W2(i, j) = 0.5 * (vmatvmatT(i, j) + vmatvmatT(j, i)) /
                 (4. * EM.eigQ().first(n - 1));
  Ellipsoid Einter = hpintersection(EM, HP);
  double h         = 2. * std::sqrt(EM.eigQ().first(n - 1));
  arma::vec c2     = Einter._c + h * v;
  Ellipsoid EM2(W2, c2);
  EM2.regQ();
  W2 = EM2._Q;
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "W2:\n" << W2 << std::endl;
  std::cout << "c2:\n" << c2 << std::endl;
#endif

  const arma::vec& c1 = EM._c;
  EM.regQ();
  const arma::mat& W1 = EM.invQ();
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "E (reg):\n" << EM << std::endl;
  std::cout << "rank(E._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "inv(E._Q):\n" << EM.invQ() << std::endl;
#endif

  double a = Ellipsoid::_ell_get_lambda(c1, W1, c2, W2, true);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "a (ell_get_lambda):" << a << std::endl;
#endif
  arma::mat W = a * W1 + (1 - a) * W2;
  arma::mat Winv(n, n, arma::fill::zeros);
  Ellipsoid::_inv(W, Winv);
  double k =
      1. - a * (1 - a) * arma::dot((c2 - c1), (W2 * Winv * W1 * (c2 - c1)));
  arma::vec q = Winv * (a * W1 * c1 + (1 - a) * W2 * c2);
  arma::mat Q = (1 + TOL) * k * Winv;
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "W:\n" << W << std::endl;
  std::cout << "Winv:\n" << Winv << std::endl;
  std::cout << "k:" << k << std::endl;
  std::cout << "q:" << q << std::endl;
  std::cout << "Q:\n" << Q << std::endl;
#endif
  return Ellipsoid(Q, q);
}

inline double
Ellipsoid::_ell_get_lambda(const arma::vec& c1, const arma::mat& W1,
                           const arma::vec& c2, const arma::mat& W2,
                           const bool isHP)
{
  assert(c1.n_elem == W1.n_cols && c2.n_elem == W2.n_cols &&
         c1.n_elem == c2.n_elem);
  const unsigned n = c1.n_elem;

  double a;
  if (options.ROOTSECANT)
  {
    try
    {
      a = _secant(.9, 1., -.01, 1.01, _ell_fusionlambda, c1, W1, c2, W2, n);
    }
    catch (Ellipsoid::Exceptions)
    {
      a = _goldsect(-.01, 1.01, _ell_fusionlambda, c1, W1, c2, W2, n);
    }
  }
  else
    a = _goldsect(-.1, 1.01, _ell_fusionlambda, c1, W1, c2, W2, n);

  if (a >= 0. && a <= 1.) return a;
  if (isHP) return 1.;
  arma::vec eigW1;
  _eigen(W1, eigW1);
  double detW1 = Ellipsoid::_det(eigW1);
  arma::vec eigW2;
  _eigen(W2, eigW2);
  double detW2 = Ellipsoid::_det(eigW2);
  if (detW1 > detW2) return 1.;
  return 0.;
}

inline double
Ellipsoid::_secant(const double x0, const double x1, const double xL,
                   const double xU, puniv f, const arma::vec& c1,
                   const arma::mat& W1, const arma::vec& c2,
                   const arma::mat& W2, const unsigned n)
{
  double xkm = std::max(xL, std::min(xU, x0));
  double fkm = f(xkm, c1, W1, c2, W2, n);
  double xk  = std::max(xL, std::min(xU, x1));

  for (unsigned it = 0; !options.ROOTMAXIT || it < options.ROOTMAXIT; it++)
  {
    double fk = f(xk, c1, W1, c2, W2, n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
    std::cout << "xk (fk):" << xk << "  (" << fk << ")" << std::endl;
#endif
    if (std::fabs(fk) < options.ROOTTOL) return xk;
    double Bk = (fk - fkm) / (xk - xkm);
    if (Bk == 0) throw Exceptions(Exceptions::ROOT);
    if (isequal(xk, xL) && fk / Bk > 0) return xk;
    if (isequal(xk, xU) && fk / Bk < 0) return xk;
    xkm = xk;
    fkm = fk;
    xk  = std::max(xL, std::min(xU, xk - fk / Bk));
  }

  throw Exceptions(Exceptions::ROOT);
}

inline double
Ellipsoid::_goldsect(const double xL, const double xU, puniv f,
                     const arma::vec& c1, const arma::mat& W1,
                     const arma::vec& c2, const arma::mat& W2, const unsigned n)
{
  const double phi = 2. - (1. + std::sqrt(5.)) / 2.;
  const double fL = f(xL, c1, W1, c2, W2, n), fU = f(xU, c1, W1, c2, W2, n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "x (f):" << xL << "  (" << fL << ")" << std::endl;
  std::cout << "x (f):" << xU << "  (" << fU << ")" << std::endl;
#endif
  if (fL * fU > 0) throw Exceptions(Exceptions::ROOT);
  const double xm = xU - phi * (xU - xL), fm = f(xm, c1, W1, c2, W2, n);
  return _goldsect_iter(true, xL, fL, xm, fm, xU, fU, f, c1, W1, c2, W2, n);
}

inline double
Ellipsoid::_goldsect_iter(const bool init, const double a, const double fa,
                          const double b, const double fb, const double c,
                          const double fc, puniv f, const arma::vec& c1,
                          const arma::mat& W1, const arma::vec& c2,
                          const arma::mat& W2, const unsigned n)
// a and c are the current bounds; the minimum is between them.
// b is a center point
{
  thread_local static unsigned iter;
  iter             = (init ? 1 : iter + 1);
  const double phi = 2. - (1. + std::sqrt(5.)) / 2.;
  bool b_then_x    = (c - b > b - a);
  double x         = (b_then_x ? b + phi * (c - b) : b - phi * (b - a));
  // if( std::fabs(c-a) < options.ROOTTOL*(std::fabs(b)+std::fabs(x))
  //  || iter > options.ROOTMAXIT ) return (c+a)/2.;
  double fx = f(x, c1, W1, c2, W2, n);
#ifdef MC__ELLIPSOID_DEBUG_INTERSECTION_EA
  std::cout << "x (f):" << x << "  (" << fx << ")" << std::endl;
#endif
  if (std::fabs(fx) < options.ROOTTOL ||
      (options.ROOTMAXIT && iter > options.ROOTMAXIT))
    return x;
  if (b_then_x)
    return (
        fa * fx < 0
            ? _goldsect_iter(false, a, fa, b, fb, x, fx, f, c1, W1, c2, W2, n)
            : _goldsect_iter(false, b, fb, x, fx, c, fc, f, c1, W1, c2, W2, n));
  return (
      fa * fb < 0
          ? _goldsect_iter(false, a, fa, x, fx, b, fb, f, c1, W1, c2, W2, n)
          : _goldsect_iter(false, x, fx, b, fb, c, fc, f, c1, W1, c2, W2, n));
}

inline double
Ellipsoid::_ell_fusionlambda(const double a, const arma::vec& c1,
                             const arma::mat& W1, const arma::vec& c2,
                             const arma::mat& W2, const unsigned n)
{
  arma::mat W = a * W1 + (1 - a) * W2;
  arma::mat Winv(n, n, arma::fill::zeros);
  Ellipsoid::_inv(W, Winv);
  double k =
      1. - a * (1 - a) * arma::dot((c2 - c1), (W2 * Winv * W1 * (c2 - c1)));
  arma::vec c = Winv * (a * W1 * c1 + (1 - a) * W2 * c2);

  arma::vec eigW;
  _eigen(W, eigW);
  double detW = Ellipsoid::_det(eigW);
  double traW = Ellipsoid::_trace(detW * Winv * (W1 - W2));
  return k * detW * traW -
         n * (detW * detW) *
             (arma::dot(c, (2 * W1 * c1 - 2 * W2 * c2 + (W2 - W1) * c)) -
              arma::dot(c1, (W1 * c1)) + arma::dot(c2, (W2 * c2)));
}

inline Ellipsoid
hpintersection(const Ellipsoid& E, const std::pair<arma::vec, double>& HP)
{
  if (E._Q.n_cols != HP.first.n_elem)
  {
    std::cerr << "mc::intersection_ea: ellipsoid and hyperplane must be of the "
                 "same dimension.\n";
    return Ellipsoid();
  }
  const unsigned n = E._Q.n_cols;
  if (!n) return Ellipsoid();

#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "dist(E,HP): " << dist(E, HP) << std::endl;
#endif
  if (dist(E, HP) > 0) return Ellipsoid();
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "E:\n" << E << std::endl;
  std::cout << "HP:\n" << HP.first << "; " << HP.second << std::endl;
#endif

  arma::vec e1(n);
  e1.zeros();
  e1(0)       = 1.;
  arma::mat T = E.align(e1, HP.first);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "T:\n" << T << std::endl;
  std::cout << "T*d:\n" << T * HP.first << std::endl;
#endif
  arma::vec f = (HP.second / (arma::dot(HP.first, HP.first))) * T * HP.first;
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "f: " << f << std::endl;
#endif
  Ellipsoid EM = mtimes(E, T, -f);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "EM:\n" << EM << std::endl;
  std::cout << "rank(EM._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "svd(EM._Q):\n" << EM.svdQ().first << std::endl;
#endif
  EM.regQ();
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "EM (reg):\n" << EM << std::endl;
  std::cout << "rank(EM._Q):\n" << EM.rankQ() << std::endl;
  std::cout << "svd(EM._Q):\n" << EM.svdQ().first << std::endl;
  std::cout << "inv(EM._Q):\n" << EM.invQ() << std::endl;
#endif

  arma::mat W(n - 1, n - 1, arma::fill::zeros);
  arma::vec w(n - 1);
  for (unsigned i = 1; i < n; i++)
  {
    for (unsigned j = i; j < n; j++) W(i - 1, j - 1) = EM.invQ()(i, j);
    w(i - 1) = EM.invQ()(i, 0);
  }
  W                = arma::symmatu(W);
  const double w11 = EM.invQ()(0, 0);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "W:\n" << W << std::endl;
  std::cout << "w:\n" << w << std::endl;
  std::cout << "w11:\n" << w11 << std::endl;
#endif
  arma::mat Winv(n - 1, n - 1, arma::fill::zeros);
  Ellipsoid::_inv(W, Winv);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "Winv:\n" << Winv << std::endl;
#endif

  arma::vec Winv_w(Winv * w);
  const double h = EM.c(0) * EM.c(0) * (w11 - arma::dot(w, Winv_w));
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "h: " << h << std::endl;
#endif

  arma::mat Z(n, n, arma::fill::zeros);
  arma::vec z(EM._c);
  z(0) = 0.;
  for (unsigned i = 1; i < n; i++)
  {
    for (unsigned j = i; j < n; j++) Z(i, j) = (1. - h) * Winv(i - 1, j - 1);
    z(i) += EM._c(0) * Winv_w(i - 1);
  }
  Z = arma::symmatu(Z);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "Z:\n" << Z << std::endl;
  std::cout << "z:\n" << z << std::endl;
#endif
  Ellipsoid I(Z, z);
#ifdef MC__ELLIPSOID_DEBUG_HPINTERSECTION
  std::cout << "I: " << I;
  std::cout << "eig(I): " << I.eigQ().first;
#endif
  I += f;
  return mtimes(I, T.t());
}

inline Ellipsoid
hpintersection(const Ellipsoid& E,
               const std::vector<std::pair<arma::vec, double> >& HP)
{
  Ellipsoid Ei(E);
  typename std::vector<std::pair<arma::vec, double> >::const_iterator HPi =
      HP.begin();
  for (; HPi != HP.end(); ++HPi) Ei = hpintersection(Ei, *HPi);
  return Ei;
}

inline Ellipsoid
ell_unitball(const unsigned n)
{
  Ellipsoid Eunit;
  Eunit._c.resize(n);
  Eunit._c.zeros();
  Eunit._Q.set_size(n, n);
  Eunit._Q.eye();
  return Eunit;
}

inline std::ostream&
operator<<(std::ostream& out, const Ellipsoid& E)
{
  if (!E.n()) return out;
  const int iprec = 5;
  out << std::scientific << std::setprecision(iprec);
  out << "\ncenter:\n" << E._c;
  out << "shape:\n" << E._Q;
  return out;
}

inline Ellipsoid
ellintersection_ia(const std::vector<std::pair<arma::vec, double> >& HP,
                   const double& TOL, const unsigned& MAXIT)
{
  if (!HP.size()) return Ellipsoid();
#ifdef MC__ELLIPSOID_USE_SDPA
  const unsigned m = HP.size();
  const unsigned n = HP.at(0).first.n_elem;
#ifdef MC__ELLIPSOID_DEBUG_ELLINTERSECTION_IA
  for (unsigned i = 0; i < m; i++)
  {
    std::cout << "Hyperplane #" << i << std::endl;
    std::cout << HP.at(i).first << HP.at(i).second << std::endl;
  }
#endif

  SDPA MAXDET;
#ifdef MC__ELLIPSOID_DEBUG_ELLINTERSECTION_IA
  SDPA::printSDPAVersion(stdout);
  MAXDET.setDisplay(stdout);
#endif

  // parameters
  MAXDET.setParameterMaxIteration(MAXIT);
  MAXDET.setParameterEpsilonStar(TOL);
  MAXDET.setParameterLambdaStar(1.0e+2);
  MAXDET.setParameterOmegaStar(2.0);
  MAXDET.setParameterLowerBound(-1.0e+5);
  MAXDET.setParameterUpperBound(1.0e+5);
  MAXDET.setParameterBetaStar(0.04);
  MAXDET.setParameterBetaBar(0.08);
  MAXDET.setParameterGammaStar(0.98);
  MAXDET.setParameterEpsilonDash(TOL);
  // MAXDET.setParameterPrintXVec((char*)"%+8.3e" );
  // MAXDET.setParameterPrintXMat((char*)"%+8.3e" );
  // MAXDET.setParameterPrintYMat((char*)"%+8.3e" );
  // MAXDET.setParameterPrintInformation((char*)"%+10.16e");
  // MAXDET.printParameters(stdout);

  // Intialize Problem 1: mDim, nblock, bLOCK_trust
  const unsigned l = std::ceil(std::log(n) / std::log(2.));
  unsigned sum     = 0.;
  for (unsigned i = 0; i < l; i++) sum += std::pow(2, i);
  const unsigned mDim = n * (n + 2) + 1 + sum, nblock = m + 4 + sum;
  MAXDET.inputConstraintNumber(mDim);
  MAXDET.inputBlockNumber(nblock);

  // Size blocks and block types (SDP-> Matrix),(LP->Inequality)
  for (unsigned i = 1; i <= m; i++)
  {
    MAXDET.inputBlockSize(i, n + 1);
    MAXDET.inputBlockType(i, SDPA::SDP);
  }
  MAXDET.inputBlockSize(m + 1, n);
  MAXDET.inputBlockType(m + 1, SDPA::SDP);
  MAXDET.inputBlockSize(m + 2, 2 * n);
  MAXDET.inputBlockType(m + 2, SDPA::SDP);
  for (unsigned i = m + 3; i <= nblock - 2; i++)
  {
    MAXDET.inputBlockSize(i, 2);
    MAXDET.inputBlockType(i, SDPA::SDP);
  }
  MAXDET.inputBlockSize(nblock - 1, -1);
  MAXDET.inputBlockType(nblock - 1, SDPA::LP);
  MAXDET.inputBlockSize(nblock, -1);
  MAXDET.inputBlockType(nblock, SDPA::LP);

  // Define coefficients in objective function
  MAXDET.initializeUpperTriangleSpace();
  for (unsigned i = 1; i <= mDim; i++)  // t is at index mDim
    MAXDET.inputCVec(i, i == mDim ? -1. : 0.);

  /* Assignment of non-zero entries in the upper triangle of Fi blocks
     non-assigned entries are taken as zero by default */

  // b and variable d (the center of ellipsoid) in the first m blocks
  for (unsigned i = 0; i <= n; i++)
    for (unsigned j = 1; j <= m; j++)
      for (unsigned k = 1; k <= n + 1; k++)
        MAXDET.inputElement(i, j, k, k,
                            i ? -HP[j - 1].first(i - 1) : -HP[j - 1].second);

  // matrix variables Z in first m blocks
  for (unsigned i = 1; i <= n; i++)
    for (unsigned j = 1; j <= i; j++)
      for (unsigned k = 1; k <= m; k++)
      {
        if (i == j)
          MAXDET.inputElement(i * (i + 1) / 2 + n, k, i, n + 1,
                              HP[k - 1].first(j - 1));
        else
        {
          MAXDET.inputElement(i * (i - 1) / 2 + j + n, k, i, n + 1,
                              HP[k - 1].first(j - 1));
          MAXDET.inputElement(i * (i - 1) / 2 + j + n, k, j, n + 1,
                              HP[k - 1].first(i - 1));
        }
      }

  // Z in block m+1
  for (unsigned i = 1; i <= n; i++)
    for (unsigned j = 1; j <= i; j++)
      MAXDET.inputElement(n + i * (i - 1) / 2 + j, m + 1, i, j, 1.);

  // D (as in the report) LMI (m+2) block
  for (unsigned i = 1; i <= n; i++)
    for (unsigned j = 1; j <= i; j++)
    {
      if (i == j)
      {
        MAXDET.inputElement(i * (i + 1) / 2 + n, m + 2, i, j, 1.);
        MAXDET.inputElement(n * (n + 3) / 2 + i * (i + 1) / 2, m + 2, j, i + n,
                            1.);
        MAXDET.inputElement(n * (n + 3) / 2 + i * (i + 1) / 2, m + 2, j + n,
                            i + n, 1.);
      }
      else
      {
        MAXDET.inputElement(i * (i - 1) / 2 + j + n, m + 2, i, j, 1.);
        MAXDET.inputElement(n * (n + 3) / 2 + i * (i - 1) / 2 + j, m + 2, i,
                            j + n, 1.);
      }
    }

  // S LMI (the rest of the blocks)
  const unsigned p = n / 2;
  const unsigned q = std::pow(2, l - 1);

  // diagonal entries of delta
  for (unsigned i = 1; i <= p; i++)
  {
    MAXDET.inputElement(n * (n + 3) / 2 + i * (2 * i - 1), m + 2 + i, 1, 1, 1.);
    MAXDET.inputElement(n * (n + 3) / 2 + i * (2 * i + 1), m + 2 + i, 2, 2, 1.);
  }
  unsigned odd = 0;
  if (n % 2)
  {
    MAXDET.inputElement(n * (n + 2), m + 3 + p, 1, 1, 1.);
    MAXDET.inputElement(mDim, m + 3 + p, 2, 2, 1.);
    odd = 1;
  }

  // t in the rest of q blocks
  for (unsigned i = p + 1 + odd; i <= q; i++)
    for (unsigned j = 1; j <= 2; j++)
      MAXDET.inputElement(mDim, m + 2 + i, j, j, 1.);

  // the adittional variables u
  for (unsigned i = 1; i <= sum; i++)
    MAXDET.inputElement(n * (n + 2) + i, m + 2 + i, 1, 2, 1.);
  for (unsigned i = 1; i <= sum - q; i++)
  {
    MAXDET.inputElement(n * (n + 2) + 2 * i - 1, m + 2 + q + i, 1, 1, 1.);
    MAXDET.inputElement(n * (n + 2) + 2 * i, m + 2 + q + i, 2, 2, 1.);
  }

  // u(sum) - t
  MAXDET.inputElement(n * (n + 2) + sum, nblock - 1, 1, 1, 1.);
  MAXDET.inputElement(mDim, nblock - 1, 1, 1, -1.);  // n*(n+2)+sum+1 = mDim
  MAXDET.inputElement(mDim, nblock, 1, 1, 1.);

  // solve SDP
  MAXDET.initializeUpperTriangle();
  MAXDET.initializeSolve();
  MAXDET.solve();

  // process results
  double* result = MAXDET.getResultXVec();
  if (MAXDET.getPhaseValue() != SDPA::pdOPT)
    throw std::runtime_error(
        "mc::ellintersection_ia\t SDPA failed to reach primal-dual optimum");

  arma::mat Z(n, n, arma::fill::zeros);
  arma::vec c(n);
  for (unsigned i = 0; i < n; i++)
  {
    c(i) = result[i];
    for (unsigned j = 0; j <= i; j++) Z(i, j) = result[n + i * (i + 1) / 2 + j];
  }
  Z           = arma::symmatl(Z);
  arma::mat Q = Z * Z;
  Q           = 0.5 * (Q + Q.t());
#ifdef MC__ELLIPSOID_DEBUG_ELLINTERSECTION_IA
  std::cout << "Inner ellipsoid:" << Ellipsoid(Q, c);
#endif

  return Ellipsoid(Q, c);

#else
  throw std::runtime_error("mc::ellintersection_ia\t SDPA solver not selected");
#endif
}

}  // namespace mc

#endif
