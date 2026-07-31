#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef MC__USE_PROFIL
#include "mcprofil.hpp"
typedef INTERVAL I;
#else
#ifdef MC__USE_FILIB
#include "mcfilib.hpp"
typedef filib::interval<double, filib::native_switched, filib::i_mode_extended>
    I;
#else
#ifdef MC__USE_BOOST
#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<
    boost::numeric::interval_lib::rounded_transc_opp<double> >
    T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round, T_boost_check>
    T_boost_policy;
typedef boost::numeric::interval<double, T_boost_policy> I;
#else
#include "interval.hpp"
typedef mc::Interval I;
#endif
#endif
#endif

#include "ellimage.hpp"

namespace py = pybind11;

namespace
{
using ELL = mc::Ellipsoid;
using EI  = mc::EllImg<I>;
using EV  = mc::EllVar<I>;
using HP  = std::pair<arma::vec, double>;

arma::vec
_to_arma_vec(py::handle obj)
{
  auto a =
      py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(
          obj);
  if (!a) throw py::type_error("expected a one-dimensional array-like object");
  py::buffer_info info = a.request();

  arma::uword n = 0;
  if (info.ndim == 1)
  {
    n = static_cast<arma::uword>(info.shape[0]);
  }
  else if (info.ndim == 2 && (info.shape[0] == 1 || info.shape[1] == 1))
  {
    n = static_cast<arma::uword>(info.shape[0] * info.shape[1]);
  }
  else
  {
    throw py::value_error("expected a vector with shape (n,), (n,1), or (1,n)");
  }

  arma::vec v(n);
  const double* data = static_cast<const double*>(info.ptr);
  for (arma::uword i = 0; i < n; ++i) v(i) = data[i];
  return v;
}

arma::mat
_to_arma_mat(py::handle obj)
{
  auto a =
      py::array_t<double, py::array::c_style | py::array::forcecast>::ensure(
          obj);
  if (!a) throw py::type_error("expected a two-dimensional array-like object");
  py::buffer_info info = a.request();

  if (info.ndim != 2)
    throw py::value_error("expected a matrix with shape (m,n)");
  const arma::uword nr = static_cast<arma::uword>(info.shape[0]);
  const arma::uword nc = static_cast<arma::uword>(info.shape[1]);

  arma::mat M(nr, nc);
  const double* data = static_cast<const double*>(info.ptr);
  for (arma::uword i = 0; i < nr; ++i)
    for (arma::uword j = 0; j < nc; ++j) M(i, j) = data[i * nc + j];
  return M;
}

py::array_t<double>
_from_arma_vec(arma::vec const& v)
{
  py::array_t<double> out({static_cast<py::ssize_t>(v.n_elem)});
  auto buf = out.mutable_unchecked<1>();
  for (arma::uword i = 0; i < v.n_elem; ++i)
    buf(static_cast<py::ssize_t>(i)) = v(i);
  return out;
}

py::array_t<double>
_from_arma_mat(arma::mat const& M)
{
  py::array_t<double> out(
      {static_cast<py::ssize_t>(M.n_rows), static_cast<py::ssize_t>(M.n_cols)});
  auto buf = out.mutable_unchecked<2>();
  for (arma::uword i = 0; i < M.n_rows; ++i)
    for (arma::uword j = 0; j < M.n_cols; ++j)
      buf(static_cast<py::ssize_t>(i), static_cast<py::ssize_t>(j)) = M(i, j);
  return out;
}

py::array_t<double>
_from_arma_spmat(arma::sp_mat const& S)
{
  return _from_arma_mat(arma::mat(S));
}

HP
_to_hp(py::handle obj)
{
  py::tuple t = py::reinterpret_borrow<py::tuple>(obj);
  if (!py::isinstance<py::tuple>(obj) || t.size() != 2)
    throw py::type_error("expected a halfspace/hyperplane as (normal, rhs)");
  return HP(_to_arma_vec(t[0]), py::cast<double>(t[1]));
}

std::vector<HP>
_to_hps(py::iterable obj)
{
  std::vector<HP> out;
  for (py::handle item : obj) out.push_back(_to_hp(item));
  return out;
}

std::vector<EV>
_to_ev_vector(py::iterable obj)
{
  std::vector<EV> out;
  for (py::handle item : obj) out.push_back(py::cast<EV>(item));
  return out;
}

std::vector<arma::vec>
_to_vec_vector(py::iterable obj)
{
  std::vector<arma::vec> out;
  for (py::handle item : obj) out.push_back(_to_arma_vec(item));
  return out;
}

std::string
_to_string(ELL const& E)
{
  std::ostringstream os;
  os << E;
  return os.str();
}

std::string
_to_string(EI const& E)
{
  std::ostringstream os;
  os << E;
  return os.str();
}

std::string
_to_string(EV const& V)
{
  std::ostringstream os;
  os << V;
  return os.str();
}
}  // namespace

void
mc_ellimage(py::module& m)
{
  // Surface the C++ ellipsoid/ellimage exceptions -- which are NOT derived from
  // std::exception -- as informative Python exceptions, instead of the opaque
  // "Caught an unknown exception!" that pybind's default handler would raise.
  static py::exception<ELL::Exceptions> pyEllipsoidError(m, "EllipsoidError");
  static py::exception<EI::Exceptions> pyEllImgError(m, "EllImgError");
  py::register_exception_translator(
      [](std::exception_ptr p)
      {
        try
        {
          if (p) std::rethrow_exception(p);
        }
        catch (EI::Exceptions& e)
        {
          py::set_error(pyEllImgError, e.what().c_str());
        }
        catch (ELL::Exceptions& e)
        {
          py::set_error(pyEllipsoidError, e.what().c_str());
        }
      });

  py::class_<ELL> pyEllipsoid(m, "Ellipsoid", R"doc(
Ellipsoid in R^n, defined by a center and a shape matrix.

The ellipsoid with center ``c`` (length-n vector) and symmetric
positive-semidefinite shape matrix ``Q`` (n-by-n) is the set::

    E(Q, c) = { c + sqrtm(Q) @ v : ||v||_2 <= 1 }

Vectors are passed as one-dimensional array-likes (lists or numpy
arrays); matrices as two-dimensional array-likes of shape (n, n),
e.g. nested lists or numpy arrays in row-major layout. Matrices are
symmetrized as ``(Q + Q.T)/2`` on input. The accessors `c` and `Q`
return numpy arrays (copies).

Ellipsoidal calculus is provided by module-level functions:
`mtimes` (affine map), `minksum_ea` / `minksum_box` /
`minksum_interval` (external approximations of Minkowski sums),
`intersection_ea`, `hpintersection`, `ellintersection_ia`, `dist`,
`inv` and `ell_unitball`. The in-place operators ``+=`` and ``-=``
with a vector translate the center.

Failures raise `pymcpp.EllipsoidError`: a non-PSD shape matrix (only
checked when ``Ellipsoid.options.PSDCHK`` is enabled), a LAPACK
failure, or a root-finding failure.

Examples
--------
>>> import pymcpp
>>> E = pymcpp.Ellipsoid([[5., 4.], [4., 5.]], [3., 4.])
>>> E.n
2
>>> E.trQ()
10.0
>>> E.c
array([3., 4.])
)doc");

  pyEllipsoid.def(py::init<>(), R"doc(
Construct an empty (0-dimensional) ellipsoid.
)doc")
      .def(py::init(
               [](py::handle Q, py::object c) {
                 return ELL(_to_arma_mat(Q),
                            c.is_none() ? arma::vec() : _to_arma_vec(c));
               }),
           py::arg("Q"), py::arg("c") = py::none(), R"doc(
Construct an ellipsoid from a dense shape matrix and optional center.

Parameters
----------
Q : array-like of shape (n, n)
    Shape matrix; symmetrized as (Q + Q.T)/2 internally. Raises
    EllipsoidError if not positive semidefinite and
    ``Ellipsoid.options.PSDCHK`` is True.
c : array-like of shape (n,), optional
    Center vector. Default is the origin.
)doc")
      .def(
          py::init(
              [](unsigned const n, std::vector<double> const& Qlt, py::object c)
              {
                if (Qlt.size() != static_cast<size_t>(n * (n + 1) / 2))
                  throw py::value_error(
                      "Qlt must contain n*(n+1)/2 lower-triangular entries");
                std::vector<double> cv;
                const double* cp = nullptr;
                if (!c.is_none())
                {
                  arma::vec cav = _to_arma_vec(c);
                  if (cav.n_elem != n)
                    throw py::value_error("center c must have length n");
                  cv.resize(n);
                  for (unsigned i = 0; i < n; ++i) cv[i] = cav(i);
                  cp = cv.data();
                }
                return ELL(n, Qlt.data(), cp);
              }),
          py::arg("n"), py::arg("Qlt"), py::arg("c") = py::none(), R"doc(
Construct an ellipsoid from the lower triangle of its shape matrix.

Parameters
----------
n : int
    Dimension of the ellipsoid.
Qlt : list of float
    The n*(n+1)/2 lower-triangular entries of the shape matrix,
    stored column-wise: [Q00, Q10, ..., Q(n-1)0, Q11, Q21, ...].
c : array-like of shape (n,), optional
    Center vector. Default is the origin.

Raises
------
ValueError
    If ``Qlt`` does not have n*(n+1)/2 entries, or ``c`` does not
    have length n.
)doc")
      .def_static(
          "from_radius",
          [](py::handle r, py::object c) {
            return ELL(_to_arma_vec(r),
                       c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("r"), py::arg("c") = py::none(), R"doc(
Construct an axis-aligned ellipsoid enclosing the box [c-r, c+r].

The diagonal shape matrix is chosen as ``Q[i,i] = r[i]*sum(r)``
(minimum-trace axis-aligned enclosure); the box corners lie on the
boundary of the ellipsoid.

Parameters
----------
r : array-like of shape (n,)
    Nonnegative box half-widths (radii), one per coordinate.
c : array-like of shape (n,), optional
    Box midpoint. Default is the origin.

Returns
-------
Ellipsoid
    A new axis-aligned ellipsoid enclosing the box.
)doc")
      .def_static(
          "from_box", [](std::vector<I> const& box)
          { return ELL(static_cast<unsigned>(box.size()), box.data()); },
          py::arg("box"), R"doc(
Construct an ellipsoid enclosing an interval box.

Equivalent to `from_radius` with the interval midpoints as center
and half-diameters as radii.

Parameters
----------
box : list of Interval
    Interval bounds, one per coordinate.

Returns
-------
Ellipsoid
    A new axis-aligned ellipsoid enclosing the interval box.
)doc")
      .def(py::init<ELL const&>(), R"doc(
Copy constructor.
)doc")
      .def_readwrite_static("options", &ELL::options, R"doc(
Ellipsoid.Options: Class-wide options for ellipsoidal calculus,
shared by all Ellipsoid instances.
)doc")
      .def("unitball", &ELL::unitball, py::arg("n"),
           py::return_value_policy::reference_internal, R"doc(
Reset this ellipsoid to the n-dimensional unit ball.

Sets the shape matrix to the identity and the center to the origin.

Parameters
----------
n : int
    Dimension.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "set",
          [](ELL& self, py::handle Q, py::object c) -> ELL&
          {
            return self.set(_to_arma_mat(Q),
                            c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("Q"), py::arg("c") = py::none(),
          py::return_value_policy::reference_internal, R"doc(
Redefine the ellipsoid from a dense shape matrix and optional center.

Same semantics as the ``Ellipsoid(Q, c)`` constructor; cached
decompositions are cleared.

Parameters
----------
Q : array-like of shape (n, n)
    Shape matrix; symmetrized internally.
c : array-like of shape (n,), optional
    Center vector. Default is the origin.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "set_lower",
          [](ELL& self, unsigned const n, std::vector<double> const& Qlt,
             py::object c) -> ELL&
          {
            if (Qlt.size() != static_cast<size_t>(n * (n + 1) / 2))
              throw py::value_error(
                  "Qlt must contain n*(n+1)/2 lower-triangular entries");
            std::vector<double> cv;
            const double* cp = nullptr;
            if (!c.is_none())
            {
              arma::vec cav = _to_arma_vec(c);
              if (cav.n_elem != n)
                throw py::value_error("center c must have length n");
              cv.resize(n);
              for (unsigned i = 0; i < n; ++i) cv[i] = cav(i);
              cp = cv.data();
            }
            return self.set(n, Qlt.data(), cp);
          },
          py::arg("n"), py::arg("Qlt"), py::arg("c") = py::none(),
          py::return_value_policy::reference_internal, R"doc(
Redefine the ellipsoid from the lower triangle of its shape matrix.

Same semantics as the ``Ellipsoid(n, Qlt, c)`` constructor.

Parameters
----------
n : int
    Dimension of the ellipsoid.
Qlt : list of float
    The n*(n+1)/2 lower-triangular shape entries, column-wise.
c : array-like of shape (n,), optional
    Center vector. Default is the origin.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "set_radius",
          [](ELL& self, py::handle r, py::object c) -> ELL&
          {
            return self.set(_to_arma_vec(r),
                            c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("r"), py::arg("c") = py::none(),
          py::return_value_policy::reference_internal, R"doc(
Redefine as an axis-aligned ellipsoid enclosing the box [c-r, c+r].

Same semantics as `from_radius`.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "set_box", [](ELL& self, std::vector<I> const& box) -> ELL&
          { return self.set(static_cast<unsigned>(box.size()), box.data()); },
          py::arg("box"), py::return_value_policy::reference_internal, R"doc(
Redefine as an axis-aligned ellipsoid enclosing an interval box.

Same semantics as `from_box`.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def("reset", &ELL::reset, py::return_value_policy::reference_internal,
           R"doc(
Reset to an empty (0-dimensional) ellipsoid.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def("reset_aux", &ELL::reset_aux,
           py::return_value_policy::reference_internal, R"doc(
Clear the cached decompositions of the shape matrix.

The eigendecomposition, square root, SVD and inverse of the shape
matrix are computed lazily and cached; call this after modifying the
shape matrix through means the class cannot track.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "extend", [](ELL& self, py::handle Qi, double const& ci) -> ELL&
          { return self.extend(arma::rowvec(_to_arma_vec(Qi).t()), ci); },
          py::arg("Qi"), py::arg("ci") = 0.,
          py::return_value_policy::reference_internal, R"doc(
Append one dimension to the ellipsoid.

The shape matrix grows from (n, n) to (n+1, n+1) by appending ``Qi``
as its new last row (and, by symmetry, last column), and ``ci`` is
appended to the center.

Parameters
----------
Qi : array-like of shape (n+1,)
    New last row of the shape matrix: the first n entries are the
    covariances with the existing dimensions, the last entry is the
    new diagonal element.
ci : float, optional
    New center entry. Default is 0.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def("O", &ELL::O, R"doc(
Return a new ellipsoid with the same shape matrix centered at the
origin.

Returns
-------
Ellipsoid
    Centered copy; this ellipsoid is not modified.
)doc")
      .def_property_readonly("n", &ELL::n, R"doc(
int: Dimension of the ellipsoid.
)doc")
      .def_property_readonly(
          "c", [](ELL const& self) { return _from_arma_vec(self.c()); },
          R"doc(
numpy.ndarray of shape (n,): Center vector (a copy).
)doc")
      .def_property_readonly(
          "Q", [](ELL const& self) { return _from_arma_mat(self.Q()); },
          R"doc(
numpy.ndarray of shape (n, n): Shape matrix (a copy).
)doc")
      .def("ci", py::overload_cast<unsigned>(&ELL::c, py::const_), py::arg("i"),
           R"doc(
Return center coefficient ``c[i]``.

Parameters
----------
i : int
    Coordinate index in {0, ..., n-1}.

Returns
-------
float
    The center entry.
)doc")
      .def("Qij", py::overload_cast<unsigned, unsigned>(&ELL::Q, py::const_),
           py::arg("i"), py::arg("j"), R"doc(
Return shape matrix coefficient ``Q[i,j]``.

Parameters
----------
i : int
    Row index in {0, ..., n-1}.
j : int
    Column index in {0, ..., n-1}.

Returns
-------
float
    The shape matrix entry.
)doc")
      .def(
          "set_ci",
          [](ELL& self, unsigned const i, double const& v) -> ELL&
          {
            if (i >= self.n())
              throw py::index_error("center index out of range");
            self.c()(i) = v;
            return self;
          },
          py::arg("i"), py::arg("value"),
          py::return_value_policy::reference_internal, R"doc(
Set center coefficient ``c[i]``.

Parameters
----------
i : int
    Coordinate index in {0, ..., n-1}.
value : float
    New center entry.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).

Raises
------
IndexError
    If ``i`` is out of range.
)doc")
      .def(
          "set_Qij",
          [](ELL& self, unsigned const i, unsigned const j,
             double const& v) -> ELL&
          {
            if (i >= self.n() || j >= self.n())
              throw py::index_error("shape index out of range");
            self.Q(i, j) = v;
            self.Q(j, i) = v;
            self.reset_aux();
            return self;
          },
          py::arg("i"), py::arg("j"), py::arg("value"),
          py::return_value_policy::reference_internal, R"doc(
Set shape matrix coefficient ``Q[i,j]`` (and ``Q[j,i]``, keeping the
matrix symmetric), and clear the cached decompositions.

Parameters
----------
i : int
    Row index in {0, ..., n-1}.
j : int
    Column index in {0, ..., n-1}.
value : float
    New shape matrix entry.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).

Raises
------
IndexError
    If ``i`` or ``j`` is out of range.
)doc")
      .def("trQ", &ELL::trQ, R"doc(
Return the trace of the shape matrix.

Returns
-------
float
    Sum of the diagonal entries of Q (0 for an empty ellipsoid).
)doc")
      .def("psdQ", &ELL::psdQ, R"doc(
Test whether the shape matrix is positive semidefinite.

The test allows eigenvalues down to ``-Ellipsoid.options.PSDTOL``.

Returns
-------
bool
    True if the smallest eigenvalue of Q is above the tolerance.
)doc")
      .def("rankQ", &ELL::rankQ, R"doc(
Return the numerical rank of the shape matrix.

Singular values below ``Ellipsoid.options.RKTOLA`` (absolute) or
``Ellipsoid.options.RKTOLR`` times the largest singular value
(relative) are treated as zero.

Returns
-------
int
    Numerical rank of Q.
)doc")
      .def(
          "eigQ",
          [](ELL& self)
          {
            auto const& e = self.eigQ();
            return py::make_tuple(_from_arma_vec(e.first),
                                  _from_arma_mat(e.second));
          },
          R"doc(
Return the eigendecomposition of the shape matrix.

The decomposition is cached; repeated calls are cheap.

Returns
-------
eigenvalues : numpy.ndarray of shape (n,)
    Eigenvalues of Q in ascending order.
eigenvectors : numpy.ndarray of shape (n, n)
    Matrix whose column k is the eigenvector for eigenvalues[k].

Raises
------
EllipsoidError
    If the LAPACK eigensolver fails.
)doc")
      .def(
          "sqrtQ", [](ELL& self, bool const complete)
          { return _from_arma_mat(self.sqrtQ(complete)); },
          py::arg("complete") = false, R"doc(
Return a matrix square root S of the shape matrix, with S @ S.T = Q.

Computed from the eigendecomposition (negative eigenvalues are
truncated at zero) and cached.

Parameters
----------
complete : bool, optional
    If True, symmetrize the cached square root as (S + S.T)/2 before
    returning it. Default is False.

Returns
-------
numpy.ndarray of shape (n, n)
    Matrix square root of Q.
)doc")
      .def(
          "svdQ",
          [](ELL& self)
          {
            auto const& s = self.svdQ();
            return py::make_tuple(_from_arma_vec(s.first),
                                  _from_arma_mat(s.second.first),
                                  _from_arma_mat(s.second.second));
          },
          R"doc(
Return the singular value decomposition of the shape matrix.

The decomposition Q = U @ diag(s) @ Vt is cached.

Returns
-------
s : numpy.ndarray of shape (n,)
    Singular values in descending order.
U : numpy.ndarray of shape (n, n)
    Left singular vectors (as columns).
Vt : numpy.ndarray of shape (n, n)
    Transpose of the matrix of right singular vectors.

Raises
------
EllipsoidError
    If the LAPACK SVD routine fails.
)doc")
      .def(
          "regQ", [](ELL& self) { return _from_arma_mat(self.regQ()); },
          R"doc(
Regularize the shape matrix in place and return it.

If Q is numerically rank-deficient (see `rankQ`), a small multiple
of the identity is added on its null space so that it becomes
invertible; otherwise Q is returned unchanged. Cached decompositions
are cleared when Q is modified.

Returns
-------
numpy.ndarray of shape (n, n)
    The (possibly regularized) shape matrix.
)doc")
      .def(
          "invQ", [](ELL& self) { return _from_arma_mat(self.invQ()); },
          R"doc(
Return the inverse (or pseudo-inverse) of the shape matrix.

The result is cached. For a singular Q the Moore-Penrose
pseudo-inverse is returned.

Returns
-------
numpy.ndarray of shape (n, n)
    Inverse or pseudo-inverse of Q.
)doc")
      .def(
          "invQij", [](ELL& self, unsigned const i, unsigned const j)
          { return self.invQ(i, j); }, py::arg("i"), py::arg("j"),
          R"doc(
Return coefficient (i, j) of the inverse shape matrix.

Parameters
----------
i : int
    Row index in {0, ..., n-1}.
j : int
    Column index in {0, ..., n-1}.

Returns
-------
float
    Entry of the (pseudo-)inverse of Q.
)doc")
      .def(
          "align",
          [](ELL const& self, py::handle v, py::handle x) {
            return _from_arma_mat(self.align(_to_arma_vec(v), _to_arma_vec(x)));
          },
          py::arg("v"), py::arg("x"), R"doc(
Return an orthogonal matrix T such that T @ x is parallel to v.

Parameters
----------
v : array-like of shape (n,)
    Target direction.
x : array-like of shape (n,)
    Vector to rotate.

Returns
-------
numpy.ndarray of shape (n, n)
    Orthogonal rotation matrix (empty array if the lengths differ).
)doc")
      .def("l", &ELL::l, py::arg("i"), R"doc(
Return the lower bound of the ellipsoid along coordinate i.

Equals ``c[i] - sqrt(Q[i,i])``.

Parameters
----------
i : int
    Coordinate index in {0, ..., n-1}.

Returns
-------
float
    Coordinate lower bound.
)doc")
      .def("u", &ELL::u, py::arg("i"), R"doc(
Return the upper bound of the ellipsoid along coordinate i.

Equals ``c[i] + sqrt(Q[i,i])``.

Parameters
----------
i : int
    Coordinate index in {0, ..., n-1}.

Returns
-------
float
    Coordinate upper bound.
)doc")
      .def("r", &ELL::r, py::arg("i"), R"doc(
Return the radius of the ellipsoid along coordinate i.

Equals ``sqrt(Q[i,i])``.

Parameters
----------
i : int
    Coordinate index in {0, ..., n-1}.

Returns
-------
float
    Coordinate radius.
)doc")
      .def(
          "translate",
          [](ELL& self, py::handle d) -> ELL&
          {
            self += _to_arma_vec(d);
            return self;
          },
          py::arg("d"), py::return_value_policy::reference_internal, R"doc(
Translate the center by the vector d (in place).

Equivalent to ``E += d``.

Parameters
----------
d : array-like of shape (n,)
    Translation vector.

Returns
-------
Ellipsoid
    This ellipsoid (for chaining).
)doc")
      .def(
          "__iadd__",
          [](ELL& self, py::handle d) -> ELL&
          {
            self += _to_arma_vec(d);
            return self;
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def(
          "__isub__",
          [](ELL& self, py::handle d) -> ELL&
          {
            self -= _to_arma_vec(d);
            return self;
          },
          py::is_operator(), py::return_value_policy::reference_internal)
      .def("__str__", static_cast<std::string (*)(ELL const&)>(&_to_string))
      .def("__repr__", static_cast<std::string (*)(ELL const&)>(&_to_string));

  // Nested Ellipsoid::Options
  py::class_<ELL::Options> pyEllipsoidOptions(pyEllipsoid, "Options", R"doc(
Options for ellipsoidal calculus with `Ellipsoid`.

The active option set is the class attribute ``Ellipsoid.options``,
shared by all instances.
)doc");
  pyEllipsoidOptions.def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def(py::init<ELL::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](ELL::Options& self) { self = ELL::Options(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("PSDCHK", &ELL::Options::PSDCHK, R"doc(
bool: Whether to check positive semidefiniteness of shape matrices
on construction and decomposition, raising EllipsoidError on
failure. Default is False.
)doc")
      .def_readwrite(
          "PSDTOL", &ELL::Options::PSDTOL, R"doc(
float: Absolute tolerance for the positive-semidefiniteness check
(most negative eigenvalue allowed). Default is 100 times the machine
precision (about 2.2e-14).
)doc")
      .def_readwrite(
          "RKTOLA", &ELL::Options::RKTOLA, R"doc(
float: Absolute tolerance below which singular values are treated as
zero in `Ellipsoid.rankQ` and `Ellipsoid.regQ`. Default is the
machine precision (about 2.2e-16).
)doc")
      .def_readwrite(
          "RKTOLR", &ELL::Options::RKTOLR, R"doc(
float: Relative tolerance (times the largest singular value) below
which singular values are treated as zero in `Ellipsoid.rankQ` and
`Ellipsoid.regQ`. Default is the machine precision (about 2.2e-16).
)doc")
      .def_readwrite(
          "ROOTTOL", &ELL::Options::ROOTTOL, R"doc(
float: Absolute stopping tolerance for the root-finding method used
in ellipsoid intersection computations. Default is 1e-10.
)doc")
      .def_readwrite("ROOTSECANT", &ELL::Options::ROOTSECANT, R"doc(
bool: Whether to try the secant method for root finding before
falling back to golden-section search. Default is False.
)doc")
      .def_readwrite(
          "ROOTMAXIT", &ELL::Options::ROOTMAXIT, R"doc(
int: Maximum number of root-finding iterations; 0 means no maximum.
Default is 0.
)doc");

  // Nested Ellipsoid::Exceptions
  py::class_<ELL::Exceptions> pyEllipsoidExceptions(pyEllipsoid, "Exceptions",
                                                    R"doc(
Error payload thrown by `Ellipsoid` operations.

In Python these errors surface as `pymcpp.EllipsoidError`; this
class only serves to inspect or construct the underlying error
codes.
)doc");
  py::enum_<ELL::Exceptions::TYPE>(pyEllipsoidExceptions, "TYPE", R"doc(
Error codes of `Ellipsoid` operations.
)doc")
      .value("NONPSD", ELL::Exceptions::TYPE::NONPSD,
             "non-positive-semidefinite shape matrix")
      .value("LAPACK", ELL::Exceptions::TYPE::LAPACK,
             "linear algebra routine failed")
      .value("ROOT", ELL::Exceptions::TYPE::ROOT, "root-finding routine failed")
      .export_values();

  pyEllipsoidExceptions.def(py::init<ELL::Exceptions::TYPE>(), R"doc(
Construct an exception payload with the given error code.
)doc")
      .def("ierr", &ELL::Exceptions::ierr, R"doc(
Return the numeric error flag (see `Ellipsoid.Exceptions.TYPE`).
)doc")
      .def("what", &ELL::Exceptions::what, R"doc(
Return the human-readable error description.
)doc");

  py::class_<EI, ELL> pyEllImg(m, "EllImg", R"doc(
Ellipsoidal image: propagate an ellipsoid through a factorable
function.

`EllImg` derives from `Ellipsoid` and computes an ellipsoidal
enclosure of the image ``{ f(x) : x in E(Qx, cx) }`` of an ellipsoid
under a factorable function ``f`` (the exact image is generally not
an ellipsoid). Internally, the ellipsoid is *lifted*: each
elementary operation of ``f`` appends one dimension to the
ellipsoid. Nonlinear univariate operations are enclosed by a linear
approximation (a degree-1 Remez minimax line by default, or a
secant line) plus an interval remainder, which is added to the
lifted ellipsoid as a Minkowski sum; bilinear products are enclosed
via a midpoint linearization or, optionally, a difference-of-convex
decomposition (see `EllImg.Options`).

Workflow: construct the `EllImg` from the host ellipsoid of the
independent variables, attach one `EllVar` per coordinate, evaluate
the function on these variables using overloaded arithmetic and the
`pymcpp` math functions, then project the lifted ellipsoid back
onto the dependents with `get`, which returns a plain `Ellipsoid`.

The propagation is not verified: rounding errors are not accounted
for. Domain violations (e.g. ``log`` of a variable whose range
contains non-positive values) raise `pymcpp.EllImgError`.

Examples
--------
>>> import pymcpp
>>> Ex = pymcpp.EllImg([[5., 4.], [4., 5.]], [3., 4.])
>>> X1 = pymcpp.EllVar(Ex, 0)
>>> X2 = pymcpp.EllVar(Ex, 1)
>>> F = [pymcpp.log(X1) + pymcpp.sqr(X2),
...      pymcpp.sin(X1) - pymcpp.cos(X2)]
>>> Ex.qdim  # 2 variables + 6 lifted operations
8
>>> Ef = Ex.get(F)  # plain Ellipsoid enclosing the image of F
>>> Ef.n
2
)doc");

  pyEllImg.def(py::init<>(), R"doc(
Construct an empty (0-dimensional) ellipsoidal image.
)doc")
      .def(py::init(
               [](py::handle Q, py::object c) {
                 return EI(_to_arma_mat(Q),
                           c.is_none() ? arma::vec() : _to_arma_vec(c));
               }),
           py::arg("Q"), py::arg("c") = py::none(), R"doc(
Construct an ellipsoidal image from a dense shape matrix and
optional center.

Parameters
----------
Q : array-like of shape (n, n)
    Shape matrix of the host ellipsoid; symmetrized internally.
c : array-like of shape (n,), optional
    Center vector. Default is the origin.
)doc")
      .def(
          py::init(
              [](unsigned const n, std::vector<double> const& Qlt, py::object c)
              {
                if (Qlt.size() != static_cast<size_t>(n * (n + 1) / 2))
                  throw py::value_error(
                      "Qlt must contain n*(n+1)/2 lower-triangular entries");
                std::vector<double> cv;
                const double* cp = nullptr;
                if (!c.is_none())
                {
                  arma::vec cav = _to_arma_vec(c);
                  if (cav.n_elem != n)
                    throw py::value_error("center c must have length n");
                  cv.resize(n);
                  for (unsigned i = 0; i < n; ++i) cv[i] = cav(i);
                  cp = cv.data();
                }
                return EI(n, Qlt.data(), cp);
              }),
          py::arg("n"), py::arg("Qlt"), py::arg("c") = py::none(), R"doc(
Construct an ellipsoidal image from the lower triangle of its shape
matrix.

Parameters
----------
n : int
    Dimension of the host ellipsoid.
Qlt : list of float
    The n*(n+1)/2 lower-triangular shape entries, stored
    column-wise: [Q00, Q10, ..., Q(n-1)0, Q11, Q21, ...].
c : array-like of shape (n,), optional
    Center vector. Default is the origin.

Raises
------
ValueError
    If ``Qlt`` does not have n*(n+1)/2 entries, or ``c`` does not
    have length n.
)doc")
      .def_static(
          "from_radius",
          [](py::handle r, py::object c) {
            return EI(_to_arma_vec(r),
                      c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("r"), py::arg("c") = py::none(), R"doc(
Construct an axis-aligned ellipsoidal image enclosing the box
[c-r, c+r].

Same enclosure as `Ellipsoid.from_radius`.

Parameters
----------
r : array-like of shape (n,)
    Nonnegative box half-widths (radii), one per coordinate.
c : array-like of shape (n,), optional
    Box midpoint. Default is the origin.

Returns
-------
EllImg
    A new ellipsoidal image.
)doc")
      .def(py::init<EI const&>(), R"doc(
Copy constructor.
)doc")
      .def_readwrite_static("options", &EI::options, R"doc(
EllImg.Options: Class-wide options for ellipsoidal image
propagation, shared by all EllImg instances. Distinct from
``Ellipsoid.options``.
)doc")
      .def_property_readonly(
          "Q_lift", [](EI& self) { return _from_arma_spmat(self.Q_lift()); },
          R"doc(
numpy.ndarray of shape (qdim, qdim): Dense copy of the lifted shape
matrix, covering the original variables and one row/column per
lifted operation.
)doc")
      .def_property_readonly(
          "c_lift", [](EI& self) { return _from_arma_vec(self.c_lift()); },
          R"doc(
numpy.ndarray of shape (qdim,): Center of the lifted ellipsoid.
)doc")
      .def_property_readonly(
          "qdim",
          [](EI& self) { return static_cast<long>(self.c_lift().n_elem); },
          R"doc(
int: Current dimension of the lifted ellipsoid (number of original
variables plus lifted operations).
)doc")
      .def(
          "set",
          [](EI& self, py::handle Q, py::object c) -> EI&
          {
            return self.set(_to_arma_mat(Q),
                            c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("Q"), py::arg("c") = py::none(),
          py::return_value_policy::reference_internal, R"doc(
Redefine the host ellipsoid from a dense shape matrix and optional
center, and reset the lifting.

Parameters
----------
Q : array-like of shape (n, n)
    Shape matrix; symmetrized internally.
c : array-like of shape (n,), optional
    Center vector. Default is the origin.

Returns
-------
EllImg
    This image (for chaining).
)doc")
      .def(
          "set_radius",
          [](EI& self, py::handle r, py::object c) -> EI&
          {
            return self.set(_to_arma_vec(r),
                            c.is_none() ? arma::vec() : _to_arma_vec(c));
          },
          py::arg("r"), py::arg("c") = py::none(),
          py::return_value_policy::reference_internal, R"doc(
Redefine as an axis-aligned image enclosing the box [c-r, c+r], and
reset the lifting.

Returns
-------
EllImg
    This image (for chaining).
)doc")
      .def(
          "reset", [](EI& self) -> EI& { return self.reset(); },
          py::return_value_policy::reference_internal, R"doc(
Discard all lifted dimensions, resetting the image to the host
ellipsoid of the original variables.

Returns
-------
EllImg
    This image (for chaining).
)doc")
      .def(
          "get",
          [](EI& self, py::iterable vars)
          {
            std::vector<EV> v = _to_ev_vector(vars);
            return self.get(static_cast<unsigned>(v.size()), v.data());
          },
          py::arg("vars"), R"doc(
Project the lifted ellipsoid onto the given variables.

Extracts the sub-ellipsoid corresponding to the rows of the lifted
center/shape occupied by ``vars``, typically the dependent
expressions computed in `EllVar` arithmetic.

Parameters
----------
vars : iterable of EllVar
    Variables/expressions to project onto; all must belong to this
    image.

Returns
-------
Ellipsoid
    Plain ellipsoid of dimension ``len(vars)`` enclosing the joint
    range of the given variables.
)doc")
      .def("__str__", static_cast<std::string (*)(EI const&)>(&_to_string))
      .def("__repr__", static_cast<std::string (*)(EI const&)>(&_to_string));

  // Nested EllImg::Options
  py::class_<EI::Options> pyEllImgOptions(pyEllImg, "Options", R"doc(
Options for ellipsoidal image propagation with `EllImg`.

The active option set is the class attribute ``EllImg.options``,
shared by all instances (and distinct from ``Ellipsoid.options``).
)doc");
  pyEllImgOptions.def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def(py::init<EI::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](EI::Options& self) { self = EI::Options(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("PREALLOC", &EI::Options::PREALLOC, R"doc(
int: Number of rows to preallocate in the lifted center when
lifting; 0 uses geometric (doubling) growth. Default is 0.
)doc")
      .def_readwrite("MINK_TOL", &EI::Options::MINK_TOL, R"doc(
float: Tolerance used in the Minkowski-sum step that adds the
linearization remainder to the lifted ellipsoid. Default is 1e-10.
)doc")
      .def_readwrite("REMEZ_USE", &EI::Options::REMEZ_USE, R"doc(
bool: Whether to linearize nonlinear univariate operations with a
degree-1 Remez minimax line; if False, a secant-line approximation
is used instead. Default is True.
)doc")
      .def_readwrite("REMEZ_MAXIT", &EI::Options::REMEZ_MAXIT, R"doc(
int: Maximum number of iterations of the Remez exchange algorithm.
Default is 5.
)doc")
      .def_readwrite("REMEZ_TOL", &EI::Options::REMEZ_TOL, R"doc(
float: Stopping tolerance of the Remez exchange algorithm.
Default is 1e-5.
)doc")
      .def_readwrite(
          "REMEZ_MIG", &EI::Options::REMEZ_MIG, R"doc(
float: Minimum domain diameter for invoking the Remez algorithm;
smaller domains fall back to the secant-line approximation.
Default is 1e-10.
)doc")
      .def_readwrite("DCPROD_USE", &EI::Options::DCPROD_USE, R"doc(
bool: Whether to lift bilinear products via a difference-of-convex
decomposition instead of the midpoint linearization.
Default is False.
)doc");

  // Nested EllImg::Exceptions
  py::class_<EI::Exceptions> pyEllImgExceptions(pyEllImg, "Exceptions", R"doc(
Error payload thrown during ellipsoidal image propagation.

In Python these errors surface as `pymcpp.EllImgError`; this class
only serves to inspect or construct the underlying error codes.
)doc");
  py::enum_<EI::Exceptions::TYPE>(pyEllImgExceptions, "TYPE", R"doc(
Error codes of ellipsoidal image propagation.
)doc")
      .value("DIV", EI::Exceptions::TYPE::DIV, "division by zero scalar")
      .value("INV", EI::Exceptions::TYPE::INV,
             "inverse operation with zero in domain")
      .value("LOG", EI::Exceptions::TYPE::LOG,
             "log operation with non-positive numbers in domain")
      .value("SQRT", EI::Exceptions::TYPE::SQRT,
             "square-root operation with negative numbers in domain")
      .value("TAN", EI::Exceptions::TYPE::TAN,
             "tangent operation with zero in cosine domain")
      .value("ACOS", EI::Exceptions::TYPE::ACOS,
             "inverse cosine operation with domain outside [-1,1]")
      .value("ASIN", EI::Exceptions::TYPE::ASIN,
             "inverse sine operation with domain outside [-1,1]")
      .value("INIT", EI::Exceptions::TYPE::INIT,
             "failed to construct ellipsoidal variable")
      .value("EIMG", EI::Exceptions::TYPE::EIMG,
             "variables belong to different ellipsoidal images")
      .value("UNDEF", EI::Exceptions::TYPE::UNDEF,
             "feature not yet implemented")
      .export_values();

  pyEllImgExceptions.def(py::init<EI::Exceptions::TYPE>(), R"doc(
Construct an exception payload with the given error code.
)doc")
      .def("ierr", &EI::Exceptions::ierr, R"doc(
Return the numeric error flag (see `EllImg.Exceptions.TYPE`).
)doc")
      .def("what", &EI::Exceptions::what, R"doc(
Return the human-readable error description.
)doc");

  py::class_<EV> pyEllVar(m, "EllVar", R"doc(
Scalar variable or expression in ellipsoidal arithmetic.

An `EllVar` is either attached to one coordinate of an `EllImg`
(an *image variable*), or represents a plain constant or interval
(not tied to any image). Arithmetic between `EllVar` objects is
overloaded: ``+``, ``-``, ``*``, ``/``, ``**`` (integer exponents)
and the `pymcpp` math functions (`exp`, `log`, `sqr`, `sqrt`,
`sin`, ...) each append a new dimension to the lifted ellipsoid of
the associated image and return a new `EllVar` indexing it, so that
evaluating a factorable function on image variables propagates the
ellipsoidal enclosure operation by operation. In-place operators
(``+=`` etc.) are supported as well.

Combining variables attached to different images raises
`pymcpp.EllImgError`; so do operations whose domain is violated by
the operand range (e.g. `log` of a range containing non-positive
values).

Examples
--------
>>> import pymcpp
>>> Ex = pymcpp.EllImg([[5., 4.], [4., 5.]], [3., 4.])
>>> X1 = pymcpp.EllVar(Ex, 0)
>>> X2 = pymcpp.EllVar(Ex, 1)
>>> Y = pymcpp.log(X1) + pymcpp.sqr(X2)
>>> Y.index  # row of Y in the lifted ellipsoid
4
)doc");

  pyEllVar.def(py::init<>(), R"doc(
Construct an uninitialized variable (not attached to any image).
)doc")
      .def(py::init<EV const&>(), R"doc(
Copy constructor.
)doc")
      .def(py::init<double const&>(), py::arg("cst"), R"doc(
Construct a constant (degenerate) variable.

Parameters
----------
cst : float
    Constant value.
)doc")
      .def(py::init<I const&>(), py::arg("range"), R"doc(
Construct a variable holding an interval range, not attached to any
image.

Parameters
----------
range : Interval
    Range of the variable.
)doc")
      .def(py::init<double const&, double const&>(), py::arg("l"), py::arg("u"),
           R"doc(
Construct a variable holding the interval range [l, u].

Parameters
----------
l : float
    Lower bound.
u : float
    Upper bound.
)doc")
      .def(py::init<EI&, unsigned const>(), py::arg("image"), py::arg("index"),
           py::keep_alive<1, 2>(), R"doc(
Construct an image variable attached to a coordinate of an `EllImg`.

The variable range is deduced from the ellipsoid as
``c[index] +/- sqrt(Q[index, index])``. The image is kept alive as
long as the variable exists.

Parameters
----------
image : EllImg
    Ellipsoidal image environment.
index : int
    Coordinate index in {0, ..., image.n - 1}.

Raises
------
EllImgError
    If the variable cannot be constructed (e.g. index out of range).
)doc")
      .def(py::init<EI&, unsigned const, I const&>(), py::arg("image"),
           py::arg("index"), py::arg("range"), py::keep_alive<1, 2>(), R"doc(
Construct an image variable with a user-supplied range.

Same as above, but ``range`` overrides the interval deduced from the
ellipsoid, which can tighten subsequent operation enclosures.

Parameters
----------
image : EllImg
    Ellipsoidal image environment.
index : int
    Coordinate index in {0, ..., image.n - 1}.
range : Interval
    Tailored range of the variable.
)doc")
      .def(
          "set", [](EV& self, EI& img, unsigned const i) -> EV&
          { return self.set(img, i); }, py::arg("image"), py::arg("index"),
          py::return_value_policy::reference_internal, py::keep_alive<1, 2>(),
          R"doc(
Re-attach this variable to coordinate ``index`` of image ``image``.

Equivalent to the ``EllVar(image, index)`` constructor.

Returns
-------
EllVar
    This variable (for chaining).
)doc")
      .def(
          "set", [](EV& self, EI& img, unsigned const i, I const& range) -> EV&
          { return self.set(img, i, range); }, py::arg("image"),
          py::arg("index"), py::arg("range"),
          py::return_value_policy::reference_internal, py::keep_alive<1, 2>(),
          R"doc(
Re-attach this variable with a user-supplied range.

Equivalent to the ``EllVar(image, index, range)`` constructor.

Returns
-------
EllVar
    This variable (for chaining).
)doc")
      .def_property_readonly("range", &EV::range, R"doc(
Interval: Interval enclosure of the variable range.
)doc")
      .def_property_readonly(
          "image", &EV::image, py::return_value_policy::reference_internal,
          R"doc(
EllImg: Ellipsoidal image the variable belongs to, or None for
constants and plain ranges.
)doc")
      .def_property_readonly(
          "index", &EV::index, R"doc(
int: Row index of the variable in the lifted ellipsoid, or -1 for
constants and plain ranges.
)doc")
      .def("__str__", static_cast<std::string (*)(EV const&)>(&_to_string))
      .def("__repr__", static_cast<std::string (*)(EV const&)>(&_to_string))
      // compound assignment: in place, returns self -> no new image reference
      .def(py::self += double())
      .def(py::self += I())
      .def(py::self += py::self)
      .def(py::self -= double())
      .def(py::self -= I())
      .def(py::self -= py::self)
      .def(py::self *= double())
      .def(py::self *= I())
      .def(py::self *= py::self)
      .def(py::self /= double())
      .def(py::self /= I())
      .def(py::self /= py::self)
      // unary and binary operators returning a NEW EllVar: tie the result's
      // lifetime to its EllVar operand(s) so the underlying EllImg stays alive
      .def(
          "__pos__", [](EV const& a) { return +a; }, py::keep_alive<0, 1>())
      .def(
          "__neg__", [](EV const& a) { return -a; }, py::keep_alive<0, 1>())
      .def(
          "__add__", [](EV const& a, EV const& b) { return a + b; },
          py::keep_alive<0, 1>(), py::keep_alive<0, 2>())
      .def(
          "__add__", [](EV const& a, double b) { return a + b; },
          py::keep_alive<0, 1>())
      .def(
          "__add__", [](EV const& a, I const& b) { return a + b; },
          py::keep_alive<0, 1>())
      .def(
          "__radd__", [](EV const& a, double b) { return b + a; },
          py::keep_alive<0, 1>())
      .def(
          "__radd__", [](EV const& a, I const& b) { return b + a; },
          py::keep_alive<0, 1>())
      .def(
          "__sub__", [](EV const& a, EV const& b) { return a - b; },
          py::keep_alive<0, 1>(), py::keep_alive<0, 2>())
      .def(
          "__sub__", [](EV const& a, double b) { return a - b; },
          py::keep_alive<0, 1>())
      .def(
          "__sub__", [](EV const& a, I const& b) { return a - b; },
          py::keep_alive<0, 1>())
      .def(
          "__rsub__", [](EV const& a, double b) { return b - a; },
          py::keep_alive<0, 1>())
      .def(
          "__rsub__", [](EV const& a, I const& b) { return b - a; },
          py::keep_alive<0, 1>())
      .def(
          "__mul__", [](EV const& a, EV const& b) { return a * b; },
          py::keep_alive<0, 1>(), py::keep_alive<0, 2>())
      .def(
          "__mul__", [](EV const& a, double b) { return a * b; },
          py::keep_alive<0, 1>())
      .def(
          "__mul__", [](EV const& a, I const& b) { return a * b; },
          py::keep_alive<0, 1>())
      .def(
          "__rmul__", [](EV const& a, double b) { return b * a; },
          py::keep_alive<0, 1>())
      .def(
          "__rmul__", [](EV const& a, I const& b) { return b * a; },
          py::keep_alive<0, 1>())
      .def(
          "__truediv__", [](EV const& a, EV const& b) { return a / b; },
          py::keep_alive<0, 1>(), py::keep_alive<0, 2>())
      .def(
          "__truediv__", [](EV const& a, double b) { return a / b; },
          py::keep_alive<0, 1>())
      .def(
          "__truediv__", [](EV const& a, I const& b) { return a / b; },
          py::keep_alive<0, 1>())
      .def(
          "__rtruediv__", [](EV const& a, double b) { return b / a; },
          py::keep_alive<0, 1>())
      .def(
          "__rtruediv__", [](EV const& a, I const& b) { return b / a; },
          py::keep_alive<0, 1>())
      .def(
          "__pow__", [](EV const& x, int const n) { return mc::pow(x, n); },
          py::keep_alive<0, 1>());

  m.def(
      "inv", [](EV const& x) { return mc::inv(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift 1/x into the ellipsoidal image; raises "
      "EllImgError if the range of x contains 0.");
  m.def(
      "sqr", [](EV const& x) { return mc::sqr(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift x**2 into the ellipsoidal image.");
  m.def(
      "sqrt", [](EV const& x) { return mc::sqrt(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift sqrt(x) into the ellipsoidal image; raises "
      "EllImgError if the range of x contains negative values.");
  m.def(
      "exp", [](EV const& x) { return mc::exp(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift exp(x) into the ellipsoidal image.");
  m.def(
      "log", [](EV const& x) { return mc::log(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift log(x) into the ellipsoidal image; raises "
      "EllImgError if the range of x contains non-positive values.");
  m.def(
      "xlog", [](EV const& x) { return mc::xlog(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift x*log(x) into the ellipsoidal image.");
  m.def(
      "cos", [](EV const& x) { return mc::cos(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift cos(x) into the ellipsoidal image.");
  m.def(
      "sin", [](EV const& x) { return mc::sin(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift sin(x) into the ellipsoidal image.");
  m.def(
      "tan", [](EV const& x) { return mc::tan(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift tan(x) into the ellipsoidal image; raises "
      "EllImgError if cos(x) can vanish on the range of x.");
  m.def(
      "acos", [](EV const& x) { return mc::acos(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift acos(x) into the ellipsoidal image; raises "
      "EllImgError if the range of x exceeds [-1,1].");
  m.def(
      "asin", [](EV const& x) { return mc::asin(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift asin(x) into the ellipsoidal image; raises "
      "EllImgError if the range of x exceeds [-1,1].");
  m.def(
      "atan", [](EV const& x) { return mc::atan(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift atan(x) into the ellipsoidal image.");
  m.def(
      "cosh", [](EV const& x) { return mc::cosh(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift cosh(x) into the ellipsoidal image.");
  m.def(
      "sinh", [](EV const& x) { return mc::sinh(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift sinh(x) into the ellipsoidal image.");
  m.def(
      "tanh", [](EV const& x) { return mc::tanh(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift tanh(x) into the ellipsoidal image.");
  m.def(
      "erf", [](EV const& x) { return mc::erf(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift erf(x) into the ellipsoidal image.");
  m.def(
      "erfc", [](EV const& x) { return mc::erfc(x); }, py::arg("x"),
      py::keep_alive<0, 1>(),
      "EllVar overload: lift erfc(x) into the ellipsoidal image.");
  m.def(
      "pow", [](EV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"), py::keep_alive<0, 1>(),
      "EllVar overload: lift x**n (integer n) into the ellipsoidal image.");
  m.def(
      "cheb", [](EV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"), py::keep_alive<0, 1>(),
      "EllVar overload: lift the Chebyshev polynomial T_n(x) into the "
      "ellipsoidal image.");

  m.def(
      "ell_unitball", [](unsigned const n) { return mc::ell_unitball(n); },
      py::arg("n"), R"doc(
Return the n-dimensional unit ball as an `Ellipsoid`.

The center is the origin and the shape matrix is the identity.

Parameters
----------
n : int
    Dimension.

Returns
-------
Ellipsoid
    The unit ball in R^n.
)doc");
  m.def(
      "mtimes",
      [](ELL const& E, py::handle A, py::object b)
      {
        return mc::mtimes(E, _to_arma_mat(A),
                          b.is_none() ? arma::vec() : _to_arma_vec(b));
      },
      py::arg("E"), py::arg("A"), py::arg("b") = py::none(), R"doc(
Apply the affine map x -> A @ x + b to an ellipsoid.

The image of an ellipsoid under an affine map is again an ellipsoid,
with center ``A @ E.c + b`` and shape matrix ``A @ E.Q @ A.T``; this
is exact (no over-approximation).

Parameters
----------
E : Ellipsoid
    Ellipsoid of dimension n.
A : array-like of shape (m, n)
    Linear map.
b : array-like of shape (m,), optional
    Translation vector. Default is zero.

Returns
-------
Ellipsoid
    The mapped ellipsoid, of dimension m.
)doc");
  m.def(
      "minksum_ea", [](ELL const& E1, ELL const& E2, double const eps)
      { return mc::minksum_ea(E1, E2, eps); }, py::arg("E1"), py::arg("E2"),
      py::arg("eps") = mc::machprec(), R"doc(
External ellipsoidal approximation of the Minkowski (geometric) sum
of two ellipsoids.

Returns the minimum-trace outer ellipsoid of ``E1 (+) E2``: the
center is ``E1.c + E2.c`` and the shape matrix is
``(s1 + s2)*(E1.Q/s1 + E2.Q/s2)`` with
``si = sqrt(trace(Ei.Q)) + eps``.

Parameters
----------
E1 : Ellipsoid
    First summand.
E2 : Ellipsoid
    Second summand (same dimension as ``E1``).
eps : float, optional
    Regularization added to each trace factor, guarding against
    degenerate (zero-trace) summands. Default is the machine
    precision.

Returns
-------
Ellipsoid
    Outer approximation of the Minkowski sum.
)doc");
  m.def(
      "minksum_ea",
      [](std::vector<ELL> const& E, py::object D) -> py::object
      {
        auto a = py::array_t<double, py::array::c_style |
                                         py::array::forcecast>::ensure(D);
        if (a && a.request().ndim == 1)
          return py::cast(mc::minksum_ea(E, _to_arma_vec(D)));
        return py::cast(mc::minksum_ea(E, _to_vec_vector(D)));
      },
      py::arg("ellipsoids"), py::arg("direction_or_directions"), R"doc(
External ellipsoidal approximation of the Minkowski sum of several
ellipsoids, tight along one or several given directions.

For a direction ``l``, the returned outer ellipsoid touches the
exact Minkowski sum in direction ``l`` (support function equality).

Parameters
----------
ellipsoids : list of Ellipsoid
    Summands, all of the same dimension n.
direction_or_directions : array-like of shape (n,), or iterable of
    such arrays. Direction(s) along which the approximation is
    tight.

Returns
-------
Ellipsoid or list of Ellipsoid
    A single outer approximation if one direction is given; a list
    with one outer approximation per direction otherwise.
)doc");
  m.def(
      "minksum_box",
      [](ELL const& E, py::handle r, py::object c, double const tol,
         double const eps)
      {
        return mc::minksum_ea(
            E,
            std::make_pair(_to_arma_vec(r),
                           c.is_none() ? arma::vec() : _to_arma_vec(c)),
            tol, eps);
      },
      py::arg("E"), py::arg("r"), py::arg("c") = py::none(),
      py::arg("tol") = 1e-10, py::arg("eps") = mc::machprec(), R"doc(
External ellipsoidal approximation of the Minkowski sum of an
ellipsoid and the interval box [c-r, c+r].

Parameters
----------
E : Ellipsoid
    Ellipsoid summand, of dimension n.
r : array-like of shape (n,)
    Nonnegative box half-widths (radii).
c : array-like of shape (n,), optional
    Box midpoint. Default is the origin.
tol : float, optional
    Tolerance of the trace-based approximation. Default is 1e-10.
eps : float, optional
    Regularization guarding against degenerate terms. Default is
    the machine precision.

Returns
-------
Ellipsoid
    Outer approximation of ``E (+) [c-r, c+r]``.
)doc");
  m.def(
      "minksum_interval",
      [](ELL const& E, std::vector<I> const& box, double const tol,
         double const eps) { return mc::minksum_ea(E, box.data(), tol, eps); },
      py::arg("E"), py::arg("box"), py::arg("tol") = 1e-10,
      py::arg("eps") = mc::machprec(), R"doc(
External ellipsoidal approximation of the Minkowski sum of an
ellipsoid and an interval box.

Same as `minksum_box` with the box given as a list of `Interval`
(midpoints as center, half-diameters as radii). The list length must
equal the dimension of ``E``.

Parameters
----------
E : Ellipsoid
    Ellipsoid summand, of dimension n.
box : list of Interval
    Interval bounds, one per coordinate.
tol : float, optional
    Tolerance of the trace-based approximation. Default is 1e-10.
eps : float, optional
    Regularization guarding against degenerate terms. Default is
    the machine precision.

Returns
-------
Ellipsoid
    Outer approximation of the Minkowski sum.
)doc");
  m.def(
      "inv", [](ELL const& E) { return mc::inv(E); }, py::arg("E"), R"doc(
Ellipsoid overload: return the ellipsoid with the same center and
the (pseudo-)inverse of the regularized shape matrix of E as shape.

Used internally by the ellipsoid intersection routines; this is not
the set-theoretic inverse image.

Parameters
----------
E : Ellipsoid
    Input ellipsoid (not modified).

Returns
-------
Ellipsoid
    Ellipsoid with inverted shape matrix.
)doc");
  m.def(
      "inv", [](std::vector<ELL> const& E) { return mc::inv(E); },
      py::arg("E"), R"doc(
Ellipsoid overload: apply the shape-matrix inversion elementwise to
a list of ellipsoids.

Parameters
----------
E : list of Ellipsoid
    Input ellipsoids.

Returns
-------
list of Ellipsoid
    New list with one inverted-shape ellipsoid per input.
)doc");
  m.def(
      "dist", [](ELL const& E, py::tuple hp)
      { return mc::dist(E, _to_hp(hp)); }, py::arg("E"), py::arg("hp"), R"doc(
Signed distance between an ellipsoid and a hyperplane.

The hyperplane is ``{ x : a @ x = b }``, given as the pair
``(a, b)``. The returned value is
``(|b - a @ E.c| - sqrt(a @ E.Q @ a)) / ||a||``: positive if the
hyperplane does not intersect the ellipsoid (then equal to the
Euclidean distance), negative if it cuts it.

Parameters
----------
E : Ellipsoid
    Ellipsoid of dimension n.
hp : tuple of (array-like of shape (n,), float)
    Hyperplane normal ``a`` and offset ``b``.

Returns
-------
float
    Signed distance (NaN if the dimensions disagree).
)doc");
  m.def(
      "dist", [](ELL const& E1, ELL const& E2) { return mc::dist(E1, E2); },
      py::arg("E1"), py::arg("E2"), R"doc(
Separation measure between two ellipsoids.

Returns a value that is positive if and only if the two ellipsoids
are disjoint, and non-positive if they overlap. The magnitude is a
scale-free separation indicator, not a Euclidean distance.

Parameters
----------
E1 : Ellipsoid
    First ellipsoid.
E2 : Ellipsoid
    Second ellipsoid (same dimension as ``E1``).

Returns
-------
float
    Positive if separated, non-positive if intersecting (NaN if the
    dimensions disagree).
)doc");
  m.def(
      "hpintersection", [](ELL const& E, py::tuple hp)
      { return mc::hpintersection(E, _to_hp(hp)); }, py::arg("E"),
      py::arg("hp"), R"doc(
Exact intersection of an ellipsoid with a hyperplane.

The intersection of ``E`` with the hyperplane ``{ x : a @ x = b }``
is itself an (degenerate, rank-deficient) ellipsoid lying inside
the hyperplane, which is returned exactly.

Parameters
----------
E : Ellipsoid
    Ellipsoid of dimension n.
hp : tuple of (array-like of shape (n,), float)
    Hyperplane normal ``a`` and offset ``b``.

Returns
-------
Ellipsoid
    The intersection; an empty (0-dimensional) ellipsoid if the
    hyperplane does not intersect ``E``.
)doc");
  m.def(
      "hpintersection", [](ELL const& E, py::iterable hp)
      { return mc::hpintersection(E, _to_hps(hp)); }, py::arg("E"),
      py::arg("hp"), R"doc(
Successive intersection of an ellipsoid with several hyperplanes.

Applies the single-hyperplane intersection for each element of
``hp`` in turn.

Parameters
----------
E : Ellipsoid
    Ellipsoid of dimension n.
hp : iterable of tuple of (array-like of shape (n,), float)
    Hyperplanes as (normal, offset) pairs.

Returns
-------
Ellipsoid
    The resulting (possibly empty) ellipsoid.
)doc");
  m.def(
      "intersection_ea", [](ELL const& E1, ELL const& E2, double const tol)
      { return mc::intersection_ea(E1, E2, tol); }, py::arg("E1"),
      py::arg("E2"), py::arg("tol") = mc::machprec(), R"doc(
External ellipsoidal approximation of the intersection of two
ellipsoids.

Computes a (small) outer ellipsoid containing ``E1 & E2`` using the
fusion parametrization of Ros et al. / the ellipsoidal toolbox.

Parameters
----------
E1 : Ellipsoid
    First ellipsoid.
E2 : Ellipsoid
    Second ellipsoid (same dimension as ``E1``).
tol : float, optional
    Relative inflation applied to the resulting shape matrix.
    Default is the machine precision.

Returns
-------
Ellipsoid
    Outer approximation of the intersection; an empty
    (0-dimensional) ellipsoid if ``E1`` and ``E2`` are disjoint.
)doc");
  m.def(
      "intersection_ea", [](ELL const& E, py::tuple hp, double const tol)
      { return mc::intersection_ea(E, _to_hp(hp), tol); }, py::arg("E"),
      py::arg("hp"), py::arg("tol") = mc::machprec(), R"doc(
External ellipsoidal approximation of the intersection of an
ellipsoid with a halfspace.

The halfspace is ``{ x : a @ x <= b }``, given as the pair
``(a, b)``. If ``E`` lies entirely inside the halfspace, ``E``
itself is returned; if they are disjoint, an empty (0-dimensional)
ellipsoid is returned.

Parameters
----------
E : Ellipsoid
    Ellipsoid of dimension n.
hp : tuple of (array-like of shape (n,), float)
    Halfspace normal ``a`` and offset ``b``.
tol : float, optional
    Relative inflation applied to the resulting shape matrix.
    Default is the machine precision.

Returns
-------
Ellipsoid
    Outer approximation of the intersection.
)doc");
  m.def(
      "ellintersection_ia",
      [](py::iterable hp, double const tol, unsigned const maxit)
      { return mc::ellintersection_ia(_to_hps(hp), tol, maxit); },
      py::arg("hp"), py::arg("tol") = 1e-4, py::arg("maxit") = 100, R"doc(
Maximum-volume ellipsoid inscribed in a polytope.

The polytope is the intersection of the halfspaces
``{ x : a_i @ x <= b_i }``. The inner ellipsoid is computed by
solving a semidefinite program (max-det problem) with the SDPA
solver.

Parameters
----------
hp : iterable of tuple of (array-like of shape (n,), float)
    Halfspaces as (normal, offset) pairs.
tol : float, optional
    SDP duality-gap tolerance. Default is 1e-4.
maxit : int, optional
    Maximum number of SDP iterations. Default is 100.

Returns
-------
Ellipsoid
    Maximum-volume inscribed ellipsoid.

Raises
------
RuntimeError
    If pymcpp was built without SDPA support, or the SDP solver
    fails to converge.
)doc");
}
