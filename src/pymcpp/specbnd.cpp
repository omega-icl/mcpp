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

#include "specbnd.hpp"

typedef mc::Specbnd<I> SB;

namespace py = pybind11;

void
mc_specbnd(py::module& m)
{
  py::class_<SB> pySpecbnd(m, "Specbnd", R"doc(
Spectral bound arithmetic for the Hessian of a factorable function.

A Specbnd variable propagates, through a factorable expression f(x_1,...,
x_n), an interval enclosure of the eigenvalue range of the Hessian matrix
of f at any point of the variable domain -- without ever forming the
Hessian. It simultaneously carries an enclosure of the function range and
of its gradient. The propagation follows Moennigmann's eigenvalue
arithmetic, using the operator overloading mechanism: build the expression
from Specbnd variables using the standard arithmetic operators (+, -, *, /,
**) and the module-level intrinsic functions (exp, log, sqr, sin, ...).

Typical use is convexity certification in global optimization: if the
spectral range ``F.SI`` of an expression F satisfies ``F.SI.l >= 0``, then f
is convex on the variable domain; if ``F.SI.u <= 0``, f is concave.

Static methods are also provided to bound the spectrum of explicitly given
symmetric (interval) Hessian matrices with Gershgorin's circle criterion,
Rohn's method or Hertz's method, e.g. for comparison with the eigenvalue
arithmetic. Depending on the function and domain, either approach may give
the tighter bounds.

Examples
--------
>>> import pymcpp
>>> N = 2
>>> x = pymcpp.Specbnd(pymcpp.Interval(-2, 1), 0, N)
>>> y = pymcpp.Specbnd(pymcpp.Interval(-1, 2), 1, N)
>>> f = x * (pymcpp.exp(x) - y)**2
>>> f.SI.l <= f.SI.u  # enclosure of the Hessian eigenvalue range
True
)doc");

  pySpecbnd
      // Constructors
      .def(py::init<>(), R"doc(
Construct an uninitialized spectral bound variable (all fields zero).
)doc")
      .def(py::init<double const&>(), py::arg("c"), R"doc(
Construct a constant with value `c`.

A constant has zero gradient and zero Hessian, hence a zero spectral
range.

Parameters
----------
c : float
    Constant value.
)doc")
      .def(py::init<I const&>(), py::arg("B"), R"doc(
Construct a constant enclosed by the interval `B`.

The variable does not depend on any independent variable; its gradient
and spectral ranges are zero.

Parameters
----------
B : Interval
    Enclosure of the constant value.
)doc")
      .def(py::init<I const&, size_t, size_t>(), py::arg("B"), py::arg("i"),
           py::arg("n"), R"doc(
Construct independent variable `i` of `n` with range `B`.

Parameters
----------
B : Interval
    Range of the independent variable.
i : int
    Index of the variable, in ``range(n)``.
n : int
    Total number of independent variables in the expression.
)doc")
      .def(py::init<SB const&>(), R"doc(
Construct a copy of another spectral bound variable.
)doc")
      // Modifiers
      .def(
          "set", [](SB& self, I const& B, size_t i, size_t n) -> SB&
          { return self.set(B, i, n); }, py::arg("B"), py::arg("i"),
          py::arg("n"), py::return_value_policy::reference_internal, R"doc(
Redefine this object as independent variable `i` of `n` with range `B`.

Same semantics as the constructor ``Specbnd(B, i, n)``; the object is
modified in place.

Parameters
----------
B : Interval
    Range of the independent variable.
i : int
    Index of the variable, in ``range(n)``.
n : int
    Total number of independent variables in the expression.

Returns
-------
var : Specbnd
    This object, updated in place.

Raises
------
RuntimeError
    If `i` is not smaller than `n`.
)doc")
      .def(
          "dep", [](SB& self, size_t i, size_t n) -> SB&
          { return self.dep(i, n); }, py::arg("i"), py::arg("n"),
          py::return_value_policy::reference_internal, R"doc(
Mark this object as independent variable `i` of `n`, keeping its range.

Unlike ``set``, the current function range is preserved; only the
dependency information and gradient seeding are (re)initialized.

Parameters
----------
i : int
    Index of the variable, in ``range(n)``.
n : int
    Total number of independent variables in the expression.

Returns
-------
var : Specbnd
    This object, updated in place.
)doc")
      // Accessors
      .def_property_readonly(
          "n", [](SB const& self) { return self.n(); }, R"doc(
Number of independent variables in the expression.
)doc")
      .def_property_readonly(
          "I", [](SB const& self) { return self.I(); }, R"doc(
Interval enclosure of the function range on the variable domain.
)doc")
      .def_property_readonly(
          "FI", [](SB const& self) { return self.FI(); }, R"doc(
Interval enclosure of the gradient. Only the first gradient component
(the partial derivative with respect to variable 0) is exposed; None if
the variable carries no gradient information.
)doc")
      .def_property_readonly(
          "SI", [](SB const& self) { return self.SI(); }, R"doc(
Interval enclosure of the spectrum (eigenvalue range) of the Hessian
matrix on the variable domain. If some dependency is known to be linear,
the enclosure is widened to contain zero.
)doc")
      .def_static(
          "spectrum",
          [](std::vector<double> const& hess)
          {
            size_t N = std::sqrt(hess.size());
            assert(N * N == hess.size());
            return SB::spectrum(N, hess.data());
          },
          py::arg("hess"), R"doc(
Compute the eigenvalue extremes of a dense symmetric real matrix.

Parameters
----------
hess : list of float
    Matrix entries in row-major order; the length must be a perfect
    square N*N. The matrix is symmetrized from its upper triangle.

Returns
-------
lmin : float
    Smallest eigenvalue.
lmax : float
    Largest eigenvalue.
)doc")
      .def_static(
          "spectrum",
          [](size_t ndim, std::vector<double> const& hess,
             std::vector<unsigned> const& row, std::vector<unsigned> const& col)
          {
            size_t NNZ = row.size();
            assert(!NNZ || NNZ == hess.size());
            return SB::spectrum(ndim, hess.data(), NNZ,
                                row.empty() ? nullptr : row.data(),
                                col.empty() ? nullptr : col.data());
          },
          py::arg("ndim"), py::arg("hess"),
          py::arg("row") = std::vector<unsigned>(),
          py::arg("col") = std::vector<unsigned>(), R"doc(
Compute the eigenvalue extremes of a sparse symmetric real matrix.

The matrix is given in coordinate format, as produced e.g. by
``FFGraph.bdiff`` applied twice to obtain a sparse Hessian.

Parameters
----------
ndim : int
    Dimension N of the matrix.
hess : list of float
    Nonzero entries; entry k is placed at position (row[k], col[k]).
    The matrix is symmetrized from its upper triangle.
row : list of int, optional
    Row indices of the nonzero entries. If empty (default), `hess` is
    interpreted as a dense row-major matrix of size N*N.
col : list of int, optional
    Column indices of the nonzero entries.

Returns
-------
lmin : float
    Smallest eigenvalue.
lmax : float
    Largest eigenvalue.

Examples
--------
>>> import pymcpp
>>> lmin, lmax = pymcpp.Specbnd.spectrum(2, [2.0, 1.0, 1.0, 2.0])
>>> (lmin, lmax)
(1.0, 3.0)
)doc")
      .def_static(
          "spectral_bound_gershgorin",
          [](std::vector<I> const& hess, std::vector<double> const& scal)
          {
            size_t N = std::sqrt(hess.size());
            assert(N * N == hess.size());
            return SB::spectral_bound_gershgorin(
                N, hess.data(), 0, nullptr, nullptr,
                scal.empty() ? nullptr : scal.data());
          },
          py::arg("hess"), py::arg("scal") = std::vector<double>(), R"doc(
Bound the spectrum of a dense symmetric interval matrix (Gershgorin).

Applies the interval variant of Gershgorin's circle criterion to an
interval Hessian matrix, i.e. an elementwise interval enclosure of all
possible Hessian matrices on a domain.

Parameters
----------
hess : list of Interval
    Matrix entries in row-major order; the length must be a perfect
    square N*N.
scal : list of float, optional
    Positive scaling factors, one per row/column, used to sharpen the
    Gershgorin discs. Default is no scaling.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def_static(
          "spectral_bound_gershgorin",
          [](size_t ndim, std::vector<I> const& hess,
             std::vector<int> const& row, std::vector<int> const& col,
             std::vector<double> const& scal)
          {
            size_t NNZ = row.size();
            assert(!NNZ || NNZ == hess.size());
            return SB::spectral_bound_gershgorin(
                ndim, hess.data(), NNZ,
                row.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(row.data()),
                col.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(col.data()),
                scal.empty() ? nullptr : scal.data());
          },
          py::arg("ndim"), py::arg("hess"), py::arg("row") = std::vector<int>(),
          py::arg("col")  = std::vector<int>(),
          py::arg("scal") = std::vector<double>(), R"doc(
Bound the spectrum of a sparse symmetric interval matrix (Gershgorin).

Same as the dense overload, with the interval matrix given in
coordinate format.

Parameters
----------
ndim : int
    Dimension N of the matrix.
hess : list of Interval
    Nonzero entries; entry k is placed at position (row[k], col[k]).
row : list of int, optional
    Row indices of the nonzero entries. If empty (default), `hess` is
    interpreted as a dense row-major matrix of size N*N.
col : list of int, optional
    Column indices of the nonzero entries.
scal : list of float, optional
    Positive scaling factors, one per row/column. Default is no scaling.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def_static(
          "spectral_bound_rohn",
          [](std::vector<I> const& hess)
          {
            size_t N = std::sqrt(hess.size());
            assert(N * N == hess.size());
            return SB::spectral_bound_rohn(N, hess.data());
          },
          py::arg("hess"), R"doc(
Bound the spectrum of a dense symmetric interval matrix (Rohn).

Rohn's method bounds the eigenvalues of all symmetric matrices contained
in an interval matrix from the eigenvalues of its midpoint and radius
matrices; it is often tighter than Gershgorin's criterion.

Parameters
----------
hess : list of Interval
    Matrix entries in row-major order; the length must be a perfect
    square N*N.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def_static(
          "spectral_bound_rohn",
          [](size_t ndim, std::vector<I> const& hess,
             std::vector<int> const& row, std::vector<int> const& col)
          {
            size_t NNZ = row.size();
            assert(!NNZ || NNZ == hess.size());
            return SB::spectral_bound_rohn(
                ndim, hess.data(), NNZ,
                row.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(row.data()),
                col.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(col.data()));
          },
          py::arg("ndim"), py::arg("hess"), py::arg("row") = std::vector<int>(),
          py::arg("col") = std::vector<int>(), R"doc(
Bound the spectrum of a sparse symmetric interval matrix (Rohn).

Same as the dense overload, with the interval matrix given in
coordinate format.

Parameters
----------
ndim : int
    Matrix dimension N.
hess : list of Interval
    Nonzero interval entries.
row : list of int, optional
    Row indices of the entries; if empty (default), `hess` is
    interpreted as a dense row-major matrix of size N*N.
col : list of int, optional
    Column indices of the entries.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def_static(
          "spectral_bound_hertz",
          [](std::vector<I> const& hess)
          {
            size_t N = std::sqrt(hess.size());
            assert(N * N == hess.size());
            return SB::spectral_bound_hertz(N, hess.data());
          },
          py::arg("hess"), R"doc(
Bound the spectrum of a dense symmetric interval matrix (Hertz).

Hertz's method computes the exact eigenvalue extremes over the interval
matrix by examining its 2**(N-1) vertex matrices; it is the tightest of
the provided methods but has exponential cost in the dimension N.

Parameters
----------
hess : list of Interval
    Matrix entries in row-major order; the length must be a perfect
    square N*N.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def_static(
          "spectral_bound_hertz",
          [](size_t ndim, std::vector<I> const& hess,
             std::vector<int> const& row, std::vector<int> const& col)
          {
            size_t NNZ = row.size();
            assert(!NNZ || NNZ == hess.size());
            return SB::spectral_bound_hertz(
                ndim, hess.data(), NNZ,
                row.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(row.data()),
                col.empty() ? nullptr
                            : reinterpret_cast<unsigned const*>(col.data()));
          },
          py::arg("ndim"), py::arg("hess"), py::arg("row") = std::vector<int>(),
          py::arg("col") = std::vector<int>(), R"doc(
Bound the spectrum of a sparse symmetric interval matrix (Hertz).

Same as the dense overload, with the interval matrix given in
coordinate format.

Parameters
----------
ndim : int
    Matrix dimension N.
hess : list of Interval
    Nonzero interval entries.
row : list of int, optional
    Row indices of the entries; if empty (default), `hess` is
    interpreted as a dense row-major matrix of size N*N.
col : list of int, optional
    Column indices of the entries.

Returns
-------
lmin : float
    Guaranteed lower bound on the smallest eigenvalue.
lmax : float
    Guaranteed upper bound on the largest eigenvalue.
)doc")
      .def("__str__",
           [](SB const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](SB const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def(+py::self)
      .def(py::self += double())
      .def(py::self += py::self)
      .def(double() + py::self)
      .def(py::self + double())
      .def(py::self + py::self)
      .def(-py::self)
      .def(py::self -= double())
      .def(py::self -= py::self)
      .def(double() - py::self)
      .def(py::self - double())
      .def(py::self - py::self)
      .def(py::self *= double())
      .def(py::self *= py::self)
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self * py::self)
      .def(double() / py::self)
      .def(py::self / double())
      .def(py::self / py::self)
      .def("__abs__", [](SB const& m) { return mc::Op<SB>::abs(m); })
      .def("__pow__", [](SB const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](SB const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__", [](SB const& m, SB const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, SB const& m) { return mc::pow(r, m); });

  m.def(
      "inv", [](SB const& x) { return mc::inv(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "1/x.");
  m.def(
      "sqr", [](SB const& x) { return mc::sqr(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "x**2.");
  m.def(
      "sqrt", [](SB const& x) { return mc::sqrt(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "sqrt(x).");
  m.def(
      "exp", [](SB const& x) { return mc::exp(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "exp(x).");
  m.def(
      "log", [](SB const& x) { return mc::log(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "log(x).");
  m.def(
      "xlog", [](SB const& x) { return mc::xlog(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "x*log(x).");
  m.def(
      "cos", [](SB const& x) { return mc::cos(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "cos(x).");
  m.def(
      "sin", [](SB const& x) { return mc::sin(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "sin(x).");
  m.def(
      "tan", [](SB const& x) { return mc::tan(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "tan(x).");
  m.def(
      "acos", [](SB const& x) { return mc::acos(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "acos(x).");
  m.def(
      "asin", [](SB const& x) { return mc::asin(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "asin(x).");
  m.def(
      "atan", [](SB const& x) { return mc::atan(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "atan(x).");
  m.def(
      "cosh", [](SB const& x) { return mc::cosh(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "cosh(x).");
  m.def(
      "sinh", [](SB const& x) { return mc::sinh(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "sinh(x).");
  m.def(
      "tanh", [](SB const& x) { return mc::tanh(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "tanh(x).");
  m.def(
      "erf", [](SB const& x) { return mc::erf(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "erf(x).");
  m.def(
      "erfc", [](SB const& x) { return mc::erfc(x); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "erfc(x).");
  m.def(
      "pow", [](SB const& x, int const n) { return mc::pow(x, n); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "x**n for integer n.");
  m.def(
      "pow", [](SB const& x, double const& r) { return mc::pow(x, r); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "x**r for real r.");
  m.def(
      "pow", [](SB const& x, SB const& y) { return mc::pow(x, y); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "x**y for Specbnd exponent y.");
  m.def(
      "pow", [](double const& r, SB const& y) { return mc::pow(r, y); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "r**y for real base r.");
  m.def(
      "cheb", [](SB const& x, unsigned const n) { return mc::cheb(x, n); },
      "Specbnd overload: bounds on value, gradient and Hessian spectrum of "
      "the Chebyshev polynomial T_n(x).");

  // Options class
  // py::class_<SB::Options> pySpecbndOptions(pySpecbnd, "Options");
  // pySpecbndOptions
  // .def(py::init<>())
  // .def(py::init<SB::Options const &>())
  // .def_readwrite(
  //   "HESSBND",
  //   &SB::Options::HESSBND,
  //   "strategy for computing spectral bounds in interval Hessian matrix"
  // )
  //;

  // Enum for HESSBND strategy
  // py::enum_<SB::Options::HESSBND_STRATEGY>(pySpecbndOptions,
  // "HESSBND_STRATEGY")
  // .value("GERSHGORIN", SB::Options::GERSHGORIN)
  // .value("HERTZROHN",  SB::Options::HERTZROHN)
  // .export_values()
  //;
}
