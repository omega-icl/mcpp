#include "ffdep.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>

namespace py = pybind11;

void
mc_ffdep(py::module_& m)
{
  typedef mc::FFDep FD;

  py::class_<FD> pyFFDep(m, "FFDep", R"doc(
Structural dependence analysis of factorable functions.

An FFDep object records, for a scalar expression, which independent
variables participate in the expression (its sparsity pattern) and the kind
of dependence on each of them. The dependence types, from most to least
specific, are:

- ``FFDep.TYPE.L`` (linear): the variable participates in linear terms only.
- ``FFDep.TYPE.Q`` (quadratic): in terms that are at most quadratic.
- ``FFDep.TYPE.P`` (polynomial): in terms that are at most polynomial.
- ``FFDep.TYPE.R`` (rational): in terms that are at most rational.
- ``FFDep.TYPE.N`` (nonlinear): in general nonlinear terms (exp, log,
  trigonometric, fractional powers, ...).

To analyze a function, create one FFDep object per variable with
``FFDep().indep(i)``, then evaluate the function on these objects - either
by calling a Python function directly, or by passing them to
``FFGraph.eval`` to propagate through a DAG. The standard arithmetic
operators (+, -, *, /) and the module-level functions (exp, log, sqrt,
pow, ...) are overloaded for FFDep and combine the dependence information
of their operands. Query the result with dep() or worst().

Examples
--------
>>> import pymcpp
>>> X = [pymcpp.FFDep().indep(i) for i in range(4)]
>>> F0 = X[2] * X[3] + X[0] / X[2]
>>> print(F0)
{ 0R 2R 3Q }
>>> F0.dep(2)
(True, <TYPE.R: 3>)
)doc");
  py::class_<FD::Exceptions> pyFDExceptions(pyFFDep, "Exceptions", R"doc(
Exception information for errors raised during FFDep arithmetic.

Instances carry an error code (ierr) and a human-readable description
(what).
)doc");

  py::enum_<FD::Exceptions::TYPE>(pyFDExceptions, "TYPE",
                                  "Error codes for FFDep exceptions")
      .value("UNDEF", FD::Exceptions::UNDEF,
             "Call to a feature that is unavailable in FFDep arithmetic")
      .value("INTERN", FD::Exceptions::INTERN, "Internal error")
      .export_values();

  pyFDExceptions
      .def("ierr", &FD::Exceptions::ierr, R"doc(
Return the error code as an integer (see FFDep.Exceptions.TYPE).
)doc")
      .def("what", &FD::Exceptions::what, R"doc(
Return a string describing the error.
)doc");

  py::enum_<FD::TYPE>(pyFFDep, "TYPE", R"doc(
Dependence type of an expression on a participating variable.

The values are ordered by increasing generality, L < Q < P < R < N: an
operation never lowers the recorded type, and combining terms keeps the
most general (worst) type per variable.
)doc")
      .value("L", FD::TYPE::L,
             "Linear: variable participates in linear terms only")
      .value("Q", FD::TYPE::Q,
             "Quadratic: variable participates in terms that are at most "
             "quadratic")
      .value("P", FD::TYPE::P,
             "Polynomial: variable participates in terms that are at most "
             "polynomial")
      .value("R", FD::TYPE::R,
             "Rational: variable participates in terms that are at most "
             "rational (quotients of polynomials)")
      .value("N", FD::TYPE::N,
             "General nonlinear: variable participates in general nonlinear "
             "terms (exp, log, trigonometric, fractional powers, ...)")
      .export_values();

  // --- Main Class ---
  pyFFDep
      // Constructors
      .def(py::init<>(), R"doc(
Construct an FFDep object with an empty dependence set (a constant).
)doc")
      .def(py::init<double const>(), py::arg("val"), R"doc(
Construct an FFDep object from a real constant.

The value itself is discarded; the resulting object has an empty
dependence set, like the default constructor.

Parameters
----------
val : float
    Constant value (ignored for the purpose of dependence analysis).
)doc")
      .def(py::init<FD const&>(), R"doc(
Construct a copy of another FFDep object.
)doc")

      // Modifiers
      .def(
          "indep", [](FD& self, int const ndx) { return self.indep(ndx); },
          py::arg("ndx"), R"doc(
Initialize as the independent variable with index `ndx`.

Any previously recorded dependencies are discarded, and the object is set
to depend linearly ('L') on the single variable index `ndx`. The object is
modified in place.

Parameters
----------
ndx : int
    Index identifying the independent variable.

Returns
-------
dep : FFDep
    The object itself, allowing the idiom
    ``X = [pymcpp.FFDep().indep(i) for i in range(NX)]``.
)doc")

      // Accessors
      .def(
          "dep", [](FD const& self, int const ndx) { return self.dep(ndx); },
          py::arg("ndx"), R"doc(
Query the dependence on the variable with index `ndx`.

Parameters
----------
ndx : int
    Index of the queried variable.

Returns
-------
depends : bool
    True if the expression depends on variable `ndx`.
dtype : FFDep.TYPE
    Dependence type on that variable. When `depends` is False, this is the
    placeholder value FFDep.TYPE.L and carries no meaning.
)doc")
      .def(
          "dep", [](FD const& self) { return self.dep(); }, R"doc(
Return the full dependence map of the expression.

Returns
-------
dep : dict of int to FFDep.TYPE
    Dictionary mapping the index of every participating variable to the
    type of dependence on that variable. Variables that do not participate
    are absent from the dictionary.
)doc")
      .def("worst", &FD::worst, R"doc(
Return the most general (worst-case) dependence type over all variables.

Returns
-------
dtype : FFDep.TYPE
    The largest dependence type among all participating variables in the
    order L < Q < P < R < N; FFDep.TYPE.L if the dependence set is empty.
)doc")
      .def("update", &FD::update, py::arg("dep"), R"doc(
Raise every recorded dependence to at least the given type, in place.

Entries whose current type is more general than `dep` are left unchanged.

Parameters
----------
dep : FFDep.TYPE
    Minimum dependence type to enforce on all participating variables.

Returns
-------
dep : FFDep
    The object itself.
)doc")

      // Operators
      .def(py::self += py::self)
      .def(py::self += double())
      .def(py::self -= py::self)
      .def(py::self -= double())
      .def(py::self *= py::self)
      .def(py::self *= double())
      .def(py::self /= py::self)
      .def(py::self /= double())

      .def(py::self + py::self)
      .def(py::self + double())
      .def(double() + py::self)
      .def(py::self - py::self)
      .def(py::self - double())
      .def(double() - py::self)
      .def(py::self * py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(py::self / py::self)
      .def(py::self / double())
      .def(double() / py::self)

      // String representation
      .def("__str__",
           [](FD const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           })
      .def("__repr__",
           [](FD const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           });

  // --- Free Functions ---
  m.def(
      "sqr", [](FD const& x) { return mc::sqr(x); }, R"doc(
FFDep overload: dependence structure of x**2; dependencies become at least
quadratic 'Q'.
)doc");
  m.def(
      "sqrt", [](FD const& x) { return mc::sqrt(x); }, R"doc(
FFDep overload: dependence structure of sqrt(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "exp", [](FD const& x) { return mc::exp(x); }, R"doc(
FFDep overload: dependence structure of exp(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "log", [](FD const& x) { return mc::log(x); }, R"doc(
FFDep overload: dependence structure of log(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "xlog", [](FD const& x) { return mc::xlog(x); }, R"doc(
FFDep overload: dependence structure of x*log(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "cos", [](FD const& x) { return mc::cos(x); }, R"doc(
FFDep overload: dependence structure of cos(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "sin", [](FD const& x) { return mc::sin(x); }, R"doc(
FFDep overload: dependence structure of sin(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "tan", [](FD const& x) { return mc::tan(x); }, R"doc(
FFDep overload: dependence structure of tan(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "acos", [](FD const& x) { return mc::acos(x); }, R"doc(
FFDep overload: dependence structure of acos(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "asin", [](FD const& x) { return mc::asin(x); }, R"doc(
FFDep overload: dependence structure of asin(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "atan", [](FD const& x) { return mc::atan(x); }, R"doc(
FFDep overload: dependence structure of atan(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "cosh", [](FD const& x) { return mc::cosh(x); }, R"doc(
FFDep overload: dependence structure of cosh(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "sinh", [](FD const& x) { return mc::sinh(x); }, R"doc(
FFDep overload: dependence structure of sinh(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "tanh", [](FD const& x) { return mc::tanh(x); }, R"doc(
FFDep overload: dependence structure of tanh(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "fabs", [](FD const& x) { return mc::fabs(x); }, R"doc(
FFDep overload: dependence structure of abs(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "erf", [](FD const& x) { return mc::erf(x); }, R"doc(
FFDep overload: dependence structure of erf(x); dependencies become
nonlinear 'N'.
)doc");
  m.def(
      "fstep", [](FD const& x) { return mc::fstep(x); }, R"doc(
FFDep overload: dependence structure of the forward step function;
dependencies become nonlinear 'N'.
)doc");
  m.def(
      "bstep", [](FD const& x) { return mc::bstep(x); }, R"doc(
FFDep overload: dependence structure of the backward step function;
dependencies become nonlinear 'N'.
)doc");
  m.def(
      "pow", [](FD const& x, int n) { return mc::pow(x, n); }, R"doc(
FFDep overload: dependence structure of x**n for integer n; dependencies
become polynomial 'P' for n > 2, rational 'R' for n < 0 (unchanged for
n in {0, 1}, quadratic 'Q' for n == 2).
)doc");
  m.def(
      "pow", [](FD const& x, double r) { return mc::pow(x, r); }, R"doc(
FFDep overload: dependence structure of x**r for real r; dependencies
become nonlinear 'N' unless r has an integer value, in which case the
integer-power rules apply.
)doc");
  m.def(
      "pow", [](FD const& x, FD const& y) { return mc::pow(x, y); }, R"doc(
FFDep overload: dependence structure of x**y; the dependence sets of x and
y are merged and all dependencies become nonlinear 'N'.
)doc");
  m.def(
      "min", [](FD const& x, FD const& y) { return mc::min(x, y); }, R"doc(
FFDep overload: dependence structure of min(x, y); the dependence sets are
merged and all dependencies become nonlinear 'N'.
)doc");
  m.def(
      "max", [](FD const& x, FD const& y) { return mc::max(x, y); }, R"doc(
FFDep overload: dependence structure of max(x, y); the dependence sets are
merged and all dependencies become nonlinear 'N'.
)doc");
  m.def(
      "cheb", [](FD const& x, unsigned n) { return mc::cheb(x, n); }, R"doc(
FFDep overload: dependence structure of the Chebyshev polynomial T_n(x);
dependencies become polynomial 'P' for n > 2.
)doc");
  m.def(
      "prod",
      [](unsigned n, std::vector<FD> const& x)
      { return mc::prod(n, x.data()); },
      R"doc(
FFDep overload: dependence structure of the product x[0]*...*x[n-1],
applying the FFDep multiplication rules term by term.
)doc");
  m.def(
      "monom",
      [](unsigned n, std::vector<FD> const& x, std::vector<unsigned> const& k)
      {
        if (x.size() != k.size())
          throw std::invalid_argument(
              "Size mismatch between variables and exponents");
        return mc::monom(n, x.data(), k.data());
      },
      R"doc(
FFDep overload: dependence structure of the monomial
x[0]**k[0]*...*x[n-1]**k[n-1], applying the FFDep power and multiplication
rules term by term.
)doc");
}
