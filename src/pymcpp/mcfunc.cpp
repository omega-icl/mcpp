// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "mcfunc.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void
mc_mcfunc(py::module_& m)
{
  m.def(
      "inv", [](double const& x) { return 1 / x; }, py::arg("x"), R"doc(
float overload: reciprocal 1/x in plain floating-point arithmetic.
)doc");
  m.def(
      "sqr", [](double const& x) { return mc::sqr(x); }, py::arg("x"), R"doc(
float overload: square x**2 in plain floating-point arithmetic.
)doc");
  m.def(
      "sqrt", [](double const& x) { return std::sqrt(x); }, py::arg("x"),
      R"doc(
float overload: square root of x in plain floating-point arithmetic.
)doc");
  m.def(
      "exp", [](double const& x) { return std::exp(x); }, py::arg("x"),
      R"doc(
float overload: exponential exp(x) in plain floating-point arithmetic.
)doc");
  m.def(
      "log", [](double const& x) { return std::log(x); }, py::arg("x"),
      R"doc(
float overload: natural logarithm log(x) in plain floating-point
arithmetic.
)doc");
  m.def(
      "cos", [](double const& x) { return std::cos(x); }, py::arg("x"),
      R"doc(
float overload: cosine of x (radians).
)doc");
  m.def(
      "sin", [](double const& x) { return std::sin(x); }, py::arg("x"),
      R"doc(
float overload: sine of x (radians).
)doc");
  m.def(
      "tan", [](double const& x) { return std::tan(x); }, py::arg("x"),
      R"doc(
float overload: tangent of x (radians).
)doc");
  m.def(
      "acos", [](double const& x) { return std::acos(x); }, py::arg("x"),
      R"doc(
float overload: arc cosine of x, in radians.
)doc");
  m.def(
      "asin", [](double const& x) { return std::asin(x); }, py::arg("x"),
      R"doc(
float overload: arc sine of x, in radians.
)doc");
  m.def(
      "atan", [](double const& x) { return std::atan(x); }, py::arg("x"),
      R"doc(
float overload: arc tangent of x, in radians.
)doc");
  m.def(
      "cosh", [](double const& x) { return std::cosh(x); }, py::arg("x"),
      R"doc(
float overload: hyperbolic cosine of x.
)doc");
  m.def(
      "sinh", [](double const& x) { return std::sinh(x); }, py::arg("x"),
      R"doc(
float overload: hyperbolic sine of x.
)doc");
  m.def(
      "tanh", [](double const& x) { return std::tanh(x); }, py::arg("x"),
      R"doc(
float overload: hyperbolic tangent of x.
)doc");
  m.def(
      "fabs", [](double const& x) { return std::fabs(x); }, py::arg("x"),
      R"doc(
float overload: absolute value |x|.
)doc");
  m.def(
      "relu", [](double const& x) { return mc::relu(x); }, py::arg("x"),
      R"doc(
float overload: rectifier max(x, 0).
)doc");
  m.def(
      "xlog", [](double const& x) { return mc::xlog(x); }, py::arg("x"),
      R"doc(
float overload: x*log(x), with the value at x = 0 taken as 0.
)doc");
  m.def(
      "fstep", [](double const& x) { return mc::fstep(x); }, py::arg("x"),
      R"doc(
float overload: forward step function, 1.0 if x >= 0, else 0.0.
)doc");
  m.def(
      "bstep", [](double const& x) { return mc::bstep(x); }, py::arg("x"),
      R"doc(
float overload: backward step function, 0.0 if x >= 0, else 1.0.
)doc");
  m.def(
      "erf", [](double const& x) { return std::erf(x); }, py::arg("x"),
      R"doc(
float overload: error function erf(x).
)doc");
  m.def(
      "erfc", [](double const& x) { return std::erfc(x); }, py::arg("x"),
      R"doc(
float overload: complementary error function erfc(x) = 1 - erf(x).
)doc");
  m.def(
      "lmtd", [](double const& x, double const& y) { return mc::lmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
float overload: log-mean temperature difference
lmtd(x, y) = (x - y)/(log(x) - log(y)), with lmtd(x, x) = x;
requires x, y > 0.
)doc");
  m.def(
      "rlmtd",
      [](double const& x, double const& y) { return mc::rlmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
float overload: reciprocal log-mean temperature difference
rlmtd(x, y) = (log(x) - log(y))/(x - y), with rlmtd(x, x) = 1/x;
requires x, y > 0.
)doc");
  m.def(
      "pow", [](double const& x, int const n) { return std::pow(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
float overload: power x**n for integer exponent n.
)doc");
  m.def(
      "pow", [](double const& x, double const& r) { return std::pow(x, r); },
      py::arg("x"), py::arg("r"), R"doc(
float overload: power x**r for real exponent r.
)doc");
  m.def(
      "cheb", [](double const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
float overload: Chebyshev polynomial of the first kind T_n(x),
evaluated by the recurrence T_n = 2*x*T_{n-1} - T_{n-2}.
)doc");
  m.def(
      "max", [](double const& x, double const& y) { return mc::max(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
float overload: maximum of x and y.

Note: the underlying C++ function converts its arguments to unsigned
integers, so the result is truncated (e.g. max(2.7, 3.9) returns 3).
Prefer Python's built-in ``max`` for float arguments.
)doc");
  m.def(
      "min", [](double const& x, double const& y) { return mc::min(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
float overload: minimum of x and y.

Note: the underlying C++ function converts its arguments to unsigned
integers, so the result is truncated (e.g. min(2.7, 3.9) returns 2).
Prefer Python's built-in ``min`` for float arguments.
)doc");

  m.def("machprec", []() { return mc::machprec(); }, R"doc(
Return the machine precision (unit round-off) for double-precision
floating-point numbers, i.e. DBL_EPSILON = 2**-52.
)doc");
}
