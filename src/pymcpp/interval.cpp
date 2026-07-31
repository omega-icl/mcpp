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

namespace py = pybind11;

void
mc_interval(py::module_& m)
{
  py::class_<I> pyInterval(m, "Interval", R"doc(
A closed real interval [l, u] supporting interval arithmetic.

An ``Interval`` represents all real numbers between a lower bound ``l``
and an upper bound ``u``. The standard arithmetic operators (``+``,
``-``, ``*``, ``/``, ``**``) and the module-level elementary functions
(``pymcpp.exp``, ``pymcpp.log``, ``pymcpp.sin``, ...) are overloaded so
that the result encloses the exact image of the operands: if x is in X
and y is in Y, then f(x, y) is in F(X, Y) (inclusion property). This
makes it possible to compute guaranteed bounds on the range of a
factorable function simply by evaluating it in ``Interval`` arithmetic.
``Interval`` also serves as the supporting bounder for the other
arithmetics of MC++, e.g. ``McCormick``.

The bindings are compiled against one of several C++ interval backends:
the built-in ``mc::Interval`` (default), PROFIL, FILIB++, or the Boost
interval library. The Python API is identical for all backends, but the
built-in backend is NOT a verified implementation: it does not perform
outward (directed) rounding, so bounds may be inexact at the level of
floating-point round-off. For fully rigorous bounds, build pymcpp
against PROFIL, FILIB++, or Boost.

Comparison operators act setwise: ``x == y`` tests equality of the
bounds, ``x <= y`` tests whether x is contained in y, and ``x < y``
tests containment in the interior of y.

The number of digits used when printing an interval is controlled by
the static ``Interval.options.DISPLAY_DIGITS`` (built-in backend only).

Examples
--------
>>> import pymcpp
>>> x = pymcpp.Interval(-2., 1.)
>>> y = pymcpp.Interval(-1., 2.)
>>> f = x * (pymcpp.exp(x) - y)**2
>>> print(f)
[ -2.7651239511697487e+01 :  1.3825619755848743e+01 ]
>>> fl, fu = f.l, f.u
)doc");
  pyInterval
      .def(py::init<>(), R"doc(
Construct an interval with uninitialized bounds.

Set ``l`` and ``u`` before use.
)doc")
      .def(py::init<double const&>(), py::arg("c"), R"doc(
Construct the degenerate (thin) interval [c, c].
)doc")
      .def(py::init<double const&, double const&>(), py::arg("l"),
           py::arg("u"), R"doc(
Construct the interval [l, u].

If ``l > u``, the bounds are swapped.
)doc")
      .def(py::init<I const&>(), py::arg("I"), R"doc(
Construct a copy of another interval.
)doc")
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && \
    !defined(MC__USE_BOOST)
      .def_readwrite_static("options", &I::options, R"doc(
Static options shared by all Interval objects, of type
``Interval.Options``. Currently only controls display precision, e.g.
``pymcpp.Interval.options.DISPLAY_DIGITS = 7``.
)doc")
      .def_property(
          "l", [](I& self) { return mc::Op<I>::l(self); },
          [](I& self, double const& l) { self.l() = l; }, R"doc(
Lower bound of the interval (float, settable).
)doc")
      .def_property(
          "u", [](I& self) { return mc::Op<I>::u(self); },
          [](I& self, double const& u) { self.u() = u; }, R"doc(
Upper bound of the interval (float, settable).
)doc")
#else
// .def_property_readonly( "l", []( I& self ){ return mc::Op<I>::l(self); } )
// .def_property_readonly( "u", []( I& self ){ return mc::Op<I>::u(self); } )
 .def_property( "l", []( I& self ){ return mc::Op<I>::l(self); },
                     []( I& self, double const& l ){ self = I(l,mc::Op<I>::u(self)); },
                R"doc(
Lower bound of the interval (float, settable).
)doc" )
 .def_property( "u", []( I& self ){ return mc::Op<I>::u(self); },
                     []( I& self, double const& u ){ self = I(mc::Op<I>::l(self),u); },
                R"doc(
Upper bound of the interval (float, settable).
)doc" )
#endif
      .def("__str__",
           [](I const& x)
           {
             std::ostringstream Iss;
             Iss << x;
             return Iss.str();
           })
      .def("__repr__",
           [](I const& x)
           {
             std::ostringstream Iss;
             Iss << x;
             return Iss.str();
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
      .def(py::self /= double())
      .def(py::self /= py::self)
      .def(double() / py::self)
      .def(py::self / double())
      .def(py::self / py::self)
      .def(
          "__abs__", [](I const& i) { return mc::Op<I>::abs(i); }, R"doc(
Return the magnitude max(|l|, |u|) of the interval as a float.

Note: this is NOT an interval enclosure of |x|; use ``pymcpp.fabs``
for that.
)doc")
      .def(
          "__pow__",
          [](I const& i, int const n) { return mc::Op<I>::pow(i, n); }, R"doc(
Enclosure of x**n for integer exponent n.
)doc")
      .def(
          "__pow__",
          [](I const& i, double const& r) { return mc::Op<I>::pow(i, r); },
          R"doc(
Enclosure of x**r for real exponent r; requires x >= 0.
)doc")
      .def(
          "__pow__",
          [](I const& i, I const& ii) { return mc::Op<I>::pow(i, ii); },
          R"doc(
Enclosure of x**y for interval exponent y, as exp(y*log(x)).
)doc")
      .def(
          "__pow__",
          [](double const& r, I const& i) { return mc::Op<I>::pow(r, i); },
          R"doc(
Enclosure of r**x for a float base r, as exp(x*log(r)).
)doc")
      .def(
          "__le__", [](I const& x, I const& y) { return mc::Op<I>::le(x, y); },
          py::is_operator(), R"doc(
Return True if x is contained in y (set inclusion).
)doc")
      .def(
          "__lt__", [](I const& x, I const& y) { return mc::Op<I>::lt(x, y); },
          py::is_operator(), R"doc(
Return True if x is contained in the interior of y.
)doc")
      .def(
          "__ge__", [](I const& x, I const& y) { return mc::Op<I>::ge(x, y); },
          py::is_operator(), R"doc(
Return True if x contains y (set inclusion).
)doc")
      .def(
          "__gt__", [](I const& x, I const& y) { return mc::Op<I>::gt(x, y); },
          py::is_operator(), R"doc(
Return True if x contains y in its interior.
)doc")
      .def(
          "__eq__", [](I const& x, I const& y) { return mc::Op<I>::eq(x, y); },
          py::is_operator(), R"doc(
Return True if x and y have identical lower and upper bounds.
)doc")
      .def(
          "__ne__", [](I const& x, I const& y) { return mc::Op<I>::ne(x, y); },
          py::is_operator(), R"doc(
Return True if x and y differ in their lower or upper bound.
)doc");

#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && \
    !defined(MC__USE_BOOST)
  py::class_<I::Options> pyIntervalOptions(pyInterval, "Options", R"doc(
Static options of the Interval class (built-in backend).

Modify the shared instance ``pymcpp.Interval.options`` rather than
creating new Options objects.
)doc");
  pyIntervalOptions.def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def(py::init<I::Options const&>(), py::arg("options"), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &I::Options::DISPLAY_DIGITS, R"doc(
Number of significant digits used when printing intervals (int).
Default is 5.
)doc");
#endif

  m.def(
      "abs", [](I const& x) { return mc::Op<I>::abs(x); }, py::arg("x"),
      R"doc(
Interval overload: magnitude max(|x.l|, |x.u|) of x, as a float.

Not an interval enclosure of |x|; see ``fabs`` for that.
)doc");
  m.def(
      "mid", [](I const& x) { return mc::Op<I>::mid(x); }, py::arg("x"),
      R"doc(
Interval overload: midpoint (x.l + x.u)/2 of x, as a float.
)doc");
  m.def(
      "diam", [](I const& x) { return mc::Op<I>::diam(x); }, py::arg("x"),
      R"doc(
Interval overload: diameter (width) x.u - x.l of x, as a float.
)doc");
  m.def(
      "inv", [](I const& x) { return mc::Op<I>::inv(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of 1/x; raises an error if 0 is in x.
)doc");
  m.def(
      "sqr", [](I const& x) { return mc::Op<I>::sqr(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of x**2.
)doc");
  m.def(
      "sqrt", [](I const& x) { return mc::Op<I>::sqrt(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of sqrt(x); requires x >= 0.
)doc");
  m.def(
      "exp", [](I const& x) { return mc::Op<I>::exp(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of exp(x).
)doc");
  m.def(
      "log", [](I const& x) { return mc::Op<I>::log(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the natural logarithm log(x);
requires x > 0.
)doc");
  m.def(
      "cos", [](I const& x) { return mc::Op<I>::cos(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of cos(x).
)doc");
  m.def(
      "sin", [](I const& x) { return mc::Op<I>::sin(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of sin(x).
)doc");
  m.def(
      "tan", [](I const& x) { return mc::Op<I>::tan(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of tan(x); x must not contain a pole.
)doc");
  m.def(
      "acos", [](I const& x) { return mc::Op<I>::acos(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of arccos(x); requires x within [-1, 1].
)doc");
  m.def(
      "asin", [](I const& x) { return mc::Op<I>::asin(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of arcsin(x); requires x within [-1, 1].
)doc");
  m.def(
      "atan", [](I const& x) { return mc::Op<I>::atan(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of arctan(x).
)doc");
  m.def(
      "cosh", [](I const& x) { return mc::Op<I>::cosh(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of cosh(x).
)doc");
  m.def(
      "sinh", [](I const& x) { return mc::Op<I>::sinh(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of sinh(x).
)doc");
  m.def(
      "tanh", [](I const& x) { return mc::Op<I>::tanh(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of tanh(x).
)doc");
  m.def(
      "fabs", [](I const& x) { return mc::Op<I>::fabs(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the absolute value |x|.
)doc");
  m.def(
      "relu", [](I const& x) { return mc::Op<I>::max(x, 0.); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the rectifier max(x, 0).
)doc");
  m.def(
      "xlog", [](I const& x) { return mc::Op<I>::xlog(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of x*log(x); requires x >= 0
(the value at 0 is taken as 0).
)doc");
  m.def(
      "fstep", [](I const& x) { return mc::Op<I>::fstep(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the forward step function
(1 if x >= 0, else 0); returns [0, 1] if x contains 0.
)doc");
  m.def(
      "bstep", [](I const& x) { return mc::Op<I>::bstep(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the backward step function
(0 if x >= 0, else 1); returns [0, 1] if x contains 0.
)doc");
  m.def(
      "erf", [](I const& x) { return mc::Op<I>::erf(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the error function erf(x).
)doc");
  m.def(
      "erfc", [](I const& x) { return mc::Op<I>::erfc(x); }, py::arg("x"),
      R"doc(
Interval overload: enclosure of the complementary error function
erfc(x) = 1 - erf(x).
)doc");
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && \
    !defined(MC__USE_BOOST)
  m.def(
      "lmtd", [](I const& x, I const& y) { return mc::lmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: enclosure of the log-mean temperature difference
lmtd(x, y) = (x - y)/(log(x) - log(y)), with lmtd(x, x) = x;
requires x, y > 0. Built-in interval backend only.
)doc");
  m.def(
      "rlmtd", [](I const& x, I const& y) { return mc::rlmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: enclosure of the reciprocal log-mean temperature
difference rlmtd(x, y) = (log(x) - log(y))/(x - y), with
rlmtd(x, x) = 1/x; requires x, y > 0. Built-in interval backend only.
)doc");
#endif
  m.def(
      "pow", [](I const& x, int const n) { return mc::Op<I>::pow(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
Interval overload: enclosure of x**n for integer exponent n.
)doc");
  m.def(
      "pow", [](I const& x, double const& r) { return mc::Op<I>::pow(x, r); },
      py::arg("x"), py::arg("r"), R"doc(
Interval overload: enclosure of x**r for real exponent r;
requires x >= 0.
)doc");
  m.def(
      "pow", [](I const& x, I const& y) { return mc::Op<I>::pow(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: enclosure of x**y for interval exponent y,
computed as exp(y*log(x)).
)doc");
  m.def(
      "pow",
      [](double const& r, I const& y)
      { return mc::Op<I>::exp(y * std::log(r)); },
      py::arg("r"), py::arg("y"), R"doc(
Interval overload: enclosure of r**y for a float base r,
computed as exp(y*log(r)); requires r > 0.
)doc");
  m.def(
      "cheb", [](I const& x, unsigned const n) { return mc::Op<I>::cheb(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
Interval overload: enclosure of the Chebyshev polynomial of the first
kind T_n(x), intersected with [-1, 1] for x within [-1, 1].
)doc");
  m.def(
      "hull", [](I const& x, I const& y) { return mc::Op<I>::hull(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: interval hull of x and y, i.e. the smallest
interval [min(x.l, y.l), max(x.u, y.u)] containing both.
)doc");
  m.def(
      "max", [](I const& x, I const& y) { return mc::Op<I>::max(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: enclosure [max(x.l, y.l), max(x.u, y.u)] of the
elementwise maximum max(x, y).
)doc");
  m.def(
      "min", [](I const& x, I const& y) { return mc::Op<I>::min(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
Interval overload: enclosure [min(x.l, y.l), min(x.u, y.u)] of the
elementwise minimum min(x, y).
)doc");
  m.def(
      "inter",
      [](I& z, I const& x, I const& y) { return mc::Op<I>::inter(z, x, y); },
      py::arg("z"), py::arg("x"), py::arg("y"), R"doc(
Interval overload: intersect x and y, storing the result in z.

The interval ``z`` passed as first argument is modified in place: on
success it is overwritten with the intersection of x and y.

Parameters
----------
z : Interval
    Output interval, overwritten with x intersected with y if the
    intersection is nonempty; unchanged otherwise.
x : Interval
    First operand.
y : Interval
    Second operand.

Returns
-------
nonempty : bool
    True if x and y intersect (z holds the intersection), False if
    they are disjoint (z is left unchanged).
)doc");
}
