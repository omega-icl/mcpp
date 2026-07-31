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

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

namespace py = pybind11;

void
mc_mccormick(py::module& m)
{
  py::class_<MC> pyMcCormick(m, "McCormick", R"doc(
McCormick convex/concave relaxations of a factorable function.

A convex relaxation of a function f on a box is a convex function that
underestimates f on that box; a concave relaxation is a concave
function that overestimates f there. A ``McCormick`` object carries,
for one variable or intermediate expression: its interval range
(``l``, ``u``, ``I``), the value ``cv`` of a convex underestimator and
the value ``cc`` of a concave overestimator at the current point, and
optionally subgradients of these relaxations (``cvsub``, ``ccsub``).

A variable is typically constructed from its interval range and its
current point value, e.g. ``pymcpp.McCormick(pymcpp.Interval(-2, 1),
0)`` for x in [-2, 1] evaluated at x = 0. Evaluating a factorable
expression in McCormick arithmetic - using the overloaded operators
``+``, ``-``, ``*``, ``/``, ``**`` and the module-level functions
(``pymcpp.exp``, ``pymcpp.sin``, ...) - propagates the interval
bounds, the relaxation values at the current point, AND the
subgradients (if seeded via ``sub``) through every operation. The
result describes valid relaxations of the whole expression on the
variable box. To trace the relaxations over the box, repeat the
evaluation at different points (see method ``c``).

Subgradients are propagated forward like directional derivatives:
first fix the number of components with ``x.sub(n, i)`` (variable x is
component i of n) or seed custom directions with
``x.sub(cvsub, ccsub)``, then evaluate the expression. From the
resulting subgradients, affine under/overestimators can be formed with
``laff``/``uaff``, e.g. for use in cutting planes or bundle methods.

The interval bounds are computed with the interval backend the package
was compiled against (built-in mc::Interval by default; see
``Interval``). McCormick relaxation values and subgradients are
computed in ordinary floating-point arithmetic, without rounding
control. Tightness of the relaxations can be traded against speed via
the static ``McCormick.options`` (e.g. ``ENVEL_USE``, ``MVCOMP_USE``).

Comparison operators combine bound comparisons: ``x == y`` tests
equality of interval bounds and relaxation values; ``x <= y`` tests
whether x has tighter bounds than y (interval containment, larger
``cv``, and smaller ``cc``).

Examples
--------
>>> import pymcpp
>>> xMC = pymcpp.McCormick(pymcpp.Interval(-2, 1), 0)
>>> yMC = pymcpp.McCormick(pymcpp.Interval(-1, 2), 1)
>>> xMC.sub(2, 0)  # seed subgradients: x is component 0 of 2
[ -2.00000e+00 :  1.00000e+00 ] [  0.00000e+00 :  0.00000e+00 ]
[ ( 1.00000e+00, 0.00000e+00) : ( 1.00000e+00, 0.00000e+00) ]
>>> yMC.sub(2, 1)  # y is component 1 of 2
[ -1.00000e+00 :  2.00000e+00 ] [  1.00000e+00 :  1.00000e+00 ]
[ ( 0.00000e+00, 1.00000e+00) : ( 0.00000e+00, 1.00000e+00) ]
>>> fMC = xMC * (pymcpp.exp(xMC) - yMC)**2
>>> fMC.l, fMC.u    # interval bounds on the box
(-27.651239511697487, 13.825619755848743)
>>> fMC.cv, fMC.cc  # relaxation values at (0, 1)
(-13.825619755848743, 8.522454227517597)
>>> fMC.cvsub       # subgradient of the convex relaxation at (0, 1)
[13.825619755848743, 0.0]
)doc");
  pyMcCormick
      .def(py::init<>(), R"doc(
Construct an empty McCormick object with uninitialized bounds.
)doc")
      .def(py::init<double const&>(), py::arg("c"), R"doc(
Construct a McCormick constant: range, convex and concave bounds are
all set to the value c.
)doc")
      .def(py::init<I const&>(), py::arg("I"), R"doc(
Construct a McCormick object with range I; the convex and concave
bounds are set to the interval bounds (weakest valid relaxations).
)doc")
      .def(py::init<I const&, double const&>(), py::arg("I"), py::arg("c"),
           R"doc(
Construct a McCormick variable with range I and current point value c
(convex and concave bounds both equal to c). This is the usual way to
initialize an independent variable before evaluating an expression.
)doc")
      .def(py::init<I const&, double const&, double const&>(), py::arg("I"),
           py::arg("cv"), py::arg("cc"), R"doc(
Construct a McCormick object with range I, convex bound cv and
concave bound cc; the bounds are cut at the interval bounds.
)doc")
      .def(py::init<MC const&>(), py::arg("MC"), R"doc(
Construct a copy of another McCormick object (including any
subgradients).
)doc")
      .def_readwrite_static("options", &MC::options, R"doc(
Static options shared by all McCormick objects, of type
``McCormick.Options``; e.g.
``pymcpp.McCormick.options.MVCOMP_USE = True``.
)doc")
      .def_property_readonly(
          "l", [](MC& self) { return self.l(); }, R"doc(
Lower interval bound of the variable range (float, read-only).
)doc")
      .def_property_readonly(
          "u", [](MC& self) { return self.u(); }, R"doc(
Upper interval bound of the variable range (float, read-only).
)doc")
      .def_property(
          "I", [](MC& self) { return self.I(); },
          [](MC& self, I const& I) { self.I() = I; }, R"doc(
Interval range of the variable (Interval, settable).
)doc")
      .def_property(
          "cv", [](MC& self) { return self.cv(); },
          [](MC& self, double const& cv) { self.cv(cv); }, R"doc(
Value of the convex underestimator at the current point (float,
settable). Setting it resets the convex subgradient to zero.
)doc")
      .def_property(
          "cc", [](MC& self) { return self.cc(); },
          [](MC& self, double const& cc) { self.cc(cc); }, R"doc(
Value of the concave overestimator at the current point (float,
settable). Setting it resets the concave subgradient to zero.
)doc")
      .def(
          "c", [](MC& self, double const& c) { return self.c(c); },
          py::arg("c"), R"doc(
Set the current point value: convex and concave bounds are both set
to c. Use this to move the point at which relaxations of an
expression are evaluated.

Parameters
----------
c : float
    New point value; must lie inside the variable range.

Returns
-------
MC : McCormick
    This object, to allow method chaining.
)doc")
      .def_property_readonly(
          "nsub", [](MC& self) { return self.nsub(); }, R"doc(
Number of subgradient components (int, read-only); 0 if subgradients
are not being propagated. Set via method ``sub``.
)doc")
      .def_property_readonly(
          "cvsub",
          [](MC& self) {
            return std::vector<double>(self.cvsub(),
                                       self.cvsub() + self.nsub());
          },
          py::return_value_policy::take_ownership, R"doc(
Subgradient of the convex underestimator at the current point, as a
list of ``nsub`` floats (read-only copy).
)doc")
      .def_property_readonly(
          "ccsub",
          [](MC& self) {
            return std::vector<double>(self.ccsub(),
                                       self.ccsub() + self.nsub());
          },
          py::return_value_policy::take_ownership, R"doc(
Subgradient of the concave overestimator at the current point, as a
list of ``nsub`` floats (read-only copy).
)doc")
      .def("sub", static_cast<MC& (MC::*)(unsigned int const)>(&MC::sub),
           py::arg("size"), R"doc(
Set the number of subgradient components to size and reset all
components to zero (zero seed direction).

Parameters
----------
size : int
    Number of subgradient components, typically the number of
    independent variables.

Returns
-------
MC : McCormick
    This object, to allow method chaining.
)doc")
      .def("sub",
           static_cast<MC& (MC::*)(unsigned int const, unsigned int const)>(
               &MC::sub),
           py::arg("size"), py::arg("index"), R"doc(
Seed this object as independent variable number index out of size.

Sets the number of subgradient components to size and seeds the unit
direction: component index of both subgradients is set to 1, all
others to 0. Call this on every independent variable (with the same
size and distinct indices, starting at 0) before evaluating an
expression, so that subgradients are propagated.

Parameters
----------
size : int
    Number of subgradient components (number of variables).
index : int
    Component associated with this variable, in 0..size-1.

Returns
-------
MC : McCormick
    This object, to allow method chaining.

Raises
------
RuntimeError
    If index >= size.
)doc")
      .def(
          "sub",
          [](MC& self, std::vector<double> const& cvsub,
             std::vector<double> const& ccsub)
          {
            if (cvsub.size() != ccsub.size())
              throw(std::runtime_error("Inconsistent subgradient sizes"));
            return self.sub(cvsub.size(), cvsub.data(), ccsub.data());
          },
          py::return_value_policy::take_ownership, py::arg("cvsub"),
          py::arg("ccsub"), R"doc(
Seed custom subgradient values, e.g. for directional subgradients.

Sets the number of subgradient components to len(cvsub) and the
subgradients of the convex and concave relaxations to cvsub and
ccsub, respectively. Seeding a single component per variable with the
components of a direction d propagates directional subgradients along
d.

Parameters
----------
cvsub : list of float
    Seed values for the convex subgradient.
ccsub : list of float
    Seed values for the concave subgradient; same length as cvsub.

Returns
-------
MC : McCormick
    This object, to allow method chaining.

Raises
------
RuntimeError
    If cvsub and ccsub have different lengths.
)doc")
      .def("cut", static_cast<MC& (MC::*)()>(&MC::cut), R"doc(
Cut the relaxation values at the interval bounds, in place: if
cv < l, set cv = l (and zero the convex subgradient); if cc > u, set
cc = u (and zero the concave subgradient). Returns this object.
)doc")
      .def(
          "laff",
          [](MC& self, std::vector<double> const& val,
             std::vector<double> const& ref)
          {
            if (val.size() != self.nsub() || ref.size() != self.nsub())
              throw(std::runtime_error("Inconsistent sizes"));
            return self.laff(val.data(), ref.data());
          },
          py::arg("val"), py::arg("ref"), R"doc(
Evaluate the affine underestimator built from the convex subgradient.

Computes cv + cvsub . (val - ref), the value at point val of the
affine (linear) underestimator supported at the reference point ref,
where ref must be the point at which this object's relaxations were
computed.

Parameters
----------
val : list of float
    Point at which to evaluate the affine underestimator; length
    nsub.
ref : list of float
    Reference point at which the relaxation and subgradient were
    computed; length nsub.

Returns
-------
value : float
    Value of the affine underestimator at val.

Raises
------
RuntimeError
    If val or ref does not have length nsub.
)doc")
      .def(
          "laff",
          [](MC& self, std::vector<I> const& rng,
             std::vector<double> const& ref)
          {
            if (rng.size() != self.nsub() || ref.size() != self.nsub())
              throw(std::runtime_error("Inconsistent sizes"));
            return self.laff(rng.data(), ref.data());
          },
          py::arg("rng"), py::arg("ref"), R"doc(
Interval-box overload: lower bound min over the box rng of the affine
underestimator supported at ref, i.e.
cv + sum_i min(cvsub[i]*(rng[i] - ref[i])).

Parameters
----------
rng : list of Interval
    Box over which to bound the affine underestimator; length nsub.
ref : list of float
    Reference point of the relaxation; length nsub.

Returns
-------
bound : float
    Lower bound of the affine underestimator over rng.
)doc")
      .def(
          "uaff",
          [](MC& self, std::vector<double> const& val,
             std::vector<double> const& ref)
          {
            if (val.size() != self.nsub() || ref.size() != self.nsub())
              throw(std::runtime_error("Inconsistent sizes"));
            return self.uaff(val.data(), ref.data());
          },
          py::arg("val"), py::arg("ref"), R"doc(
Evaluate the affine overestimator built from the concave subgradient.

Computes cc + ccsub . (val - ref), the value at point val of the
affine overestimator supported at the reference point ref, where ref
must be the point at which this object's relaxations were computed.

Parameters
----------
val : list of float
    Point at which to evaluate the affine overestimator; length nsub.
ref : list of float
    Reference point of the relaxation; length nsub.

Returns
-------
value : float
    Value of the affine overestimator at val.

Raises
------
RuntimeError
    If val or ref does not have length nsub.
)doc")
      .def(
          "uaff",
          [](MC& self, std::vector<I> const& rng,
             std::vector<double> const& ref)
          {
            if (rng.size() != self.nsub() || ref.size() != self.nsub())
              throw(std::runtime_error("Inconsistent sizes"));
            return self.uaff(rng.data(), ref.data());
          },
          py::arg("rng"), py::arg("ref"), R"doc(
Interval-box overload: upper bound max over the box rng of the affine
overestimator supported at ref.

Parameters
----------
rng : list of Interval
    Box over which to bound the affine overestimator; length nsub.
ref : list of float
    Reference point of the relaxation; length nsub.

Returns
-------
bound : float
    Upper bound of the affine overestimator over rng.
)doc")
      //.def( "copy", []( MC const& self ){ return MC( self ); } )
      //.def( "assign", py::overload_cast<double const&>( &MC::operator= ),
      // py::arg("Value") ) .def( "assign", py::overload_cast<I const&>(
      //&MC::operator= ), py::arg("Value") ) .def( "assign",
      // py::overload_cast<MC const&>( &MC::operator= ), py::arg("Value") )
      .def("__str__",
           [](MC const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](MC const& self)
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
      .def(py::self /= double())
      .def(py::self /= py::self)
      .def(double() / py::self)
      .def(py::self / double())
      .def(py::self / py::self)
      .def(
          "__abs__", [](MC const& m) { return mc::Op<MC>::abs(m); }, R"doc(
Return the magnitude max(|l|, |u|) of the interval range as a float.

Note: this is NOT a relaxation of |x|; use ``pymcpp.fabs`` for that.
)doc")
      .def(
          "__pow__",
          [](MC const& m, int const n) { return mc::Op<MC>::pow(m, n); },
          R"doc(
Relaxation of x**n for integer exponent n.
)doc")
      .def(
          "__pow__",
          [](MC const& m, double const& r) { return mc::Op<MC>::pow(m, r); },
          R"doc(
Relaxation of x**r for real exponent r; requires x >= 0.
)doc")
      .def(
          "__pow__",
          [](MC const& m, MC const& mm) { return mc::Op<MC>::pow(m, mm); },
          R"doc(
Relaxation of x**y for a McCormick exponent y, as exp(y*log(x)).
)doc")
      .def(
          "__pow__",
          [](double const& r, MC const& m) { return mc::Op<MC>::pow(r, m); },
          R"doc(
Relaxation of r**x for a float base r, as exp(x*log(r)).
)doc")
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self <= py::self)
      .def(py::self >= py::self)
      .def(py::self < py::self)
      .def(py::self > py::self);

  py::class_<MC::Options> pyMcCormickOptions(pyMcCormick, "Options", R"doc(
Static options of the McCormick class.

Modify the shared instance ``pymcpp.McCormick.options`` rather than
creating new Options objects; the settings affect all subsequent
McCormick computations.
)doc");
  pyMcCormickOptions
      .def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def(py::init<MC::Options const&>(), py::arg("options"), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("ENVEL_USE", &MC::Options::ENVEL_USE, R"doc(
Whether to compute convex/concave envelopes for univariate functions
that are neither convex nor concave (odd powers, sin, cos, asin,
acos, tan, atan, erf, erfc). Tighter relaxations but more time
consuming. Default is True.
)doc")
      .def_readwrite("ENVEL_MAXIT", &MC::Options::ENVEL_MAXIT, R"doc(
Maximum number of iterations for locating junction points in the
convex/concave envelopes of univariate terms. Default is 100.
)doc")
      .def_readwrite("ENVEL_TOL", &MC::Options::ENVEL_TOL, R"doc(
Termination tolerance for locating junction points in the
convex/concave envelopes of univariate terms. Default is 1e-10.
)doc")
      .def_readwrite("MVCOMP_USE", &MC::Options::MVCOMP_USE, R"doc(
Whether to use Tsoukalas & Mitsos's multivariate composition result
for min/max, product and division terms. Tighter relaxations but more
time consuming. Default is False.
)doc")
      .def_readwrite("MVCOMP_TOL", &MC::Options::MVCOMP_TOL, R"doc(
Tolerance for equality tests in subgradient propagation of product
terms with the multivariate composition result. Default is 1e-10.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &MC::Options::DISPLAY_DIGITS, R"doc(
Number of significant digits used when printing McCormick objects
(int). Default is 5.
)doc");

  m.def(
      "cut", [](MC const& x) { return mc::cut(x); }, py::arg("x"), R"doc(
McCormick overload: return a copy of x with the relaxation values cut
at the interval bounds (cv raised to at least l, cc lowered to at
most u); see method ``McCormick.cut``.
)doc");
  m.def(
      "inv", [](MC const& x) { return mc::inv(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of 1/x; the range of x must not
contain 0.
)doc");
  m.def(
      "sqr", [](MC const& x) { return mc::sqr(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of x**2.
)doc");
  m.def(
      "sqrt", [](MC const& x) { return mc::sqrt(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of sqrt(x); requires x >= 0.
)doc");
  m.def(
      "exp", [](MC const& x) { return mc::exp(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of exp(x).
)doc");
  m.def(
      "log", [](MC const& x) { return mc::log(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the natural logarithm log(x);
requires x > 0.
)doc");
  m.def(
      "cos", [](MC const& x) { return mc::cos(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of cos(x).
)doc");
  m.def(
      "sin", [](MC const& x) { return mc::sin(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of sin(x).
)doc");
  m.def(
      "tan", [](MC const& x) { return mc::tan(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of tan(x); the range of x must not
contain a pole.
)doc");
  m.def(
      "acos", [](MC const& x) { return mc::acos(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of arccos(x); requires x within
[-1, 1].
)doc");
  m.def(
      "asin", [](MC const& x) { return mc::asin(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of arcsin(x); requires x within
[-1, 1].
)doc");
  m.def(
      "atan", [](MC const& x) { return mc::atan(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of arctan(x).
)doc");
  m.def(
      "cosh", [](MC const& x) { return mc::cosh(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of cosh(x).
)doc");
  m.def(
      "sinh", [](MC const& x) { return mc::sinh(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of sinh(x).
)doc");
  m.def(
      "tanh", [](MC const& x) { return mc::tanh(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of tanh(x).
)doc");
  m.def(
      "fabs", [](MC const& x) { return mc::fabs(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the absolute value |x|.
)doc");
  m.def(
      "relu", [](MC const& x) { return mc::max(x, 0.); }, py::arg("x"),
      R"doc(
McCormick overload: relaxation of the rectifier max(x, 0).
)doc");
  m.def(
      "xlog", [](MC const& x) { return mc::xlog(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of x*log(x); requires x >= 0.
)doc");
  m.def(
      "fstep", [](MC const& x) { return mc::fstep(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the forward step function
(1 if x >= 0, else 0), after Wechsung & Barton (2013).
)doc");
  m.def(
      "bstep", [](MC const& x) { return mc::bstep(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the backward step function
(0 if x >= 0, else 1), after Wechsung & Barton (2013).
)doc");
  m.def(
      "erf", [](MC const& x) { return mc::erf(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the error function erf(x).
)doc");
  m.def(
      "erfc", [](MC const& x) { return mc::erfc(x); }, py::arg("x"), R"doc(
McCormick overload: relaxation of the complementary error function
erfc(x) = 1 - erf(x).
)doc");
  m.def(
      "lmtd", [](MC const& x, MC const& y) { return mc::lmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxation of the log-mean temperature difference
lmtd(x, y) = (x - y)/(log(x) - log(y)), with lmtd(x, x) = x;
requires x, y > 0.
)doc");
  m.def(
      "rlmtd", [](MC const& x, MC const& y) { return mc::rlmtd(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxation of the reciprocal log-mean temperature
difference rlmtd(x, y) = (log(x) - log(y))/(x - y), with
rlmtd(x, x) = 1/x; requires x, y > 0.
)doc");
  m.def(
      "pow", [](MC const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
McCormick overload: relaxation of x**n for integer exponent n.
)doc");
  m.def(
      "pow", [](MC const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"), R"doc(
McCormick overload: relaxation of x**r for real exponent r;
requires x >= 0.
)doc");
  m.def(
      "pow", [](MC const& x, MC const& y) { return mc::pow(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxation of x**y for a McCormick exponent y,
computed as exp(y*log(x)).
)doc");
  m.def(
      "pow", [](double const& r, MC const& y) { return mc::pow(r, y); },
      py::arg("r"), py::arg("y"), R"doc(
McCormick overload: relaxation of r**y for a float base r, computed
as exp(y*log(r)); requires r > 0.
)doc");
  m.def(
      "cheb", [](MC const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"), R"doc(
McCormick overload: relaxation of the Chebyshev polynomial of the
first kind T_n(x) for x within [-1, 1].
)doc");
  m.def(
      "hull", [](MC const& x, MC const& y) { return mc::hull(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxations of the union of x and y (interval
hull of the ranges, with valid convex/concave bounds and subgradients
for the union).
)doc");
  m.def(
      "max", [](MC const& x, MC const& y) { return mc::max(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxation of the pointwise maximum max(x, y).
)doc");
  m.def(
      "min", [](MC const& x, MC const& y) { return mc::min(x, y); },
      py::arg("x"), py::arg("y"), R"doc(
McCormick overload: relaxation of the pointwise minimum min(x, y).
)doc");
  m.def(
      "inter",
      [](MC& z, MC const& x, MC const& y) { return mc::inter(z, x, y); },
      py::arg("z"), py::arg("x"), py::arg("y"), R"doc(
McCormick overload: intersect the enclosures x and y, storing the
result in z.

The object ``z`` passed as first argument is modified in place: on
success its interval range is set to the intersection of the ranges,
its convex bound to max(x.cv, y.cv), its concave bound to
min(x.cc, y.cc), with the corresponding subgradients.

Parameters
----------
z : McCormick
    Output object, overwritten with the intersection if it is
    nonempty.
x : McCormick
    First operand.
y : McCormick
    Second operand.

Returns
-------
nonempty : bool
    True if the intersection is nonempty (z holds the result),
    False if the enclosures are inconsistent (z's range may already
    have been partially updated in this case).

Raises
------
RuntimeError
    If x and y carry subgradients of different sizes.
)doc");
  m.def(
      "ltcond", [](I const& z, MC const& x, MC const& y)
      { return mc::ltcond(z, x, y); }, py::arg("z"), py::arg("x"),
      py::arg("y"), R"doc(
McCormick overload: relaxation of the conditional
{ x if z <= 0; y otherwise } for an Interval condition z, after
Wechsung & Barton (2013).
)doc");
  m.def(
      "ltcond", [](MC const& z, MC const& x, MC const& y)
      { return mc::ltcond(z, x, y); }, py::arg("z"), py::arg("x"),
      py::arg("y"), R"doc(
McCormick overload: relaxation of the conditional
{ x if z <= 0; y otherwise } for a McCormick condition z.
)doc");
  m.def(
      "gtcond", [](I const& z, MC const& x, MC const& y)
      { return mc::gtcond(z, x, y); }, py::arg("z"), py::arg("x"),
      py::arg("y"), R"doc(
McCormick overload: relaxation of the conditional
{ x if z >= 0; y otherwise } for an Interval condition z.
)doc");
  m.def(
      "gtcond", [](MC const& z, MC const& x, MC const& y)
      { return mc::gtcond(z, x, y); }, py::arg("z"), py::arg("x"),
      py::arg("y"), R"doc(
McCormick overload: relaxation of the conditional
{ x if z >= 0; y otherwise } for a McCormick condition z.
)doc");
}
