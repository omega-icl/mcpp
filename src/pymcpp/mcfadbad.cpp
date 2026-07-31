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

#include "badiff.h"
#include "fadbad.h"
#include "fadiff.h"
#include "mcfadbad.hpp"
#include "tadiff.h"
typedef fadbad::FTypeName<MC> FMC;
typedef fadbad::BTypeName<MC> BMC;
typedef fadbad::TTypeName<MC> TMC;

namespace py = pybind11;

void
mc_fadbad_forward(py::module& m)
{
  py::class_<FMC> pyFadbadF(m, "FadbadF", R"doc(
Forward automatic differentiation over McCormick relaxations.

FadbadF combines the FADBAD++ forward AD type F<> with the MC++ McCormick
arithmetic: every FadbadF holds a McCormick relaxation of a function value
together with McCormick relaxations of its partial derivatives. Evaluating
a factorable expression of FadbadF variables therefore yields guaranteed
bounds and convex/concave relaxations of both the function AND its
derivatives -- e.g. to propagate subgradient information or to bound the
entries of a Jacobian over a box.

Workflow: wrap each independent variable's McCormick object in a FadbadF,
seed it with ``diff(i, N)`` to declare it the i-th of N independents, then
build the expression with the overloaded arithmetic operators (+, -, *, /,
**, comparisons) and the module-level functions (exp, log, sqr, ...). On
the result, ``value`` is the McCormick relaxation of the function and
``deriv(i)`` (or ``[i]``) the McCormick relaxation of the i-th partial
derivative.

Examples
--------
>>> import pymcpp
>>> x = pymcpp.FadbadF(pymcpp.McCormick(pymcpp.Interval(1, 2), 1.5))
>>> _ = x.diff(0, 1)              # seed x as variable 0 of 1
>>> f = pymcpp.sqr(x) + x
>>> f.value.l, f.value.u          # bounds on f = x**2 + x on [1,2]
(2.0, 6.0)
>>> f.deriv(0).l, f.deriv(0).u    # bounds on df/dx = 2x + 1 on [1,2]
(3.0, 5.0)
)doc");
  pyFadbadF
      .def(py::init<MC const&>(), R"doc(
Construct from a McCormick relaxation, with no derivative information
attached yet (use ``diff`` to seed it as an independent variable).
)doc")
      .def(py::init<const double&>(), R"doc(
Construct from a real constant, with no derivative information attached.
)doc")
      .def(py::init<const FMC&>(), R"doc(
Construct a copy of another forward AD variable.
)doc")
      .def_property_readonly(
          "value", [](FMC const& self) { return self.val(); }, R"doc(
McCormick relaxation of the function value (a McCormick object).
)doc")
      .def_property_readonly(
          "size", [](FMC const& self) { return self.size(); }, R"doc(
Length of the derivative vector: the number of independent variables N
used when seeding, or 0 if the object carries no derivative
information.
)doc")
      .def(
          "__getitem__", [](FMC const& self, unsigned int i)
          { return self.operator[](i); }, R"doc(
Return the McCormick relaxation of the i-th partial derivative.

Requires ``0 <= i < size``; unlike ``deriv``, out-of-range indices are
an error.
)doc")
      .def(
          "variable", [](FMC& self) { return self.x(); }, R"doc(
Return the McCormick relaxation of the function value.

Equivalent to the ``value`` property.

Returns
-------
val : McCormick
    McCormick relaxation of the function value.
)doc")
      .def(
          "deriv", [](FMC& self, const unsigned int i)
          { return self.deriv(i); }, py::arg("i"), R"doc(
Return the McCormick relaxation of the i-th partial derivative.

Parameters
----------
i : int
    Index of the independent variable to differentiate against.

Returns
-------
d : McCormick
    McCormick relaxation of the partial derivative with respect to
    variable `i`; zero if the object carries no derivative in slot `i`.
)doc")
      .def(
          "diff", [](FMC& self, unsigned int i, unsigned int N)
          { return self.diff(i, N); }, py::arg("i"), py::arg("N"), R"doc(
Seed this object as independent variable `i` of `N`.

Allocates a derivative vector of length `N` and sets it to the i-th
unit vector (derivative 1 with respect to itself, 0 otherwise). Must be
called on each independent variable before building the expression.

Parameters
----------
i : int
    Index of this variable, in ``range(N)``.
N : int
    Total number of independent variables.

Returns
-------
d : McCormick
    The seeded derivative component (the constant 1).
)doc")
      .def(
          "depend", [](FMC const& self) { return self.depend(); }, R"doc(
Return True if the object carries derivative information (it was seeded
with ``diff`` or results from operations on seeded variables).
)doc")
      .def(
          "setDepend", [](FMC& self, FMC const& val) { self.setDepend(val); },
          py::arg("val"), R"doc(
Allocate this object's derivative vector with the same length as that
of `val`. Low-level helper; normal use is ``diff`` on the independent
variables instead.
)doc")
      .def(
          "setDepend", [](FMC& self, FMC const& val1, FMC const& val2)
          { self.setDepend(val1, val2); }, py::arg("val1"), py::arg("val2"),
          R"doc(
Allocate this object's derivative vector with the common length of
those of `val1` and `val2`, which must agree.
)doc")

      .def(+py::self)
      .def(py::self + double())
      .def(double() + py::self)
      .def(py::self + py::self)
      .def(py::self += double())
      .def(py::self += py::self)
      .def(-py::self)
      .def(py::self - double())
      .def(double() - py::self)
      .def(py::self - py::self)
      .def(py::self -= double())
      .def(py::self -= py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(py::self * py::self)
      .def(py::self *= double())
      .def(py::self *= py::self)
      .def(py::self / double())
      .def(double() / py::self)
      .def(py::self / py::self)
      .def(py::self /= double())
      .def(py::self /= py::self)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self < py::self)
      .def(py::self <= py::self)
      .def(py::self > py::self)
      .def(py::self >= py::self)
      .def(py::self == double())
      .def(double() == py::self)
      .def(py::self != double())
      .def(double() != py::self)
      .def(py::self < double())
      .def(double() < py::self)
      .def(py::self <= double())
      .def(double() <= py::self)
      .def(py::self > double())
      .def(double() > py::self)
      .def(py::self >= double())
      .def(double() >= py::self)
      .def("__pow__",
           [](double const& a, FMC const& m) { return fadbad::pow(a, m); })
      .def("__pow__",
           [](int const& a, FMC const& m) { return fadbad::pow(a, m); })
      .def("__pow__",
           [](FMC const& m, double const& a) { return fadbad::pow(m, a); })
      .def("__pow__",
           [](FMC const& m, int const& a) { return fadbad::pow(m, a); })
      .def("__pow__",
           [](FMC const& m1, FMC const& m2) { return fadbad::pow(m1, m2); })
      .def("__pos__", [](FMC const& m) { return operator+(m); })
      .def("__neg__", [](FMC const& m) { return operator-(m); });

  m.def(
      "sqr", [](FMC const& m) { return fadbad::sqr(m); },
      "FadbadF overload: McCormick relaxation of x**2 and of its "
      "derivatives (forward AD).");
  m.def(
      "exp", [](FMC const& m) { return fadbad::exp(m); },
      "FadbadF overload: McCormick relaxation of exp(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "log", [](FMC const& m) { return fadbad::log(m); },
      "FadbadF overload: McCormick relaxation of log(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "sqrt", [](FMC const& m) { return fadbad::sqrt(m); },
      "FadbadF overload: McCormick relaxation of sqrt(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "sin", [](FMC const& m) { return fadbad::sin(m); },
      "FadbadF overload: McCormick relaxation of sin(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "cos", [](FMC const& m) { return fadbad::cos(m); },
      "FadbadF overload: McCormick relaxation of cos(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "tan", [](FMC const& m) { return fadbad::tan(m); },
      "FadbadF overload: McCormick relaxation of tan(x) and of its "
      "derivatives (forward AD).");
  // m.def( "cot",   []( FMC const& m ){ return fadbad::cot(m); } );
  m.def(
      "asin", [](FMC const& m) { return fadbad::asin(m); },
      "FadbadF overload: McCormick relaxation of asin(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "acos", [](FMC const& m) { return fadbad::acos(m); },
      "FadbadF overload: McCormick relaxation of acos(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "atan", [](FMC const& m) { return fadbad::atan(m); },
      "FadbadF overload: McCormick relaxation of atan(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "sinh", [](FMC const& m) { return fadbad::sinh(m); },
      "FadbadF overload: McCormick relaxation of sinh(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "cosh", [](FMC const& m) { return fadbad::cosh(m); },
      "FadbadF overload: McCormick relaxation of cosh(x) and of its "
      "derivatives (forward AD).");
  m.def(
      "tanh", [](FMC const& m) { return fadbad::tanh(m); },
      "FadbadF overload: McCormick relaxation of tanh(x) and of its "
      "derivatives (forward AD).");
  // m.def( "coth",  []( FMC const& m ){ return fadbad::coth(m); } );
  // cot and coth are not defined for the FMC
}

void
mc_fadbad_backward(py::module& m)
{
  py::class_<BMC> pyFadbadB(m, "FadbadB", R"doc(
Backward (reverse-mode) automatic differentiation over McCormick
relaxations.

FadbadB combines the FADBAD++ backward AD type B<> with the MC++ McCormick
arithmetic. As in FadbadF, all quantities are McCormick relaxations, so the
derivatives obtained are guaranteed enclosures/relaxations of the true
partial derivatives over the variable ranges. Reverse mode computes the
derivatives with respect to ALL inputs in a single sweep, which is
efficient for functions with many inputs and few outputs.

Workflow (note the reversal compared with FadbadF): build the expression
from FadbadB variables first; then call ``diff(i, N)`` on the i-th of N
DEPENDENT (output) variables to trigger the reverse sweep; finally read the
derivatives off the INPUT variables with ``deriv(i)``, giving the McCormick
relaxation of d(output i)/d(input).

Examples
--------
>>> import pymcpp
>>> x = pymcpp.FadbadB(pymcpp.McCormick(pymcpp.Interval(1, 2), 1.5))
>>> f = pymcpp.sqr(x) + x
>>> f.diff(0, 1)                  # declare f as dependent 0 of 1
>>> x.deriv(0).l, x.deriv(0).u    # bounds on df/dx = 2x + 1 on [1,2]
(3.0, 5.0)
)doc");
  pyFadbadB
      .def(py::init<MC const&>(), R"doc(
Construct from a McCormick relaxation; the object records operations
applied to it so that a later ``diff`` call on an output can propagate
derivatives back to it.
)doc")
      .def(py::init<const double&>(), R"doc(
Construct from a real constant.
)doc")
      .def(py::init<const BMC&>(), R"doc(
Construct a copy of another backward AD variable.
)doc")
      .def_property_readonly(
          "value", [](BMC const& self) { return self.val(); }, R"doc(
McCormick relaxation of the function value (a McCormick object).
)doc")
      .def(
          "variable", [](BMC& self) { return self.x(); }, R"doc(
Return the McCormick relaxation of the function value.

Equivalent to the ``value`` property.

Returns
-------
val : McCormick
    McCormick relaxation of the function value.
)doc")
      .def(
          "deriv", [](BMC& self, const unsigned int i)
          { return self.deriv(i); }, py::arg("i"), R"doc(
Return the derivative of dependent `i` with respect to this variable.

Only valid on an INPUT variable after the reverse sweep has been
triggered by calling ``diff`` on every dependent and after all recorded
operations have been propagated; FADBAD++ raises an error if
unpropagated dependencies remain.

Parameters
----------
i : int
    Index of the dependent (output) variable.

Returns
-------
d : McCormick
    McCormick relaxation of the partial derivative of dependent `i`
    with respect to this variable.
)doc")
      .def(
          "diff", [](BMC self, const unsigned int i, const unsigned int N)
          { return self.diff(i, N); }, py::arg("i"), py::arg("N"), R"doc(
Declare this object as dependent (output) `i` of `N` and start the
reverse derivative sweep from it.

Call this on each output of the expression; afterwards, retrieve the
derivatives from the input variables with ``deriv``.

Parameters
----------
i : int
    Index of this dependent among the outputs, in ``range(N)``.
N : int
    Total number of dependents (outputs).
)doc")

      .def(py::self += double())
      .def(py::self += py::self)
      .def(py::self -= double())
      .def(py::self -= py::self)
      .def(py::self *= double())
      .def(py::self *= py::self)
      .def(py::self /= double())
      .def(py::self /= py::self)

      .def("__eq__", [](BMC const& m1, BMC const& m2) { return m1 == m2; })
      .def("__eq__", [](BMC const& m1, double const& m2) { return m1 == m2; })
      .def("__eq__", [](double const& m1, BMC const& m2) { return m1 == m2; })
      .def("__ne__", [](BMC const& m1, BMC const& m2) { return m1 != m2; })
      .def("__ne__", [](BMC const& m1, double const& m2) { return m1 != m2; })
      .def("__ne__", [](double const& m1, BMC const& m2) { return m1 != m2; })
      .def("__lt__", [](BMC const& m1, BMC const& m2) { return m1 < m2; })
      .def("__lt__", [](BMC const& m1, double const& m2) { return m1 < m2; })
      .def("__lt__", [](double const& m1, BMC const& m2) { return m1 < m2; })
      .def("__le__", [](BMC const& m1, BMC const& m2) { return m1 <= m2; })
      .def("__le__", [](BMC const& m1, double const& m2) { return m1 <= m2; })
      .def("__le__", [](double const& m1, BMC const& m2) { return m1 <= m2; })
      .def("__gt__", [](BMC const& m1, BMC const& m2) { return m1 > m2; })
      .def("__gt__", [](BMC const& m1, double const& m2) { return m1 > m2; })
      .def("__gt__", [](double const& m1, BMC const& m2) { return m1 > m2; })
      .def("__ge__", [](BMC const& m1, BMC const& m2) { return m1 >= m2; })
      .def("__ge__", [](BMC const& m1, double const& m2) { return m1 >= m2; })
      .def("__ge__", [](double const& m1, BMC const& m2) { return m1 >= m2; })
      .def("__pos__", [](BMC const& m) { return operator+(m); })
      .def("__neg__", [](BMC const& m) { return operator-(m); })
      .def("__pow__",
           [](BMC const& m1, BMC const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](BMC const& m1, double const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](double const& m1, BMC const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](BMC const& m, int const& n) { return fadbad::pow(m, n); })
      .def("__pow__",
           [](int const& n, BMC const& m) { return fadbad::pow(n, m); })

      .def(py::self + double())
      .def(double() + py::self)
      .def(py::self + py::self)
      .def(py::self - double())
      .def(double() - py::self)
      .def(py::self - py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(py::self * py::self)
      .def(py::self / double())
      .def(double() / py::self)
      .def(py::self / py::self);

  m.def(
      "sqr", [](BMC const& m) { return fadbad::sqr(m); },
      "FadbadB overload: McCormick relaxation of x**2, recorded for "
      "reverse-mode differentiation.");
  m.def(
      "sqrt", [](BMC const& m) { return fadbad::sqrt(m); },
      "FadbadB overload: McCormick relaxation of sqrt(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "exp", [](BMC const& m) { return fadbad::exp(m); },
      "FadbadB overload: McCormick relaxation of exp(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "log", [](BMC const& m) { return fadbad::log(m); },
      "FadbadB overload: McCormick relaxation of log(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "sin", [](BMC const& m) { return fadbad::sin(m); },
      "FadbadB overload: McCormick relaxation of sin(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "cos", [](BMC const& m) { return fadbad::cos(m); },
      "FadbadB overload: McCormick relaxation of cos(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "tan", [](BMC const& m) { return fadbad::tan(m); },
      "FadbadB overload: McCormick relaxation of tan(x), recorded for "
      "reverse-mode differentiation.");
  // m.def( "cot",   []( BMC const& m ){ return fadbad::cot(m); } );
  m.def(
      "asin", [](BMC const& m) { return fadbad::asin(m); },
      "FadbadB overload: McCormick relaxation of asin(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "acos", [](BMC const& m) { return fadbad::acos(m); },
      "FadbadB overload: McCormick relaxation of acos(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "atan", [](BMC const& m) { return fadbad::atan(m); },
      "FadbadB overload: McCormick relaxation of atan(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "sinh", [](BMC const& m) { return fadbad::sinh(m); },
      "FadbadB overload: McCormick relaxation of sinh(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "cosh", [](BMC const& m) { return fadbad::cosh(m); },
      "FadbadB overload: McCormick relaxation of cosh(x), recorded for "
      "reverse-mode differentiation.");
  m.def(
      "tanh", [](BMC const& m) { return fadbad::tanh(m); },
      "FadbadB overload: McCormick relaxation of tanh(x), recorded for "
      "reverse-mode differentiation.");
  // m.def( "coth",  []( BMC const& m ){ return fadbad::coth(m); } );

  ;
}

void
mc_fadbad_taylor(py::module& m)
{
  py::class_<TMC> pyFadbadT(m, "FadbadT", R"doc(
Taylor expansion arithmetic over McCormick relaxations.

FadbadT combines the FADBAD++ Taylor AD type T<> with the MC++ McCormick
arithmetic: it propagates the Taylor coefficients of a univariate expansion
of a factorable function, each coefficient being a McCormick relaxation.
This yields guaranteed bounds and convex/concave relaxations of the Taylor
coefficients (i.e. of the higher-order derivatives divided by factorials)
of the function along the expansion parameter.

Workflow: create a FadbadT for the expansion variable, set its first-order
coefficient to 1 with ``x[1] = 1.0`` to mark it as the expansion parameter,
build the expression, then call ``eval(k)`` on the result to compute its
Taylor coefficients up to order `k`; retrieve coefficient j with ``f[j]``.
Call ``reset`` before re-evaluating after changing inputs.

Examples
--------
>>> import pymcpp
>>> x = pymcpp.FadbadT(pymcpp.McCormick(pymcpp.Interval(1, 2), 1.5))
>>> x[1] = 1.0                  # mark x as the expansion variable
>>> f = pymcpp.sqr(x)           # f = x**2
>>> _ = f.eval(2)               # Taylor coefficients up to order 2
>>> f[2].l, f[2].u              # coefficient of order 2 is exactly 1
(1.0, 1.0)
)doc");
  pyFadbadT
      .def(py::init<MC const&>(), R"doc(
Construct from a McCormick relaxation used as the zeroth-order Taylor
coefficient (the expansion point value).
)doc")
      .def(py::init<const double&>(), R"doc(
Construct from a real constant used as the zeroth-order Taylor
coefficient.
)doc")
      .def(py::init<const TMC&>(), R"doc(
Construct a copy of another Taylor AD variable.
)doc")
      .def_property_readonly(
          "value", [](TMC const& self) { return self.val(); }, R"doc(
McCormick relaxation of the zeroth-order Taylor coefficient (the
function value at the expansion point).
)doc")
      .def_property_readonly(
          "size", [](TMC const& self) { return self.length(); }, R"doc(
Number of Taylor coefficients currently computed (highest computed
order plus one).
)doc")
      .def(
          "__getitem__", [](TMC const& self, unsigned int i)
          { return self[i]; }, R"doc(
Return the McCormick relaxation of the i-th Taylor coefficient.
)doc")
      .def(
          "__setitem__", [](TMC& self, unsigned int i, double const& val)
          { self[i] = val; }, R"doc(
Set the i-th Taylor coefficient to the real constant `val`; typically
used to seed the expansion variable with ``x[1] = 1.0``.
)doc")
      .def(
          "reset", [](TMC& self) { return self.reset(); }, R"doc(
Discard the computed Taylor coefficients so that the expression can be
re-evaluated (e.g. after changing the coefficients of its inputs).
)doc")
      .def(
          "eval", [](TMC& self, const unsigned int i) { return self.eval(i); },
          py::arg("i"), R"doc(
Compute the Taylor coefficients of the expression up to order `i`.

Parameters
----------
i : int
    Highest expansion order to compute.

Returns
-------
order : int
    Number of valid coefficients after the call (highest computed
    order plus one).
)doc")
      //  .def(
      //    "diff",
      //    []( const TMC& self, const int i ){ return fadbad::diff(self,i); },
      //    "differentiate f to order i"
      //  )

      .def(py::self += double())
      .def(py::self += py::self)
      .def(py::self -= double())
      .def(py::self -= py::self)
      .def(py::self *= double())
      .def(py::self *= py::self)
      .def(py::self /= double())
      .def(py::self /= py::self)

      .def("__eq__", [](TMC const& m1, TMC const& m2) { return m1 == m2; })
      .def("__eq__", [](TMC const& m1, double const& m2) { return m1 == m2; })
      .def("__eq__", [](double const& m1, TMC const& m2) { return m1 == m2; })
      .def("__ne__", [](TMC const& m1, TMC const& m2) { return m1 != m2; })
      .def("__ne__", [](TMC const& m1, double const& m2) { return m1 != m2; })
      .def("__ne__", [](double const& m1, TMC const& m2) { return m1 != m2; })
      .def("__lt__", [](TMC const& m1, TMC const& m2) { return m1 < m2; })
      .def("__lt__", [](TMC const& m1, double const& m2) { return m1 < m2; })
      .def("__lt__", [](double const& m1, TMC const& m2) { return m1 < m2; })
      .def("__le__", [](TMC const& m1, TMC const& m2) { return m1 <= m2; })
      .def("__le__", [](TMC const& m1, double const& m2) { return m1 <= m2; })
      .def("__le__", [](double const& m1, TMC const& m2) { return m1 <= m2; })
      .def("__gt__", [](TMC const& m1, TMC const& m2) { return m1 > m2; })
      .def("__gt__", [](TMC const& m1, double const& m2) { return m1 > m2; })
      .def("__gt__", [](double const& m1, TMC const& m2) { return m1 > m2; })
      .def("__ge__", [](TMC const& m1, TMC const& m2) { return m1 >= m2; })
      .def("__ge__", [](TMC const& m1, double const& m2) { return m1 >= m2; })
      .def("__ge__", [](double const& m1, TMC const& m2) { return m1 >= m2; })
      .def("__pos__", [](TMC const& m) { return operator+(m); })
      .def("__neg__", [](TMC const& m) { return operator-(m); })
      .def("__pow__",
           [](TMC const& m1, TMC const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](TMC const& m1, double const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](double const& m1, TMC const& m2) { return fadbad::pow(m1, m2); })
      .def("__pow__",
           [](TMC const& m, int const& n) { return fadbad::pow(m, n); })
      .def("__pow__",
           [](int const& m, TMC const& n) { return fadbad::pow(m, n); })

      .def(py::self + double())
      .def(double() + py::self)
      .def(py::self + py::self)
      .def(py::self - double())
      .def(double() - py::self)
      .def(py::self - py::self)
      .def(py::self * double())
      .def(double() * py::self)
      .def(py::self * py::self)
      .def(py::self / double())
      .def(double() / py::self)
      .def(py::self / py::self);

  m.def(
      "sqr", [](TMC const& m) { return fadbad::sqr(m); },
      "FadbadT overload: Taylor coefficients of x**2 in McCormick "
      "arithmetic.");
  m.def(
      "sqrt", [](TMC const& m) { return fadbad::sqrt(m); },
      "FadbadT overload: Taylor coefficients of sqrt(x) in McCormick "
      "arithmetic.");
  m.def(
      "exp", [](TMC const& m) { return fadbad::exp(m); },
      "FadbadT overload: Taylor coefficients of exp(x) in McCormick "
      "arithmetic.");
  m.def(
      "log", [](TMC const& m) { return fadbad::log(m); },
      "FadbadT overload: Taylor coefficients of log(x) in McCormick "
      "arithmetic.");
  m.def(
      "sin", [](TMC const& m) { return fadbad::sin(m); },
      "FadbadT overload: Taylor coefficients of sin(x) in McCormick "
      "arithmetic.");
  m.def(
      "cos", [](TMC const& m) { return fadbad::cos(m); },
      "FadbadT overload: Taylor coefficients of cos(x) in McCormick "
      "arithmetic.");
  m.def(
      "tan", [](TMC const& m) { return fadbad::tan(m); },
      "FadbadT overload: Taylor coefficients of tan(x) in McCormick "
      "arithmetic.");
  // m.def( "cot",   []( TMC const& m ){ return fadbad::cot(m); } );
  m.def(
      "asin", [](TMC const& m) { return fadbad::asin(m); },
      "FadbadT overload: Taylor coefficients of asin(x) in McCormick "
      "arithmetic.");
  m.def(
      "acos", [](TMC const& m) { return fadbad::acos(m); },
      "FadbadT overload: Taylor coefficients of acos(x) in McCormick "
      "arithmetic.");
  m.def(
      "atan", [](TMC const& m) { return fadbad::atan(m); },
      "FadbadT overload: Taylor coefficients of atan(x) in McCormick "
      "arithmetic.");
  m.def(
      "sinh", [](TMC const& m) { return fadbad::sinh(m); },
      "FadbadT overload: Taylor coefficients of sinh(x) in McCormick "
      "arithmetic.");
  m.def(
      "cosh", [](TMC const& m) { return fadbad::cosh(m); },
      "FadbadT overload: Taylor coefficients of cosh(x) in McCormick "
      "arithmetic.");
  m.def(
      "tanh", [](TMC const& m) { return fadbad::tanh(m); },
      "FadbadT overload: Taylor coefficients of tanh(x) in McCormick "
      "arithmetic.");
  // m.def( "coth",  []( TMC const& m ){ return fadbad::coth(m); } );

  ;
}
