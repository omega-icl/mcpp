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

#include "pwlu.hpp"
#include "supmodel.hpp"
typedef mc::PWLU PWLU;
typedef mc::SupModel<mc::PWLU> PWLSM;
typedef mc::SupVar<mc::PWLU> PWLSV;

#include "pwcu.hpp"
typedef mc::PWCU PWCU;
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;

namespace py = pybind11;

void
mc_supmodel(py::module& m)
{
  py::class_<PWLSM> pyPWLSModel(m, "PWLSModel", R"doc(
Environment for superposition relaxations with piecewise-linear estimators.

A superposition relaxation brackets a multivariate function f(x_1,...,x_n)
between two separable functions: an underestimator sum_i fu_i(x_i) and an
overestimator sum_i fo_i(x_i), where each fu_i and fo_i is a univariate
function of a single variable. PWLSModel is the environment ("model") in
which such relaxations are propagated with piecewise-linear univariate
estimators on adaptive partitions (class PWLU) as the summands.

Typical workflow: create a PWLSModel with the number of participating
variables, create one PWLSVar per variable with its index and range, then
evaluate a factorable expression of these variables using the overloaded
arithmetic operators and the module-level functions (exp, log, sqr, ...).
The result is a PWLSVar holding guaranteed under- and overestimators of the
expression on the variable ranges.

The companion class PWCSModel implements the same arithmetic with
piecewise-constant univariate estimators (class PWCU) on fixed grid
partitions instead.

Examples
--------
>>> import pymcpp
>>> mod = pymcpp.PWLSModel(2)
>>> x = pymcpp.PWLSVar(mod, 0, pymcpp.Interval(1, 2), 8)
>>> y = pymcpp.PWLSVar(mod, 1, pymcpp.Interval(0, 1), 8)
>>> f = x * pymcpp.exp(x + y**2) - y**2
>>> f.l() <= f.u()  # guaranteed range of f on [1,2]x[0,1]
True
)doc");
  py::class_<PWLSM::Options> pyPWLSModelOptions(pyPWLSModel, "Options", R"doc(
Option set for PWLSModel.

Accessed through the ``options`` attribute of a PWLSModel instance; fields
may be assigned directly, e.g. ``mod.options.PROD_METH = mod.options.FULL``.
)doc");
  py::class_<PWLSV> pyPWLSVar(m, "PWLSVar", R"doc(
Superposition relaxation variable with piecewise-linear estimators.

A PWLSVar represents either a constant or an intermediate/final result of a
superposition relaxation computed in a PWLSModel environment. It stores, for
every participating variable index i, a piecewise-linear univariate
underestimator and overestimator (class PWLU) such that the sums
sum_i uest[i](x_i) and sum_i oest[i](x_i) bracket the represented function
on the variable ranges.

Standard arithmetic operators (+, -, *, /, **) are overloaded and propagate
the relaxation; module-level functions (exp, log, sqr, sqrt, sin, ...)
provide the intrinsic operations. Bounds are obtained with ``l()``/``u()``,
and the under-/overestimators can be evaluated at a point with
``uval()``/``oval()``.

When option ``USE_SHADOW`` is enabled in the model, an auxiliary "shadow"
under-/overestimator pair may be maintained alongside the primal pair to
tighten the relaxation; the ``opt`` argument of several methods selects
between primal (0), shadow (1) and best-of-both (2).

Examples
--------
>>> import pymcpp
>>> mod = pymcpp.PWLSModel(2)
>>> x = pymcpp.PWLSVar(mod, 0, pymcpp.Interval(-3, 3), 4)
>>> y = pymcpp.PWLSVar(mod, 1, pymcpp.Interval(-3, 3), 4)
>>> z = x * (pymcpp.exp(x) - y)**2
>>> z.uval({0: 0.5, 1: 0.5}) <= z.oval({0: 0.5, 1: 0.5})
True
)doc");
  py::class_<PWLU> pyPWLU(m, "PWLU", R"doc(
Piecewise-linear univariate function on an adaptive partition.

PWLU implements the univariate summands used by PWLSModel/PWLSVar
superposition relaxations. A PWLU stores a continuous piecewise-linear
function on an interval [xL, xU] as the initial point (``xL``, ``yL``)
together with the list of segment widths ``dx`` and segment slopes ``dy``.
Breakpoints are inserted adaptively as the relaxation is propagated.

Linear operations (+, -, scalar *, scalar /) are overloaded; ``min``/``max``
clip the function at a constant level; ``clean``, ``merge`` and ``reduce``
simplify the partition. Bounds over the whole range are available through
the properties ``l``, ``u`` and ``w``, and pointwise evaluation through
``lval``/``uval``. Since a PWLU holds a single piecewise-linear function,
``lval`` and ``uval`` coincide; both exist for interface symmetry with PWCU.

Examples
--------
>>> import pymcpp
>>> est = pymcpp.PWLU(0.0, 1.0, 4)  # identity function on [0,1], 4 segments
>>> est.l, est.u
(0.0, 1.0)
)doc");
  py::class_<PWLU::Options> pyPWLUOptions(pyPWLU, "Options", R"doc(
Static option set for PWLU.

Accessed through the class-level attribute ``PWLU.options`` and shared by
all PWLU instances.
)doc");

  pyPWLSModel
      .def(py::init<size_t const&>(), py::arg("nvar"), R"doc(
Construct a superposition model environment for `nvar` variables.

Parameters
----------
nvar : int
    Number of independent variables participating in the model. Each
    PWLSVar variable subsequently attached to this model must have an
    index in ``range(nvar)``.
)doc")
      .def_readwrite("options", &PWLSM::options, R"doc(
Option set of this model (a ``PWLSModel.Options`` instance).
)doc")
      .def_property_readonly(
          "nvar", [](PWLSM& self) { return self.nvar(); }, R"doc(
Number of independent variables in the model.
)doc")
      .def_property_readonly(
          "lbdvar", [](PWLSM& self) { return self.lbdvar(); }, R"doc(
Lower bounds of the independent variables, as a list of float of length
``nvar``. Entries are meaningful only for variables already defined in
the model.
)doc")
      .def_property_readonly(
          "ubdvar", [](PWLSM& self) { return self.ubdvar(); }, R"doc(
Upper bounds of the independent variables, as a list of float of length
``nvar``. Entries are meaningful only for variables already defined in
the model.
)doc")
      .def(
          "min", [](PWLSM& self, PWLSV& var, double const& minval)
          { return self.min(var, minval); }, py::arg("var"),
          py::arg("minval"), R"doc(
Intersect a superposition relaxation with the lower bound `minval`.

Cuts the relaxation stored in `var` off from below at level `minval`,
i.e. exploits the knowledge that the represented function is at least
`minval`. `var` is modified in place.

Parameters
----------
var : PWLSVar
    Variable whose relaxation is refined in place.
minval : float
    Known (natural) lower bound on the represented function.
)doc")
      .def(
          "max", [](PWLSM& self, PWLSV& var, double const& maxval)
          { return self.max(var, maxval); }, py::arg("var"),
          py::arg("maxval"), R"doc(
Intersect a superposition relaxation with the upper bound `maxval`.

Cuts the relaxation stored in `var` off from above at level `maxval`,
i.e. exploits the knowledge that the represented function is at most
`maxval`. `var` is modified in place.

Parameters
----------
var : PWLSVar
    Variable whose relaxation is refined in place.
maxval : float
    Known (natural) upper bound on the represented function.
)doc")
      .def(
          "uref", [](PWLSM& self, PWLSV& var, double const& lbd)
          { return self.min(var, lbd); }, py::arg("var"), py::arg("lbd"), R"doc(
Refine the superposition underestimator with a natural lower bound.

Tightens the underestimator of `var` so that it is consistent with the
known lower bound `lbd` of the represented function. `var` is modified
in place.

Parameters
----------
var : PWLSVar
    Variable whose underestimator is refined in place.
lbd : float
    Known (natural) lower bound on the represented function.
)doc")
      .def(
          "oref", [](PWLSM& self, PWLSV& var, double const& ubd)
          { return self.oref(var, ubd); }, py::arg("var"), py::arg("ubd"),
          R"doc(
Refine the superposition overestimator with a natural upper bound.

Tightens the overestimator of `var` so that it is consistent with the
known upper bound `ubd` of the represented function. `var` is modified
in place.

Parameters
----------
var : PWLSVar
    Variable whose overestimator is refined in place.
ubd : float
    Known (natural) upper bound on the represented function.
)doc");

  pyPWLSModelOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<PWLSM::Options const&>(), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("PROD_METH", &PWLSM::Options::PROD_METH, R"doc(
Reformulation method used for product terms, as a value of the
``PWLSModel.Options.PROD_REF`` enumeration (NONE, PARTIAL, FULL or LOG).
Default is PARTIAL (DC decomposition with range rescaling).
)doc")
      .def_readwrite("PROD_CUT", &PWLSM::Options::PROD_CUT, R"doc(
Whether to cut superposition relaxations of product terms at the interval
product bounds. Default is False.
)doc")
      .def_readwrite("SUM_TOL", &PWLSM::Options::SUM_TOL, R"doc(
Tolerance on range for univariate estimator propagation. Default is
1e2 * DBL_EPSILON (about 2.2e-14).
)doc")
      .def_readwrite("REF_WEIGHT", &PWLSM::Options::REF_WEIGHT, R"doc(
Weight in [0,1] used in the binary functions ``min`` and ``max``: 0
selects the first operand as reference, values >= 1 the second operand,
and intermediate values a weighted combination of both. Default is 0.5.
)doc")
      .def_readwrite("MAX_SUBDIV", &PWLSM::Options::MAX_SUBDIV, R"doc(
Maximal number of subdivisions for univariate estimators on adaptive
partitions (piecewise-linear estimators only). Default is 0 (no
restriction).
)doc")
      .def_readwrite("USE_SHADOW", &PWLSM::Options::USE_SHADOW, R"doc(
Whether to maintain auxiliary "shadow" estimators alongside the primal
estimators to tighten the relaxation. Default is False.
)doc")
      .def_readwrite("DISPLAY_SHADOW", &PWLSM::Options::DISPLAY_SHADOW, R"doc(
Whether to include shadow estimators when printing a variable. Default
is True.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &PWLSM::Options::DISPLAY_DIGITS, R"doc(
Number of significant digits used when printing a variable. Default is 5.
)doc")
      .def(
          "reset", [](PWLSM::Options& self) { return self.reset(); }, R"doc(
Reset all options to their default values.
)doc");

  py::enum_<PWLSM::Options::PROD_REF>(pyPWLSModelOptions, "PROD_REF", R"doc(
Reformulation methods for relaxing product terms in a superposition model.
)doc")
      .value("NONE", PWLSM::Options::PROD_REF::NONE,
             "DC (difference-of-convex) decomposition without rescaling.")
      .value("PARTIAL", PWLSM::Options::PROD_REF::PARTIAL,
             "DC decomposition with range rescaling.")
      .value("FULL", PWLSM::Options::PROD_REF::FULL,
             "DC decomposition with range and midpoint rescaling.")
      .value("LOG", PWLSM::Options::PROD_REF::LOG,
             "Log-transform with range and midpoint rescaling.")
      .export_values();

  pyPWLSVar
      .def(py::init<double const&>(), py::arg("cst") = 0., R"doc(
Construct a constant, not attached to any model.

Parameters
----------
cst : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<PWLSM&>(), py::arg("mod"), R"doc(
Construct a zero-valued variable attached to model `mod`.

Parameters
----------
mod : PWLSModel
    Superposition model environment.
)doc")
      .def(py::init<PWLSM&, unsigned int, I const&, size_t const>(),
           py::arg("mod"), py::arg("ndx"), py::arg("bnd"), py::arg("ndiv") = 1,
           R"doc(
Construct an independent variable of a superposition model.

Parameters
----------
mod : PWLSModel
    Superposition model environment.
ndx : int
    Index of the variable in the model, in ``range(mod.nvar)``.
bnd : Interval
    Range of the variable.
ndiv : int, optional
    Initial number of equal-size partition segments. Default is 1.

Raises
------
RuntimeError
    If `ndx` is not smaller than ``mod.nvar``.
)doc")
      .def(py::init<PWLSM&, unsigned int, PWLU const&>(), py::arg("mod"),
           py::arg("ndx"), py::arg("est"), R"doc(
Construct an independent variable from an existing univariate estimator.

Parameters
----------
mod : PWLSModel
    Superposition model environment.
ndx : int
    Index of the variable in the model, in ``range(mod.nvar)``.
est : PWLU
    Univariate estimator used for both the under- and overestimator of
    the variable; its abscissa range defines the variable range.
)doc")
      .def(py::init<PWLSV const&>(), R"doc(
Construct a copy of another variable.
)doc")
      .def(
          "set", [](PWLSV& self, PWLSM& mod) { return self.set(mod); },
          py::arg("mod"), R"doc(
Attach the variable to model `mod` and reset it to zero.

Parameters
----------
mod : PWLSModel
    Superposition model environment.

Returns
-------
var : PWLSVar
    The updated variable (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWLSV& self, PWLSM& mod, unsigned int ndx, I const& bnd,
             size_t const ndiv) { return self.set(mod, ndx, bnd, ndiv); },
          py::arg("mod"), py::arg("ndx"), py::arg("bnd"), py::arg("ndiv") = 1,
          R"doc(
Redefine the variable as independent variable `ndx` with range `bnd`.

Same semantics as the constructor ``PWLSVar(mod, ndx, bnd, ndiv)``; the
object is modified in place.

Returns
-------
var : PWLSVar
    The updated variable.
)doc")
      .def(
          "set", [](PWLSV& self, PWLSM& mod, unsigned int ndx, PWLU const& est)
          { return self.set(mod, ndx, est); }, py::arg("mod"), py::arg("ndx"),
          py::arg("est"), R"doc(
Redefine the variable as independent variable `ndx` with estimator `est`.

Same semantics as the constructor ``PWLSVar(mod, ndx, est)``; the object
is modified in place.

Returns
-------
var : PWLSVar
    The updated variable.
)doc")
      .def(
          "set", [](PWLSV& self, double const& cst) { return self.set(cst); },
          py::arg("cst"), R"doc(
Redefine the variable as the constant `cst`, detaching it from any model.

Returns
-------
var : PWLSVar
    The updated variable.
)doc")
      .def_property_readonly(
          "ndep", [](PWLSV& self) { return self.ndep(); }, R"doc(
Number of independent variables the relaxation depends on.
)doc")
      .def_property_readonly(
          "sdep", [](PWLSV& self) { return self.sdep(); }, R"doc(
Set of indices of the independent variables the relaxation depends on.
)doc")
      .def_property_readonly(
          "cst", [](PWLSV& self) { return self.cst(); }, R"doc(
Constant value of the variable; only meaningful when the variable does
not depend on any independent variable (``ndep == 0``).
)doc")
      .def(
          "uest", [](PWLSV& self, unsigned int opt) { return self.uest(opt); },
          py::arg("opt") = 0, R"doc(
Return the univariate summands of the superposition underestimator.

Parameters
----------
opt : int, optional
    0 for the primal underestimator (default); any nonzero value for the
    shadow underestimator (empty list unless shadow estimators are
    enabled and populated).

Returns
-------
uest : list of PWLU
    New list of the univariate underestimator summands, indexed by
    variable index; only entries whose index is in ``sdep`` are
    meaningful.
)doc")
      .def(
          "oest", [](PWLSV& self, unsigned int opt) { return self.oest(opt); },
          py::arg("opt") = 0, R"doc(
Return the univariate summands of the superposition overestimator.

Parameters
----------
opt : int, optional
    0 for the primal overestimator (default); any nonzero value for the
    shadow overestimator (empty list unless shadow estimators are
    enabled and populated).

Returns
-------
oest : list of PWLU
    New list of the univariate overestimator summands, indexed by
    variable index; only entries whose index is in ``sdep`` are
    meaningful.
)doc")
      .def(
          "l", [](PWLSV& self, unsigned int opt) { return self.l(opt); },
          py::arg("opt") = 2, R"doc(
Return a guaranteed lower bound of the superposition relaxation.

The bound is the minimum of the underestimator over the variable ranges.
For a constant variable, returns the constant value.

Parameters
----------
opt : int, optional
    0: primal underestimator only; 1: shadow underestimator only
    (-inf if absent); 2: best of primal and shadow (default); >= 3:
    best of both, recomputed from scratch.

Returns
-------
lb : float
    Lower bound on the represented function over the variable ranges.
)doc")
      .def(
          "u", [](PWLSV& self, unsigned int opt) { return self.u(opt); },
          py::arg("opt") = 2, R"doc(
Return a guaranteed upper bound of the superposition relaxation.

The bound is the maximum of the overestimator over the variable ranges.
For a constant variable, returns the constant value.

Parameters
----------
opt : int, optional
    0: primal overestimator only; 1: shadow overestimator only (+inf if
    absent); 2: best of primal and shadow (default); >= 3: best of
    both, recomputed from scratch.

Returns
-------
ub : float
    Upper bound on the represented function over the variable ranges.
)doc")
      .def(
          "uval",
          [](PWLSV& self, std::map<unsigned int, double> const& x,
             unsigned int opt) { return self.uval(x, opt); },
          py::arg("x"), py::arg("opt") = 2, R"doc(
Evaluate the superposition underestimator at a point.

Parameters
----------
x : dict of int to float
    Coordinates of the evaluation point, keyed by variable index; an
    entry is required for every index in ``sdep``.
opt : int, optional
    0: primal underestimator only; 1: shadow underestimator only; >= 2:
    best (largest) of primal and shadow (default).

Returns
-------
uval : float
    Value of the underestimator at `x`; guaranteed not to exceed the
    represented function value at `x`.
)doc")
      .def(
          "oval",
          [](PWLSV& self, std::map<unsigned int, double> const& x,
             unsigned int opt) { return self.oval(x, opt); },
          py::arg("x"), py::arg("opt") = 2, R"doc(
Evaluate the superposition overestimator at a point.

Parameters
----------
x : dict of int to float
    Coordinates of the evaluation point, keyed by variable index; an
    entry is required for every index in ``sdep``.
opt : int, optional
    0: primal overestimator only; 1: shadow overestimator only; >= 2:
    best (smallest) of primal and shadow (default).

Returns
-------
oval : float
    Value of the overestimator at `x`; guaranteed not to be below the
    represented function value at `x`.
)doc")
      //.def( "copy", []( PWLSV const& self ){ return PWLSV( self ); } )
      //.def( "assign", py::overload_cast<double const&>( &PWLSV::operator= ),
      // py::arg("val") ) .def( "assign", py::overload_cast<PWLSV const&>(
      //&PWLSV::operator= ), py::arg("var") )
      .def("__str__",
           [](PWLSV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](PWLSV const& self)
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
      // .def( "__abs__", []( PWLSV const& m ){ return mc::Op<PWLSV>::abs(m); }
      // )
      .def("__pow__", [](PWLSV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](PWLSV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](PWLSV const& m, PWLSV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, PWLSV const& m) { return mc::pow(r, m); })
      // .def( py::self == py::self )
      // .def( py::self != py::self )
      // .def( py::self <= py::self )
      // .def( py::self >= py::self )
      // .def( py::self < py::self )
      // .def( py::self > py::self )
      ;

  m.def(
      "inv", [](PWLSV const& x) { return mc::inv(x); },
      "PWLSVar overload: superposition relaxation of 1/x with "
      "piecewise-linear estimators.");
  m.def(
      "sqr", [](PWLSV const& x) { return mc::sqr(x); },
      "PWLSVar overload: superposition relaxation of x**2 with "
      "piecewise-linear estimators.");
  m.def(
      "sqrt", [](PWLSV const& x) { return mc::sqrt(x); },
      "PWLSVar overload: superposition relaxation of sqrt(x) with "
      "piecewise-linear estimators.");
  m.def(
      "exp", [](PWLSV const& x) { return mc::exp(x); },
      "PWLSVar overload: superposition relaxation of exp(x) with "
      "piecewise-linear estimators.");
  m.def(
      "log", [](PWLSV const& x) { return mc::log(x); },
      "PWLSVar overload: superposition relaxation of log(x) with "
      "piecewise-linear estimators.");
  m.def(
      "cos", [](PWLSV const& x) { return mc::cos(x); },
      "PWLSVar overload: superposition relaxation of cos(x) with "
      "piecewise-linear estimators.");
  m.def(
      "sin", [](PWLSV const& x) { return mc::sin(x); },
      "PWLSVar overload: superposition relaxation of sin(x) with "
      "piecewise-linear estimators.");
  m.def(
      "tan", [](PWLSV const& x) { return mc::tan(x); },
      "PWLSVar overload: superposition relaxation of tan(x) with "
      "piecewise-linear estimators.");
  m.def(
      "acos", [](PWLSV const& x) { return mc::acos(x); },
      "PWLSVar overload: superposition relaxation of acos(x) with "
      "piecewise-linear estimators.");
  m.def(
      "asin", [](PWLSV const& x) { return mc::asin(x); },
      "PWLSVar overload: superposition relaxation of asin(x) with "
      "piecewise-linear estimators.");
  m.def(
      "atan", [](PWLSV const& x) { return mc::atan(x); },
      "PWLSVar overload: superposition relaxation of atan(x) with "
      "piecewise-linear estimators.");
  m.def(
      "cosh", [](PWLSV const& x) { return mc::cosh(x); },
      "PWLSVar overload: superposition relaxation of cosh(x) with "
      "piecewise-linear estimators.");
  m.def(
      "sinh", [](PWLSV const& x) { return mc::sinh(x); },
      "PWLSVar overload: superposition relaxation of sinh(x) with "
      "piecewise-linear estimators.");
  m.def(
      "tanh", [](PWLSV const& x) { return mc::tanh(x); },
      "PWLSVar overload: superposition relaxation of tanh(x) with "
      "piecewise-linear estimators.");
  m.def(
      "fabs", [](PWLSV const& x) { return mc::fabs(x); },
      "PWLSVar overload: superposition relaxation of abs(x) with "
      "piecewise-linear estimators.");
  m.def(
      "relu", [](PWLSV const& x) { return mc::relu(x); },
      "PWLSVar overload: superposition relaxation of max(x, 0) with "
      "piecewise-linear estimators.");
  m.def(
      "xlog", [](PWLSV const& x) { return mc::xlog(x); },
      "PWLSVar overload: superposition relaxation of x*log(x) with "
      "piecewise-linear estimators.");
  m.def(
      "erf", [](PWLSV const& x) { return mc::erf(x); },
      "PWLSVar overload: superposition relaxation of erf(x) with "
      "piecewise-linear estimators.");
  m.def(
      "erfc", [](PWLSV const& x) { return mc::erfc(x); },
      "PWLSVar overload: superposition relaxation of erfc(x) with "
      "piecewise-linear estimators.");
  m.def(
      "pow", [](PWLSV const& x, int const n) { return mc::pow(x, n); },
      "PWLSVar overload: superposition relaxation of x**n for integer n.");
  m.def(
      "pow", [](PWLSV const& x, double const& r) { return mc::pow(x, r); },
      "PWLSVar overload: superposition relaxation of x**r for real r.");
  m.def(
      "pow", [](PWLSV const& x, PWLSV const& y) { return mc::pow(x, y); },
      "PWLSVar overload: superposition relaxation of x**y for PWLSVar "
      "exponent y.");
  m.def(
      "pow", [](double const& r, PWLSV const& y) { return mc::pow(r, y); },
      "PWLSVar overload: superposition relaxation of r**y for real base r.");
  m.def(
      "cheb", [](PWLSV const& x, unsigned const n) { return mc::cheb(x, n); },
      "PWLSVar overload: superposition relaxation of the Chebyshev "
      "polynomial T_n(x).");
  m.def(
      "max", [](PWLSV const& x, double const& y) { return mc::max(x, y); },
      "PWLSVar overload: superposition relaxation of max(x, y) with scalar "
      "y.");
  m.def(
      "min", [](PWLSV const& x, double const& y) { return mc::min(x, y); },
      "PWLSVar overload: superposition relaxation of min(x, y) with scalar "
      "y.");
  m.def(
      "max", [](PWLSV const& x, PWLSV const& y) { return mc::max(x, y); },
      "PWLSVar overload: superposition relaxation of max(x, y); the "
      "reference operand is selected via option REF_WEIGHT.");
  m.def(
      "min", [](PWLSV const& x, PWLSV const& y) { return mc::min(x, y); },
      "PWLSVar overload: superposition relaxation of min(x, y); the "
      "reference operand is selected via option REF_WEIGHT.");

  pyPWLU
      .def(py::init<>(), R"doc(
Construct an empty estimator (no range, no segments).
)doc")
      .def(py::init<double const&, double const&, double const&>(),
           py::arg("xL"), py::arg("xU"), py::arg("y"), R"doc(
Construct the constant function y on the interval [xL, xU].

Parameters
----------
xL : float
    Lower abscissa of the range.
xU : float
    Upper abscissa of the range.
y : float
    Constant function value.
)doc")
      .def(py::init<double const&, double const&, size_t const>(),
           py::arg("xL"), py::arg("xU"), py::arg("N") = 1, R"doc(
Construct the identity function x on [xL, xU] with N equal segments.

This is the usual starting estimator for an independent variable.

Parameters
----------
xL : float
    Lower abscissa of the range.
xU : float
    Upper abscissa of the range.
N : int, optional
    Number of equal-width partition segments. Default is 1.
)doc")
      .def(py::init<double const&, std::vector<double> const&>(), py::arg("xL"),
           py::arg("dx"), R"doc(
Construct the identity function x on a custom partition.

Parameters
----------
xL : float
    Lower abscissa of the range.
dx : list of float
    Positive widths of the partition segments; the range upper bound is
    ``xL + sum(dx)``.
)doc")
      .def(py::init<PWLU const&>(), R"doc(
Construct a copy of another estimator.
)doc")
      .def_readwrite_static("options", &PWLU::options, R"doc(
Class-level option set (a ``PWLU.Options`` instance) shared by all PWLU
instances.
)doc")
      .def(
          "set", [](PWLU& self, double const& y) { return self.set(y); },
          py::arg("y"), R"doc(
Reset the estimator to the constant function y on its current range.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWLU& self, double const& xL, double const& xU, double const& y)
          { return self.set(xL, xU, y); }, py::arg("xL"), py::arg("xU"),
          py::arg("y"), R"doc(
Reset the estimator to the constant function y on [xL, xU].

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWLU& self, double const& xL, double const& xU, size_t const N)
          { return self.set(xL, xU, N); }, py::arg("xL"), py::arg("xU"),
          py::arg("N") = 1, R"doc(
Reset the estimator to the identity on [xL, xU] with N equal segments.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "set", [](PWLU& self, double const& xL, std::vector<double> const& dx)
          { return self.set(xL, dx); }, py::arg("xL"), py::arg("dx"), R"doc(
Reset the estimator to the identity on a custom partition.

Parameters
----------
xL : float
    Lower abscissa of the range.
dx : list of float
    Positive widths of the partition segments.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "insert", [](PWLU& self, std::vector<double> const& dx)
          { return self.insert(dx); }, py::arg("dx"), R"doc(
Insert several breakpoints into the partition.

Breakpoints closer than the tolerances ``options.BKPTATOL`` /
``options.BKPTRTOL`` to an existing breakpoint are skipped. The function
values are unchanged.

Parameters
----------
dx : list of float
    Abscissae of the breakpoints to insert.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "insert", [](PWLU& self, double const& x) { return self.insert(x); },
          py::arg("x"), R"doc(
Insert a single breakpoint at abscissa x, if not already present.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def_property_readonly(
          "l", [](PWLU& self) { return self.l(); }, R"doc(
Minimum of the piecewise-linear function over its range.
)doc")
      .def_property_readonly(
          "u", [](PWLU& self) { return self.u(); }, R"doc(
Maximum of the piecewise-linear function over its range.
)doc")
      .def_property_readonly(
          "w", [](PWLU& self) { return self.w(); }, R"doc(
Width of the function range, equal to ``u - l``.
)doc")
      .def(
          "lval", [](PWLU& self, double const& x) { return self.l(x); },
          py::arg("x"), R"doc(
Evaluate the piecewise-linear function at abscissa x.

Equivalent to ``uval`` since a PWLU stores a single piecewise-linear
function.

Parameters
----------
x : float
    Evaluation point; must lie within the estimator range.

Returns
-------
y : float
    Function value at x.

Raises
------
RuntimeError
    If x lies outside the estimator range.
)doc")
      .def(
          "uval", [](PWLU& self, double const& x) { return self.l(x); },
          py::arg("x"), R"doc(
Evaluate the piecewise-linear function at abscissa x.

Equivalent to ``lval`` since a PWLU stores a single piecewise-linear
function.

Parameters
----------
x : float
    Evaluation point; must lie within the estimator range.

Returns
-------
y : float
    Function value at x.
)doc")
      .def_property(
          "xL", [](PWLU& self) { return self.xL(); },
          [](PWLU& self, double const& x) { self.xL() = x; }, R"doc(
Lower abscissa of the estimator range (read/write).
)doc")
      .def_property(
          "yL", [](PWLU& self) { return self.yL(); },
          [](PWLU& self, double const& y) { self.yL() = y; }, R"doc(
Function value at the lower abscissa ``xL`` (read/write).
)doc")
      .def_property(
          "dx", [](PWLU& self) { return self.dx(); },
          [](PWLU& self, std::vector<double> const& dx) { self.dx() = dx; },
          R"doc(
Widths of the partition segments, as a list of float (read/write).
)doc")
      .def_property(
          "dy", [](PWLU& self) { return self.dy(); },
          [](PWLU& self, std::vector<double> const& dy) { self.dy() = dy; },
          R"doc(
Slopes of the function on the partition segments, as a list of float
(read/write); same length as ``dx``.
)doc")
      .def("__str__",
           [](PWLU const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](PWLU const& self)
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
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self /= double())
      .def(py::self / double())
      .def(
          "min", [](PWLU& self, double const& cutoff)
          { return self.min(cutoff); }, py::arg("cutoff"), R"doc(
Clip the function from above at level `cutoff` (pointwise minimum).

Replaces the function f(x) by min(f(x), cutoff), inserting breakpoints
where the function crosses the cutoff level.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "max", [](PWLU& self, double const& cutoff)
          { return self.max(cutoff); }, py::arg("cutoff"), R"doc(
Clip the function from below at level `cutoff` (pointwise maximum).

Replaces the function f(x) by max(f(x), cutoff), inserting breakpoints
where the function crosses the cutoff level.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "clean",
          [](PWLU& self, bool const under) { return self.clean(under); },
          py::arg("under"), R"doc(
Remove partition segments narrower than the breakpoint tolerances.

Segments with width below ``options.BKPTATOL`` / ``options.BKPTRTOL``
are merged into their neighbors, rounding in the direction that keeps
the function a valid estimator.

Parameters
----------
under : bool
    True if the function is used as an underestimator (merging rounds
    downward), False for an overestimator (merging rounds upward).

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "merge",
          [](PWLU& self, bool const under) { return self.clean(under); },
          py::arg("under"), R"doc(
Merge adjacent partition segments with identical slope.

Slopes equal within the tolerances ``options.BKPTATOL`` /
``options.BKPTRTOL`` are considered identical; merging rounds in the
direction that keeps the function a valid estimator.

Parameters
----------
under : bool
    True if the function is used as an underestimator, False for an
    overestimator.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "reduce", [](PWLU& self, bool const under, size_t const nseg)
          { return self.clean(under); }, py::arg("under"), py::arg("nseg"),
          R"doc(
Reduce the partition to at most `nseg` segments by successive relaxation.

Breakpoints whose removal loses the least tightness are eliminated
first, so the result remains a valid (weaker) estimator.

Parameters
----------
under : bool
    True if the function is used as an underestimator, False for an
    overestimator.
nseg : int
    Target maximal number of partition segments; 0 means no reduction.

Returns
-------
est : PWLU
    The updated estimator (the object itself is modified in place).
)doc");

  pyPWLUOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<PWLU::Options const&>(), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("BKPTATOL", &PWLU::Options::BKPTATOL, R"doc(
Absolute tolerance for managing (inserting, merging) breakpoints.
Default is 1e4 * DBL_EPSILON (about 2.2e-12).
)doc")
      .def_readwrite("BKPTRTOL", &PWLU::Options::BKPTRTOL, R"doc(
Relative tolerance for managing (inserting, merging) breakpoints.
Default is 1e4 * DBL_EPSILON (about 2.2e-12).
)doc")
      .def_readwrite("REDUCEMETH", &PWLU::Options::REDUCEMETH, R"doc(
Method for breakpoint reduction: <= 0 uses the standard-library
implementation (default 0); > 0 uses a tailored heap implementation
whenever the number of segments exceeds the target by this threshold.
)doc")
      .def_readwrite("DISPNUM", &PWLU::Options::DISPNUM, R"doc(
Number of significant digits used when printing an estimator. Default
is 5.
)doc")
      .def(
          "reset", [](PWLU::Options& self) { return self.reset(); }, R"doc(
Reset all options to their default values.
)doc");

  py::class_<PWCSM> pyPWCSModel(m, "PWCSModel", R"doc(
Environment for superposition relaxations with piecewise-constant
estimators.

Same superposition relaxation arithmetic as PWLSModel, but the univariate
under- and overestimator summands are piecewise-constant functions on fixed
grid partitions (class PWCU) instead of piecewise-linear functions on
adaptive partitions (class PWLU).

See Also
--------
PWLSModel : piecewise-linear variant, including a workflow description.

Examples
--------
>>> import pymcpp
>>> mod = pymcpp.PWCSModel(2)
>>> x = pymcpp.PWCSVar(mod, 0, pymcpp.Interval(1, 2), 8)
>>> y = pymcpp.PWCSVar(mod, 1, pymcpp.Interval(0, 1), 8)
>>> f = x * pymcpp.exp(x + y**2) - y**2
>>> f.l() <= f.u()
True
)doc");
  py::class_<PWCSM::Options> pyPWCSModelOptions(pyPWCSModel, "Options", R"doc(
Option set for PWCSModel.

Accessed through the ``options`` attribute of a PWCSModel instance.
)doc");
  py::class_<PWCSV> pyPWCSVar(m, "PWCSVar", R"doc(
Superposition relaxation variable with piecewise-constant estimators.

Same role and interface as PWLSVar, but the univariate summands are
piecewise-constant estimators (class PWCU) on fixed grid partitions.
Standard arithmetic operators (+, -, *, /, **) are overloaded and propagate
the relaxation; module-level functions (exp, log, sqr, ...) provide the
intrinsic operations.

See Also
--------
PWLSVar : piecewise-linear variant, with a full description of the
    accessors ``l``, ``u``, ``uval``, ``oval``, ``uest``, ``oest``.
)doc");
  py::class_<PWCU> pyPWCU(m, "PWCU", R"doc(
Piecewise-constant univariate estimator pair on a fixed grid partition.

PWCU implements the univariate summands used by PWCSModel/PWCSVar
superposition relaxations. On an interval [xL, xU] divided into N segments
of equal width, a PWCU stores for each segment a lower level (``yL``) and
an upper level (``yU``), so that a single PWCU represents a
piecewise-constant under/over bracket of a univariate function. Slopes may
additionally be propagated internally to sharpen the bracket.

Linear operations (+, -, scalar *, scalar /) are overloaded; ``min``/``max``
clip the estimator at a constant level. Bounds over the whole range are
available through the properties ``l``, ``u`` and ``w``, and pointwise
evaluation through ``lval``/``uval``.
)doc");
  py::class_<PWCU::Options> pyPWCUOptions(pyPWCU, "Options", R"doc(
Static option set for PWCU.

Accessed through the class-level attribute ``PWCU.options`` and shared by
all PWCU instances.
)doc");

  pyPWCSModel
      .def(py::init<size_t const&>(), py::arg("nvar"), R"doc(
Construct a superposition model environment for `nvar` variables.

Parameters
----------
nvar : int
    Number of independent variables participating in the model. Each
    PWCSVar variable subsequently attached to this model must have an
    index in ``range(nvar)``.
)doc")
      .def_readwrite("options", &PWCSM::options, R"doc(
Option set of this model (a ``PWCSModel.Options`` instance).
)doc")
      .def_property_readonly(
          "nvar", [](PWCSM& self) { return self.nvar(); }, R"doc(
Number of independent variables in the model.
)doc")
      .def_property_readonly(
          "lbdvar", [](PWCSM& self) { return self.lbdvar(); }, R"doc(
Lower bounds of the independent variables, as a list of float of length
``nvar``. Entries are meaningful only for variables already defined in
the model.
)doc")
      .def_property_readonly(
          "ubdvar", [](PWCSM& self) { return self.ubdvar(); }, R"doc(
Upper bounds of the independent variables, as a list of float of length
``nvar``. Entries are meaningful only for variables already defined in
the model.
)doc")
      .def(
          "min", [](PWCSM& self, PWCSV& var, double const& minval)
          { return self.min(var, minval); }, py::arg("var"),
          py::arg("minval"), R"doc(
Intersect a superposition relaxation with the lower bound `minval`.

Cuts the relaxation stored in `var` off from below at level `minval`.
`var` is modified in place.

Parameters
----------
var : PWCSVar
    Variable whose relaxation is refined in place.
minval : float
    Known (natural) lower bound on the represented function.
)doc")
      .def(
          "max", [](PWCSM& self, PWCSV& var, double const& maxval)
          { return self.max(var, maxval); }, py::arg("var"),
          py::arg("maxval"), R"doc(
Intersect a superposition relaxation with the upper bound `maxval`.

Cuts the relaxation stored in `var` off from above at level `maxval`.
`var` is modified in place.

Parameters
----------
var : PWCSVar
    Variable whose relaxation is refined in place.
maxval : float
    Known (natural) upper bound on the represented function.
)doc")
      .def(
          "uref", [](PWCSM& self, PWCSV& var, double const& lbd)
          { return self.min(var, lbd); }, py::arg("var"), py::arg("lbd"), R"doc(
Refine the superposition underestimator with a natural lower bound.

Tightens the underestimator of `var` so that it is consistent with the
known lower bound `lbd` of the represented function. `var` is modified
in place.

Parameters
----------
var : PWCSVar
    Variable whose underestimator is refined in place.
lbd : float
    Known (natural) lower bound on the represented function.
)doc")
      .def(
          "oref", [](PWCSM& self, PWCSV& var, double const& ubd)
          { return self.oref(var, ubd); }, py::arg("var"), py::arg("ubd"),
          R"doc(
Refine the superposition overestimator with a natural upper bound.

Tightens the overestimator of `var` so that it is consistent with the
known upper bound `ubd` of the represented function. `var` is modified
in place.

Parameters
----------
var : PWCSVar
    Variable whose overestimator is refined in place.
ubd : float
    Known (natural) upper bound on the represented function.
)doc");

  pyPWCSModelOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<PWCSM::Options const&>(), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("PROD_METH", &PWCSM::Options::PROD_METH, R"doc(
Reformulation method used for product terms, as a value of the
``PWCSModel.Options.PROD_REF`` enumeration (NONE, PARTIAL, FULL or LOG).
Default is PARTIAL (DC decomposition with range rescaling).
)doc")
      .def_readwrite("PROD_CUT", &PWCSM::Options::PROD_CUT, R"doc(
Whether to cut superposition relaxations of product terms at the interval
product bounds. Default is False.
)doc")
      .def_readwrite("SUM_TOL", &PWCSM::Options::SUM_TOL, R"doc(
Tolerance on range for univariate estimator propagation. Default is
1e2 * DBL_EPSILON (about 2.2e-14).
)doc")
      .def_readwrite("REF_WEIGHT", &PWCSM::Options::REF_WEIGHT, R"doc(
Weight in [0,1] used in the binary functions ``min`` and ``max``: 0
selects the first operand as reference, values >= 1 the second operand,
and intermediate values a weighted combination of both. Default is 0.5.
)doc")
      .def_readwrite("MAX_SUBDIV", &PWCSM::Options::MAX_SUBDIV, R"doc(
Maximal number of subdivisions for univariate estimators on adaptive
partitions; not relevant for the fixed grid partitions of PWCU. Default
is 0 (no restriction).
)doc")
      .def_readwrite("USE_SHADOW", &PWCSM::Options::USE_SHADOW, R"doc(
Whether to maintain auxiliary "shadow" estimators alongside the primal
estimators to tighten the relaxation. Default is False.
)doc")
      .def_readwrite("DISPLAY_SHADOW", &PWCSM::Options::DISPLAY_SHADOW, R"doc(
Whether to include shadow estimators when printing a variable. Default
is True.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &PWCSM::Options::DISPLAY_DIGITS, R"doc(
Number of significant digits used when printing a variable. Default is 5.
)doc")
      .def(
          "reset", [](PWCSM::Options& self) { return self.reset(); }, R"doc(
Reset all options to their default values.
)doc");

  py::enum_<PWCSM::Options::PROD_REF>(pyPWCSModelOptions, "PROD_REF", R"doc(
Reformulation methods for relaxing product terms in a superposition model.
)doc")
      .value("NONE", PWCSM::Options::PROD_REF::NONE,
             "DC (difference-of-convex) decomposition without rescaling.")
      .value("PARTIAL", PWCSM::Options::PROD_REF::PARTIAL,
             "DC decomposition with range rescaling.")
      .value("FULL", PWCSM::Options::PROD_REF::FULL,
             "DC decomposition with range and midpoint rescaling.")
      .value("LOG", PWCSM::Options::PROD_REF::LOG,
             "Log-transform with range and midpoint rescaling.")
      .export_values();

  pyPWCSVar
      .def(py::init<double const&>(), py::arg("cst") = 0., R"doc(
Construct a constant, not attached to any model.

Parameters
----------
cst : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<PWCSM&>(), py::arg("mod"), R"doc(
Construct a zero-valued variable attached to model `mod`.

Parameters
----------
mod : PWCSModel
    Superposition model environment.
)doc")
      .def(py::init<PWCSM&, unsigned int, I const&, size_t const>(),
           py::arg("mod"), py::arg("ndx"), py::arg("bnd"), py::arg("ndiv") = 1,
           R"doc(
Construct an independent variable of a superposition model.

Parameters
----------
mod : PWCSModel
    Superposition model environment.
ndx : int
    Index of the variable in the model, in ``range(mod.nvar)``.
bnd : Interval
    Range of the variable.
ndiv : int, optional
    Number of equal-size partition segments. Default is 1.

Raises
------
RuntimeError
    If `ndx` is not smaller than ``mod.nvar``.
)doc")
      .def(py::init<PWCSM&, unsigned int, PWCU const&>(), py::arg("mod"),
           py::arg("ndx"), py::arg("est"), R"doc(
Construct an independent variable from an existing univariate estimator.

Parameters
----------
mod : PWCSModel
    Superposition model environment.
ndx : int
    Index of the variable in the model, in ``range(mod.nvar)``.
est : PWCU
    Univariate estimator used for both the under- and overestimator of
    the variable; its abscissa range defines the variable range.
)doc")
      .def(py::init<PWCSV const&>(), R"doc(
Construct a copy of another variable.
)doc")
      .def(
          "set", [](PWCSV& self, PWCSM& mod) { return self.set(mod); },
          py::arg("mod"), R"doc(
Attach the variable to model `mod` and reset it to zero.

Parameters
----------
mod : PWCSModel
    Superposition model environment.

Returns
-------
var : PWCSVar
    The updated variable (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWCSV& self, PWCSM& mod, unsigned int ndx, I const& bnd,
             size_t const ndiv) { return self.set(mod, ndx, bnd, ndiv); },
          py::arg("mod"), py::arg("ndx"), py::arg("bnd"), py::arg("ndiv") = 1,
          R"doc(
Redefine the variable as independent variable `ndx` with range `bnd`.

Same semantics as the constructor ``PWCSVar(mod, ndx, bnd, ndiv)``; the
object is modified in place.

Returns
-------
var : PWCSVar
    The updated variable.
)doc")
      .def(
          "set", [](PWCSV& self, PWCSM& mod, unsigned int ndx, PWCU const& est)
          { return self.set(mod, ndx, est); }, py::arg("mod"), py::arg("ndx"),
          py::arg("est"), R"doc(
Redefine the variable as independent variable `ndx` with estimator `est`.

Same semantics as the constructor ``PWCSVar(mod, ndx, est)``; the object
is modified in place.

Returns
-------
var : PWCSVar
    The updated variable.
)doc")
      .def(
          "set", [](PWCSV& self, double const& cst) { return self.set(cst); },
          py::arg("cst"), R"doc(
Redefine the variable as the constant `cst`, detaching it from any model.

Returns
-------
var : PWCSVar
    The updated variable.
)doc")
      .def_property_readonly(
          "ndep", [](PWCSV& self) { return self.ndep(); }, R"doc(
Number of independent variables the relaxation depends on.
)doc")
      .def_property_readonly(
          "sdep", [](PWCSV& self) { return self.sdep(); }, R"doc(
Set of indices of the independent variables the relaxation depends on.
)doc")
      .def_property_readonly(
          "cst", [](PWCSV& self) { return self.cst(); }, R"doc(
Constant value of the variable; only meaningful when the variable does
not depend on any independent variable (``ndep == 0``).
)doc")
      .def(
          "uest", [](PWCSV& self, unsigned int opt) { return self.uest(opt); },
          py::arg("opt") = 0, R"doc(
Return the univariate summands of the superposition underestimator.

Parameters
----------
opt : int, optional
    0 for the primal underestimator (default); any nonzero value for
    the shadow underestimator (empty list unless shadow estimators are
    enabled and populated).

Returns
-------
uest : list of PWCU
    New list of the univariate underestimator summands, indexed by
    variable index; only entries whose index is in ``sdep`` are
    meaningful.
)doc")
      .def(
          "oest", [](PWCSV& self, unsigned int opt) { return self.oest(opt); },
          py::arg("opt") = 0, R"doc(
Return the univariate summands of the superposition overestimator.

Parameters
----------
opt : int, optional
    0 for the primal overestimator (default); any nonzero value for the
    shadow overestimator (empty list unless shadow estimators are
    enabled and populated).

Returns
-------
oest : list of PWCU
    New list of the univariate overestimator summands, indexed by
    variable index; only entries whose index is in ``sdep`` are
    meaningful.
)doc")
      .def(
          "l", [](PWCSV& self, unsigned int opt) { return self.l(opt); },
          py::arg("opt") = 2, R"doc(
Return a guaranteed lower bound of the superposition relaxation.

Parameters
----------
opt : int, optional
    0: primal underestimator only; 1: shadow underestimator only
    (-inf if absent); 2: best of primal and shadow (default); >= 3:
    best of both, recomputed from scratch.

Returns
-------
lb : float
    Lower bound on the represented function over the variable ranges.
)doc")
      .def(
          "u", [](PWCSV& self, unsigned int opt) { return self.u(opt); },
          py::arg("opt") = 2, R"doc(
Return a guaranteed upper bound of the superposition relaxation.

Parameters
----------
opt : int, optional
    0: primal overestimator only; 1: shadow overestimator only (+inf if
    absent); 2: best of primal and shadow (default); >= 3: best of
    both, recomputed from scratch.

Returns
-------
ub : float
    Upper bound on the represented function over the variable ranges.
)doc")
      .def(
          "uval",
          [](PWCSV& self, std::map<unsigned int, double> const& x,
             unsigned int opt) { return self.uval(x, opt); },
          py::arg("x"), py::arg("opt") = 2, R"doc(
Evaluate the superposition underestimator at a point.

Parameters
----------
x : dict of int to float
    Coordinates of the evaluation point, keyed by variable index; an
    entry is required for every index in ``sdep``.
opt : int, optional
    0: primal underestimator only; 1: shadow underestimator only; >= 2:
    best (largest) of primal and shadow (default).

Returns
-------
uval : float
    Value of the underestimator at `x`; guaranteed not to exceed the
    represented function value at `x`.
)doc")
      .def(
          "oval",
          [](PWCSV& self, std::map<unsigned int, double> const& x,
             unsigned int opt) { return self.oval(x, opt); },
          py::arg("x"), py::arg("opt") = 2, R"doc(
Evaluate the superposition overestimator at a point.

Parameters
----------
x : dict of int to float
    Coordinates of the evaluation point, keyed by variable index; an
    entry is required for every index in ``sdep``.
opt : int, optional
    0: primal overestimator only; 1: shadow overestimator only; >= 2:
    best (smallest) of primal and shadow (default).

Returns
-------
oval : float
    Value of the overestimator at `x`; guaranteed not to be below the
    represented function value at `x`.
)doc")
      //.def( "copy", []( PWCSV const& self ){ return PWCSV( self ); } )
      //.def( "assign", py::overload_cast<double const&>( &PWCSV::operator= ),
      // py::arg("val") ) .def( "assign", py::overload_cast<PWCSV const&>(
      //&PWCSV::operator= ), py::arg("var") )
      .def("__str__",
           [](PWCSV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](PWCSV const& self)
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
      // .def( "__abs__", []( PWCSV const& m ){ return mc::Op<PWCSV>::abs(m); }
      // )
      .def("__pow__", [](PWCSV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](PWCSV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](PWCSV const& m, PWCSV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, PWCSV const& m) { return mc::pow(r, m); })
      // .def( py::self == py::self )
      // .def( py::self != py::self )
      // .def( py::self <= py::self )
      // .def( py::self >= py::self )
      // .def( py::self < py::self )
      // .def( py::self > py::self )
      ;

  m.def(
      "inv", [](PWCSV const& x) { return mc::inv(x); },
      "PWCSVar overload: superposition relaxation of 1/x with "
      "piecewise-constant estimators.");
  m.def(
      "sqr", [](PWCSV const& x) { return mc::sqr(x); },
      "PWCSVar overload: superposition relaxation of x**2 with "
      "piecewise-constant estimators.");
  m.def(
      "sqrt", [](PWCSV const& x) { return mc::sqrt(x); },
      "PWCSVar overload: superposition relaxation of sqrt(x) with "
      "piecewise-constant estimators.");
  m.def(
      "exp", [](PWCSV const& x) { return mc::exp(x); },
      "PWCSVar overload: superposition relaxation of exp(x) with "
      "piecewise-constant estimators.");
  m.def(
      "log", [](PWCSV const& x) { return mc::log(x); },
      "PWCSVar overload: superposition relaxation of log(x) with "
      "piecewise-constant estimators.");
  m.def(
      "cos", [](PWCSV const& x) { return mc::cos(x); },
      "PWCSVar overload: superposition relaxation of cos(x) with "
      "piecewise-constant estimators.");
  m.def(
      "sin", [](PWCSV const& x) { return mc::sin(x); },
      "PWCSVar overload: superposition relaxation of sin(x) with "
      "piecewise-constant estimators.");
  m.def(
      "tan", [](PWCSV const& x) { return mc::tan(x); },
      "PWCSVar overload: superposition relaxation of tan(x) with "
      "piecewise-constant estimators.");
  m.def(
      "acos", [](PWCSV const& x) { return mc::acos(x); },
      "PWCSVar overload: superposition relaxation of acos(x) with "
      "piecewise-constant estimators.");
  m.def(
      "asin", [](PWCSV const& x) { return mc::asin(x); },
      "PWCSVar overload: superposition relaxation of asin(x) with "
      "piecewise-constant estimators.");
  m.def(
      "atan", [](PWCSV const& x) { return mc::atan(x); },
      "PWCSVar overload: superposition relaxation of atan(x) with "
      "piecewise-constant estimators.");
  m.def(
      "cosh", [](PWCSV const& x) { return mc::cosh(x); },
      "PWCSVar overload: superposition relaxation of cosh(x) with "
      "piecewise-constant estimators.");
  m.def(
      "sinh", [](PWCSV const& x) { return mc::sinh(x); },
      "PWCSVar overload: superposition relaxation of sinh(x) with "
      "piecewise-constant estimators.");
  m.def(
      "tanh", [](PWCSV const& x) { return mc::tanh(x); },
      "PWCSVar overload: superposition relaxation of tanh(x) with "
      "piecewise-constant estimators.");
  m.def(
      "fabs", [](PWCSV const& x) { return mc::fabs(x); },
      "PWCSVar overload: superposition relaxation of abs(x) with "
      "piecewise-constant estimators.");
  m.def(
      "relu", [](PWCSV const& x) { return mc::relu(x); },
      "PWCSVar overload: superposition relaxation of max(x, 0) with "
      "piecewise-constant estimators.");
  m.def(
      "xlog", [](PWCSV const& x) { return mc::xlog(x); },
      "PWCSVar overload: superposition relaxation of x*log(x) with "
      "piecewise-constant estimators.");
  m.def(
      "erf", [](PWCSV const& x) { return mc::erf(x); },
      "PWCSVar overload: superposition relaxation of erf(x) with "
      "piecewise-constant estimators.");
  m.def(
      "erfc", [](PWCSV const& x) { return mc::erfc(x); },
      "PWCSVar overload: superposition relaxation of erfc(x) with "
      "piecewise-constant estimators.");
  m.def(
      "pow", [](PWCSV const& x, int const n) { return mc::pow(x, n); },
      "PWCSVar overload: superposition relaxation of x**n for integer n.");
  m.def(
      "pow", [](PWCSV const& x, double const& r) { return mc::pow(x, r); },
      "PWCSVar overload: superposition relaxation of x**r for real r.");
  m.def(
      "pow", [](PWCSV const& x, PWCSV const& y) { return mc::pow(x, y); },
      "PWCSVar overload: superposition relaxation of x**y for PWCSVar "
      "exponent y.");
  m.def(
      "pow", [](double const& r, PWCSV const& y) { return mc::pow(r, y); },
      "PWCSVar overload: superposition relaxation of r**y for real base r.");
  m.def(
      "cheb", [](PWCSV const& x, unsigned const n) { return mc::cheb(x, n); },
      "PWCSVar overload: superposition relaxation of the Chebyshev "
      "polynomial T_n(x).");
  m.def(
      "max", [](PWCSV const& x, double const& y) { return mc::max(x, y); },
      "PWCSVar overload: superposition relaxation of max(x, y) with scalar "
      "y.");
  m.def(
      "min", [](PWCSV const& x, double const& y) { return mc::min(x, y); },
      "PWCSVar overload: superposition relaxation of min(x, y) with scalar "
      "y.");
  m.def(
      "max", [](PWCSV const& x, PWCSV const& y) { return mc::max(x, y); },
      "PWCSVar overload: superposition relaxation of max(x, y); the "
      "reference operand is selected via option REF_WEIGHT.");
  m.def(
      "min", [](PWCSV const& x, PWCSV const& y) { return mc::min(x, y); },
      "PWCSVar overload: superposition relaxation of min(x, y); the "
      "reference operand is selected via option REF_WEIGHT.");

  pyPWCU
      .def(py::init<>(), R"doc(
Construct an empty estimator (no range, no segments).
)doc")
      .def(py::init<double const&, double const&, double const&,
                    size_t const&>(),
           py::arg("xL"), py::arg("xU"), py::arg("y"), py::arg("N"), R"doc(
Construct the constant function y on [xL, xU] with N grid segments.

Parameters
----------
xL : float
    Lower abscissa of the range.
xU : float
    Upper abscissa of the range.
y : float
    Constant function value.
N : int
    Number of equal-width partition segments.
)doc")
      .def(py::init<double const&, double const&, size_t const>(),
           py::arg("xL"), py::arg("xU"), py::arg("N"), R"doc(
Construct the identity function x on [xL, xU] with N grid segments.

This is the usual starting estimator for an independent variable; on
each segment, the lower level is the segment's left endpoint and the
upper level its right endpoint.

Parameters
----------
xL : float
    Lower abscissa of the range.
xU : float
    Upper abscissa of the range.
N : int
    Number of equal-width partition segments.
)doc")
      .def(py::init<PWCU const&>(), R"doc(
Construct a copy of another estimator.
)doc")
      .def_readwrite_static("options", &PWCU::options, R"doc(
Class-level option set (a ``PWCU.Options`` instance) shared by all PWCU
instances.
)doc")
      .def(
          "set", [](PWCU& self, double const& y) { return self.set(y); },
          py::arg("y"), R"doc(
Reset the estimator to the constant function y on its current range and
partition.

Returns
-------
est : PWCU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWCU& self, double const& xL, double const& xU, double const& y,
             size_t const& N) { return self.set(xL, xU, y, N); },
          py::arg("xL"), py::arg("xU"), py::arg("y"), py::arg("N"), R"doc(
Reset the estimator to the constant function y on [xL, xU] with N grid
segments.

Returns
-------
est : PWCU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "set",
          [](PWCU& self, double const& xL, double const& xU, size_t const& N)
          { return self.set(xL, xU, N); }, py::arg("xL"), py::arg("xU"),
          py::arg("N"), R"doc(
Reset the estimator to the identity on [xL, xU] with N grid segments.

Returns
-------
est : PWCU
    The updated estimator (the object itself is modified in place).
)doc")
      .def_property_readonly(
          "l", [](PWCU& self) { return self.l(); }, R"doc(
Lower bound of the estimator over its range (minimum of the per-segment
lower levels).
)doc")
      .def_property_readonly(
          "u", [](PWCU& self) { return self.u(); }, R"doc(
Upper bound of the estimator over its range (maximum of the per-segment
upper levels).
)doc")
      .def_property_readonly(
          "w", [](PWCU& self) { return self.w(); }, R"doc(
Width of the estimator range, equal to ``u - l``.
)doc")
      .def(
          "lval", [](PWCU& self, double const& x) { return self.l(x); },
          py::arg("x"), R"doc(
Evaluate the piecewise-constant underestimator at abscissa x.

Parameters
----------
x : float
    Evaluation point; must lie within the estimator range.

Returns
-------
y : float
    Lower level of the segment containing x (tightened using slope
    information when available).
)doc")
      .def(
          "uval", [](PWCU& self, double const& x) { return self.l(x); },
          py::arg("x"), R"doc(
Evaluate the piecewise-constant overestimator at abscissa x.

Parameters
----------
x : float
    Evaluation point; must lie within the estimator range.

Returns
-------
y : float
    Upper level of the segment containing x (tightened using slope
    information when available).
)doc")
      .def_property_readonly(
          "n", [](PWCU& self) { return self.size(); }, R"doc(
Number of segments in the grid partition.
)doc")
      .def_property(
          "xL", [](PWCU& self) { return self.xL(); },
          [](PWCU& self, double const& x) { self.xL() = x; }, R"doc(
Lower abscissa of the estimator range (read/write).
)doc")
      .def_property(
          "xU", [](PWCU& self) { return self.xU(); },
          [](PWCU& self, double const& x) { self.xU() = x; }, R"doc(
Upper abscissa of the estimator range (read/write).
)doc")
      .def_property(
          "yL", [](PWCU& self) { return self.yL(); },
          [](PWCU& self, std::vector<double> const& y) { self.yL() = y; },
          R"doc(
Per-segment lower levels of the estimator, as a list of float of length
``n`` (read/write).
)doc")
      .def_property(
          "yU", [](PWCU& self) { return self.yU(); },
          [](PWCU& self, std::vector<double> const& y) { self.yU() = y; },
          R"doc(
Per-segment upper levels of the estimator, as a list of float of length
``n`` (read/write).
)doc")
      .def("__str__",
           [](PWCU const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](PWCU const& self)
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
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self /= double())
      .def(py::self / double())
      .def(
          "min", [](PWCU& self, double const& cutoff)
          { return self.min(cutoff); }, py::arg("cutoff"), R"doc(
Clip the estimator from above at level `cutoff` (pointwise minimum).

Returns
-------
est : PWCU
    The updated estimator (the object itself is modified in place).
)doc")
      .def(
          "max", [](PWCU& self, double const& cutoff)
          { return self.max(cutoff); }, py::arg("cutoff"), R"doc(
Clip the estimator from below at level `cutoff` (pointwise maximum).

Returns
-------
est : PWCU
    The updated estimator (the object itself is modified in place).
)doc");

  pyPWCUOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<PWCU::Options const&>(), R"doc(
Construct a copy of another option set.
)doc")
      .def_readwrite("BKPTATOL", &PWCU::Options::BKPTATOL, R"doc(
Absolute tolerance for managing breakpoints. Default is
1e2 * DBL_EPSILON (about 2.2e-14).
)doc")
      .def_readwrite("BKPTRTOL", &PWCU::Options::BKPTRTOL, R"doc(
Relative tolerance for managing breakpoints. Default is
1e2 * DBL_EPSILON (about 2.2e-14).
)doc")
      .def_readwrite("DISPNUM", &PWCU::Options::DISPNUM, R"doc(
Number of significant digits used when printing an estimator. Default
is 5.
)doc")
      .def(
          "reset", [](PWCU::Options& self) { return self.reset(); }, R"doc(
Reset all options to their default values.
)doc");
}
