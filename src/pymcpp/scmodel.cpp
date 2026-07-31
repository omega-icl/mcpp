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

#include "scmodel.hpp"

namespace py = pybind11;

void
mc_scmodel(py::module& m)
{
  typedef mc::SCModel<I> SCM;
  typedef mc::SCVar<I> SCV;

  py::class_<SCM> pySCModel(m, "SCModel", R"doc(
Sparse Chebyshev model environment.

A q-th order Chebyshev model of a function f over an interval domain D
consists of a q-th order multivariate polynomial P in Chebyshev basis,
plus an interval remainder bound R such that f(x) belongs to P(x) + R
for every x in D. In the sparse implementation only the nonzero
monomial terms are stored, as a dictionary mapping monomials (``SMon``)
to float coefficients, so the number of variables need not be fixed
upfront: variables are registered on the fly, when an ``SCVar`` is
constructed with this environment, an index and a domain interval. The
polynomial coefficients are propagated as floating-point numbers, so
the enclosure is guaranteed up to rounding in these coefficients (the
implementation is not fully verified).

The environment fixes the maximal polynomial order and stores the
propagation options (attribute ``options``). Expressions of ``SCVar``
objects are evaluated in sparse Chebyshev model arithmetic via the
overloaded operators and the module-level functions (``exp``, ``log``,
``sin``, ...). The environment can also append auxiliary variables,
used by ``SCVar.lift`` to trade remainder width for extra dimensions.

Related environments: ``TModel``/``CModel`` (dense Taylor/Chebyshev
storage with a fixed number of variables), ``SICModel`` (sparse
Chebyshev with interval coefficients).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.SCModel(2)            # maximal order 2
>>> x = pymcpp.SCVar(env, 0, pymcpp.Interval(0, 2))
>>> f = x * pymcpp.exp(-x)             # sparse Chebyshev model
>>> f.coefmon                          # dict of SMon -> float
{1: 0.2578713111329687, [0]: 0.10094058590842268, [0]^2: -0.10597552347094877}
>>> f.B                                # enclosure of the range of f
[ -3.3121008850814523e-02 :  4.5994115272719505e-01 ]
)doc");

  py::class_<SCM::Options> pySCModelOptions(pySCModel, "Options", R"doc(
Options of a sparse Chebyshev model environment.

Modify the fields of ``SCModel.options`` in place to control range
bounding and sparse Chebyshev model propagation, e.g.
``env.options.BOUNDER_TYPE = pymcpp.SCModel.Options.BOUNDER.NAIVE``.
)doc");

  py::enum_<SCM::Options::BOUNDER>(pySCModelOptions, "BOUNDER")
      .value("NAIVE", SCM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", SCM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .export_values();

  pySCModel
      // Constructors
      .def(py::init<unsigned const>(), py::arg("maxord") = 3, R"doc(
Construct a sparse Chebyshev model environment.

Parameters
----------
maxord : int, optional
    Maximal order of the polynomial inclusion (maximal total degree
    of the polynomial part). Default is 3.
)doc")
      // Accessors
      .def_readwrite("options", &SCM::options, R"doc(
Propagation and bounding options of this environment (in-place
modifiable ``SCModel.Options`` instance).
)doc")
      .def_property_readonly(
          "maxord", [](SCM const& self) { return self.maxord(); },
          R"doc(
Maximal order of the model environment, i.e. maximal total degree of
the polynomial part (int).
)doc")
      .def_property_readonly(
          "setvar", [](SCM const& self) { return self.setvar(); },
          R"doc(
Indices of the variables registered in this environment (set of int).
)doc")
      .def_property_readonly(
          "bndvar", [](SCM const& self) { return self.bndvar(); },
          R"doc(
Domains of the registered variables, as a dict mapping variable index
(int) to domain (Interval).
)doc")
      .def_property_readonly(
          "refvar", [](SCM const& self) { return self.refvar(); },
          R"doc(
Reference points of the registered variables, as a dict mapping
variable index (int) to the domain midpoint (float) used for scaling
to [-1,1].
)doc")
      .def_property_readonly(
          "scalvar", [](SCM const& self) { return self.scalvar(); },
          R"doc(
Scaling factors of the registered variables, as a dict mapping
variable index (int) to the domain half-diameter (float) used for
scaling to [-1,1].
)doc")
      .def_property_readonly(
          "setaux", [](SCM const& self) { return self.setvar(); },
          R"doc(
Currently an alias for ``setvar``: returns the indices of the
registered model variables (set of int).
)doc")
      .def(
          "append_aux", [](SCM& self) { return self.append_aux(); },
          R"doc(
Append a new auxiliary variable to the environment, with domain
[-1,1], and return its index (int). Auxiliary variables are used by
``SCVar.lift`` to absorb remainder terms.
)doc")
      .def(
          "reset_aux", [](SCM& self) { self.reset_aux(); },
          R"doc(
Remove all auxiliary variables from the environment.
)doc");

  py::class_<SCV> pySCVar(m, "SCVar", R"doc(
Sparse Chebyshev model variable.

An ``SCVar`` holds a sparse multivariate polynomial, stored as a
dictionary mapping monomials (``SMon``) to float coefficients, together
with an interval remainder bound, valid over the variable domains
recorded in its ``SCModel`` environment. With the default option
``BASIS = CHEB``, a monomial entry (k, e) denotes the Chebyshev
polynomial T_e of variable k scaled to [-1,1]; with ``BASIS = MONOM``
it denotes the plain power.

Participating variables are created by attaching an environment, a
variable index and a domain interval; composite models are obtained via
the overloaded arithmetic operators (+, -, *, /, **, unary -, abs) and
the module-level functions (``exp``, ``log``, ``sin``, ...), which all
propagate both the polynomial part and the remainder bound. Operations
between ``SCVar`` objects linked to different environments raise an
error (``SCModel.Exceptions``).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.SCModel(4)
>>> x = pymcpp.SCVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.SCVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2
>>> f.B                      # range enclosure
[ -1.2600326699615843e+01 :  1.3998742094975464e+01 ]
>>> f.R                      # remainder bound
[ -6.8562421377154859e-01 :  6.8562421377154859e-01 ]
>>> f.P({0: 0.5, 1: 1.5})    # polynomial part at a point
-0.05811645661952752
)doc");

  pySCVar
      // Constructors
      .def(py::init<double const&, SCM*>(), py::arg("cst"),
           py::arg("env") = nullptr, R"doc(
Construct a constant sparse Chebyshev variable from a scalar.

Parameters
----------
cst : float
    Constant value.
env : SCModel, optional
    Model environment to attach to. Default is None.
)doc")
      .def(py::init<I const&, SCM*>(), py::arg("bnd"), py::arg("env") = nullptr,
           R"doc(
Construct a constant sparse Chebyshev variable from an interval.

The interval becomes the remainder bound; the polynomial part is zero.

Parameters
----------
bnd : Interval
    Constant enclosure.
env : SCModel, optional
    Model environment to attach to. Default is None.
)doc")
      .def(py::init<SCM*, unsigned const&, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"), R"doc(
Construct the sparse Chebyshev variable of index ``ndx`` in
environment ``env``, registering the variable in the environment.

Parameters
----------
env : SCModel
    Sparse Chebyshev model environment to attach to.
ndx : int
    Variable index; registered in the environment if new.
rng : Interval
    Domain (range) of the variable.
)doc")
      .def(py::init<SCV const&>(), R"doc(
Copy constructor.
)doc")
      // Modifiers
      .def(
          "set", [](SCV& self, SCM* env) -> SCV& { return self.set(env); },
          py::arg("env"), py::return_value_policy::reference_internal, R"doc(
Attach this variable to environment ``env``. Returns this object, so
that calls can be chained.
)doc")
      .def(
          "set", [](SCV& self, I const& rem) -> SCV& { return self.set(rem); },
          py::arg("rem"), py::return_value_policy::reference_internal, R"doc(
Set the remainder bound of this variable to ``rem`` (Interval),
leaving the polynomial part unchanged. Returns this object.
)doc")
      .def(
          "set", [](SCV& self, SCV::t_poly const& coefmon, I const& rem) -> SCV&
          { return self.set(coefmon); }, py::arg("coefmon") = 0.,
          py::arg("rem") = 0., py::return_value_policy::reference_internal,
          R"doc(
Set the sparse polynomial part of this variable from a monomial
coefficient dictionary. Returns this object.

Parameters
----------
coefmon : dict of SMon to float
    New monomial-coefficient map of the polynomial part.
rem : Interval, optional
    Currently ignored: the remainder bound is reset to [0,0]
    regardless of this argument.
)doc")
      .def(
          "set",
          [](SCV& self, SCM* env, unsigned const& ndx, I const& rng) -> SCV&
          { return self.set(env, ndx, rng); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as the variable of index ``ndx`` in
environment ``env``, with domain ``rng``. Returns this object.
)doc")
      // Accessors
      .def("env", &SCV::env, R"doc(
Return the model environment this variable is attached to
(``SCModel``, or None for a constant).
)doc")
      .def_property_readonly(
          "nord", [](SCV const& self) { return self.nord(); }, R"doc(
Maximal total order of the monomial terms present in this variable
(int).
)doc")
      .def_property_readonly(
          "nvar", [](SCV const& self) { return self.nvar(); },
          R"doc(
Number of variables participating in this variable (int).
)doc")
      .def_property_readonly(
          "nmon", [](SCV const& self) { return self.nmon(); },
          R"doc(
Number of nonzero monomial terms in the polynomial part (int).
)doc")
      .def_property_readonly(
          "coefmon", py::overload_cast<>(&SCV::coefmon, py::const_),
          R"doc(
Monomial-coefficient map of the polynomial part: dict mapping ``SMon``
monomials to float coefficients. With the default option
``BASIS = CHEB``, an entry (k, e) in a monomial denotes the Chebyshev
polynomial T_e of variable k scaled to [-1,1].
)doc")
      .def_property_readonly("ndxvar",
                             py::overload_cast<>(&SCV::ndxvar, py::const_),
                             //[]( SCV const& self )
                             //  { return self.ndxvar(); },
                             R"doc(
Indices of the variables participating in this variable (set of int).
)doc")
      .def_property_readonly(
          "bound", [](SCV const& self) { return self.bound(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus
the remainder bound. Same as ``B``.
)doc")
      .def_property_readonly(
          "B", [](SCV const& self) { return self.B(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus
the remainder bound. Same as ``bound``.
)doc")
      .def(
          "bndpol", [](SCV const& self) { return self.bndpol(); },
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using the bounder selected in option ``BOUNDER_TYPE``. The remainder
bound is not added.
)doc")
      .def(
          "bndpol", [](SCV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using bounder ``type`` (a ``SCModel.Options.BOUNDER`` value).
)doc")
      .def("bndord", &SCV::bndord, py::arg("minord"), R"doc(
Return a bound on the sum of all monomial terms whose total order is
``minord`` or higher (Interval).
)doc")
      .def(
          "P", [](SCV const& self, std::map<unsigned, double> const& x)
          { return self.P(x); }, py::arg("x"),
          R"doc(
Evaluate the polynomial part at point ``x``.

Parameters
----------
x : dict of int to float
    Point in the original (unscaled) variable domain, mapping the
    index of each participating variable to its value.

Returns
-------
val : float
    Value of the polynomial part at ``x`` (remainder not included).
)doc")
      .def(
          "P", [](SCV const& self) { return self.P(); },
          R"doc(
Return a new sparse Chebyshev variable with the same polynomial part
and a zero remainder bound. This variable is not modified.
)doc")
      .def_property_readonly("R", &SCV::R, R"doc(
Remainder bound of the sparse Chebyshev model (Interval).
)doc")
      .def("center", &SCV::center, R"doc(
Center the remainder term in place: the midpoint of the remainder
bound is added to the constant polynomial coefficient and subtracted
from the remainder. Returns this object.
)doc")
      .def("C", &SCV::C, R"doc(
Center the remainder term in place; shortcut for ``center``.
)doc")
      .def("constant", &SCV::constant, py::arg("reset") = false, R"doc(
Return the coefficient of the constant term (float).

Parameters
----------
reset : bool, optional
    If True, also remove the constant term from this variable.
    Default is False.
)doc")
      .def("linear", &SCV::linear, py::arg("id"), py::arg("reset") = false,
           R"doc(
Return the coefficient of the linear term in variable ``id`` (float),
expressed in the scaled variable. Returns 0.0 if no such term is
present.

Parameters
----------
id : int
    Variable index.
reset : bool, optional
    If True, also remove this linear term from the variable. Default
    is False.
)doc")
      .def(
          "lift",
          [](SCV& self, SCM& scm, double const& atol, double const& rtol)
          { return self.lift(&scm, atol, rtol); }, py::arg("scm"),
          py::arg("atol"), py::arg("rtol"),
          R"doc(
Absorb the remainder term into a new auxiliary variable.

If the remainder radius exceeds ``atol`` plus ``rtol`` times half the
polynomial range diameter, the model is centered, a new auxiliary
variable is appended to environment ``scm``, the remainder radius
becomes the coefficient of that auxiliary variable, and the remainder
bound is reset to zero. Otherwise the variable is left unchanged.
This variable is modified in place and also returned.

Parameters
----------
scm : SCModel
    Environment receiving the new auxiliary variable.
atol : float
    Absolute threshold on the remainder radius.
rtol : float
    Relative threshold on the remainder radius.

Returns
-------
var : SCVar
    The lifted variable.
)doc")
      .def(
          "project", [](SCV& self, bool const reset)
          { return self.project(reset); }, py::arg("reset") = false,
          R"doc(
Remove all monomials in which auxiliary variables participate.

The removed contributions are enclosed back into the remainder bound
(and partly reallocated to reduced monomials). This variable is
modified in place and also returned.

Parameters
----------
reset : bool, optional
    If True, also remove all auxiliary variables from the
    environment afterwards. Default is False.
)doc")
      .def(
          "project", [](SCV& self, unsigned const& id)
          { return self.project(id); }, py::arg("id"),
          R"doc(
Remove all monomials in which variable ``id`` participates, enclosing
the removed contributions into the remainder bound. This variable is
modified in place and also returned.
)doc")
      .def(
          "scale", [](SCV& self, unsigned const& id, I const& dom)
          { return self.scale(id, dom); }, py::arg("id"), py::arg("dom"),
          R"doc(
Rescale the polynomial coefficients so that the model is expressed
over the modified domain ``dom`` of variable ``id``. This variable is
modified in place and also returned.
)doc")
      .def(
          "rescale", [](SCV& self, std::map<unsigned, I> const& dom)
          { return self.scale(dom); }, py::arg("dom"),
          R"doc(
Rescale the polynomial coefficients so that the model is expressed
over the modified variable domains ``dom`` (dict mapping variable
index to new domain interval). This variable is modified in place and
also returned.
)doc")
      .def("unscale", &SCV::unscale, R"doc(
Return the monomial-coefficient map of the polynomial part expressed
in the original (unscaled) variables, as a dict mapping ``SMon``
monomials to float coefficients. This variable is not modified.
)doc")
      .def("simplify", &SCV::simplify, py::arg("atol") = 0e0,
           py::arg("rtol") = 0e0, py::arg("tord") = -1, R"doc(
Simplify the polynomial part in place.

Monomial terms whose coefficient magnitude is below ``atol`` plus
``rtol`` times half the polynomial range diameter, or whose total
order exceeds ``tord`` (if ``tord`` is nonnegative), are removed and
their range contribution is added to the remainder bound. Returns
this variable.

Parameters
----------
atol : float, optional
    Absolute coefficient threshold. Default is 0.
rtol : float, optional
    Relative coefficient threshold. Default is 0.
tord : int, optional
    Maximal total order kept; -1 keeps all orders. Default is -1.
)doc")
      .def(
          "to_monomial", [](SCV const& self, bool const scaled)
          { return self.to_monomial(scaled); }, py::arg("scaled") = false,
          R"doc(
Convert the polynomial part to monomial (power) basis.

Parameters
----------
scaled : bool, optional
    If True, express the coefficients in the variables scaled to
    [-1,1]; if False, in the original variables. Default is False.

Returns
-------
coefmon : dict of SMon to float
    New monomial-coefficient map in monomial basis. This variable is
    not modified.
)doc")
      .def(
          "to_monomial",
          [](SCV const& self, bool const scaled, double const& atol,
             double const& rtol, int const tord)
          { return self.to_monomial(scaled, atol, rtol, tord); },
          py::arg("scaled"), py::arg("atol"), py::arg("rtol"),
          py::arg("tord") = -1,
          R"doc(
Convert the polynomial part to monomial (power) basis, discarding
small or high-order terms.

Parameters
----------
scaled : bool
    If True, express the coefficients in the variables scaled to
    [-1,1]; if False, in the original variables.
atol : float
    Absolute coefficient threshold below which terms are removed.
rtol : float
    Relative coefficient threshold below which terms are removed.
tord : int, optional
    Terms of total order greater than ``tord`` are removed if
    ``tord`` is nonnegative; -1 keeps all orders. Default is -1.

Returns
-------
coefmon : dict of SMon to float
    New monomial-coefficient map in monomial basis, without the
    removed terms.
rem : Interval
    Bound on the sum of the removed terms.
)doc")
      .def("__str__",
           [](SCV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](SCV const& self)
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
      .def("__abs__", [](SCV const& m) { return mc::Op<SCV>::abs(m); })
      .def("__pow__", [](SCV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](SCV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](SCV const& m, SCV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, SCV const& m) { return mc::pow(r, m); });

  m.def(
      "inv", [](SCV const& x) { return mc::inv(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of 1/x.");
  m.def(
      "sqr", [](SCV const& x) { return mc::sqr(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of x**2.");
  m.def(
      "sqrt", [](SCV const& x) { return mc::sqrt(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of sqrt(x).");
  m.def(
      "exp", [](SCV const& x) { return mc::exp(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of exp(x).");
  m.def(
      "log", [](SCV const& x) { return mc::log(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of log(x).");
  m.def(
      "xlog", [](SCV const& x) { return mc::xlog(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of x*log(x).");
  m.def(
      "cos", [](SCV const& x) { return mc::cos(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of cos(x).");
  m.def(
      "sin", [](SCV const& x) { return mc::sin(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of sin(x).");
  m.def(
      "tan", [](SCV const& x) { return mc::tan(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of tan(x).");
  m.def(
      "acos", [](SCV const& x) { return mc::acos(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of acos(x).");
  m.def(
      "asin", [](SCV const& x) { return mc::asin(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of asin(x).");
  m.def(
      "atan", [](SCV const& x) { return mc::atan(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of atan(x).");
  m.def(
      "cosh", [](SCV const& x) { return mc::cosh(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of cosh(x).");
  m.def(
      "sinh", [](SCV const& x) { return mc::sinh(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of sinh(x).");
  m.def(
      "tanh", [](SCV const& x) { return mc::tanh(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of tanh(x).");
  m.def(
      "erf", [](SCV const& x) { return mc::erf(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of erf(x).");
  m.def(
      "erfc", [](SCV const& x) { return mc::erfc(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of erfc(x).");
  m.def(
      "pow", [](SCV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"),
      "Sparse Chebyshev model overload: Chebyshev model of x**n for integer "
      "n.");
  m.def(
      "pow", [](SCV const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"),
      "Sparse Chebyshev model overload: Chebyshev model of x**r for real r.");
  m.def(
      "pow", [](SCV const& x, SCV const& y) { return mc::pow(x, y); },
      py::arg("x"), py::arg("y"),
      "Sparse Chebyshev model overload: Chebyshev model of x**y.");
  m.def(
      "pow", [](double const& r, SCV const& y) { return mc::pow(r, y); },
      py::arg("r"), py::arg("y"),
      "Sparse Chebyshev model overload: Chebyshev model of r**y for real r.");
  m.def(
      "cheb", [](SCV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"),
      "Sparse Chebyshev model overload: Chebyshev model of the Chebyshev "
      "polynomial T_n(x).");
  m.def(
      "fabs", [](SCV const& x) { return mc::fabs(x); }, py::arg("x"),
      "Sparse Chebyshev model overload: Chebyshev model of abs(x).");
  m.def(
      "hull", [](SCV const& x, SCV const& y) { return mc::hull(x, y); },
      py::arg("x"), py::arg("y"),
      "Sparse Chebyshev model overload: Chebyshev model enclosing the union "
      "of the enclosures x and y (see option REF_POLY for the polynomial "
      "part).");
  m.def(
      "inter",
      [](SCV& xy, SCV const& x, SCV const& y) { return mc::inter(xy, x, y); },
      py::arg("xy"), py::arg("x"), py::arg("y"),
      "Sparse Chebyshev model overload: overwrite xy with a Chebyshev model "
      "of the intersection of the enclosures x and y; returns False if the "
      "intersection is empty.");

  // Nested class Options

  pySCModelOptions.def(py::init<>(), "Construct options with default values.")
      .def(py::init<SCM::Options const&>(), "Copy constructor.")
      .def(
          "reset", [](SCM::Options& self) { self = SCM::Options(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("BASIS", &SCM::Options::BASIS, R"doc(
Basis representation of the monomials, as a
``SCModel.Options.MONBASIS`` value: 0 (MONOM) for monomial basis, 1
(CHEB) for Chebyshev basis. Default is 1 (CHEB).
)doc")
      .def_readwrite("LIFT_USE", &SCM::Options::LIFT_USE, R"doc(
Whether to lift the remainder term in nonlinear operations by
introducing auxiliary variables in the model (bool). Default is
False.
)doc")
      .def_readwrite("LIFT_ATOL", &SCM::Options::LIFT_ATOL, R"doc(
Absolute tolerance for lifting the remainder term; only used if
``LIFT_USE`` is True. Default is 1e-10.
)doc")
      .def_readwrite("LIFT_RTOL", &SCM::Options::LIFT_RTOL, R"doc(
Relative tolerance for lifting the remainder term; only used if
``LIFT_USE`` is True. Default is 1e-3.
)doc")
      .def_readwrite("REMEZ_USE", &SCM::Options::REMEZ_USE, R"doc(
Whether to use the Remez algorithm for computing a minimax
approximation of univariate terms (bool). Default is True.
)doc")
      .def_readwrite("REMEZ_MAXIT", &SCM::Options::REMEZ_MAXIT, R"doc(
Maximal number of iterations in the Remez algorithm for computing a
minimax approximation of univariate terms (int). Default is 10.
)doc")
      .def_readwrite("REMEZ_TOL", &SCM::Options::REMEZ_TOL, R"doc(
Stopping tolerance in the Remez algorithm for computing a minimax
approximation of univariate terms (float). Default is 1e-5.
)doc")
      .def_readwrite("REMEZ_MIG", &SCM::Options::REMEZ_MIG, R"doc(
Threshold on the interval width below which the Remez algorithm is
not used (float). Default is 1e-10.
)doc")
      .def_readwrite("INTERP_EXTRA", &SCM::Options::INTERP_EXTRA, R"doc(
Number of extra terms in the Chebyshev interpolation of univariate
functions (int): 0 uses a Chebyshev interpolation of the model order;
extra terms allow approximating the truncated Chebyshev series.
Default is 0.
)doc")
      .def_readwrite("INTERP_THRES", &SCM::Options::INTERP_THRES, R"doc(
Threshold below which coefficients in the Chebyshev expansion of
transcendental univariates are ignored during bounding (float).
Default is 1e2 times the machine epsilon.
)doc")
      .def_readwrite("BOUNDER_TYPE", &SCM::Options::BOUNDER_TYPE, R"doc(
Chebyshev model range bounder, as a ``SCModel.Options.BOUNDER``
value: NAIVE or LSB. Default is LSB.
)doc")
      .def_readwrite("MIG_USE", &SCM::Options::MIG_USE, R"doc(
Whether to simplify monomial terms with small magnitude during
propagation (bool). Default is False.
)doc")
      .def_readwrite("MIG_ATOL", &SCM::Options::MIG_ATOL, R"doc(
Absolute tolerance for simplifying monomial terms; only used if
``MIG_USE`` is True. Default is 0.
)doc")
      .def_readwrite("MIG_RTOL", &SCM::Options::MIG_RTOL, R"doc(
Relative tolerance for simplifying monomial terms; only used if
``MIG_USE`` is True. Default is the machine epsilon.
)doc")
      .def_readwrite("MIXED_IA", &SCM::Options::MIXED_IA, R"doc(
Whether to intersect the internal model bounds with bounds propagated
in the underlying (interval) arithmetic (bool). Default is False.
)doc")
      .def_readwrite("REF_POLY", &SCM::Options::REF_POLY, R"doc(
Scalar in [0,1] selecting the polynomial part in the functions
``inter`` and ``hull``: 0 selects the polynomial part of the left
operand, 1 that of the right operand, intermediate values
interpolate. Default is 0.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &SCM::Options::DISPLAY_DIGITS, R"doc(
Number of digits used when printing Chebyshev model coefficients
(int). Default is 7.
)doc");


  py::enum_<SCM::Options::MONBASIS>(pySCModelOptions, "MONBASIS")
      .value("MONOM", SCM::Options::MONBASIS::MONOM, "Monomial basis")
      .value("CHEB", SCM::Options::MONBASIS::CHEB, "Chebyshev basis")
      .export_values();

  // Nested class Exceptions
  py::class_<SCM::Exceptions> pySCModelExceptions(pySCModel, "Exceptions",
                                                  R"doc(
Exception data thrown by sparse Chebyshev model arithmetic.

Carries an error code (``ierr``) and a description (``what``); see the
``SCModel.Exceptions.TYPE`` enumeration for the possible errors.
)doc");

  py::enum_<SCM::Exceptions::TYPE>(pySCModelExceptions, "TYPE")
      .value("DIV", SCM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", SCM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", SCM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", SCM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("DPOW", SCM::Exceptions::TYPE::DPOW,
             "Real power operation with negative numbers in range")
      .value("TAN", SCM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", SCM::Exceptions::TYPE::ACOS,
             "Cosine inverse operation with range outside [-1,1]")
      .value("ASIN", SCM::Exceptions::TYPE::ASIN,
             "Sine inverse operation with range outside [-1,1]")
      .value("COMPOSE", SCM::Exceptions::TYPE::COMPOSE,
             "Failed to compose Chebyshev variable")
      .value("INIT", SCM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", SCM::Exceptions::TYPE::INCON,
             "Inconsistent bounds with template parameter arithmetic")
      .value("MODEL", SCM::Exceptions::TYPE::MODEL,
             "Operation between variables linked to different models")
      .value("INTERNAL", SCM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", SCM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented")
      .export_values();

  pySCModelExceptions.def(py::init<SCM::Exceptions::TYPE>())
      .def("ierr", &SCM::Exceptions::ierr, "Return the error code (int).")
      .def("what", &SCM::Exceptions::what, "Return the error description (str).");
}
