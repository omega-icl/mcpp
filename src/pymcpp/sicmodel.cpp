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

#include "sicmodel.hpp"

namespace py = pybind11;

void
mc_sicmodel(py::module& m)
{
  typedef mc::SICModel<I> SICM;
  typedef mc::SICVar<I> SICV;

  py::class_<SICM> pySICModel(m, "SICModel", R"doc(
Sparse interval Chebyshev model environment.

A q-th order sparse interval Chebyshev model of a function f over an
interval domain D is a q-th order multivariate polynomial P in
Chebyshev basis with INTERVAL coefficients, such that f(x) belongs to
P(x) for every x in D. Unlike ``SCModel``, there is no separate
remainder term: the distinction between polynomial part and remainder
is blurred by distributing the approximation error over the interval
coefficients of the participating monomials, which in general enables
tighter enclosures (Zha & Chachuat, 2021). Only the nonzero monomial
terms are stored, as a dictionary mapping monomials (``SMon``) to
Interval coefficients, and the number of variables need not be fixed
upfront: variables are registered on the fly, when an ``SICVar`` is
constructed with this environment, an index and a domain interval.

The environment fixes the maximal polynomial order and stores the
propagation options (attribute ``options``). Expressions of ``SICVar``
objects are evaluated in sparse interval Chebyshev model arithmetic via
the overloaded operators and the module-level functions (``exp``,
``log``, ``sin``, ...). The environment can also append auxiliary
variables, used by ``SICVar.lift`` to trade coefficient uncertainty for
extra dimensions.

Related environments: ``TModel``/``CModel`` (dense Taylor/Chebyshev
storage), ``SCModel`` (sparse Chebyshev with scalar coefficients and a
separate remainder bound).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.SICModel(2)           # maximal order 2
>>> x = pymcpp.SICVar(env, 0, pymcpp.Interval(0, 2))
>>> f = x * pymcpp.exp(-x)             # sparse interval Chebyshev model
>>> f.coefmon                          # dict of SMon -> Interval
{1: [  1.9035560495396409e-01 :  3.2538701731197323e-01 ],
 [0]: [  8.4380081483015332e-02 :  1.1750109033382998e-01 ],
 [0]^2: [ -1.0597552347094877e-01 : -1.0597552347094877e-01 ]}
>>> f.B                                # enclosure of the range of f
[ -3.3121008850814815e-02 :  4.5994115272719521e-01 ]
)doc");

  py::class_<SICM::Options> pySICModelOptions(pySICModel, "Options", R"doc(
Options of a sparse interval Chebyshev model environment.

Modify the fields of ``SICModel.options`` in place to control range
bounding and model propagation, e.g.
``env.options.BOUNDER_TYPE = pymcpp.SICModel.Options.BOUNDER.NAIVE``.
)doc");

  py::enum_<SICM::Options::BOUNDER>(pySICModelOptions, "BOUNDER")
      .value("NAIVE", SICM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", SICM::Options::BOUNDER::LSB,
             "Lin & Stadtherr range bounder")
      .value("BERNSTEIN", SICM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .export_values();

  pySICModel
      // Constructors
      .def(py::init<unsigned const>(), py::arg("maxord") = 3, R"doc(
Construct a sparse interval Chebyshev model environment.

Parameters
----------
maxord : int, optional
    Maximal order of the polynomial inclusion (maximal total degree
    of the polynomial part). Default is 3.
)doc")
      // Accessors
      .def_readwrite("options", &SICM::options, R"doc(
Propagation and bounding options of this environment (in-place
modifiable ``SICModel.Options`` instance).
)doc")
      .def_property_readonly(
          "maxord", [](SICM const& self) { return self.maxord(); },
          R"doc(
Maximal order of the model environment, i.e. maximal total degree of
the polynomial part (int).
)doc")
      .def_property_readonly(
          "setvar", [](SICM const& self) { return self.setvar(); },
          R"doc(
Indices of the variables registered in this environment (set of int).
)doc")
      .def_property_readonly(
          "bndvar", [](SICM const& self) { return self.bndvar(); },
          R"doc(
Domains of the registered variables, as a dict mapping variable index
(int) to domain (Interval).
)doc")
      .def_property_readonly(
          "refvar", [](SICM const& self) { return self.refvar(); },
          R"doc(
Reference points of the registered variables, as a dict mapping
variable index (int) to the domain midpoint (float) used for scaling
to [-1,1].
)doc")
      .def_property_readonly(
          "scalvar", [](SICM const& self) { return self.scalvar(); },
          R"doc(
Scaling factors of the registered variables, as a dict mapping
variable index (int) to the domain half-diameter (float) used for
scaling to [-1,1].
)doc")
      .def(
          "append_aux", [](SICM& self) { return self.append_aux(); },
          R"doc(
Append a new auxiliary variable to the environment, with domain
[-1,1], and return its index (int). Auxiliary variables are used by
``SICVar.lift`` to absorb coefficient uncertainty.
)doc")
      .def(
          "reset_aux", [](SICM& self) { self.reset_aux(); },
          R"doc(
Remove all auxiliary variables from the environment.
)doc");

  py::class_<SICV> pySICVar(m, "SICVar", R"doc(
Sparse interval Chebyshev model variable.

An ``SICVar`` holds a sparse multivariate polynomial with INTERVAL
coefficients, stored as a dictionary mapping monomials (``SMon``) to
Interval coefficients, valid over the variable domains recorded in its
``SICModel`` environment. There is no separate remainder bound: the
approximation error is carried by the widths of the interval
coefficients. With the default option ``BASIS = CHEB``, a monomial
entry (k, e) denotes the Chebyshev polynomial T_e of variable k scaled
to [-1,1]; with ``BASIS = MONOM`` it denotes the plain power.

Participating variables are created by attaching an environment, a
variable index and a domain interval; composite models are obtained via
the overloaded arithmetic operators (+, -, *, /, **, unary -, abs) and
the module-level functions (``exp``, ``log``, ``sin``, ...). Operations
between ``SICVar`` objects linked to different environments raise an
error (``SICModel.Exceptions``).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.SICModel(4)
>>> x = pymcpp.SICVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.SICVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2
>>> f.B                       # range enclosure
[ -1.2600326699615895e+01 :  1.3998742094975512e+01 ]
>>> f.IP({0: 0.5, 1: 1.5})    # interval polynomial value at a point
[ -6.2901096554861147e-01 :  5.1277805230955553e-01 ]
)doc");

  pySICVar
      // Constructors
      .def(py::init<double const&, SICM*>(), py::arg("cst"),
           py::arg("env") = nullptr, R"doc(
Construct a constant sparse interval Chebyshev variable from a scalar.

Parameters
----------
cst : float
    Constant value.
env : SICModel, optional
    Model environment to attach to. Default is None.
)doc")
      .def(py::init<I const&, SICM*>(), py::arg("bnd"),
           py::arg("env") = nullptr, R"doc(
Construct a constant sparse interval Chebyshev variable from an
interval, stored as the constant coefficient.

Parameters
----------
bnd : Interval
    Constant enclosure.
env : SICModel, optional
    Model environment to attach to. Default is None.
)doc")
      .def(py::init<SICM*, unsigned const&, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"), R"doc(
Construct the sparse interval Chebyshev variable of index ``ndx`` in
environment ``env``, registering the variable in the environment.

Parameters
----------
env : SICModel
    Sparse interval Chebyshev model environment to attach to.
ndx : int
    Variable index; registered in the environment if new.
rng : Interval
    Domain (range) of the variable.
)doc")
      .def(py::init<SICV const&>(), R"doc(
Copy constructor.
)doc")
      // Modifiers
      .def(
          "set", [](SICV& self, SICM* env) -> SICV& { return self.set(env); },
          py::arg("env"), py::return_value_policy::reference_internal, R"doc(
Attach this variable to environment ``env``. Returns this object, so
that calls can be chained.
)doc")
      .def(
          "set",
          [](SICV& self, I const& bnd) -> SICV& { return self.set(bnd); },
          py::arg("bnd"), py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as an interval-valued constant: the
polynomial part reduces to the constant coefficient ``bnd``. Returns
this object.
)doc")
      .def(
          "set",
          [](SICV& self, SICM* env, unsigned const& ndx, I const& rng) -> SICV&
          { return self.set(env, ndx, rng); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as the variable of index ``ndx`` in
environment ``env``, with domain ``rng``. Returns this object.
)doc")
      // Accessors
      .def("env", &SICV::env, R"doc(
Return the model environment this variable is attached to
(``SICModel``, or None for a constant).
)doc")
      .def_property_readonly(
          "nord", [](SICV const& self) { return self.nord(); },
          R"doc(
Maximal total order of the monomial terms present in this variable
(int).
)doc")
      .def_property_readonly(
          "nvar", [](SICV const& self) { return self.nvar(); },
          R"doc(
Number of variables participating in this variable (int).
)doc")
      .def_property_readonly(
          "nmon", [](SICV const& self) { return self.nmon(); },
          R"doc(
Number of nonzero monomial terms in the polynomial part (int).
)doc")
      .def_property_readonly(
          "coefmon", py::overload_cast<>(&SICV::coefmon, py::const_),
          R"doc(
Monomial-coefficient map of the interval polynomial: dict mapping
``SMon`` monomials to Interval coefficients. With the default option
``BASIS = CHEB``, an entry (k, e) in a monomial denotes the Chebyshev
polynomial T_e of variable k scaled to [-1,1].
)doc")
      .def_property_readonly("ndxvar",
                             py::overload_cast<>(&SICV::ndxvar, py::const_),
                             R"doc(
Indices of the variables participating in this variable (set of int).
)doc")
      .def_property_readonly(
          "bound", [](SICV const& self) { return self.bound(); },
          R"doc(
Enclosure of the variable range (Interval), computed from the interval
polynomial with the bounder selected in option ``BOUNDER_TYPE``. Same
as ``B``.
)doc")
      .def_property_readonly(
          "B", [](SICV const& self) { return self.B(); },
          R"doc(
Enclosure of the variable range (Interval), computed from the interval
polynomial with the bounder selected in option ``BOUNDER_TYPE``. Same
as ``bound``.
)doc")
      .def(
          "bndpol", [](SICV const& self) { return self.bndpol(); },
          R"doc(
Return a bound on the interval-valued polynomial (Interval), using the
bounder selected in option ``BOUNDER_TYPE``.
)doc")
      .def("bndord", &SICV::bndord, py::arg("minord"), R"doc(
Return a bound on the sum of all monomial terms whose total order is
``minord`` or higher (Interval).
)doc")
      .def(
          "P", [](SICV const& self, std::map<unsigned, double> const& x)
          { return self.P(x); }, py::arg("x"),
          R"doc(
Evaluate the mid-point polynomial at point ``x``: each interval
coefficient is replaced by its midpoint.

Parameters
----------
x : dict of int to float
    Point in the original (unscaled) variable domain, mapping the
    index of each participating variable to its value.

Returns
-------
val : float
    Value of the mid-point polynomial at ``x``.
)doc")
      .def(
          "IP", [](SICV const& self, std::map<unsigned, double> const& x)
          { return self.IP(x); }, py::arg("x"),
          R"doc(
Evaluate the interval-valued polynomial at point ``x``.

Parameters
----------
x : dict of int to float
    Point in the original (unscaled) variable domain, mapping the
    index of each participating variable to its value.

Returns
-------
val : Interval
    Enclosure of f(x): the interval polynomial evaluated at ``x``.
)doc")
      .def(
          "P", [](SICV const& self) { return self.P(); },
          R"doc(
Return a new variable in which every non-constant coefficient is
replaced by its midpoint and the associated uncertainty is aggregated
into the constant coefficient. This variable is not modified.
)doc")
      .def("constant", &SICV::constant, py::arg("reset") = false, R"doc(
Return the interval coefficient of the constant term (Interval).

Parameters
----------
reset : bool, optional
    If True, also remove the constant term from this variable.
    Default is False.
)doc")
      .def("linear", &SICV::linear, py::arg("id"), py::arg("reset") = false,
           R"doc(
Return the interval coefficient of the linear term in variable ``id``
(Interval), expressed in the scaled variable. Returns [0,0] if no such
term is present.

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
          [](SICV& self, SICM& sicm, double const& atol, double const& rtol)
          { return self.lift(&sicm, atol, rtol); }, py::arg("sicm"),
          py::arg("atol"), py::arg("rtol"),
          R"doc(
Absorb the coefficient uncertainty into a new auxiliary variable.

If the aggregated uncertainty of the interval coefficients exceeds the
``atol``/``rtol`` thresholds, a new auxiliary variable is appended to
environment ``sicm`` and receives that uncertainty as its coefficient,
tightening the remaining coefficients to scalars. Otherwise the
variable is left unchanged. This variable is modified in place and
also returned.

Parameters
----------
sicm : SICModel
    Environment receiving the new auxiliary variable.
atol : float
    Absolute threshold on the uncertainty magnitude.
rtol : float
    Relative threshold on the uncertainty magnitude.

Returns
-------
var : SICVar
    The lifted variable.
)doc")
      .def(
          "project", [](SICV& self, bool const reset)
          { return self.project(reset); }, py::arg("reset") = false,
          R"doc(
Remove all monomials in which auxiliary variables participate.

The removed contributions are enclosed back into the interval
coefficients of the remaining monomials. This variable is modified in
place and also returned.

Parameters
----------
reset : bool, optional
    If True, also remove all auxiliary variables from the
    environment afterwards. Default is False.
)doc")
      .def(
          "project", [](SICV& self, unsigned const& id)
          { return self.project(id); }, py::arg("id"),
          R"doc(
Remove all monomials in which variable ``id`` participates, enclosing
the removed contributions into the remaining coefficients. This
variable is modified in place and also returned.
)doc")
      .def(
          "scale", [](SICV& self, unsigned const& id, I const& dom)
          { return self.scale(id, dom); }, py::arg("id"), py::arg("dom"),
          R"doc(
Rescale the polynomial coefficients so that the model is expressed
over the modified domain ``dom`` of variable ``id``. This variable is
modified in place and also returned.
)doc")
      .def(
          "rescale", [](SICV& self, std::map<unsigned, I> const& dom)
          { return self.scale(dom); }, py::arg("dom"),
          R"doc(
Rescale the polynomial coefficients so that the model is expressed
over the modified variable domains ``dom`` (dict mapping variable
index to new domain interval). This variable is modified in place and
also returned.
)doc")
      .def("unscale", &SICV::unscale, R"doc(
Return the monomial-coefficient map of the interval polynomial
expressed in the original (unscaled) variables, as a dict mapping
``SMon`` monomials to Interval coefficients. This variable is not
modified.
)doc")
      .def("simplify", &SICV::simplify, py::arg("atol") = 0e0,
           py::arg("rtol") = 0e0, py::arg("tord") = -1, R"doc(
Simplify the interval polynomial in place.

Monomial terms whose coefficient magnitude is below ``atol`` plus
``rtol`` times half the polynomial range diameter, or whose total
order exceeds ``tord`` (if ``tord`` is nonnegative), are removed and
their range contribution is enclosed into the constant coefficient.
Returns this variable.

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
          "to_monomial", [](SICV const& self, bool const scaled)
          { return self.to_monomial(scaled); }, py::arg("scaled") = false,
          R"doc(
Convert the interval polynomial to monomial (power) basis.

Parameters
----------
scaled : bool, optional
    If True, express the coefficients in the variables scaled to
    [-1,1]; if False, in the original variables. Default is False.

Returns
-------
coefmon : dict of SMon to Interval
    New monomial-coefficient map in monomial basis. This variable is
    not modified.
)doc")
      .def(
          "to_monomial",
          [](SICV const& self, bool const scaled, double const& atol,
             double const& rtol, int const tord)
          { return self.to_monomial(scaled, atol, rtol, tord); },
          py::arg("scaled"), py::arg("atol"), py::arg("rtol"),
          py::arg("tord") = -1,
          R"doc(
Convert the interval polynomial to monomial (power) basis, discarding
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
coefmon : dict of SMon to Interval
    New monomial-coefficient map in monomial basis, without the
    removed terms.
rem : Interval
    Bound on the sum of the removed terms.
)doc")
      .def("__str__",
           [](SICV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](SICV const& self)
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
      .def("__abs__", [](SICV const& m) { return mc::Op<SICV>::abs(m); })
      .def("__pow__", [](SICV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](SICV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](SICV const& m, SICV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, SICV const& m) { return mc::pow(r, m); });

  m.def(
      "inv", [](SICV const& x) { return mc::inv(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of 1/x.");
  m.def(
      "sqr", [](SICV const& x) { return mc::sqr(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of x**2.");
  m.def(
      "sqrt", [](SICV const& x) { return mc::sqrt(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of sqrt(x).");
  m.def(
      "exp", [](SICV const& x) { return mc::exp(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of exp(x).");
  m.def(
      "log", [](SICV const& x) { return mc::log(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of log(x).");
  m.def(
      "xlog", [](SICV const& x) { return mc::xlog(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of x*log(x).");
  m.def(
      "cos", [](SICV const& x) { return mc::cos(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of cos(x).");
  m.def(
      "sin", [](SICV const& x) { return mc::sin(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of sin(x).");
  m.def(
      "tan", [](SICV const& x) { return mc::tan(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of tan(x).");
  m.def(
      "acos", [](SICV const& x) { return mc::acos(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of acos(x).");
  m.def(
      "asin", [](SICV const& x) { return mc::asin(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of asin(x).");
  m.def(
      "atan", [](SICV const& x) { return mc::atan(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of atan(x).");
  m.def(
      "cosh", [](SICV const& x) { return mc::cosh(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of cosh(x).");
  m.def(
      "sinh", [](SICV const& x) { return mc::sinh(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of sinh(x).");
  m.def(
      "tanh", [](SICV const& x) { return mc::tanh(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of tanh(x).");
  m.def(
      "erf", [](SICV const& x) { return mc::erf(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of erf(x).");
  m.def(
      "erfc", [](SICV const& x) { return mc::erfc(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of erfc(x).");
  m.def(
      "pow", [](SICV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"),
      "Sparse interval Chebyshev model overload: model of x**n for integer "
      "n.");
  m.def(
      "pow", [](SICV const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"),
      "Sparse interval Chebyshev model overload: model of x**r for real r.");
  m.def(
      "pow", [](SICV const& x, SICV const& y) { return mc::pow(x, y); },
      py::arg("x"), py::arg("y"),
      "Sparse interval Chebyshev model overload: model of x**y.");
  m.def(
      "pow", [](double const& r, SICV const& y) { return mc::pow(r, y); },
      py::arg("r"), py::arg("y"),
      "Sparse interval Chebyshev model overload: model of r**y for real r.");
  m.def(
      "cheb", [](SICV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"),
      "Sparse interval Chebyshev model overload: model of the Chebyshev "
      "polynomial T_n(x).");
  m.def(
      "fabs", [](SICV const& x) { return mc::fabs(x); }, py::arg("x"),
      "Sparse interval Chebyshev model overload: model of abs(x).");
  m.def(
      "hull", [](SICV const& x, SICV const& y) { return mc::hull(x, y); },
      py::arg("x"), py::arg("y"),
      "Sparse interval Chebyshev model overload: model enclosing the union "
      "of the enclosures x and y (see option REF_POLY for the polynomial "
      "part).");
  m.def(
      "inter",
      [](SICV& xy, SICV const& x, SICV const& y)
      { return mc::inter(xy, x, y); },
      py::arg("xy"), py::arg("x"), py::arg("y"),
      "Sparse interval Chebyshev model overload: overwrite xy with a model "
      "of the intersection of the enclosures x and y; returns False if the "
      "intersection is empty.");

  // Nested class Options

  pySICModelOptions.def(py::init<>(), "Construct options with default values.")
      .def(py::init<SICM::Options const&>(), "Copy constructor.")
      .def(
          "reset", [](SICM::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("BASIS", &SICM::Options::BASIS, R"doc(
Basis representation of the monomials, as a
``SICModel.Options.MONBASIS`` value: 0 (MONOM) for monomial basis, 1
(CHEB) for Chebyshev basis. Default is 1 (CHEB).
)doc")
      .def_readwrite("HOT_SPLIT", &SICM::Options::HOT_SPLIT, R"doc(
Strategy for allocating higher-order terms (HOT) arising in sparse
products to the interval coefficients of the monomials, as a
``SICModel.Options.ALLOCATION`` value: NONE, SIMPLE or FULL. Default
is FULL.
)doc")
      .def_readwrite("PRODBND_SPLIT", &SICM::Options::PRODBND_SPLIT, R"doc(
Whether to distribute the uncertainty over the monomial coefficients
when multiplying a model by an interval bound (bool). Default is
False.
)doc")
      .def_readwrite("LIFT_USE", &SICM::Options::LIFT_USE, R"doc(
Whether to lift the uncertainty in nonlinear operations by
introducing auxiliary variables in the model (bool). Default is
False.
)doc")
      .def_readwrite("LIFT_ATOL", &SICM::Options::LIFT_ATOL, R"doc(
Absolute tolerance for lifting the uncertainty; only used if
``LIFT_USE`` is True. Default is 1e-10.
)doc")
      .def_readwrite("LIFT_RTOL", &SICM::Options::LIFT_RTOL, R"doc(
Relative tolerance for lifting the uncertainty; only used if
``LIFT_USE`` is True. Default is 1e-3.
)doc")
      .def_readwrite("REMEZ_USE", &SICM::Options::REMEZ_USE, R"doc(
Whether to use the Remez algorithm for computing a minimax
approximation of univariate terms (bool). Default is True.
)doc")
      .def_readwrite("REMEZ_MAXIT", &SICM::Options::REMEZ_MAXIT, R"doc(
Maximal number of iterations in the Remez algorithm for computing a
minimax approximation of univariate terms (int). Default is 10.
)doc")
      .def_readwrite("REMEZ_TOL", &SICM::Options::REMEZ_TOL, R"doc(
Stopping tolerance in the Remez algorithm for computing a minimax
approximation of univariate terms (float). Default is 1e-5.
)doc")
      .def_readwrite("REMEZ_MIG", &SICM::Options::REMEZ_MIG, R"doc(
Threshold on the interval width below which the Remez algorithm is
not used (float). Default is 1e-10.
)doc")
      .def_readwrite("INTERP_EXTRA", &SICM::Options::INTERP_EXTRA, R"doc(
Number of extra terms in the Chebyshev interpolation of univariate
functions (int): 0 uses a Chebyshev interpolation of the model order;
extra terms allow approximating the truncated Chebyshev series.
Default is 0.
)doc")
      .def_readwrite("INTERP_THRES", &SICM::Options::INTERP_THRES, R"doc(
Threshold below which coefficients in the Chebyshev interpolation of
univariates are ignored (float). Default is 1e2 times the machine
epsilon.
)doc")
      .def_readwrite("BOUNDER_TYPE", &SICM::Options::BOUNDER_TYPE, R"doc(
Chebyshev model range bounder, as a ``SICModel.Options.BOUNDER``
value: NAIVE, LSB or BERNSTEIN. Default is LSB.
)doc")
      .def_readwrite("BERNSTEIN_ORDER", &SICM::Options::BERNSTEIN_ORDER, R"doc(
Degree of the Bernstein basis for range bounding when
``BOUNDER_TYPE`` is BERNSTEIN. Default is 0 (same as the model
order).
)doc")
      .def_readwrite("MIG_USE", &SICM::Options::MIG_USE, R"doc(
Whether to simplify monomial terms with small magnitude during
propagation (bool). Default is False.
)doc")
      .def_readwrite("MIG_ATOL", &SICM::Options::MIG_ATOL, R"doc(
Absolute tolerance for simplifying monomial terms; only used if
``MIG_USE`` is True. Default is 0.
)doc")
      .def_readwrite("MIG_RTOL", &SICM::Options::MIG_RTOL, R"doc(
Relative tolerance for simplifying monomial terms; only used if
``MIG_USE`` is True. Default is the machine epsilon.
)doc")
      .def_readwrite("MIXED_IA", &SICM::Options::MIXED_IA, R"doc(
Whether to intersect the internal model bounds with bounds propagated
in the underlying (interval) arithmetic (bool). Default is False.
)doc")
      .def_readwrite("REF_POLY", &SICM::Options::REF_POLY, R"doc(
Scalar in [0,1] selecting the polynomial part in the functions
``inter`` and ``hull``: 0 selects the polynomial part of the left
operand, 1 that of the right operand, intermediate values
interpolate. Default is 0.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &SICM::Options::DISPLAY_DIGITS, R"doc(
Number of digits used when printing Chebyshev model coefficients
(int). Default is 7.
)doc");


  py::enum_<SICM::Options::ALLOCATION>(pySICModelOptions, "ALLOCATION")
      .value("NONE", SICM::Options::ALLOCATION::NONE,
             "No split of HOT; allocate to constant coefficient")
      .value("SIMPLE", SICM::Options::ALLOCATION::SIMPLE,
             "Split HOT between current term and new variable")
      .value("FULL", SICM::Options::ALLOCATION::FULL,
             "Split HOT among all variables")
      .export_values();

  py::enum_<SICM::Options::MONBASIS>(pySICModelOptions, "MONBASIS")
      .value("MONOM", SICM::Options::MONBASIS::MONOM, "Monomial basis")
      .value("CHEB", SICM::Options::MONBASIS::CHEB, "Chebyshev basis")
      .export_values();

  // Nested class Exceptions
  py::class_<SICM::Exceptions> pySICModelExceptions(pySICModel, "Exceptions",
                                                    R"doc(
Exception data thrown by sparse interval Chebyshev model arithmetic.

Carries an error code (``ierr``) and a description (``what``); see the
``SICModel.Exceptions.TYPE`` enumeration for the possible errors.
)doc");

  py::enum_<SICM::Exceptions::TYPE>(pySICModelExceptions, "TYPE")
      .value("DIV", SICM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", SICM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", SICM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", SICM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("DPOW", SICM::Exceptions::TYPE::DPOW,
             "Real power operation with negative numbers in range")
      .value("TAN", SICM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", SICM::Exceptions::TYPE::ACOS,
             "Cosine inverse operation with range outside [-1,1]")
      .value("ASIN", SICM::Exceptions::TYPE::ASIN,
             "Sine inverse operation with range outside [-1,1]")
      .value("COMPOSE", SICM::Exceptions::TYPE::COMPOSE,
             "Failed to compose Chebyshev variable")
      .value("INIT", SICM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", SICM::Exceptions::TYPE::INCON,
             "Inconsistent bounds with template parameter arithmetic")
      .value("MODEL", SICM::Exceptions::TYPE::MODEL,
             "Operation between variables linked to different models")
      .value("INTERNAL", SICM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", SICM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented")
      .export_values();

  pySICModelExceptions.def(py::init<SICM::Exceptions::TYPE>())
      .def("ierr", &SICM::Exceptions::ierr, "Return the error code (int).")
      .def("what", &SICM::Exceptions::what,
           "Return the error description (str).");
}
