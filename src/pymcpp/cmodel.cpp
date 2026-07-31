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

#include "cmodel.hpp"

namespace py = pybind11;

void
mc_cmodel(py::module& m)
{
  typedef mc::CModel<I> CM;
  typedef mc::CVar<I> CV;

  py::class_<CM> pyCModel(m, "CModel", R"doc(
Chebyshev model environment (dense storage).

A q-th order Chebyshev model of a function f over an interval domain D
consists of a q-th order multivariate polynomial P expressed in the
Chebyshev basis (products of Chebyshev polynomials of the variables
scaled to [-1,1]), plus an interval remainder bound R such that f(x)
belongs to P(x) + R for every x in D. Compared with a Taylor model of
the same order (``TModel``), the Chebyshev basis typically yields much
tighter remainder bounds over wide variable domains. The polynomial
coefficients are propagated as floating-point numbers, so the enclosure
is guaranteed up to rounding in these coefficients (the implementation
is not fully verified).

A ``CModel`` fixes the number of variables and the expansion order, and
stores the propagation options (attribute ``options``). Variables of
type ``CVar`` are created by attaching them to an environment together
with a variable index and a domain interval; expressions of these
variables are then evaluated in Chebyshev model arithmetic via the
overloaded operators and the module-level functions (``exp``, ``log``,
``sin``, ...).

All coefficients up to the model order are stored (dense storage), in
graded lexicographic order: monomial terms are grouped by increasing
total degree, and sorted lexicographically by exponent tuple within
each degree.

Related environments: ``TModel`` (dense Taylor expansion), ``SCModel``
(sparse Chebyshev), ``SICModel`` (sparse Chebyshev with interval
coefficients).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.CModel(2, 4)          # 2 variables, order 4
>>> x = pymcpp.CVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.CVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2   # Chebyshev model of f(x,y)
>>> B = f.B                            # enclosure of the range of f
>>> R = f.R                            # remainder bound
)doc");

  py::class_<CM::Options> pyCModelOptions(pyCModel, "Options", R"doc(
Options of a Chebyshev model environment.

Modify the fields of ``CModel.options`` in place to control range
bounding and Chebyshev model propagation, e.g.
``env.options.BOUNDER_TYPE = pymcpp.CModel.Options.BOUNDER.NAIVE``.
)doc");

  py::enum_<CM::Options::BOUNDER>(pyCModelOptions, "BOUNDER")
      .value("NAIVE", CM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", CM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .value("EIGEN", CM::Options::BOUNDER::EIGEN,
             "Eigenvalue decomposition-based bounder")
      .value("BERNSTEIN", CM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .value("HYBRID", CM::Options::BOUNDER::HYBRID,
             "Hybrid LSB + EIGEN range bounder")
      .export_values();

  pyCModel
      // Constructors
      .def(py::init<unsigned const, unsigned const, bool const>(),
           py::arg("nvar"), py::arg("maxord"), py::arg("sparse") = false,
           R"doc(
Construct a Chebyshev model environment.

Parameters
----------
nvar : int
    Number of variables in the model.
maxord : int
    Order of the polynomial inclusion (maximal total degree of the
    polynomial part).
sparse : bool, optional
    Whether to enable sparse operations internally. Default is False.
)doc")
      // Accessors
      .def_readwrite("options", &CM::options, R"doc(
Propagation and bounding options of this environment (in-place
modifiable ``CModel.Options`` instance).
)doc")
      // Inherited from PolyModel
      .def_property_readonly("sparse", &CM::sparse, R"doc(
Whether sparse operations are enabled in this environment (bool).
)doc")
      .def_property_readonly("nvar", &CM::nvar, R"doc(
Number of variables in the model environment (int).
)doc")
      .def_property_readonly("nord", &CM::nord, R"doc(
Order of the model environment, i.e. maximal total degree of the
polynomial part (int).
)doc")
      .def_property_readonly("nmon", &CM::nmon, R"doc(
Total number of monomial terms in the polynomial part (int), equal to
binomial(nvar + nord, nord).
)doc")
      .def(
          "expmon",
          [](CM const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), R"doc(
Return the variable exponents of monomial term ``i``.

Monomial terms are indexed 0 to ``nmon``-1 in graded lexicographic
order: grouped by increasing total degree, and sorted lexicographically
by exponent tuple within each degree. An exponent ``e`` for variable
``k`` denotes the Chebyshev polynomial T_e of the scaled variable
``k``.

Parameters
----------
i : int
    Monomial index, in 0 to ``nmon``-1.

Returns
-------
iexp : list of int
    New list of length ``nvar``; ``iexp[k]`` is the Chebyshev degree
    of variable ``k`` in monomial term ``i``.

Raises
------
IndexError
    If ``i`` is not smaller than ``nmon``.
)doc")
      .def(
          "loc_expmon",
          [](CM const& self, std::vector<unsigned> const& iexp)
          {
            if (iexp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.loc_expmon(iexp.data());
          },
          py::arg("iexp"), R"doc(
Return the index of the monomial term with given variable exponents.

Inverse of ``expmon``: ``loc_expmon(expmon(i)) == i``.

Parameters
----------
iexp : list of int
    Variable exponents, one per variable (length ``nvar``).

Returns
-------
i : int
    Index of the corresponding monomial term.

Raises
------
ValueError
    If ``iexp`` does not have length ``nvar``.
)doc")
      .def_property_readonly(
          "posord",
          [](CM const& self)
          {
            const unsigned* p = self.posord();
            return std::vector<unsigned>(p, p + self.nord() + 2);
          },
          R"doc(
Start index of each total-degree group of monomial terms.

Returns a new list of length ``nord``+2 whose entry ``i`` (for i in
0 to ``nord``) is the index of the first monomial term of total order
``i``; the last entry equals ``nmon``. Hence the terms of total order
``i`` occupy indices ``posord[i]`` to ``posord[i+1]``-1.
)doc")
      .def("get_binom", &CM::get_binom, py::arg("n"), py::arg("k"), R"doc(
Return the binomial coefficient n-choose-k (int).
)doc")
      // CModel specific
      .def("bndvar", &CM::bndvar, py::arg("i"), R"doc(
Return the domain interval of model variable ``i`` (Interval).
)doc")
      .def("refvar", &CM::refvar, py::arg("i"), R"doc(
Return the reference point of model variable ``i`` (float), i.e. the
midpoint of its domain used for scaling to [-1,1].
)doc")
      .def("scalvar", &CM::scalvar, py::arg("i"), R"doc(
Return the scaling factor of model variable ``i`` (float), i.e. the
half-diameter of its domain used for scaling to [-1,1].
)doc")
      //.def(
      //  "reset",
      //  &CM::reset,
      //  "reset the bounds on Chebyshev basis functions"
      //)
      .def(
          "lift", [](CM const& self, unsigned const n, std::vector<I> const& X)
          { return self.lift(n, X.data()); }, py::arg("n"), py::arg("X"),
          py::return_value_policy::take_ownership,
          R"doc(
Return a new Chebyshev model environment with ``n`` extra variables.

The new environment has the same order, the first ``nvar`` variables
keep their current domains, and the appended variables take the domains
in ``X``. This environment is not modified.

Parameters
----------
n : int
    Number of extra variables.
X : list of Interval
    Domains of the appended variables.

Returns
-------
env : CModel
    New environment with ``nvar + n`` variables.
)doc");

  py::class_<CV> pyCVar(m, "CVar", R"doc(
Chebyshev model variable (dense storage).

A ``CVar`` holds a multivariate polynomial in the Chebyshev basis
(dense coefficient array over products of Chebyshev polynomials of the
variables scaled to [-1,1]) together with an interval remainder bound,
valid over the variable domains recorded in its ``CModel`` environment.
Participating variables are created by attaching an environment, a
variable index and a domain interval; composite models are obtained via
the overloaded arithmetic operators (+, -, *, /, **, unary -, abs) and
the module-level functions (``exp``, ``log``, ``sin``, ...), which all
propagate both the polynomial part and the remainder bound.

Operations between ``CVar`` objects linked to different environments
raise an error (``CModel.Exceptions``).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.CModel(2, 4)
>>> x = pymcpp.CVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.CVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2
>>> f.B          # range enclosure using the selected bounder
[ -1.2603052508837877e+01 :  1.4000743049036982e+01 ]
>>> f.R          # remainder bound
[ -6.8904566365688225e-01 :  6.8904566365688225e-01 ]
>>> f.P([0.5, 1.5])  # value of the polynomial part at a point
-0.05770280538083117
)doc");

  pyCVar
      // Constructors
      .def(py::init<double const>(), py::arg("cst") = 0., R"doc(
Construct a constant Chebyshev variable from a scalar.

Parameters
----------
cst : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<I const&>(), py::arg("bnd"), R"doc(
Construct a constant Chebyshev variable from an interval.

The interval becomes the remainder bound; the polynomial part is zero.

Parameters
----------
bnd : Interval
    Constant enclosure.
)doc")
      .def(py::init<CM*, unsigned const, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"), R"doc(
Construct the Chebyshev variable of index ``ndx`` in environment
``env``.

Parameters
----------
env : CModel
    Chebyshev model environment to attach to.
ndx : int
    Variable index in the environment, in 0 to ``env.nvar``-1.
rng : Interval
    Domain (range) of the variable.
)doc")
      .def(py::init<CV const&>(), R"doc(
Copy constructor.
)doc")
      // Modifiers
      .def(
          "set", [](CV& self, CM* env, unsigned const ix, I const& X) -> CV&
          { return self.set(env, ix, X); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as the variable of index ``ndx`` in
environment ``env``, with domain ``rng``. Returns this object, so that
calls can be chained.
)doc")
      .def(
          "set", [](CV& self, CM* env, bool const reset) -> CV&
          { return self.set(env, reset); }, py::arg("env"),
          py::arg("reset") = false, py::return_value_policy::reference_internal,
          R"doc(
Attach this variable to environment ``env``.

Parameters
----------
env : CModel
    New model environment.
reset : bool, optional
    If True, reset all polynomial coefficients and the remainder to
    zero. Default is False.
)doc")
      // Accessors
      .def("env", &CV::env, R"doc(
Return the model environment this variable is attached to
(``CModel``, or None for a constant).
)doc")
      // Inherited from PolyVar
      .def_property_readonly("nord", &CV::nord, R"doc(
Order of the model environment (int); 0 for a constant.
)doc")
      .def_property_readonly("nvar", &CV::nvar, R"doc(
Number of variables in the model environment (int); 0 for a constant.
)doc")
      .def_property_readonly("nmon", &CV::nmon, R"doc(
Total number of monomial terms in the polynomial part (int).
)doc")
      .def(
          "coefmon",
          [](CV const& self)
          {
            // Return all coefficients as a list/vector
            // Note: dense CVar stores them contiguously
            return std::vector<double>(self.coefmon().second,
                                       self.coefmon().second + self.nmon());
          },
          R"doc(
Return all monomial coefficients of the polynomial part.

Returns
-------
coef : list of float
    New list of length ``nmon``; ``coef[i]`` is the coefficient of
    Chebyshev monomial term ``i`` in graded lexicographic order (the
    exponents of term ``i`` are given by ``expmon(i)``). Term ``i``
    is the product over all variables ``k`` of the Chebyshev
    polynomial of degree ``expmon(i)[k]`` in variable ``k`` scaled to
    [-1,1].
)doc")
      .def(
          "coefmon",
          [](CV const& self, std::vector<unsigned> const& exp)
          {
            if (exp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.coefmon(exp.data());
          },
          py::arg("exp"), R"doc(
Return the coefficient of the Chebyshev monomial with variable
exponents ``exp``.

Parameters
----------
exp : list of int
    Chebyshev degrees, one per variable (length ``nvar``).

Returns
-------
coef : float
    Coefficient of that monomial term (0.0 if absent).

Raises
------
ValueError
    If ``exp`` does not have length ``nvar``.
)doc")
      .def(
          "expmon",
          [](CV const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), R"doc(
Return the variable exponents of monomial term ``i``.

Returns
-------
iexp : list of int
    New list of length ``nvar``; ``iexp[k]`` is the Chebyshev degree
    of variable ``k`` in monomial term ``i`` (graded lexicographic
    ordering, see ``CModel.expmon``).

Raises
------
IndexError
    If ``i`` is not smaller than ``nmon``.
)doc")
      .def_property_readonly("ndxmon", &CV::ndxmon, R"doc(
Indices of the nonzero monomial terms in sparse representation
(set of int). Empty when the environment uses dense storage.
)doc")
      .def_property_readonly(
          "bound", [](CV const& self) { return self.bound(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus the
remainder bound. Same as ``B``.
)doc")
      .def_property_readonly(
          "B", [](CV const& self) { return self.B(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus the
remainder bound. Same as ``bound``.
)doc")
      .def(
          "bndpol", [](CV const& self) { return self.bndpol(); },
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using the bounder selected in option ``BOUNDER_TYPE``. The remainder
bound is not added.
)doc")
      .def(
          "bndpol", [](CV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using bounder ``type`` (a ``CModel.Options.BOUNDER`` value).
)doc")
      .def("bndord", &CV::bndord, py::arg("minord"), R"doc(
Return a bound on the sum of all monomial terms whose total order
equals ``minord`` (Interval). Returns [0,0] if ``minord`` exceeds the
model order.
)doc")
      .def(
          "P", [](CV const& self, std::vector<double> const& x)
          { return self.P(x.data()); }, py::arg("x"),
          R"doc(
Evaluate the polynomial part at point ``x``.

Parameters
----------
x : list of float
    Coordinates in the original (unscaled) variable domain, one per
    variable (length ``nvar``), ordered by variable index.

Returns
-------
val : float
    Value of the polynomial part at ``x`` (remainder not included).
)doc")
      .def(
          "P", [](CV const& self) { return self.P(); },
          R"doc(
Return a new Chebyshev variable with the same polynomial part and a
zero remainder bound. This variable is not modified.
)doc")
      .def_property_readonly(
          "R", [](CV const& self) { return self.R(); }, R"doc(
Remainder bound of the Chebyshev model (Interval).
)doc")
      .def("center", &CV::center, R"doc(
Center the remainder term in place: the midpoint of the remainder bound
is added to the constant polynomial coefficient and subtracted from the
remainder. Returns this object.
)doc")
      .def("C", &CV::C, R"doc(
Center the remainder term in place; shortcut for ``center``.
)doc")
      .def("constant", &CV::constant, py::arg("reset") = false, R"doc(
Return the coefficient of the constant term (float).

Parameters
----------
reset : bool, optional
    If True, also reset the constant term to zero in this variable.
    Default is False.
)doc")
      .def(
          "linear", [](CV& self, unsigned const ivar, bool const reset)
          { return self.linear(ivar, reset); }, py::arg("id"),
          py::arg("reset") = false, R"doc(
Return the coefficient of the linear term in variable ``id`` (float),
expressed in the original (unscaled) variable. Returns 0.0 if ``id`` is
out of range or the model order is 0.

Parameters
----------
id : int
    Variable index.
reset : bool, optional
    If True, also reset this linear coefficient to zero. Default is
    False.
)doc")
      //.def(
      //  "scale",
      //  []( CV const& self, std::vector<I> const& X )
      //    { return self.scale( X.data() ); },
      //  py::arg("dom"),
      //  "Scale coefficients in Chebyshev variable for the reduced variable
      //  dom"
      //)
      .def("__str__",
           [](CV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](CV const& self)
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
      .def("__abs__", [](CV const& m) { return mc::Op<CV>::abs(m); })
      .def("__pow__",
           [](CV const& m, int const n) { return mc::Op<CV>::pow(m, n); })
      .def("__pow__",
           [](CV const& m, double const& r) { return mc::Op<CV>::pow(m, r); })
      .def("__pow__",
           [](CV const& m, CV const& mm) { return mc::Op<CV>::pow(m, mm); })
      .def("__pow__",
           [](double const& r, CV const& m) { return mc::Op<CV>::pow(r, m); });

  m.def(
      "inv", [](CV const& x) { return mc::inv(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of 1/x.");
  m.def(
      "sqr", [](CV const& x) { return mc::sqr(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of x**2.");
  m.def(
      "sqrt", [](CV const& x) { return mc::sqrt(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of sqrt(x).");
  m.def(
      "exp", [](CV const& x) { return mc::exp(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of exp(x).");
  m.def(
      "log", [](CV const& x) { return mc::log(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of log(x).");
  m.def(
      "xlog", [](CV const& x) { return mc::xlog(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of x*log(x).");
  m.def(
      "cos", [](CV const& x) { return mc::cos(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of cos(x).");
  m.def(
      "sin", [](CV const& x) { return mc::sin(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of sin(x).");
  m.def(
      "tan", [](CV const& x) { return mc::tan(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of tan(x).");
  m.def(
      "acos", [](CV const& x) { return mc::acos(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of acos(x).");
  m.def(
      "asin", [](CV const& x) { return mc::asin(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of asin(x).");
  m.def(
      "atan", [](CV const& x) { return mc::atan(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of atan(x).");
  m.def(
      "cosh", [](CV const& x) { return mc::cosh(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of cosh(x).");
  m.def(
      "sinh", [](CV const& x) { return mc::sinh(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of sinh(x).");
  m.def(
      "tanh", [](CV const& x) { return mc::tanh(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of tanh(x).");
  m.def(
      "pow", [](CV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"),
      "Dense Chebyshev model overload: Chebyshev model of x**n for integer "
      "n.");
  m.def(
      "pow", [](CV const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"),
      "Dense Chebyshev model overload: Chebyshev model of x**r for real r.");
  m.def(
      "pow", [](CV const& x, CV const& y) { return mc::pow(x, y); },
      py::arg("x"), py::arg("y"),
      "Dense Chebyshev model overload: Chebyshev model of x**y.");
  m.def(
      "pow", [](double const& r, CV const& y) { return mc::pow(r, y); },
      py::arg("r"), py::arg("y"),
      "Dense Chebyshev model overload: Chebyshev model of r**y for real r.");
  m.def(
      "cheb", [](CV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"),
      "Dense Chebyshev model overload: Chebyshev model of the Chebyshev "
      "polynomial T_n(x).");
  m.def(
      "fabs", [](CV const& x) { return mc::fabs(x); }, py::arg("x"),
      "Dense Chebyshev model overload: Chebyshev model of abs(x).");
  m.def(
      "hull", [](CV const& x, CV const& y) { return mc::hull(x, y); },
      py::arg("x"), py::arg("y"),
      "Dense Chebyshev model overload: Chebyshev model enclosing the union "
      "of the enclosures x and y (see option REF_POLY for the polynomial "
      "part).");
  m.def(
      "inter",
      [](CV& xy, CV const& x, CV const& y) { return mc::inter(xy, x, y); },
      py::arg("xy"), py::arg("x"), py::arg("y"),
      "Dense Chebyshev model overload: overwrite xy with a Chebyshev model "
      "of the intersection of the enclosures x and y; returns False if the "
      "intersection is empty.");

  // Nested class Options


  pyCModelOptions.def(py::init<>(), "Construct options with default values.")
      .def(py::init<CM::Options const&>(), "Copy constructor.")
      .def("reset", &CM::Options::reset, R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("INTERP_EXTRA", &CM::Options::INTERP_EXTRA, R"doc(
Number of extra terms in the Chebyshev interpolation of univariate
functions (int): 0 uses a Chebyshev interpolation of the model order;
extra terms allow approximating the truncated Chebyshev series.
Default is 0.
)doc")
      .def_readwrite("INTERP_THRES", &CM::Options::INTERP_THRES, R"doc(
Threshold below which coefficients in the Chebyshev expansion of
transcendental univariates are ignored during bounding (float).
Default is 1e2 times the machine epsilon.
)doc")
      .def_readwrite("BOUNDER_TYPE", &CM::Options::BOUNDER_TYPE, R"doc(
Chebyshev model range bounder, as a ``CModel.Options.BOUNDER`` value:
NAIVE, LSB, EIGEN, BERNSTEIN or HYBRID. Default is LSB.
)doc")
      .def_readwrite("BERNSTEIN_ORDER", &CM::Options::BERNSTEIN_ORDER, R"doc(
Degree of the Bernstein basis for range bounding when ``BOUNDER_TYPE``
is BERNSTEIN. Default is 0 (same as the model order).
)doc")
      .def_readwrite("MIXED_IA", &CM::Options::MIXED_IA, R"doc(
Whether to intersect the internal model bounds with bounds propagated
in the underlying (interval) arithmetic (bool). Default is False.
)doc")
      .def_readwrite("REF_POLY", &CM::Options::REF_POLY, R"doc(
Scalar in [0,1] selecting the polynomial part in the functions
``inter`` and ``hull``: 0 selects the polynomial part of the left
operand, 1 that of the right operand, intermediate values interpolate.
Default is 0.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &CM::Options::DISPLAY_DIGITS, R"doc(
Number of digits used when printing Chebyshev model coefficients (int).
Default is 7.
)doc");


  // Nested class Exceptions
  py::class_<CM::Exceptions> pyCModelExceptions(pyCModel, "Exceptions", R"doc(
Exception data thrown by dense Chebyshev model arithmetic.

Carries an error code (``ierr``) and a description (``what``); see the
``CModel.Exceptions.TYPE`` enumeration for the possible errors.
)doc");

  py::enum_<CM::Exceptions::TYPE>(pyCModelExceptions, "TYPE")
      .value("DIV", CM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", CM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", CM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", CM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("TAN", CM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", CM::Exceptions::TYPE::ACOS,
             "Sine/Cosine inverse operation with range outside [-1,1]")
      .value("EIGEN", CM::Exceptions::TYPE::EIGEN,
             "Failed to compute eigenvalue decomposition in range bounder "
             "CModel::Options::EIGEN")
      .value("INIT", CM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", CM::Exceptions::TYPE::INCON,
             "Chebyshev model bound does not intersect with bound in template "
             "parameter arithmetic")
      .value("CMODEL", CM::Exceptions::TYPE::CMODEL,
             "Operation between Chebyshev variables linked to different "
             "Chebyshev models")
      .value("INTERNAL", CM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", CM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented in module")
      .export_values();

  pyCModelExceptions.def(py::init<CM::Exceptions::TYPE>())
      .def("ierr", &CM::Exceptions::ierr, "Return the error code (int).")
      .def("what", &CM::Exceptions::what, "Return the error description (str).");
}
