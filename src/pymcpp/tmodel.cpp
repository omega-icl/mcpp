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

#include "tmodel.hpp"

namespace py = pybind11;

void
mc_tmodel(py::module& m)
{
  typedef mc::TModel<I> TM;
  typedef mc::TVar<I> TV;

  py::class_<TM> pyTModel(m, "TModel", R"doc(
Taylor model environment (dense storage).

A q-th order Taylor model of a function f over an interval domain D
consists of a q-th order multivariate Taylor polynomial P, expanded
around a reference point of D, plus an interval remainder bound R such
that f(x) belongs to P(x) + R for every x in D. The polynomial part is
propagated symbolically through arithmetic operations, while the
remainder is bounded using interval arithmetic; the enclosure is
guaranteed up to floating-point rounding in the polynomial coefficients
(the implementation is not fully verified).

A ``TModel`` fixes the number of variables and the expansion order, and
stores the propagation options (attribute ``options``). Variables of
type ``TVar`` are created by attaching them to an environment together
with a variable index and a domain interval; expressions of these
variables are then evaluated in Taylor model arithmetic simply by using
the overloaded operators and the module-level functions (``exp``,
``log``, ``sin``, ...).

All coefficients up to the model order are stored (dense storage), in
graded lexicographic order: monomial terms are grouped by increasing
total degree, and sorted lexicographically by exponent tuple within
each degree.

Related environments: ``CModel`` (dense Chebyshev basis instead of
Taylor expansion), ``SCModel`` (sparse Chebyshev), ``SICModel`` (sparse
Chebyshev with interval coefficients).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.TModel(2, 4)          # 2 variables, order 4
>>> x = pymcpp.TVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.TVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2   # Taylor model of f(x,y)
>>> B = f.B                            # enclosure of the range of f
>>> R = f.R                            # remainder bound
)doc");

  py::class_<TM::Options> pyTModelOptions(pyTModel, "Options", R"doc(
Options of a Taylor model environment.

Modify the fields of ``TModel.options`` in place to control range
bounding and Taylor model propagation, e.g.
``env.options.BOUNDER_TYPE = pymcpp.TModel.Options.BOUNDER.NAIVE``.
)doc");

  py::enum_<TM::Options::BOUNDER>(pyTModelOptions, "BOUNDER")
      .value("NAIVE", TM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", TM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .value("EIGEN", TM::Options::BOUNDER::EIGEN,
             "Eigenvalue decomposition-based bounder")
      .value("BERNSTEIN", TM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .value("HYBRID", TM::Options::BOUNDER::HYBRID,
             "Hybrid LSB + EIGEN range bounder")
      .export_values();

  pyTModel
      // Constructors
      .def(py::init<unsigned const, unsigned const>(), py::arg("nvar"),
           py::arg("maxord"), R"doc(
Construct a Taylor model environment.

Parameters
----------
nvar : int
    Number of variables in the model.
maxord : int
    Order of the Taylor expansion (maximal total degree of the
    polynomial part).
)doc")
      // Accessors
      .def_readwrite("options", &TM::options, R"doc(
Propagation and bounding options of this environment (in-place
modifiable ``TModel.Options`` instance).
)doc")
      // Inherited from PolyModel
      .def_property_readonly("sparse", &TM::sparse, R"doc(
Whether sparse operations are enabled in this environment (bool).
)doc")
      .def_property_readonly("nvar", &TM::nvar, R"doc(
Number of variables in the model environment (int).
)doc")
      .def_property_readonly("nord", &TM::nord, R"doc(
Order of the model environment, i.e. maximal total degree of the
polynomial part (int).
)doc")
      .def_property_readonly("nmon", &TM::nmon, R"doc(
Total number of monomial terms in the polynomial part (int), equal to
binomial(nvar + nord, nord).
)doc")
      .def(
          "expmon",
          [](TM const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), R"doc(
Return the variable exponents of monomial term ``i``.

Monomial terms are indexed 0 to ``nmon``-1 in graded lexicographic
order: grouped by increasing total degree, and sorted lexicographically
by exponent tuple within each degree.

Parameters
----------
i : int
    Monomial index, in 0 to ``nmon``-1.

Returns
-------
iexp : list of int
    New list of length ``nvar``; ``iexp[k]`` is the exponent of
    variable ``k`` in monomial term ``i``.

Raises
------
IndexError
    If ``i`` is not smaller than ``nmon``.
)doc")
      .def(
          "loc_expmon",
          [](TM const& self, std::vector<unsigned> const& iexp)
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
          [](TM const& self)
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
      .def("get_binom", &TM::get_binom, py::arg("n"), py::arg("k"), R"doc(
Return the binomial coefficient n-choose-k (int).
)doc")
      // TModel specific and inherited from PolyModel
      .def(
          "bndvar",
          [](TM const& self, unsigned const i)
          {
            if (i >= self.nvar()) throw py::index_error();
            return self.bndvar()[i];
          },
          py::arg("i"), R"doc(
Return the domain interval of model variable ``i`` (Interval).

Raises
------
IndexError
    If ``i`` is not smaller than ``nvar``.
)doc")
      .def("refvar", &TM::refvar, py::arg("i"), R"doc(
Return the reference (expansion) point of model variable ``i`` (float).
)doc")
      .def("scalvar", &TM::scalvar, py::arg("i"), R"doc(
Return the scaling factor of model variable ``i`` (float).

Equals half the domain diameter when option ``SCALE_VARIABLES`` is
True, and 1 otherwise.
)doc")
      .def("reset", &TM::reset, R"doc(
Reset the internally stored bounds on powers of the variable ranges.
)doc")
      .def(
          "lift", [](TM const& self, unsigned const n, std::vector<I> const& X)
          { return self.lift(n, X.data()); }, py::arg("n"), py::arg("X"),
          py::return_value_policy::take_ownership,
          R"doc(
Return a new Taylor model environment with ``n`` extra variables.

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
env : TModel
    New environment with ``nvar + n`` variables.
)doc");

  py::class_<TV> pyTVar(m, "TVar", R"doc(
Taylor model variable.

A ``TVar`` holds a multivariate Taylor polynomial (dense coefficient
array) together with an interval remainder bound, valid over the
variable domains recorded in its ``TModel`` environment. Participating
variables are created by attaching an environment, a variable index and
a domain interval; composite models are obtained via the overloaded
arithmetic operators (+, -, *, /, **, unary -, abs) and the
module-level functions (``exp``, ``log``, ``sin``, ...), which all
propagate both the polynomial part and the remainder bound.

Operations between ``TVar`` objects linked to different environments
raise an error (``TModel.Exceptions``).

Examples
--------
>>> import pymcpp
>>> env = pymcpp.TModel(2, 4)
>>> x = pymcpp.TVar(env, 0, pymcpp.Interval(-2, 1))
>>> y = pymcpp.TVar(env, 1, pymcpp.Interval(-1, 2))
>>> f = x * (pymcpp.exp(x) - y) ** 2
>>> f.B          # range enclosure using the selected bounder
[ -1.9163792263857939e+01 :  1.8327270962414772e+01 ]
>>> f.R          # remainder bound
[ -6.9054396969907295e+00 :  6.2123279370094693e+00 ]
>>> f.P([0.5, 1.5])  # value of the polynomial part at a point
-0.04398461354988953
)doc");

  pyTVar
      // Constructors
      .def(py::init<double const>(), py::arg("cst") = 0., R"doc(
Construct a constant Taylor variable from a scalar.

Parameters
----------
cst : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<I const&>(), py::arg("bnd"), R"doc(
Construct a constant Taylor variable from an interval.

The interval becomes the remainder bound; the polynomial part is zero.

Parameters
----------
bnd : Interval
    Constant enclosure.
)doc")
      .def(py::init<TM*, unsigned const, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"), R"doc(
Construct the Taylor variable of index ``ndx`` in environment ``env``.

The expansion reference point is the midpoint of ``rng``.

Parameters
----------
env : TModel
    Taylor model environment to attach to.
ndx : int
    Variable index in the environment, in 0 to ``env.nvar``-1.
rng : Interval
    Domain (range) of the variable.
)doc")
      .def(py::init<TM*, unsigned const, I const&, double const>(),
           py::arg("env"), py::arg("ndx"), py::arg("rng"), py::arg("ref"),
           R"doc(
Construct the Taylor variable of index ``ndx`` in environment ``env``
with expansion reference point ``ref``.

Parameters
----------
env : TModel
    Taylor model environment to attach to.
ndx : int
    Variable index in the environment, in 0 to ``env.nvar``-1.
rng : Interval
    Domain (range) of the variable.
ref : float
    Reference point of the Taylor expansion; must lie in ``rng``.
)doc")
      .def(py::init<TV const&>(), R"doc(
Copy constructor.
)doc")
      // Modifiers
      .def(
          "set",
          [](TV& self, TM* env, unsigned const ix, I const& X,
             double const& ref) -> TV& { return self.set(env, ix, X, ref); },
          py::arg("env"), py::arg("ndx"), py::arg("rng"), py::arg("ref"),
          py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as the variable of index ``ndx`` in
environment ``env``, with domain ``rng`` and reference point ``ref``.

Returns this object, so that calls can be chained.
)doc")
      .def(
          "set", [](TV& self, TM* env, unsigned const ix, I const& X) -> TV&
          { return self.set(env, ix, X); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal, R"doc(
Redefine this object in place as the variable of index ``ndx`` in
environment ``env``, with domain ``rng`` (reference point at the
midpoint of ``rng``). Returns this object.
)doc")
      .def(
          "set", [](TV& self, TM* env, bool const reset) -> TV&
          { return self.set(env, reset); }, py::arg("env"),
          py::arg("reset") = false, py::return_value_policy::reference_internal,
          R"doc(
Attach this variable to environment ``env``.

Parameters
----------
env : TModel
    New model environment.
reset : bool, optional
    If True, reset all polynomial coefficients and the remainder to
    zero. Default is False.
)doc")
      // Accessors
      .def("env", &TV::env, R"doc(
Return the model environment this variable is attached to
(``TModel``, or None for a constant).
)doc")
      // Inherited from PolyVar
      .def_property_readonly("nord", &TV::nord, R"doc(
Order of the model environment (int); 0 for a constant.
)doc")
      .def_property_readonly("nvar", &TV::nvar, R"doc(
Number of variables in the model environment (int); 0 for a constant.
)doc")
      .def_property_readonly("nmon", &TV::nmon, R"doc(
Total number of monomial terms in the polynomial part (int).
)doc")
      .def(
          "coefmon",
          [](TV const& self)
          {
            // Return all coefficients as a list/vector
            // Note: dense TVar stores them contiguously
            return std::vector<double>(self.coefmon().second,
                                       self.coefmon().second + self.nmon());
          },
          R"doc(
Return all monomial coefficients of the polynomial part.

Returns
-------
coef : list of float
    New list of length ``nmon``; ``coef[i]`` is the coefficient of
    monomial term ``i`` in graded lexicographic order (the exponents
    of term ``i`` are given by ``expmon(i)``). The monomials are in
    the shifted variables x_k - refvar(k), divided by scalvar(k) when
    option ``SCALE_VARIABLES`` is enabled.
)doc")
      .def(
          "coefmon",
          [](TV const& self, std::vector<unsigned> const& exp)
          {
            if (exp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.coefmon(exp.data());
          },
          py::arg("exp"), R"doc(
Return the coefficient of the monomial with variable exponents ``exp``.

Parameters
----------
exp : list of int
    Variable exponents, one per variable (length ``nvar``).

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
          [](TV const& self, unsigned const i)
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
    New list of length ``nvar``; ``iexp[k]`` is the exponent of
    variable ``k`` in monomial term ``i`` (graded lexicographic
    ordering, see ``TModel.expmon``).

Raises
------
IndexError
    If ``i`` is not smaller than ``nmon``.
)doc")
      .def_property_readonly("ndxmon", &TV::ndxmon, R"doc(
Indices of the nonzero monomial terms in sparse representation
(set of int). Empty when the environment uses dense storage.
)doc")
      .def_property_readonly(
          "bound", [](TV const& self) { return self.bound(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus the
remainder bound. Same as ``B``.
)doc")
      .def_property_readonly(
          "B", [](TV const& self) { return self.B(); },
          R"doc(
Enclosure of the variable range (Interval): polynomial range bound,
computed with the bounder selected in option ``BOUNDER_TYPE``, plus the
remainder bound. Same as ``bound``.
)doc")
      .def(
          "bndpol", [](TV const& self) { return self.bndpol(); },
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using the bounder selected in option ``BOUNDER_TYPE``. The remainder
bound is not added.
)doc")
      .def(
          "bndpol", [](TV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          R"doc(
Return a bound on the multivariate polynomial part only (Interval),
using bounder ``type`` (a ``TModel.Options.BOUNDER`` value).
)doc")
      .def("bndord", &TV::bndord, py::arg("minord"), R"doc(
Return a bound on the sum of all monomial terms whose total order
equals ``minord`` (Interval). Returns [0,0] if ``minord`` exceeds the
model order.
)doc")
      .def(
          "P", [](TV const& self, std::vector<double> const& x)
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
          "P", [](TV const& self) { return self.P(); },
          R"doc(
Return a new Taylor variable with the same polynomial part and a zero
remainder bound. This variable is not modified.
)doc")
      .def_property_readonly(
          "R", [](TV const& self) { return self.R(); }, R"doc(
Remainder bound of the Taylor model (Interval).
)doc")
      .def("center", &TV::center, R"doc(
Center the remainder term in place: the midpoint of the remainder bound
is added to the constant polynomial coefficient and subtracted from the
remainder. Returns this object.
)doc")
      .def("C", &TV::C, R"doc(
Center the remainder term in place; shortcut for ``center``.
)doc")
      .def("constant", &TV::constant, py::arg("reset") = false, R"doc(
Return the coefficient of the constant term (float).

Parameters
----------
reset : bool, optional
    If True, also reset the constant term to zero in this variable.
    Default is False.
)doc")
      .def(
          "linear", [](TV& self, unsigned const ivar, bool const reset)
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
      .def(
          "linear",
          [](TV const& self)
          {
            double* p = self.linear();
            std::vector<double> res(p, p + self.nvar());
            delete[] p;
            return res;
          },
          R"doc(
Return the coefficients of all linear terms.

Returns
-------
lin : list of float
    New list of length ``nvar``; ``lin[k]`` is the coefficient of the
    linear term in variable ``k``, expressed in the original
    (unscaled) variables.
)doc")
      .def("__str__",
           [](TV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](TV const& self)
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
      .def("__abs__", [](TV const& m) { return mc::Op<TV>::abs(m); })
      .def("__pow__", [](TV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](TV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__", [](TV const& m, TV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, TV const& m) { return mc::pow(r, m); });

  m.def(
      "inv", [](TV const& x) { return mc::inv(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of 1/x.");
  m.def(
      "sqr", [](TV const& x) { return mc::sqr(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of x**2.");
  m.def(
      "sqrt", [](TV const& x) { return mc::sqrt(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of sqrt(x).");
  m.def(
      "exp", [](TV const& x) { return mc::exp(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of exp(x).");
  m.def(
      "log", [](TV const& x) { return mc::log(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of log(x).");
  m.def(
      "xlog", [](TV const& x) { return mc::xlog(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of x*log(x).");
  m.def(
      "cos", [](TV const& x) { return mc::cos(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of cos(x).");
  m.def(
      "sin", [](TV const& x) { return mc::sin(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of sin(x).");
  m.def(
      "tan", [](TV const& x) { return mc::tan(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of tan(x).");
  m.def(
      "acos", [](TV const& x) { return mc::acos(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of acos(x).");
  m.def(
      "asin", [](TV const& x) { return mc::asin(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of asin(x).");
  m.def(
      "atan", [](TV const& x) { return mc::atan(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of atan(x).");
  m.def(
      "cosh", [](TV const& x) { return mc::cosh(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of cosh(x).");
  m.def(
      "sinh", [](TV const& x) { return mc::sinh(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of sinh(x).");
  m.def(
      "tanh", [](TV const& x) { return mc::tanh(x); }, py::arg("x"),
      "Taylor model overload: Taylor model of tanh(x).");
  m.def(
      "pow", [](TV const& x, int const n) { return mc::pow(x, n); },
      py::arg("x"), py::arg("n"),
      "Taylor model overload: Taylor model of x**n for integer n.");
  m.def(
      "pow", [](TV const& x, double const& r) { return mc::pow(x, r); },
      py::arg("x"), py::arg("r"),
      "Taylor model overload: Taylor model of x**r for real r.");
  m.def(
      "pow", [](TV const& x, TV const& y) { return mc::pow(x, y); },
      py::arg("x"), py::arg("y"),
      "Taylor model overload: Taylor model of x**y.");
  m.def(
      "pow", [](double const& r, TV const& y) { return mc::pow(r, y); },
      py::arg("r"), py::arg("y"),
      "Taylor model overload: Taylor model of r**y for real r.");
  m.def(
      "cheb", [](TV const& x, unsigned const n) { return mc::cheb(x, n); },
      py::arg("x"), py::arg("n"),
      "Taylor model overload: Taylor model of the Chebyshev polynomial "
      "T_n(x).");
  m.def(
      "hull", [](TV const& x, TV const& y) { return mc::hull(x, y); },
      py::arg("x"), py::arg("y"),
      "Taylor model overload: Taylor model enclosing the union of the "
      "enclosures x and y (see option REF_POLY for the polynomial part).");
  m.def(
      "inter",
      [](TV& xy, TV const& x, TV const& y) { return mc::inter(xy, x, y); },
      py::arg("xy"), py::arg("x"), py::arg("y"),
      "Taylor model overload: overwrite xy with a Taylor model of the "
      "intersection of the enclosures x and y; returns False if the "
      "intersection is empty.");


  pyTModelOptions.def(py::init<>(), "Construct options with default values.")
      .def(py::init<TM::Options const&>(), "Copy constructor.")
      .def("reset", &TM::Options::reset, R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("BOUNDER_TYPE", &TM::Options::BOUNDER_TYPE, R"doc(
Taylor model range bounder, as a ``TModel.Options.BOUNDER`` value:
NAIVE, LSB, EIGEN, BERNSTEIN or HYBRID. Default is LSB.
)doc")
      .def_readwrite("BERNSTEIN_ORDER", &TM::Options::BERNSTEIN_ORDER, R"doc(
Order of the Bernstein polynomial for range bounding when
``BOUNDER_TYPE`` is BERNSTEIN; must be no less than the model order.
Default is 0 (same as the model order).
)doc")
      .def_readwrite("SCALE_VARIABLES", &TM::Options::SCALE_VARIABLES, R"doc(
Whether to scale the variable ranges to [-1,1] internally (bool).
Default is False.
)doc")
      .def_readwrite("CENTER_REMAINDER", &TM::Options::CENTER_REMAINDER, R"doc(
Whether to center the remainder term during Taylor model propagation
(bool). Default is False.
)doc")
      .def_readwrite("REF_MIDPOINT", &TM::Options::REF_MIDPOINT, R"doc(
Whether to take the midpoint of the inner range as the reference point
in the composition with an outer univariate function (True), as opposed
to the constant coefficient of the inner Taylor model (False). Default
is False.
)doc")
      .def_readwrite("REF_POLY", &TM::Options::REF_POLY, R"doc(
Scalar in [0,1] selecting the polynomial part in the functions
``inter`` and ``hull``: 0 selects the polynomial part of the left
operand, 1 that of the right operand, intermediate values interpolate.
Default is 0.
)doc")
      .def_readwrite("BERNSTEIN_USE", &TM::Options::BERNSTEIN_USE, R"doc(
Whether to compute a Bernstein model of the outer function in a
univariate composition and use it instead of the Taylor model when its
remainder is smaller (bool). Tighter on wide ranges, at extra cost.
Default is False.
)doc")
      .def_readwrite("BERNSTEIN_OPT", &TM::Options::BERNSTEIN_OPT, R"doc(
Whether to compute exact remainder bounds for Bernstein models of the
convex/concave univariates exp, log, inv and sqrt (bool), using golden
section search. Default is True.
)doc")
      .def_readwrite("BERNSTEIN_MAXIT", &TM::Options::BERNSTEIN_MAXIT, R"doc(
Maximum number of iterations for the determination of the exact
remainder bounds in a Bernstein model (int). Default is 100.
)doc")
      .def_readwrite("BERNSTEIN_TOL", &TM::Options::BERNSTEIN_TOL, R"doc(
Termination tolerance for the determination of the exact remainder
bounds in a Bernstein model (float). Default is 1e-10.
)doc")
      .def_readwrite("DISPLAY_DIGITS", &TM::Options::DISPLAY_DIGITS, R"doc(
Number of digits used when printing Taylor model coefficients (int).
Default is 7.
)doc");


  // Nested class Exceptions
  py::class_<TM::Exceptions> pyTModelExceptions(pyTModel, "Exceptions", R"doc(
Exception data thrown by Taylor model arithmetic.

Carries an error code (``ierr``) and a description (``what``); see the
``TModel.Exceptions.TYPE`` enumeration for the possible errors.
)doc");

  py::enum_<TM::Exceptions::TYPE>(pyTModelExceptions, "TYPE")
      .value("DIV", TM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", TM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", TM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", TM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("ASIN", TM::Exceptions::TYPE::ASIN,
             "Sine/Cosine inverse operation with range outside [-1,1]")
      .value("EIGEN", TM::Exceptions::TYPE::EIGEN,
             "Failed to compute eigenvalue decomposition in range bounder "
             "TModel::Options::EIGEN")
      .value("BERNSTEIN", TM::Exceptions::TYPE::BERNSTEIN,
             "Failed to compute the maximum gap between a univariate term and "
             "its Bernstein model")
      .value("INIT", TM::Exceptions::TYPE::INIT,
             "Failed to construct Taylor variable")
      .value("INCON", TM::Exceptions::TYPE::INCON,
             "Taylor model bound does not intersect with bound in template "
             "parameter arithmetic")
      .value("TMODEL", TM::Exceptions::TYPE::TMODEL,
             "Operation between Taylor variables linked to different Taylor "
             "models")
      .value("UNDEF", TM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented in mc::TModel")
      .export_values();

  pyTModelExceptions.def(py::init<TM::Exceptions::TYPE>())
      .def("ierr", &TM::Exceptions::ierr, "Return the error code (int).")
      .def("what", &TM::Exceptions::what, "Return the error description (str).");
}
