#include "spoly.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ffunc.hpp"

namespace py = pybind11;

void
mc_spoly(py::module_& m)
{
  typedef unsigned KEY;
  typedef std::less<unsigned> COMP;
  typedef mc::SPoly<KEY, COMP> SP;
  typedef mc::SMon<KEY, COMP> SM;

  py::class_<SP> pySPoly(m, "SPoly", R"doc(
Sparse multivariate polynomial with integer-indexed variables.

An SPoly stores a polynomial as a sparse coefficient map from monomials
(SMon objects) to float coefficients, ordered in graded lexicographic
order. Variables are identified by nonnegative integer indices. The
coefficient map is exposed through the ``coefmon`` property, and the set
of participating variable indices through ``setvar``.

The monomial keys are interpreted either in the monomial (power) basis or
in the Chebyshev basis, as selected by the static option
``SPoly.options.BASIS``; this choice affects arithmetic, display,
differentiation and evaluation. Use convert() to change the basis of a
given polynomial.

Polynomials are built from constants, variables, monomials or coefficient
dictionaries via the constructors, and combined with the overloaded
arithmetic operators +, -, *, / (division only by a scalar or by a
polynomial that divides exactly, otherwise SPoly.Exceptions is raised) and
** with a nonnegative integer exponent.

Examples
--------
>>> import pymcpp
>>> x = pymcpp.SPoly(0)   # variable with index 0
>>> y = pymcpp.SPoly(1)   # variable with index 1
>>> p = (x + y)**2
>>> print(p.nmon)
3
>>> p.eval({0: 1.0, 1: 2.0})
9.0
)doc");

  // --- Exceptions Nested Class ---
  py::class_<SP::Exceptions> pySPExceptions(pySPoly, "Exceptions", R"doc(
Exception information for errors raised during SPoly arithmetic.

Instances carry an error code (ierr) and a human-readable description
(what).
)doc");
  py::enum_<SP::Exceptions::TYPE>(pySPExceptions, "TYPE",
                                  "Error codes for SPoly exceptions")
      .value("DIVZERO", SP::Exceptions::DIVZERO, "Scalar division by zero")
      .value("DIVPOLY", SP::Exceptions::DIVPOLY,
             "Division between two polynomials that is not exact")
      .value("INTERNAL", SP::Exceptions::INTERNAL, "Internal error")
      .export_values();

  pySPExceptions
      .def("ierr", &SP::Exceptions::ierr, R"doc(
Return the error code as an integer (see SPoly.Exceptions.TYPE).
)doc")
      .def("what", &SP::Exceptions::what, R"doc(
Return a string describing the error.
)doc");

  // --- Options Nested Struct ---
  py::class_<SP::Options> pySPOptions(pySPoly, "Options", R"doc(
Options controlling sparse polynomial arithmetic and display.

Options take effect through the static member ``SPoly.options``, shared by
all SPoly objects.
)doc");

  py::enum_<SP::Options::BASIS_TYPE>(pySPOptions, "BASIS_TYPE",
                                     "Available polynomial bases")
      .value("MONOM", SP::Options::MONOM,
             "Monomial (power) basis: monomial keys stand for products "
             "x_i**k_i")
      .value("CHEB", SP::Options::CHEB,
             "Chebyshev basis: monomial keys stand for products T_{k_i}(x_i) "
             "of Chebyshev polynomials")
      .export_values();

  pySPOptions
      .def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def_readwrite("BASIS", &SP::Options::BASIS, R"doc(
Basis in which the monomial keys of the coefficient maps are interpreted:
0 for the monomial (power) basis, 1 for the Chebyshev basis (see
SPoly.Options.BASIS_TYPE). Default is 0.
)doc")
      .def_readwrite("REMZERO", &SP::Options::REMZERO, R"doc(
Whether to remove zero coefficient entries from the coefficient map during
arithmetic operations (bool). Default is True.
)doc")
      .def_readwrite("DISPLEN", &SP::Options::DISPLEN, R"doc(
Number of significant digits used for coefficients in string
representations (int). Default is 5.
)doc")
      .def_readwrite("DISPLINE", &SP::Options::DISPLINE, R"doc(
Whether string representations put the whole polynomial on a single line
(True) or one monomial term per line (False). Default is True.
)doc");

  // --- SPoly Main Class ---
  pySPoly
      // Static Options
      .def_readwrite_static("options", &SP::options, R"doc(
Static options shared by all SPoly objects (see SPoly.Options).
)doc")
      // Constructors
      .def(py::init<>(), R"doc(
Construct the zero polynomial (empty coefficient map).
)doc")
      .def(py::init<double const>(), py::arg("coef"), R"doc(
Construct a constant polynomial.

Parameters
----------
coef : float
    Value of the constant.
)doc")
      .def(py::init<KEY const&, double const&>(), py::arg("var"),
           py::arg("coef") = 1., R"doc(
Construct the linear polynomial coef * var.

Parameters
----------
var : int
    Index of the participating variable.
coef : float, optional
    Coefficient multiplying the variable. Default is 1.0.
)doc")
      .def(
          // py::init<std::pair<SM, double> const&>(),
          py::init([](SM const& mon, double const& coef)
                   { return new SP({mon, coef}); }),
          py::arg("mon"), py::arg("coef") = 1., R"doc(
Construct the single-term polynomial coef * mon.

Parameters
----------
mon : SMon
    Monomial making up the single term.
coef : float, optional
    Coefficient multiplying the monomial. Default is 1.0.
)doc")
      .def(py::init<SP::t_poly const&>(),
           // py::init( []( std::vector< std::pair< SM, double > > const&
           // coefmon )
           //           { return new SP( SP::t_poly( coefmon.cbegin(),
           //           coefmon.cend() ) ); } ),
           py::arg("coefmon"), R"doc(
Construct a polynomial from a coefficient dictionary.

Parameters
----------
coefmon : dict of SMon to float
    Dictionary mapping each monomial to its coefficient, e.g.
    ``{mon1: coef1, mon2: coef2}``.
)doc")
      .def(py::init<SP const&>(), R"doc(
Construct a copy of another sparse polynomial.
)doc")
      // Accessors
      .def("maxord", &SP::maxord, R"doc(
Return the maximal total order among all monomial terms (0 for the zero or
constant polynomial).
)doc")
      .def("minord", py::overload_cast<>(&SP::minord, py::const_), R"doc(
Return the minimal total order among all monomial terms (0 if a constant
term is present or the polynomial is zero).
)doc")
      .def("minord", py::overload_cast<KEY const&>(&SP::minord, py::const_),
           py::arg("x"), R"doc(
Return the minimal exponent of variable `x` over all monomial terms.

Parameters
----------
x : int
    Index of the queried variable.

Returns
-------
ord : int
    Smallest exponent of `x` among the monomial terms; 0 if any term is
    free of `x` or the polynomial has a constant term.
)doc")
      .def_property_readonly("nmon", &SP::nmon, R"doc(
Number of monomial terms in the polynomial (int).
)doc")
      .def("coef", &SP::coef, py::arg("mon") = SM(), R"doc(
Return the coefficient of a given monomial.

Parameters
----------
mon : SMon, optional
    Queried monomial. Default is the constant monomial, so ``coef()``
    returns the constant term.

Returns
-------
coef : float
    Coefficient of `mon` in the polynomial; 0.0 if the monomial is not
    present.
)doc")
      .def("nvar", &SP::nvar, R"doc(
Return the number of participating variables (int).
)doc")
      // Properties
      .def_property_readonly(
          "coefmon", py::overload_cast<>(&SP::mapmon, py::const_), R"doc(
Coefficient map of the polynomial, as a dictionary.

Returns a dict mapping each monomial (SMon) to its float coefficient. The
monomial keys are interpreted in the basis given by
``SPoly.options.BASIS``.
)doc")
      .def_property_readonly("setvar",
                             py::overload_cast<>(&SP::setvar, py::const_),
                             R"doc(
Set of the indices (int) of all participating variables.
)doc")
      // Methods
      .def("swap", &SP::swap, py::arg("spoly"), R"doc(
Exchange the contents of this polynomial with those of `spoly`, in place.

Parameters
----------
spoly : SPoly
    Polynomial to swap contents with; it is modified as well.
)doc")
      .def("clean", &SP::clean, py::arg("tol") = 0., R"doc(
Remove all monomial terms whose coefficient magnitude is at most `tol`,
in place.

Parameters
----------
tol : float, optional
    Threshold below which (in absolute value) coefficients are discarded.
    Default is 0.0, removing exact zeros only.
)doc")
      .def("convert", &SP::convert, py::arg("basis"), R"doc(
Convert the coefficient map to another basis, in place.

The current coefficients are interpreted in the basis given by
``SPoly.options.BASIS`` and converted to `basis`; no conversion is
performed if the two coincide. Note that ``SPoly.options.BASIS`` itself is
not modified.

Parameters
----------
basis : int
    Target basis: 0 for the monomial (power) basis, 1 for the Chebyshev
    basis.

Returns
-------
spoly : SPoly
    The converted polynomial.
)doc")
      .def("ismultiple", &SP::ismultiple, py::arg("spoly"), R"doc(
Test whether this polynomial is a scalar multiple of `spoly`.

Parameters
----------
spoly : SPoly
    Polynomial to compare against.

Returns
-------
res : bool
    True if the two polynomials have identical monomials and
    proportional coefficients.
mult : float
    Proportionality factor such that ``self == mult * spoly``; 0.0 when
    `res` is False.
)doc")
      .def("factor", &SP::factor, py::arg("x"), R"doc(
Collect the polynomial by powers of variable `x`.

Returns the decomposition ``p = sum_k b_k(x) * q_k`` where ``b_k`` is the
k-th basis function in `x` (``x**k`` in the monomial basis, ``T_k(x)`` in
the Chebyshev basis) and the coefficient polynomials ``q_k`` do not
involve `x`.

Parameters
----------
x : int
    Index of the variable to factor with respect to.

Returns
-------
factors : dict of int to SPoly
    Dictionary mapping each participating order k of `x` to the
    corresponding coefficient polynomial ``q_k``. Empty if `x` does not
    participate in the polynomial.
)doc")
      .def("diff", &SP::diff, py::arg("x"), R"doc(
Differentiate the polynomial with respect to variable `x`.

The differentiation rules follow the current basis
``SPoly.options.BASIS``.

Parameters
----------
x : int
    Index of the variable to differentiate with respect to.

Returns
-------
der : SPoly
    New polynomial holding the partial derivative; the zero polynomial if
    `x` does not participate.
)doc")
      .def("var", &SP::var, py::arg("x"), R"doc(
Reset the polynomial to equal the variable `x`, in place.

Parameters
----------
x : int
    Index of the variable.

Returns
-------
spoly : SPoly
    The polynomial itself after being set to `x`.
)doc")
      .def_static(
          "display",
          [](std::vector<std::pair<SM, double>> const& coefmon, int basis,
             int len, bool line)
          {
            SP::t_poly internal_map;
            for (auto const& [k, v] : coefmon) internal_map.insert({k, v});
            return SP::display(internal_map, basis, len, line);
          },
          py::arg("coefmon"), py::arg("BASIS") = SP::options.BASIS,
          py::arg("DISPLEN") = SP::options.DISPLEN,
          py::arg("ONELINE") = SP::options.DISPLINE, R"doc(
Format a coefficient map as a polynomial expression string.

Parameters
----------
coefmon : list of tuple of (SMon, float)
    Coefficient map given as (monomial, coefficient) pairs, e.g. the items
    of the ``coefmon`` property of a polynomial.
BASIS : int, optional
    Basis used to render the monomials: 0 monomial, 1 Chebyshev. Default
    is the current ``SPoly.options.BASIS``.
DISPLEN : int, optional
    Number of significant digits for the coefficients. Default is the
    current ``SPoly.options.DISPLEN``.
ONELINE : bool, optional
    Whether to place all terms on a single line rather than one term per
    line. Default is the current ``SPoly.options.DISPLINE``.

Returns
-------
s : str
    Formatted polynomial expression.
)doc")
      .def(
          "eval", [](SP const& self, std::map<KEY, double, COMP> const& val)
          { return self.eval(val); }, py::arg("x"), R"doc(
Evaluate the polynomial at a point using float arithmetic.

The monomials are interpreted in the basis given by
``SPoly.options.BASIS``.

Parameters
----------
x : dict of int to float
    Value of each variable, keyed by variable index. Every participating
    variable must be present, otherwise IndexError is raised.

Returns
-------
val : float
    Value of the polynomial at the given point.
)doc")
      // Operators
      .def(py::self += py::self)
      .def(py::self += double())
      .def(py::self += std::pair<SM, double>())
      //.def( "add_mon", []( SP& self, SM const& m, double c ){ self +=
      // std::make_pair(m, c); return self; } )

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

      .def(
          "__eq__", [](SP const& x, SP const& y) { return x == y; },
          py::is_operator())
      .def(
          "__ne__", [](SP const& x, SP const& y) { return !(x == y); },
          py::is_operator())
      .def(
          "__pow__", [](SP const& s, unsigned n) { return mc::pow(s, n); },
          py::is_operator())

      // String Representation
      .def("__str__",
           [](SP const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           })
      .def("__repr__",
           [](SP const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           });

  m.def(
      "sqr", [](SP const& s) { return mc::sqr(s); }, R"doc(
SPoly overload: square of a sparse polynomial, computed in the basis given
by SPoly.options.BASIS.
)doc");
  m.def(
      "pow", [](SP const& s, unsigned n) { return mc::pow(s, n); }, R"doc(
SPoly overload: n-th power of a sparse polynomial for a nonnegative
integer n, computed in the basis given by SPoly.options.BASIS.
)doc");
  m.def(
      "cheb", [](SP const& s, unsigned n) { return mc::cheb(s, n); }, R"doc(
SPoly overload: Chebyshev polynomial of the first kind T_n applied to a
sparse polynomial, using the recurrence T_n = 2*s*T_{n-1} - T_{n-2}.
)doc");
  m.def(
      "prod",
      [](unsigned n, std::vector<SP> const& v)
      { return mc::prod(n, v.data()); },
      R"doc(
SPoly overload: product v[0]*...*v[n-1] of sparse polynomials.
)doc");
}
