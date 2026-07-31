#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ffunc.hpp"
#include "smon.hpp"
#include "spoly.hpp"

namespace py = pybind11;

void
mc_ffmon(py::module_& m)
{
  typedef mc::FFVar const* KEY;
  typedef mc::lt_FFVar COMP;
  typedef mc::SMon<KEY, COMP> FFMon;
  typedef mc::SPoly<KEY, COMP> FFPoly;

  py::class_<FFMon> pyFFMon(m, "FFMon", R"doc(
Sparse multivariate monomial over DAG variables (FFVar).

FFMon is the same sparse monomial representation as SMon, except that the
participating variables are FFVar nodes of a DAG (FFGraph) instead of
plain integer indices: the attribute ``expr`` holds (FFVar, exponent)
pairs and ``tord`` the total order (sum of all exponents). The empty
monomial represents the constant 1.

Depending on the basis chosen when displaying or evaluating (see
FFPoly.Options.BASIS), a monomial with exponents k_i is interpreted either
as the power product ``prod_i x_i**k_i`` (monomial basis) or as the
Chebyshev product ``prod_i T_{k_i}(x_i)`` (Chebyshev basis).

FFMon objects are primarily used as the keys of the coefficient map of
sparse polynomials over DAG variables (FFPoly); they are hashable and
comparable for equality. The overloaded operators act on exponents:
``mon1 + mon2`` merges two monomials by adding exponents (the product of
the power products), ``mon1 - mon2`` subtracts exponents (requires mon2 to
be a subset), ``mon * k`` and ``mon / k`` multiply and divide every
exponent by the integer k.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(2, "x")
>>> mon = pymcpp.FFMon({X[0]: 1, X[1]: 2})
>>> mon.tord
3
>>> print(mon)
[x0]·[x1]^2
)doc");

  // --- Exceptions Nested Class ---
  py::class_<FFMon::Exceptions> pyFFMonExceptions(pyFFMon, "Exceptions",
                                                  R"doc(
Exception information for errors raised during FFMon operations.

Instances carry an error code (ierr) and a human-readable description
(what).
)doc");
  py::enum_<FFMon::Exceptions::TYPE>(pyFFMonExceptions, "TYPE",
                                     "Error codes for FFMon exceptions")
      .value("SUB", FFMon::Exceptions::SUB,
             "Subtraction of a monomial that is not a proper subset")
      .value("DIV", FFMon::Exceptions::DIV,
             "Division by a factor greater than the greatest common exponent")
      .value("CONV", FFMon::Exceptions::CONV,
             "Conversion with different variable indexing failed")
      .export_values();

  pyFFMonExceptions
      .def("ierr", &FFMon::Exceptions::ierr, R"doc(
Return the error code as an integer (see FFMon.Exceptions.TYPE).
)doc")
      .def("what", &FFMon::Exceptions::what, R"doc(
Return a string describing the error.
)doc");

  // --- FFMon Main Class ---
  pyFFMon
      // Constructors
      .def(py::init<>(), R"doc(
Construct the constant monomial 1 (no variables, total order 0).
)doc")
      .def(py::init<KEY const&, unsigned const>(), py::arg("var"),
           py::arg("ord") = 1, R"doc(
Construct the univariate monomial var**ord.

Parameters
----------
var : FFVar
    Participating DAG variable.
ord : int, optional
    Exponent of the variable. Default is 1. If 0, the constant monomial 1
    is constructed.
)doc")
      .def(py::init<std::map<KEY, unsigned, COMP> const&>(), py::arg("expr"),
           R"doc(
Construct a monomial from a dictionary of variables and exponents.

Parameters
----------
expr : dict of FFVar to int
    Dictionary mapping each participating DAG variable to its exponent.
    Entries with zero exponent are discarded.
)doc")
      .def(py::init<FFMon const&>(), R"doc(
Construct a copy of another monomial.
)doc")
      // Accessors
      .def_readonly("tord", &FFMon::tord, R"doc(
Total order of the monomial: the sum of all variable exponents (int).
)doc")
      .def_property_readonly(
          "expr",
          //&FFMon::expr,
          [](FFMon& self)
          {
            auto const& expr = self.expr;
            return std::vector<std::pair<KEY, unsigned> >(expr.cbegin(),
                                                          expr.cend());
          },
          R"doc(
Participating variables and exponents, as a new list of tuples.

Each element is a tuple ``(var, ord)`` of the DAG variable (FFVar) and its
positive exponent. Empty for the constant monomial 1.
)doc")
      // Methods
      .def("exp", &FFMon::exp, py::arg("var"), R"doc(
Return the exponent of variable `var` in the monomial.

Parameters
----------
var : FFVar
    Queried DAG variable.

Returns
-------
ord : int
    Exponent of `var`; 0 if the variable does not participate.
)doc")
      .def("gcexp", &FFMon::gcexp, R"doc(
Return the greatest common divisor of all variable exponents.

Returns
-------
gce : int
    Greatest common exponent; equals `tord` for monomials with at most one
    participating variable.
)doc")
      .def("lexp", &FFMon::lexp, R"doc(
Return the least exponent among all participating variables.

Returns
-------
le : int
    Smallest exponent; equals `tord` for monomials with at most one
    participating variable.
)doc")
      .def("gexp", &FFMon::gexp, R"doc(
Return the greatest exponent among all participating variables.

Returns
-------
ge : int
    Largest exponent; equals `tord` for monomials with at most one
    participating variable.
)doc")
      .def("hull", &FFMon::hull, py::arg("mon"), R"doc(
Merge with another monomial by taking the largest exponent per variable.

The monomial is modified in place: for every variable in `mon`, its
exponent is raised to that of `mon` if larger, and `tord` is updated
accordingly. The argument `mon` is not modified.

Parameters
----------
mon : FFMon
    Monomial to take the union (elementwise maximum of exponents) with.
)doc")
      .def("inter", &FFMon::inter, py::arg("mon"), R"doc(
Test whether the two monomials share at least one variable.

Parameters
----------
mon : FFMon
    Monomial to test against.

Returns
-------
res : bool
    True if any variable participates in both monomials, regardless of
    exponents.
)doc")
      .def("subset", &FFMon::subset, py::arg("mon"), R"doc(
Test whether this monomial is a proper subset (strict divisor) of `mon`.

Parameters
----------
mon : FFMon
    Monomial to test against.

Returns
-------
res : bool
    True if every variable exponent in this monomial is no larger than in
    `mon` and the total order is strictly smaller.
)doc")
      .def("subseteq", &FFMon::subseteq, py::arg("mon"), R"doc(
Test whether this monomial is a subset (divisor) of `mon`.

Parameters
----------
mon : FFMon
    Monomial to test against.

Returns
-------
res : bool
    True if every variable exponent in this monomial is no larger than in
    `mon`; equal monomials also qualify.
)doc")
      .def("display", &FFMon::display, py::arg("basis") = 0, R"doc(
Return a string representation of the monomial in the given basis.

Parameters
----------
basis : int, optional
    0 for the monomial basis (e.g. ``[x0]·[x1]^2``), 1 for the Chebyshev
    basis (e.g. ``T1[x0]·T2[x1]``). Default is 0.

Returns
-------
s : str
    String representation; "1" for the constant monomial.
)doc")
      // Operators
      .def(
          "__hash__",
          [](FFMon const& self)
          {
            // Hash is computed based on 'tord' and the 'expr' map content.
            // We convert the map to a tuple of items (which is sorted by key in
            // C++ map) to ensure consistent hashing in Python.
            py::list expr_items;
            for (auto const& [var, ord] : self.expr)
            {
              expr_items.append(py::make_tuple(var, ord));
            }
            // Return hash((tord, tuple(items)))
            return py::hash(py::make_tuple(self.tord, py::tuple(expr_items)));
          },
          R"doc(
Return a hash based on the total order and the variable-exponent pairs, so
that FFMon objects can be used as dictionary keys.
)doc")
      .def(
          "__getitem__", [](FFMon const& self, KEY const& var)
          { return self[var]; }, R"doc(
Extract the univariate sub-monomial var**exponent for variable `var`;
returns the constant monomial 1 if `var` does not participate.
)doc")
      .def(py::self += py::self)
      .def(py::self -= py::self)
      .def(py::self + py::self)
      .def(py::self - py::self)
      .def(py::self == py::self)
      .def(py::self != py::self)
      .def(py::self *= unsigned())
      .def(py::self /= unsigned())
      .def(py::self * unsigned())
      .def(py::self / unsigned())

      // Comparison operators
      .def(
          "__eq__", [](FFMon const& Mon1, FFMon const& Mon2)
          { return Mon1 == Mon2; }, py::is_operator())
      .def(
          "__ne__", [](FFMon const& Mon1, FFMon const& Mon2)
          { return !(Mon1 == Mon2); }, py::is_operator())

      // String representation
      .def("__str__", [](FFMon const& self)
           { return self.display(FFPoly::options.BASIS); })
      .def("__repr__", [](FFMon const& self)
           { return self.display(FFPoly::options.BASIS); });
}
