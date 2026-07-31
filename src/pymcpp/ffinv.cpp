// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include "ffinv.hpp"

#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include <sstream>

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFInv::Options::NLINV>);

void
mc_ffinv(py::module_& m)
{
  typedef mc::FFInv FI;

  py::class_<FI> pyFFInv(m, "FFInv", R"doc(
Invertible structure detection in factorable functions.

An FFInv object records, for a scalar expression, which independent
variables participate in the expression and whether the expression could
be inverted with respect to each of them - that is, whether an equation
``f(x) = c`` could be solved symbolically for that variable. The possible
outcomes are:

- ``FFInv.TYPE.L`` (linear): the variable participates linearly and is
  invertible.
- ``FFInv.TYPE.S`` (separably linear): the variable participates in a
  single linear or multilinear term.
- ``FFInv.TYPE.N`` (separably nonlinear): the variable participates in a
  single nonlinear term whose operations are all invertible, as specified
  by the option set ``FFInv.options.INVOP``.
- ``FFInv.TYPE.U`` (undetermined): invertibility could not be established,
  for instance because the variable occurs multiple times or participates
  in a non-invertible operation.

To analyze a function, create one FFInv object per variable with
``FFInv().indep(i)``, then evaluate the function on these objects - either
by calling a Python function directly, or by passing them to
``FFGraph.eval`` to propagate through a DAG. The standard arithmetic
operators (+, -, *, /) and the module-level functions (exp, log, sqrt,
pow, ...) are overloaded for FFInv. Query the result with inv().

Which nonlinear operations count as invertible is controlled by the static
member ``FFInv.options``; it applies to all FFInv computations. This class
is typically used alongside FFDep, which detects the type of dependence
instead of invertibility.

Examples
--------
>>> import pymcpp
>>> X = [pymcpp.FFInv().indep(i) for i in range(4)]
>>> F0 = X[2] * X[3] + X[0] / X[2]
>>> print(F0)
{ 0S 2U 3S }
>>> F0.inv(2)
(True, <TYPE.U: 3>)
)doc");
  py::class_<FI::Exceptions> pyFIExceptions(pyFFInv, "Exceptions", R"doc(
Exception information for errors raised during FFInv arithmetic.

Instances carry an error code (ierr) and a human-readable description
(what).
)doc");
  py::class_<FI::Options> pyFIOptions(pyFFInv, "Options", R"doc(
Options controlling invertible structure detection with FFInv.

The only option is INVOP, the set of nonlinear operations that are
considered invertible. Options take effect through the static member
``FFInv.options``, shared by all FFInv objects.
)doc");

  py::enum_<FI::Exceptions::TYPE>(pyFIExceptions, "TYPE",
                                  "Error codes for FFInv exceptions")
      .value("UNDEF", FI::Exceptions::UNDEF,
             "Call to a feature that is unavailable in FFInv arithmetic")
      .value("INTERN", FI::Exceptions::INTERN, "Internal error")
      .export_values();

  pyFIExceptions
      .def("ierr", &FI::Exceptions::ierr, R"doc(
Return the error code as an integer (see FFInv.Exceptions.TYPE).
)doc")
      .def("what", &FI::Exceptions::what, R"doc(
Return a string describing the error.
)doc");

  py::enum_<FI::Options::NLINV>(pyFIOptions, "NLINV", R"doc(
Nonlinear operations that may be declared invertible.

Members of this enumeration are inserted into the option set
``FFInv.options.INVOP`` to declare the corresponding operation invertible.
An expression applying an operation from INVOP to an invertible
subexpression is tagged separably nonlinear 'N'; applying an operation not
in INVOP yields undetermined 'U'.
)doc")
      .value("INV", FI::Options::INV, "Reciprocal, 1/x")
      .value("SQRT", FI::Options::SQRT, "Square root, sqrt(x)")
      .value("EXP", FI::Options::EXP, "Exponential, exp(x)")
      .value("LOG", FI::Options::LOG, "Natural logarithm, log(x)")
      .value("IPOW", FI::Options::IPOW, "Integer power, x**n")
      .value("RPOW", FI::Options::RPOW, "Real power, x**r")
      .export_values();

  py::class_<std::set<FI::Options::NLINV>>(m, "NLINVSet", R"doc(
Mutable set of FFInv.Options.NLINV values.

Wraps the underlying C++ set of invertible nonlinear operations used for
the option ``FFInv.options.INVOP``, so that modifications such as
``FFInv.options.INVOP.add(...)`` act directly on the options. Supports
add, remove, discard, membership testing with ``in``, ``len`` and
iteration, like a Python set.
)doc")
      .def(py::init<>(), R"doc(
Construct an empty set of invertible operations.
)doc")
      .def(
          "add",
          [](std::set<FI::Options::NLINV>& s, FI::Options::NLINV const& k)
          { s.insert(k); }, py::arg("op"), R"doc(
Add an operation to the set; no effect if already present.

Parameters
----------
op : FFInv.Options.NLINV
    Operation to declare invertible.
)doc")
      .def(
          "remove",
          [](std::set<FI::Options::NLINV>& s, FI::Options::NLINV const& k)
          {
            if (s.find(k) == s.end()) throw py::key_error();
            s.erase(k);
          },
          py::arg("op"), R"doc(
Remove an operation from the set.

Parameters
----------
op : FFInv.Options.NLINV
    Operation to remove.

Raises
------
KeyError
    If the operation is not in the set.
)doc")
      .def(
          "discard",
          [](std::set<FI::Options::NLINV>& s, FI::Options::NLINV const& k)
          { s.erase(k); }, py::arg("op"), R"doc(
Remove an operation from the set if present; no error otherwise.

Parameters
----------
op : FFInv.Options.NLINV
    Operation to remove.
)doc")
      .def("__contains__",
           [](std::set<FI::Options::NLINV>& s, FI::Options::NLINV const& k)
           { return s.find(k) != s.end(); })
      .def("__len__", [](std::set<FI::Options::NLINV>& s) { return s.size(); })
      .def(
          "__iter__", [](std::set<FI::Options::NLINV>& s)
          { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>())
      .def("__repr__",
           [](std::set<FI::Options::NLINV>& s)
           {
             std::ostringstream oss;
             oss << "{";
             bool first = true;
             for (auto const& v : s)
             {
               if (!first) oss << ", ";
               switch (v)
               {
                 case FI::Options::INV:
                   oss << "INV";
                   break;
                 case FI::Options::SQRT:
                   oss << "SQRT";
                   break;
                 case FI::Options::EXP:
                   oss << "EXP";
                   break;
                 case FI::Options::LOG:
                   oss << "LOG";
                   break;
                 case FI::Options::IPOW:
                   oss << "IPOW";
                   break;
                 case FI::Options::RPOW:
                   oss << "RPOW";
                   break;
               }
               // Casting to int for simple enum representation, or you can map
               // back to string
               // oss << v;//static_cast<int>(v);
               first = false;
             }
             oss << "}";
             return oss.str();
           });

  pyFIOptions
      .def(py::init<>(), R"doc(
Construct an option set with the default invertible operations.
)doc")
      .def(py::init<FI::Options const&>(), R"doc(
Construct a copy of another option set.
)doc")
      .def(
          "reset", [](FI::Options& self) { self = FI::Options(); }, R"doc(
Reset the options to their default values.

Restores INVOP to {INV, SQRT, EXP, LOG, RPOW}; IPOW is not included by
default.
)doc")
      //.def_readwrite( "INVOP", &FI::Options::INVOP, "Set of allowed invertible
      // operations" );
      .def_property(
          "INVOP", [](FI::Options& self) -> std::set<FI::Options::NLINV>&
          { return self.INVOP; },
          [](FI::Options& self, py::object const& val)
          {
            // 1. Optimized path: Assignment from another NLINVSet (C++ copy)
            if (py::isinstance<std::set<FI::Options::NLINV>>(val))
            {
              self.INVOP = val.cast<std::set<FI::Options::NLINV> const&>();
              return;
            }
            // 2. Flexible path: Assignment from Python set/list/tuple
            else if (py::isinstance<py::iterable>(val))
            {
              self.INVOP.clear();
              for (auto item : val)
                // Cast items to the Enum type and insert
                self.INVOP.insert(item.cast<FI::Options::NLINV>());
              return;
            }
            throw py::type_error(
                "Cannot assign to INVOP: Expected NLINVSet or iterable of "
                "NLINV enums");
          },
          py::return_value_policy::reference_internal,
          R"doc(
Set of nonlinear operations considered invertible.

Reading the attribute returns an NLINVSet that references the underlying
option, so in-place modifications like
``FFInv.options.INVOP.add(FFInv.options.IPOW)`` take effect immediately.
The attribute may be assigned from an NLINVSet or from any iterable of
FFInv.Options.NLINV values. Default is {INV, SQRT, EXP, LOG, RPOW}; IPOW
is not included by default.
)doc");

  py::enum_<FI::TYPE>(pyFFInv, "TYPE", R"doc(
Invertibility of an expression with respect to a participating variable.
)doc")
      .value("L", FI::TYPE::L,
             "Linear: variable participates linearly and is invertible")
      .value("S", FI::TYPE::S,
             "Separably linear: variable participates in a single linear or "
             "multilinear term")
      .value("N", FI::TYPE::N,
             "Separably nonlinear: variable participates in a single "
             "nonlinear term whose operations are all invertible (as "
             "specified in FFInv.options.INVOP)")
      .value("U", FI::TYPE::U,
             "Undetermined: invertibility could not be established, e.g. due "
             "to multiple occurrences or non-invertible operations")
      .export_values();

  // --- Main Class ---
  pyFFInv
      // Constructors
      .def(py::init<>(), R"doc(
Construct an FFInv object with an empty invertibility map (a constant).
)doc")
      .def(py::init<double const>(), py::arg("val"), R"doc(
Construct an FFInv object from a real constant.

The value itself is discarded; the resulting object has an empty
invertibility map, like the default constructor.

Parameters
----------
val : float
    Constant value (ignored for the purpose of invertibility analysis).
)doc")
      .def(py::init<FI const&>(), R"doc(
Construct a copy of another FFInv object.
)doc")

      // Modifiers
      .def(
          "indep", [](FI& self, int const ndx) { return self.indep(ndx); },
          py::arg("ndx"), R"doc(
Initialize as the independent variable with index `ndx`.

Any previously recorded information is discarded, and the object is set to
depend linearly ('L') on the single variable index `ndx`. The object is
modified in place.

Parameters
----------
ndx : int
    Index identifying the independent variable.

Returns
-------
inv : FFInv
    The object itself, allowing the idiom
    ``X = [pymcpp.FFInv().indep(i) for i in range(NX)]``.
)doc")

      // Accessors
      .def_readwrite_static("options", &FI::options, R"doc(
Static options shared by all FFInv computations (see FFInv.Options).
)doc")
      .def(
          "inv", [](FI const& self, int const ndx) { return self.inv(ndx); },
          py::arg("ndx"), R"doc(
Query the invertibility with respect to the variable with index `ndx`.

Parameters
----------
ndx : int
    Index of the queried variable.

Returns
-------
participates : bool
    True if the expression depends on variable `ndx`.
itype : FFInv.TYPE
    Invertibility type with respect to that variable. When `participates`
    is False, this is the placeholder value FFInv.TYPE.U and carries no
    meaning.
)doc")
      .def(
          "inv", [](FI const& self) { return self.inv(); }, R"doc(
Return the full invertibility map of the expression.

Returns
-------
inv : dict of int to FFInv.TYPE
    Dictionary mapping the index of every participating variable to the
    invertibility type with respect to that variable. Variables that do
    not participate are absent from the dictionary.
)doc")
      .def("update", &FI::update, py::arg("invmin"), R"doc(
Raise every recorded invertibility tag to at least the given type, in
place.

Entries whose current type is already more general than `invmin` are left
unchanged.

Parameters
----------
invmin : FFInv.TYPE
    Minimum invertibility type to enforce on all participating variables.

Returns
-------
inv : FFInv
    The object itself.
)doc")

      // Operators
      .def(py::self += py::self)
      .def(py::self += double())
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

      // String representation
      .def("__str__",
           [](FI const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           })
      .def("__repr__",
           [](FI const& self)
           {
             std::ostringstream oss;
             oss << self;
             return oss.str();
           });

  m.def(
      "sqr", [](FI const& x) { return mc::sqr(x); }, R"doc(
FFInv overload: invertible structure of x**2; invertible ('N') only if
IPOW is in FFInv.options.INVOP, otherwise undetermined ('U').
)doc");
  m.def(
      "sqrt", [](FI const& x) { return mc::sqrt(x); }, R"doc(
FFInv overload: invertible structure of sqrt(x); invertible ('N') only if
SQRT is in FFInv.options.INVOP, otherwise undetermined ('U').
)doc");
  m.def(
      "exp", [](FI const& x) { return mc::exp(x); }, R"doc(
FFInv overload: invertible structure of exp(x); invertible ('N') only if
EXP is in FFInv.options.INVOP, otherwise undetermined ('U').
)doc");
  m.def(
      "log", [](FI const& x) { return mc::log(x); }, R"doc(
FFInv overload: invertible structure of log(x); invertible ('N') only if
LOG is in FFInv.options.INVOP, otherwise undetermined ('U').
)doc");
  m.def(
      "xlog", [](FI const& x) { return mc::xlog(x); }, R"doc(
FFInv overload: invertible structure of x*log(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "cos", [](FI const& x) { return mc::cos(x); }, R"doc(
FFInv overload: invertible structure of cos(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "sin", [](FI const& x) { return mc::sin(x); }, R"doc(
FFInv overload: invertible structure of sin(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "tan", [](FI const& x) { return mc::tan(x); }, R"doc(
FFInv overload: invertible structure of tan(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "acos", [](FI const& x) { return mc::acos(x); }, R"doc(
FFInv overload: invertible structure of acos(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "asin", [](FI const& x) { return mc::asin(x); }, R"doc(
FFInv overload: invertible structure of asin(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "atan", [](FI const& x) { return mc::atan(x); }, R"doc(
FFInv overload: invertible structure of atan(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "cosh", [](FI const& x) { return mc::cosh(x); }, R"doc(
FFInv overload: invertible structure of cosh(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "sinh", [](FI const& x) { return mc::sinh(x); }, R"doc(
FFInv overload: invertible structure of sinh(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "tanh", [](FI const& x) { return mc::tanh(x); }, R"doc(
FFInv overload: invertible structure of tanh(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "fabs", [](FI const& x) { return mc::fabs(x); }, R"doc(
FFInv overload: invertible structure of abs(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "erf", [](FI const& x) { return mc::erf(x); }, R"doc(
FFInv overload: invertible structure of erf(x); all dependencies become
undetermined ('U').
)doc");
  m.def(
      "fstep", [](FI const& x) { return mc::fstep(x); }, R"doc(
FFInv overload: invertible structure of the forward step function; all
dependencies become undetermined ('U').
)doc");
  m.def(
      "bstep", [](FI const& x) { return mc::bstep(x); }, R"doc(
FFInv overload: invertible structure of the backward step function; all
dependencies become undetermined ('U').
)doc");
  m.def(
      "pow", [](FI const& x, int n) { return mc::pow(x, n); }, R"doc(
FFInv overload: invertible structure of x**n for integer n; the outcome
depends on whether IPOW (and INV for negative n) is in
FFInv.options.INVOP.
)doc");
  m.def(
      "pow", [](FI const& x, double r) { return mc::pow(x, r); }, R"doc(
FFInv overload: invertible structure of x**r for real r; invertible ('N')
only if RPOW is in FFInv.options.INVOP, otherwise undetermined ('U').
)doc");
  m.def(
      "pow", [](FI const& x, FI const& y) { return mc::pow(x, y); }, R"doc(
FFInv overload: invertible structure of x**y; the dependence sets of x and
y are merged and all dependencies become undetermined ('U').
)doc");
  m.def(
      "min", [](FI const& x, FI const& y) { return mc::min(x, y); }, R"doc(
FFInv overload: invertible structure of min(x, y); the dependence sets are
merged and all dependencies become undetermined ('U').
)doc");
  m.def(
      "max", [](FI const& x, FI const& y) { return mc::max(x, y); }, R"doc(
FFInv overload: invertible structure of max(x, y); the dependence sets are
merged and all dependencies become undetermined ('U').
)doc");
  m.def(
      "cheb", [](FI const& x, unsigned n) { return mc::cheb(x, n); }, R"doc(
FFInv overload: invertible structure of the Chebyshev polynomial T_n(x);
all dependencies become undetermined ('U') for n > 2.
)doc");
  m.def(
      "prod",
      [](unsigned n, std::vector<FI> const& x)
      { return mc::prod(n, x.data()); },
      R"doc(
FFInv overload: invertible structure of the product x[0]*...*x[n-1],
applying the FFInv multiplication rules term by term.
)doc");
  m.def(
      "monom",
      [](unsigned n, std::vector<FI> const& x, std::vector<unsigned> const& k)
      {
        if (x.size() != k.size())
          throw std::invalid_argument(
              "Size mismatch between variables and exponents");
        return mc::monom(n, x.data(), k.data());
      },
      R"doc(
FFInv overload: invertible structure of the monomial
x[0]**k[0]*...*x[n-1]**k[n-1], applying the FFInv power and multiplication
rules term by term.
)doc");
}
