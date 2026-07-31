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
    boost::numeric::interval_lib::rounded_transc_opp<double>>
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

#include "specbnd.hpp"
typedef mc::Specbnd<I> SB;

#include "tmodel.hpp"
typedef mc::TVar<I> T;

#include "cmodel.hpp"
typedef mc::CVar<I> C;

#include "scmodel.hpp"
typedef mc::SCVar<I> SC;

#include "pwlu.hpp"
#include "supmodel.hpp"
typedef mc::SupVar<mc::PWLU> PWLS;
#include "pwcu.hpp"
typedef mc::SupVar<mc::PWCU> PWCS;

#include "polimage.hpp"
typedef mc::PolVar<I> PV;

#include "ellimage.hpp"
typedef mc::EllVar<I> EV;

#include "ffdep.hpp"
typedef mc::FFDep FD;

#include "ffinv.hpp"
typedef mc::FFInv FI;

#include <chrono>
#include <fstream>

#include "ffexpr.hpp"
#include "ffunc.hpp"

namespace py = pybind11;

void
mc_ffunc(py::module_& m)
{
  py::class_<mc::FFNum> pyFFNum(m, "FFNum", R"doc(
Numeric field of a DAG constant, holding an integer or real value.

`FFNum` stores the value attached to constant nodes of a factorable-
function DAG (`FFVar` nodes of type `CINT` or `CREAL`). It is normally
obtained from `FFVar.num()` rather than constructed directly.
)doc");
  pyFFNum
      .def(py::init<int const>(), R"doc(
Construct an integer constant (also the default constructor).

Parameters
----------
i : int, optional
    Integer value. Default is 0.
)doc",
           py::arg("i") = 0)
      .def(py::init<double const&>(), R"doc(
Construct a real constant.

Parameters
----------
d : float
    Real value.
)doc",
           py::arg("d"))
      .def(py::init<mc::FFNum const&>(), R"doc(
Copy constructor.
)doc",
           py::arg("num"))
      .def_property_readonly("val", &mc::FFNum::val, R"doc(
Constant value held by the numeric field, as a float.
)doc")
      .def("__str__",
           [](mc::FFNum const& V)
           {
             std::ostringstream Vss;
             Vss << V.val();
             return Vss.str();
           })
      .def("__repr__",
           [](mc::FFNum const& V)
           {
             std::ostringstream Vss;
             Vss << V.val();
             return Vss.str();
           });

  py::class_<mc::FFVar> pyFFVar(m, "FFVar", R"doc(
Variable node in the DAG of a factorable function.

An `FFVar` is a handle to a node of a directed acyclic graph (DAG)
recorded in an `FFGraph` environment: an original variable, an auxiliary
variable (the result of an operation), or an integer/real constant.
Applying the overloaded arithmetic operators ``+``, ``-``, ``*``, ``/``,
``**`` or the module functions (`pymcpp.exp`, `pymcpp.sqrt`, ...) to
`FFVar` operands does not compute numerical values; it appends new
auxiliary nodes to the DAG and returns `FFVar` handles to them. Common
subexpressions are detected and reused automatically.

Variables are typically created attached to a DAG with
``pymcpp.FFVar(DAG, name)`` or via `FFGraph.add_var`. The resulting
dependent expressions can then be differentiated symbolically
(`FFGraph.fdiff`, `FFGraph.bdiff`, `FFGraph.tdiff`) and evaluated in a
variety of arithmetics (`FFGraph.eval`).

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = X[0] * pymcpp.exp(X[0] * X[1]) + X[1] ** 2  # adds nodes to DAG
>>> DAG.eval([F], X, [1.0, 2.0])  # evaluate at X0=1, X1=2
[11.38905609893065]
)doc");
  pyFFVar
      .def(py::init<int const>(), R"doc(
Construct an unattached integer constant (also the default constructor).

Parameters
----------
i : int, optional
    Integer value. Default is 0.
)doc",
           py::arg("i") = 0)
      .def(py::init<double const&>(), R"doc(
Construct an unattached real constant.

Parameters
----------
d : float
    Real value.
)doc",
           py::arg("d"))
      .def(py::init<mc::FFBase*, std::string const&>(), R"doc(
Construct a new original variable attached to a DAG.

A new variable node is appended to the DAG ``dag``.

Parameters
----------
dag : FFGraph
    DAG environment recording the factorable function.
name : str, optional
    Variable name used in string representations and DOT scripts.
    Default is an automatically generated name.
)doc",
           py::arg("dag"), py::arg("name") = "")
      .def(py::init<mc::FFVar const&>(), R"doc(
Copy constructor; the copy refers to the same DAG node.
)doc",
           py::arg("var"))
      .def("set",
           py::overload_cast<mc::FFBase*, std::string const&>(&mc::FFVar::set),
           R"doc(
Attach this variable to a DAG as a new original variable.

Parameters
----------
dag : FFGraph
    DAG environment to attach the variable to.
name : str, optional
    Variable name. Default is an automatically generated name.
)doc",
           py::arg("dag"), py::arg("name") = "")
      .def("set", py::overload_cast<int const>(&mc::FFVar::set, py::const_),
           R"doc(
Fix the variable at an integer constant value.

The node remains in the DAG; it is marked constant with value ``i``.
Use `unset` to make it a free variable again.

Parameters
----------
i : int, optional
    Constant value. Default is 0.
)doc",
           py::arg("i") = 0)
      .def("set", py::overload_cast<double const&>(&mc::FFVar::set, py::const_),
           R"doc(
Fix the variable at a real constant value.

The node remains in the DAG; it is marked constant with value ``d``.
Use `unset` to make it a free variable again.

Parameters
----------
d : float
    Constant value.
)doc",
           py::arg("d"))
      .def("set",
           py::overload_cast<std::string const&>(&mc::FFVar::set, py::const_),
           R"doc(
Set the variable name.

Parameters
----------
name : str
    New name, used in string representations and DOT scripts.
)doc",
           py::arg("name"))
      .def("unset", &mc::FFVar::unset, R"doc(
Release a constant value previously fixed with `set`, making the node a
free variable again.
)doc")
      .def("num", &mc::FFVar::num, R"doc(
Return the numeric field of the node.

Returns
-------
num : FFNum
    Numeric field holding the constant value; only meaningful for
    constant nodes (see `cst`).
)doc")
      .def("cst", &mc::FFVar::cst, R"doc(
Return whether the variable is currently a constant.

Returns
-------
cst : bool
    True if the node is a constant or was fixed with `set`.
)doc")
      .def("dag", &mc::FFVar::dag, R"doc(
Return the DAG environment the variable is attached to, or None if the
variable is unattached.
)doc")
      .def("str",
           [](mc::FFVar const& V) { return mc::FFExpr::dep(V).ostr().str(); },
           R"doc(
Return a one-line infix string expression of this node in terms of the
DAG variables (e.g. ``"X0*exp(X0*X1)+sqr(X1)"``), by traversing its
subgraph.
)doc")
      .def_property_readonly("opdef",
                             py::overload_cast<>(&mc::FFVar::opdef, py::const_),
                             R"doc(
Defining operation of this node, as a tuple ``(op, k)`` where ``op`` is
the `FFOp` producing the node and ``k`` the index of this node among the
operation's outputs. ``op`` is None for unreferenced constants.
)doc")
      .def_property_readonly("id",
                             py::overload_cast<>(&mc::FFVar::id, py::const_),
                             R"doc(
Identifier of the node, as a tuple ``(type, index)`` with ``type`` an
`FFVar.TYPE` value and ``index`` the node index within that type.
)doc")
      .def("__str__",
           [](mc::FFVar const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
           })
      .def("__repr__",
           [](mc::FFVar const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
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
      .def("__pow__",
           [](mc::FFVar const& V, int const n) { return mc::pow(V, n); })
      .def("__pow__",
           [](mc::FFVar const& V, double const& r) { return mc::pow(V, r); })
      .def("__pow__", [](mc::FFVar const& V, mc::FFVar const& W)
           { return mc::pow(V, W); });

  m.def("inv", [](mc::FFVar const& V) { return mc::inv(V); },
        "FFVar overload: add a DAG node for the reciprocal 1/x and return it.",
        py::arg("x"));
  m.def("sqr", [](mc::FFVar const& V) { return mc::sqr(V); },
        "FFVar overload: add a DAG node for the square x**2 and return it.",
        py::arg("x"));
  m.def("sqrt", [](mc::FFVar const& V) { return mc::sqrt(V); },
        "FFVar overload: add a DAG node for the square root of x and return "
        "it.",
        py::arg("x"));
  m.def("exp", [](mc::FFVar const& V) { return mc::exp(V); },
        "FFVar overload: add a DAG node for the exponential exp(x) and return "
        "it.",
        py::arg("x"));
  m.def("log", [](mc::FFVar const& V) { return mc::log(V); },
        "FFVar overload: add a DAG node for the natural logarithm log(x) and "
        "return it.",
        py::arg("x"));
  m.def("cos", [](mc::FFVar const& V) { return mc::cos(V); },
        "FFVar overload: add a DAG node for cos(x) and return it.",
        py::arg("x"));
  m.def("sin", [](mc::FFVar const& V) { return mc::sin(V); },
        "FFVar overload: add a DAG node for sin(x) and return it.",
        py::arg("x"));
  m.def("tan", [](mc::FFVar const& V) { return mc::tan(V); },
        "FFVar overload: add a DAG node for tan(x) and return it.",
        py::arg("x"));
  m.def("acos", [](mc::FFVar const& V) { return mc::acos(V); },
        "FFVar overload: add a DAG node for arccos(x) and return it.",
        py::arg("x"));
  m.def("asin", [](mc::FFVar const& V) { return mc::asin(V); },
        "FFVar overload: add a DAG node for arcsin(x) and return it.",
        py::arg("x"));
  m.def("atan", [](mc::FFVar const& V) { return mc::atan(V); },
        "FFVar overload: add a DAG node for arctan(x) and return it.",
        py::arg("x"));
  m.def("cosh", [](mc::FFVar const& V) { return mc::cosh(V); },
        "FFVar overload: add a DAG node for cosh(x) and return it.",
        py::arg("x"));
  m.def("sinh", [](mc::FFVar const& V) { return mc::sinh(V); },
        "FFVar overload: add a DAG node for sinh(x) and return it.",
        py::arg("x"));
  m.def("tanh", [](mc::FFVar const& V) { return mc::tanh(V); },
        "FFVar overload: add a DAG node for tanh(x) and return it.",
        py::arg("x"));
  m.def("fabs", [](mc::FFVar const& V) { return mc::fabs(V); },
        "FFVar overload: add a DAG node for the absolute value |x| and return "
        "it.",
        py::arg("x"));
  m.def("relu", [](mc::FFVar const& V) { return mc::max(V, 0.); },
        "FFVar overload: add a DAG node for the rectifier max(x, 0) and "
        "return it.",
        py::arg("x"));
  m.def("xlog", [](mc::FFVar const& V) { return mc::xlog(V); },
        "FFVar overload: add a DAG node for x*log(x) and return it.",
        py::arg("x"));
  m.def("fstep", [](mc::FFVar const& V) { return mc::fstep(V); },
        "FFVar overload: add a DAG node for the forward unit step at 0 "
        "(1 if x >= 0, else 0) and return it.",
        py::arg("x"));
  m.def("bstep", [](mc::FFVar const& V) { return mc::bstep(V); },
        "FFVar overload: add a DAG node for the backward unit step at 0 "
        "(0 if x >= 0, else 1) and return it.",
        py::arg("x"));
  m.def("erf", [](mc::FFVar const& V) { return mc::erf(V); },
        "FFVar overload: add a DAG node for the error function erf(x) and "
        "return it.",
        py::arg("x"));
  m.def("erfc", [](mc::FFVar const& V) { return mc::erfc(V); },
        "FFVar overload: add a DAG node for the complementary error function "
        "erfc(x) and return it.",
        py::arg("x"));
  m.def("pow", [](I const& x, int const n) { return mc::Op<I>::pow(x, n); },
        "Interval overload: enclosure of x**n (integer exponent) computed "
        "with interval arithmetic.",
        py::arg("x"), py::arg("n"));
  m.def("pow",
        [](I const& x, double const& r) { return mc::Op<I>::pow(x, r); },
        "Interval overload: enclosure of x**r (real exponent) computed with "
        "interval arithmetic.",
        py::arg("x"), py::arg("r"));
  m.def("pow", [](I const& x, I const& y) { return mc::Op<I>::pow(x, y); },
        "Interval overload: enclosure of x**y for interval exponent y.",
        py::arg("x"), py::arg("y"));
  m.def("pow",
        [](double const& r, I const& y) { return mc::Op<I>::pow(r, y); },
        "Interval overload: enclosure of r**y for real base r and interval "
        "exponent y.",
        py::arg("r"), py::arg("y"));
  m.def("pow", [](mc::FFVar const& V, int const n) { return mc::pow(V, n); },
        "FFVar overload: add a DAG node for x**n (integer exponent) and "
        "return it.",
        py::arg("x"), py::arg("n"));
  m.def("pow",
        [](mc::FFVar const& V, double const& r) { return mc::pow(V, r); },
        "FFVar overload: add a DAG node for x**r (real exponent) and return "
        "it.",
        py::arg("x"), py::arg("r"));
  m.def("pow",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::pow(V, W); },
        "FFVar overload: add a DAG node for x**y with FFVar exponent y and "
        "return it.",
        py::arg("x"), py::arg("y"));
  m.def("pow",
        [](double const& r, mc::FFVar const& W) { return mc::pow(r, W); },
        "FFVar overload: add a DAG node for r**y with real base r and FFVar "
        "exponent y, and return it.",
        py::arg("r"), py::arg("y"));
  m.def("cheb",
        [](mc::FFVar const& V, unsigned const n) { return mc::cheb(V, n); },
        "FFVar overload: add a DAG node for the Chebyshev polynomial of the "
        "first kind T_n(x) and return it.",
        py::arg("x"), py::arg("n"));
  m.def("max",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::max(V, W); },
        "FFVar overload: add a DAG node for the maximum of x and y and return "
        "it.",
        py::arg("x"), py::arg("y"));
  m.def("min",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::min(V, W); },
        "FFVar overload: add a DAG node for the minimum of x and y and return "
        "it.",
        py::arg("x"), py::arg("y"));
  m.def("inter",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::inter(V, W); },
        "FFVar overload: add a DAG node representing the intersection of x "
        "and y, which is meaningful when the DAG is evaluated in a set "
        "arithmetic (e.g. intervals), and return it.",
        py::arg("x"), py::arg("y"));

  py::enum_<mc::FFVar::TYPE>(pyFFVar, "TYPE",
                             "Kind of a DAG node (see `FFVar.id`).")
      .value("VAR", mc::FFVar::TYPE::VAR, "Original (independent) variable.")
      .value("AUX", mc::FFVar::TYPE::AUX,
             "Auxiliary variable, defined as the result of an operation.")
      .value("CINT", mc::FFVar::TYPE::CINT, "Integer constant.")
      .value("CREAL", mc::FFVar::TYPE::CREAL, "Real constant.")
      .export_values();

  py::class_<mc::FFOp> pyFFOp(m, "FFOp", R"doc(
Operation node in the DAG of a factorable function.

An `FFOp` links its input nodes (`varin`) to the output nodes it defines
(`varout`). Operations are created implicitly while building expressions
from `FFVar` operands; they are typically inspected via `FFVar.opdef` or
by traversing an `FFSubgraph`.
)doc");
  pyFFOp
      .def_readwrite("type", &mc::FFOp::type,
                     "Operation type, as an `FFOp.TYPE` enumeration value.")
      .def_readwrite("varin", &mc::FFOp::varin,
                     "Input (operand) FFVar nodes of the operation.")
      .def_readwrite("varout", &mc::FFOp::varout,
                     "Output FFVar nodes defined by the operation.")
      .def("name", &mc::FFOp::name,
           "Return the operation name as a string (e.g. 'EXP').")
      .def("__str__",
           [](mc::FFOp const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFOp const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           });

  py::enum_<mc::FFOp::TYPE>(pyFFOp, "TYPE",
                            "Type of a DAG operation (see `FFOp.type`).")
      .value("CNST", mc::FFOp::TYPE::CNST, "Constant.")
      .value("VAR", mc::FFOp::TYPE::VAR, "Original variable.")
      .value("PLUS", mc::FFOp::TYPE::PLUS, "Binary addition.")
      .value("SHIFT", mc::FFOp::TYPE::SHIFT, "Addition of a constant.")
      .value("NEG", mc::FFOp::TYPE::NEG, "Unary negation.")
      .value("MINUS", mc::FFOp::TYPE::MINUS, "Binary subtraction.")
      .value("TIMES", mc::FFOp::TYPE::TIMES, "Binary multiplication.")
      .value("SCALE", mc::FFOp::TYPE::SCALE, "Multiplication by a constant.")
      .value("DIV", mc::FFOp::TYPE::DIV, "Binary division.")
      .value("INV", mc::FFOp::TYPE::INV, "Reciprocal 1/x.")
      .value("PROD", mc::FFOp::TYPE::PROD, "N-ary product.")
      .value("IPOW", mc::FFOp::TYPE::IPOW, "Power with integer exponent.")
      .value("DPOW", mc::FFOp::TYPE::DPOW, "Power with real exponent.")
      .value("CHEB", mc::FFOp::TYPE::CHEB,
             "Chebyshev polynomial of the first kind.")
      .value("SQR", mc::FFOp::TYPE::SQR, "Square x**2.")
      .value("SQRT", mc::FFOp::TYPE::SQRT, "Square root.")
      .value("EXP", mc::FFOp::TYPE::EXP, "Exponential.")
      .value("LOG", mc::FFOp::TYPE::LOG, "Natural logarithm.")
      .value("XLOG", mc::FFOp::TYPE::XLOG, "x*log(x).")
      .value("SIN", mc::FFOp::TYPE::SIN, "Sine.")
      .value("COS", mc::FFOp::TYPE::COS, "Cosine.")
      .value("TAN", mc::FFOp::TYPE::TAN, "Tangent.")
      .value("ASIN", mc::FFOp::TYPE::ASIN, "Inverse sine.")
      .value("ACOS", mc::FFOp::TYPE::ACOS, "Inverse cosine.")
      .value("ATAN", mc::FFOp::TYPE::ATAN, "Inverse tangent.")
      .value("SINH", mc::FFOp::TYPE::SINH, "Hyperbolic sine.")
      .value("COSH", mc::FFOp::TYPE::COSH, "Hyperbolic cosine.")
      .value("TANH", mc::FFOp::TYPE::TANH, "Hyperbolic tangent.")
      .value("ERF", mc::FFOp::TYPE::ERF, "Error function.")
      .value("FABS", mc::FFOp::TYPE::FABS, "Absolute value.")
      .value("FSTEP", mc::FFOp::TYPE::FSTEP, "Forward unit step at 0.")
      .value("MINF", mc::FFOp::TYPE::MINF, "Binary minimum.")
      .value("MAXF", mc::FFOp::TYPE::MAXF, "Binary maximum.")
      .value("INTER", mc::FFOp::TYPE::INTER,
             "Intersection of operands (set arithmetics).")
      .value("EXTERN", mc::FFOp::TYPE::EXTERN,
             "External (user-defined) operation.")
      .export_values();

  py::class_<mc::FFSubgraph> pyFFSubgraph(m, "FFSubgraph", R"doc(
Subgraph of a DAG: the ordered list of operations needed to evaluate a
given subset of dependents.

Instances are created by `FFBase.subgraph` and can be passed to
`FFGraph.eval`, `FFGraph.reval` and `FFGraph.veval` to avoid
re-extracting the operation list on every evaluation of the same
dependents.
)doc");
  pyFFSubgraph.def(py::init<>(), "Construct an empty subgraph.")
      .def(py::init<mc::FFSubgraph const&>(), "Copy constructor.")
      .def("clear", &mc::FFSubgraph::clear, R"doc(
Reset to an empty subgraph.
)doc")
      .def_readonly("len_tap", &mc::FFSubgraph::len_tap,
                    "Length of the work array (evaluation tape) needed to "
                    "evaluate the subgraph.")
      .def_readonly("len_wrk", &mc::FFSubgraph::len_wrk,
                    "Length of the extra work array used for moving n-ary "
                    "operations during evaluation.");

  py::class_<mc::FFBase> pyFFBase(m, "FFBase", R"doc(
Base DAG environment of a factorable function.

`FFBase` stores the nodes (`FFVar`) and operations (`FFOp`) of the
directed acyclic graph and provides construction and inspection
facilities: adding variables, extracting subgraphs, and printing or
exporting them. Use the derived class `FFGraph` for differentiation and
evaluation capabilities.
)doc");
  pyFFBase.def(py::init<>(), "Construct an empty DAG environment.")
      .def(
          "add_var", [](mc::FFBase& G, std::string const& name)
          { return G.add_var(name); }, py::arg("name") = "",
          R"doc(
Append a new original variable to the DAG.

Parameters
----------
name : str, optional
    Variable name. Default is an automatically generated name.

Returns
-------
var : FFVar
    The new variable node.
)doc")
      .def(
          "add_vars", [](mc::FFBase& G, size_t dim, std::string const& name)
          { return G.add_vars(dim, name); }, py::arg("dim"),
          py::arg("name") = "",
          R"doc(
Append several new original variables to the DAG.

Parameters
----------
dim : int
    Number of variables to add.
name : str, optional
    Common base name for the new variables. Default is an automatically
    generated name.

Returns
-------
vars : list of FFVar
    The new variable nodes.
)doc")
      .def("clear", &mc::FFBase::clear, R"doc(
Erase all variables and operations from the DAG.

Any FFVar handle referring to this DAG becomes invalid.
)doc")
      .def(
          "subgraph", [](mc::FFBase& G, std::vector<mc::FFVar const*> const& V)
          { return G.subgraph(V); }, py::arg("vDep"),
          R"doc(
Extract the subgraph of operations needed to evaluate given dependents.

Parameters
----------
vDep : list of FFVar
    Dependent nodes whose evaluation subgraph is sought.

Returns
-------
sgDep : FFSubgraph
    Operations participating in the evaluation of ``vDep``, in order of
    appearance. Pass it to `FFGraph.eval`, `FFGraph.reval` or
    `FFGraph.veval` to avoid re-extracting it on repeated evaluations.
)doc")
      .def(
          "output",
          [](mc::FFBase& G, std::vector<mc::FFVar const*> const& V,
             std::string const& S) { mc::FFBase::output(G.subgraph(V), S); },
          py::arg("vDep"), py::arg("str") = "",
          R"doc(
Print the subgraph of the dependents ``vDep`` to standard output.

Parameters
----------
vDep : list of FFVar
    Dependents whose subgraph is displayed.
str : str, optional
    Text appended to the header lines of the printout. Default is "".
)doc")
      .def(
          "output",
          [](mc::FFBase const& G, mc::FFSubgraph const& SG,
             std::string const& S) { mc::FFBase::output(SG, S); },
          py::arg("sgDep"), py::arg("str") = "",
          R"doc(
Print a precomputed subgraph to standard output.

Parameters
----------
sgDep : FFSubgraph
    Subgraph to display, as returned by `subgraph`.
str : str, optional
    Text appended to the header lines of the printout. Default is "".
)doc")
      .def(
          "dot_script",
          [](mc::FFBase const& G, std::vector<mc::FFVar const*> const& V,
             std::string const& fname)
          {
            if (fname == "") return G.dot_script(V);
            std::ofstream ofs(fname);
            return G.dot_script(V, ofs);
          },
          py::arg("vDep"), py::arg("fname"),
          R"doc(
Generate a DOT (Graphviz) script depicting the subgraph of ``vDep``.

Parameters
----------
vDep : list of FFVar
    Dependents whose subgraph is exported.
fname : str
    Output file name. If empty, the script is printed to standard
    output instead.
)doc")
      .def("__str__",
           [](mc::FFBase const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def("__repr__",
           [](mc::FFBase const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           });

  py::class_<mc::FFGraph, mc::FFBase> pyFFGraph(m, "FFGraph", R"doc(
DAG environment for construction, differentiation and evaluation of
factorable functions.

An `FFGraph` records the directed acyclic graph (DAG) of factorable
expressions built from its `FFVar` variables. On top of the storage and
inspection facilities inherited from `FFBase`, it provides:

- symbolic differentiation: `fdiff` (forward mode), `bdiff` (reverse
  mode), both returning sparse Jacobians as DAG nodes, and `tdiff`
  (Taylor expansion of ODE solutions);
- DAG manipulation: `compose`, `insert`, `substitute`;
- evaluation of any subset of dependents in a range of arithmetics via
  `eval` (floats, `Interval`, `McCormick`, `Specbnd`, Taylor/Chebyshev
  models, superposition models, polyhedral and ellipsoidal images,
  dependency and invariant detection);
- reverse (constraint) propagation via `reval` and vectorized
  multi-scenario evaluation via `veval`.

Behavior is controlled by the `options` attribute (`FFGraph.Options`).

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = [X[0] * X[1] - 1, pymcpp.exp(X[0]) + X[1]]
>>> DAG.eval(F, X, [1.0, 2.0])
[1.0, 4.718281828459045]
>>> rows, cols, jac = DAG.bdiff(F, X)  # sparse Jacobian nodes
)doc");
  py::class_<mc::FFGraph::Options> pyFFGraphOptions(pyFFGraph, "Options",
                                                    R"doc(
Option set of an `FFGraph`, accessed via the `FFGraph.options` attribute.

Fields can be assigned directly, e.g. ``DAG.options.MAXTHREAD = 4``.
)doc");

  pyFFGraph.def(py::init<>(), "Construct an empty DAG environment.")
      .def_readwrite("options", &mc::FFGraph::options,
                     "Option set of this DAG (an `FFGraph.Options` instance).")
      .def(
          "fdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDFAD(vDep, vIndep, std::vector<mc::FFVar const*>()); },
          py::arg("vDep"), py::arg("vIndep"),
          py::return_value_policy::reference_internal,
          R"doc(
Compute the sparse Jacobian of the dependents by forward-mode automatic
differentiation.

The differentiation is symbolic: new DAG nodes representing the nonzero
partial derivatives dF[i]/dx[j] are appended to this graph; no numerical
values are computed. The result is returned in sparse coordinate
(triplet) format.

Parameters
----------
vDep : list of FFVar
    Dependents (function outputs) to differentiate.
vIndep : list of FFVar
    Independents (variables) to differentiate with respect to.

Returns
-------
rows : list of int
    Row index (position in ``vDep``) of each nonzero Jacobian entry.
cols : list of int
    Column index (position in ``vIndep``) of each nonzero entry.
jac : list of FFVar
    DAG nodes representing the nonzero entries dF[rows[k]]/dx[cols[k]],
    added to this DAG. They can be evaluated with `eval` or
    differentiated further to obtain higher-order derivatives.

See Also
--------
bdiff : reverse-mode counterpart, often more economical when there are
    fewer dependents than independents.
tdiff : Taylor expansion of ODE solutions.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = [X[0] * X[1], X[1] ** 2]
>>> rows, cols, jac = DAG.fdiff(F, X)
>>> J = [[0.0] * len(X) for _ in F]  # dense reconstruction
>>> for k in range(len(jac)):
...     J[rows[k]][cols[k]] = DAG.eval([jac[k]], X, [3.0, 2.0])[0]
>>> J
[[2.0, 3.0], [0.0, 4.0]]
)doc")
      .def(
          "fdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep,
             std::vector<mc::FFVar const*> const& vDir)
          { return G.SDFAD(vDep, vIndep, vDir); },
          py::arg("vDep"), py::arg("vIndep"), py::arg("vDir"),
          py::return_value_policy::reference_internal,
          R"doc(
Directional variant: Jacobian-vector product by forward-mode automatic
differentiation.

Appends DAG nodes for the directional derivatives
(dF/dx * vDir)[i] = sum_j dF[i]/dx[j] * vDir[j] instead of the full
Jacobian.

Parameters
----------
vDep : list of FFVar
    Dependents to differentiate.
vIndep : list of FFVar
    Independents to differentiate with respect to.
vDir : list of FFVar
    Direction, one entry per independent (same length as ``vIndep``).

Returns
-------
rows : list of int
    Row index (position in ``vDep``) of each nonzero entry.
cols : list of int
    All zero, since the Jacobian-vector product is a single column.
jvp : list of FFVar
    DAG nodes representing the nonzero entries
    (dF/dx * vDir)[rows[k]], added to this DAG.
)doc")
      .def(
          "bdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDBAD(vDep, std::vector<mc::FFVar const*>(), vIndep); },
          py::arg("vDep"), py::arg("vIndep"),
          py::return_value_policy::reference_internal,
          R"doc(
Compute the sparse Jacobian of the dependents by reverse-mode (adjoint)
automatic differentiation.

Same sparse coordinate-format result as `fdiff`, but the derivative
nodes are constructed with the reverse mode, which usually yields a more
economical DAG when there are fewer dependents than independents.

Parameters
----------
vDep : list of FFVar
    Dependents (function outputs) to differentiate.
vIndep : list of FFVar
    Independents (variables) to differentiate with respect to.

Returns
-------
rows : list of int
    Row index (position in ``vDep``) of each nonzero Jacobian entry.
cols : list of int
    Column index (position in ``vIndep``) of each nonzero entry.
jac : list of FFVar
    DAG nodes representing the nonzero entries dF[rows[k]]/dx[cols[k]],
    added to this DAG.

See Also
--------
fdiff : forward-mode counterpart.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = [X[0] * X[1] - 1, pymcpp.exp(X[0]) + X[1]]
>>> rows, cols, jac = DAG.bdiff(F, X)
>>> J = [[0.0] * len(X) for _ in F]  # dense reconstruction
>>> for k in range(len(jac)):
...     J[rows[k]][cols[k]] = DAG.eval([jac[k]], X, [3.0, 2.0])[0]
>>> J
[[2.0, 3.0], [20.085536923187668, 1.0]]
)doc")
      .def(
          "bdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vDir,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDBAD(vDep, vDir, vIndep); },
          py::arg("vDep"), py::arg("vDir"), py::arg("vIndep"),
          py::return_value_policy::reference_internal,
          R"doc(
Directional variant: vector-Jacobian product by reverse-mode automatic
differentiation.

Appends DAG nodes for the adjoint directional derivatives
(vDir . dF/dx)[j] = sum_i vDir[i] * dF[i]/dx[j] instead of the full
Jacobian.

Parameters
----------
vDep : list of FFVar
    Dependents to differentiate.
vDir : list of FFVar
    Adjoint direction, one entry per dependent (same length as
    ``vDep``).
vIndep : list of FFVar
    Independents to differentiate with respect to.

Returns
-------
rows : list of int
    All zero, since the vector-Jacobian product is a single row.
cols : list of int
    Column index (position in ``vIndep``) of each nonzero entry.
vjp : list of FFVar
    DAG nodes representing the nonzero entries
    (vDir . dF/dx)[cols[k]], added to this DAG.
)doc")
      .def(
          "tdiff",
          [](mc::FFGraph& G, unsigned int const ordermax,
             std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vVar,
             mc::FFVar const* const pIndep)
          { return G.TAD(ordermax, vDep, vVar, pIndep); },
          py::arg("ordermax"), py::arg("vDep"), py::arg("vVar"),
          py::arg("pIndep"), py::return_value_policy::reference_internal,
          R"doc(
Expand the DAG with Taylor coefficients of an ODE solution.

The dependents ``vDep`` are interpreted as the right-hand side f of the
ODE system dx/dt = f(x), with matching state variables ``vVar`` (same
length as ``vDep``). DAG nodes are appended for the Taylor coefficients
phi_q of the ODE solution in time, defined recursively by phi_0 := x and
phi_q := (1/q) * d(phi_{q-1})/dx * f(x) for q >= 1.

Parameters
----------
ordermax : int
    Maximum expansion order.
vDep : list of FFVar
    Right-hand-side expressions f, one per state variable.
vVar : list of FFVar
    State variables x; must have the same length as ``vDep``.
pIndep : FFVar or None
    Independent (time-like) variable, to account for a non-autonomous
    right-hand side; pass None for an autonomous system.

Returns
-------
coef : list of FFVar
    Taylor coefficient nodes, of length (ordermax + 1) * len(vDep),
    ordered by increasing order with all states contiguous at each
    order: ``coef[q * len(vDep) + j]`` is the coefficient phi_q of
    state j. The 0th-order block consists of the state variables
    themselves.

Examples
--------
For dx/dt = x, the Taylor coefficients of the solution are x/q!:

>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> T = pymcpp.FFVar(DAG, "T")
>>> X = pymcpp.FFVar(DAG, "X")
>>> coef = DAG.tdiff(3, [X], [X], T)
>>> DAG.eval(coef, [X], [1.0])  # [1, 1, 1/2, 1/6]
[1.0, 1.0, 0.5, 0.16666666666666666]
)doc")
      .def(
          "compose",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDepOut,
             std::vector<std::pair<mc::FFVar const*, mc::FFVar const*>> const&
                 vDepIn) { return G.compose(vDepOut, vDepIn); },
          py::arg("vDepOut"), py::arg("vDepIn"),
          py::return_value_policy::reference_internal,
          R"doc(
Compose dependents with inner expressions, replacing DAG variables.

Each pair ``(var, expr)`` in ``vDepIn`` substitutes the DAG variable
``var`` by the expression ``expr`` inside the dependents ``vDepOut``.
New DAG nodes for the composed expressions are appended to this graph.

Parameters
----------
vDepOut : list of FFVar
    Outer dependents in which the substitution is performed.
vDepIn : list of tuple of (FFVar, FFVar)
    Pairs ``(var, expr)``; ``var`` must be an original (leaf) DAG
    variable. Use `substitute` to replace auxiliary variables.

Returns
-------
vDep : list of FFVar
    New dependents representing the composed expressions, in ``vDepOut``
    order.
)doc")
      .def(
          "insert",
          [](mc::FFGraph& G, mc::FFGraph& dag,
             std::vector<mc::FFVar> const& vDepIn,
             std::vector<mc::FFVar> vDepOut)
          {
            G.insert(&dag, vDepIn, vDepOut);
            return vDepOut;
          },
          py::arg("dag"), py::arg("vDepIn"),
          py::arg("vDepOut") = std::vector<mc::FFVar>(),
          R"doc(
Insert dependents from another DAG into this graph.

The subgraphs of the dependents ``vDepIn`` of ``dag`` are copied into
this graph. Participating variables are matched by index: they share the
same indices in both DAGs.

Parameters
----------
dag : FFGraph
    Source DAG holding the dependents to insert.
vDepIn : list of FFVar
    Dependents of ``dag`` to insert.
vDepOut : list of FFVar, optional
    Initial output list; entries are overwritten with the inserted
    dependents. Default is an empty list.

Returns
-------
vDepOut : list of FFVar
    Nodes of this graph representing the inserted dependents, in
    ``vDepIn`` order.
)doc")
      .def(
          "substitute",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDepOut,
             std::vector<mc::FFVar> const& vAuxTarg,
             std::vector<mc::FFVar> const& vAuxSubst)
          { return G.substitute(vDepOut, vAuxTarg, vAuxSubst); },
          py::arg("vDepOut"), py::arg("vAuxTarg"), py::arg("vAuxSubst"),
          R"doc(
Substitute variables or auxiliaries inside dependents.

Each target node ``vAuxTarg[k]`` is replaced by the expression
``vAuxSubst[k]`` within the dependents ``vDepOut``. Unlike `compose`,
the targets may be auxiliary (non-leaf) nodes produced by intermediate
operations; the subgraph is pruned at the substitution targets so their
now-unused upstream operations are not evaluated.

Parameters
----------
vDepOut : list of FFVar
    Dependents in which the substitution is performed.
vAuxTarg : list of FFVar
    Target nodes to be replaced (variables or auxiliaries).
vAuxSubst : list of FFVar
    Replacement expressions, in the same order as ``vAuxTarg``.

Returns
-------
vDep : list of FFVar
    New dependents with the substitutions applied, in ``vDepOut`` order.
)doc")
      .def("__str__",
           [](mc::FFGraph const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def("__repr__",
           [](mc::FFGraph const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<double> const& DVar)
          {
            size_t const nDep = vDep.size();
            std::vector<double> DDep(nDep);
            G.eval(SgDep, vDep, DDep, vVar, DVar);
            return DDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("DVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in double-precision arithmetic, using a precomputed
subgraph.

The DAG operations in ``SgDep`` are propagated with the variable values
``DVar`` to produce the values of the dependents ``vDep``.

Parameters
----------
SgDep : FFSubgraph
    Subgraph of the dependents, as returned by `subgraph`. Passing it
    avoids re-extracting the operation list on every call, e.g. in
    repeated evaluations of the same dependents.
vDep : list of FFVar
    Dependent nodes (outputs) to evaluate.
vVar : list of FFVar
    Variable nodes whose values are set.
DVar : list of float
    Values of the variables, in the same order as ``vVar``.

Returns
-------
DDep : list of float
    Values of the dependents, in the same order as ``vDep``. A new list
    is returned; the inputs are not modified.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<double> const& DVar)
          {
            size_t const nDep = vDep.size();
            std::vector<double> DDep(nDep);
            G.eval(vDep, DDep, vVar, DVar);
            return DDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("DVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in double-precision arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.

Parameters
----------
vDep : list of FFVar
    Dependent nodes to evaluate.
vVar : list of FFVar
    Variable nodes whose values are set.
DVar : list of float
    Values of the variables, in ``vVar`` order.

Returns
-------
DDep : list of float
    Values of the dependents, in ``vDep`` order.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = [X[0] * X[1] - 1, pymcpp.exp(X[0]) + X[1]]
>>> DAG.eval(F, X, [1.0, 2.0])
[1.0, 4.718281828459045]
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I> const& IVar)
          {
            size_t const nDep = vDep.size();
            std::vector<I> IDep(nDep);
            G.eval(SgDep, vDep, IDep, vVar, IVar);
            return IDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("IVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in interval arithmetic, using a precomputed
subgraph.

Given `Interval` bounds ``IVar`` on the variables ``vVar``, returns
guaranteed `Interval` enclosures of the ranges of the dependents
``vDep`` (in matching order). ``SgDep`` avoids re-extracting the
subgraph on every call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I> const& IVar)
          {
            size_t const nDep = vDep.size();
            std::vector<I> IDep(nDep);
            G.eval(vDep, IDep, vVar, IVar);
            return IDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("IVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in interval arithmetic.

Given `Interval` bounds ``IVar`` on the variables ``vVar``, returns
guaranteed `Interval` enclosures of the dependents ``vDep`` (in matching
order); the subgraph is extracted internally on each call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<MC> const& MCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<MC> MCDep(nDep);
            G.eval(SgDep, vDep, MCDep, vVar, MCVar);
            return MCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("MCVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in McCormick relaxation arithmetic, using a
precomputed subgraph.

Given `McCormick` variables ``MCVar`` (interval bounds plus convex and
concave relaxation values at a point), returns `McCormick` objects for
the dependents ``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<MC> const& MCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<MC> MCDep(nDep);
            G.eval(vDep, MCDep, vVar, MCVar);
            return MCDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("MCVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in McCormick relaxation arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SB> const& SBVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SB> SBDep(nDep);
            G.eval(SgDep, vDep, SBDep, vVar, SBVar);
            return SBDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("SBVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in spectral bound arithmetic, using a precomputed
subgraph.

Given `Specbnd` variables ``SBVar``, returns `Specbnd` objects for the
dependents ``vDep`` (in matching order), enclosing the spectrum of their
Hessian matrices.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SB> const& SBVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SB> SBDep(nDep);
            G.eval(vDep, SBDep, vVar, SBVar);
            return SBDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("SBVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in spectral bound arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<T> const& TVar)
          {
            size_t const nDep = vDep.size();
            std::vector<T> TDep(nDep);
            G.eval(SgDep, vDep, TDep, vVar, TVar);
            return TDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("TVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in Taylor model arithmetic, using a precomputed
subgraph.

Given `TVar` variables ``TVar`` (multivariate Taylor polynomial plus
guaranteed remainder bound), returns `TVar` models of the dependents
``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<T> const& TVar)
          {
            size_t const nDep = vDep.size();
            std::vector<T> TDep(nDep);
            G.eval(vDep, TDep, vVar, TVar);
            return TDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("TVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in Taylor model arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<C> const& CVar)
          {
            size_t const nDep = vDep.size();
            std::vector<C> CDep(nDep);
            G.eval(SgDep, vDep, CDep, vVar, CVar);
            return CDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("CVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dense Chebyshev model arithmetic, using a
precomputed subgraph.

Given `CVar` variables ``CVar`` (dense Chebyshev polynomial expansion
plus guaranteed remainder bound), returns `CVar` models of the
dependents ``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<C> const& CVar)
          {
            size_t const nDep = vDep.size();
            std::vector<C> CDep(nDep);
            G.eval(vDep, CDep, vVar, CVar);
            return CDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("CVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dense Chebyshev model arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SC> const& SCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SC> SCDep(nDep);
            G.eval(SgDep, vDep, SCDep, vVar, SCVar);
            return SCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("SCVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in sparse Chebyshev model arithmetic, using a
precomputed subgraph.

Given `SCVar` variables ``SCVar`` (sparse Chebyshev polynomial expansion
plus guaranteed remainder bound), returns `SCVar` models of the
dependents ``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SC> const& SCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SC> SCDep(nDep);
            G.eval(vDep, SCDep, vVar, SCVar);
            return SCDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("SCVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in sparse Chebyshev model arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWLS> const& PWLSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWLS> PWLSDep(nDep);
            G.eval(SgDep, vDep, PWLSDep, vVar, PWLSVar);
            return PWLSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"),
          py::arg("PWLSVar"), py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-linear univariate estimators, using a precomputed subgraph.

Given `PWLSVar` variables ``PWLSVar``, returns `PWLSVar` relaxations of
the dependents ``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWLS> const& PWLSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWLS> PWLSDep(nDep);
            G.eval(vDep, PWLSDep, vVar, PWLSVar);
            return PWLSDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("PWLSVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-linear univariate estimators.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWCS> const& PWCSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWCS> PWCSDep(nDep);
            G.eval(SgDep, vDep, PWCSDep, vVar, PWCSVar);
            return PWCSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"),
          py::arg("PWCSVar"), py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-constant univariate estimators, using a precomputed subgraph.

Given `PWCSVar` variables ``PWCSVar``, returns `PWCSVar` relaxations of
the dependents ``vDep`` (in matching order).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWCS> const& PWCSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWCS> PWCSDep(nDep);
            G.eval(vDep, PWCSDep, vVar, PWCSVar);
            return PWCSDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("PWCSVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-constant univariate estimators.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<PV> const& PVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PV> PVDep(nDep);
            G.eval(SgDep, vDep, PVDep, vVar, PVVar);
            return PVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("PVVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in polyhedral image arithmetic, using a precomputed
subgraph.

Given `PolVar` variables ``PVVar`` attached to a `PolImg` environment,
returns `PolVar` objects for the dependents ``vDep`` (in matching
order); auxiliary variables and linear cuts describing a polyhedral
relaxation of the dependents are added to the `PolImg` environment.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<PV> const& PVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PV> PVDep(nDep);
            G.eval(vDep, PVDep, vVar, PVVar);
            return PVDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("PVVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in polyhedral image arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<EV> const& EVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<EV> EVDep(nDep);
            G.eval(SgDep, vDep, EVDep, vVar, EVVar);
            return EVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("EVVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in ellipsoidal image arithmetic, using a precomputed
subgraph.

Given `EllVar` variables ``EVVar`` attached to an `EllImg` environment,
returns `EllVar` objects for the dependents ``vDep`` (in matching
order), describing an ellipsoidal enclosure of the image of the
dependents.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<EV> const& EVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<EV> EVDep(nDep);
            G.eval(vDep, EVDep, vVar, EVVar);
            return EVDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("EVVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in ellipsoidal image arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FD> const& FDVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FD> FDDep(nDep);
            G.eval(SgDep, vDep, FDDep, vVar, FDVar);
            return FDDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("FDVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dependency (sparsity pattern) arithmetic, using a
precomputed subgraph.

Given `FFDep` variables ``FDVar``, returns `FFDep` objects for the
dependents ``vDep`` (in matching order), describing which variables each
dependent depends on and how (linearly, polynomially, ...).
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FD> const& FDVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FD> FDDep(nDep);
            G.eval(vDep, FDDep, vVar, FDVar);
            return FDDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("FDVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dependency (sparsity pattern) arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FI> const& FIVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FI> FIDep(nDep);
            G.eval(SgDep, vDep, FIDep, vVar, FIVar);
            return FIDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar"), py::arg("FIVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in invariant (invertible structure detection)
arithmetic, using a precomputed subgraph.

Given `FFInv` variables ``FIVar``, returns `FFInv` objects for the
dependents ``vDep`` (in matching order), detecting invertible structure
in the factorable expressions.
)doc")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FI> const& FIVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FI> FIDep(nDep);
            G.eval(vDep, FIDep, vVar, FIVar);
            return FIDep;
          },
          py::arg("vDep"), py::arg("vVar"), py::arg("FIVar"),
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in invariant (invertible structure detection)
arithmetic.

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "reval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep, std::vector<I>& IDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I>& IVar,
             I const& InfVal, unsigned const MaxPass, double const& ThresPass)
          {
            G.reval(SgDep, vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass);
            return std::make_pair(IVar, IDep);
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("IDep"), py::arg("vVar"),
          py::arg("IVar"), py::arg("InfVal"), py::arg("MaxPass") = 5,
          py::arg("ThresPass") = 0e0,
          // py::return_value_policy::take_ownership,
          R"doc(
Tighten variable and dependent bounds by forward/backward interval
propagation (constraint propagation), using a precomputed subgraph.

Performs up to ``MaxPass`` combined forward/backward sweeps through the
subgraph: forward interval evaluation of the dependents, intersection
with the a priori dependent enclosures ``IDep``, then reverse
propagation through the operations to contract the variable enclosures
``IVar``.

Parameters
----------
SgDep : FFSubgraph
    Subgraph of the dependents, as returned by `subgraph`; passing it
    avoids re-extracting it on every call.
vDep : list of FFVar
    Dependent nodes (e.g. constraint expressions).
IDep : list of Interval
    A priori enclosures (constraints) on the dependents, in ``vDep``
    order.
vVar : list of FFVar
    Variable nodes.
IVar : list of Interval
    Initial enclosures of the variables, in ``vVar`` order.
InfVal : Interval
    Fallback value assigned when backward propagation fails for an
    operation (e.g. unbounded inversion); typically a very large
    interval such as ``1e20 * Interval(-1, 1)``.
MaxPass : int, optional
    Maximum number of forward/backward passes. Default is 5.
ThresPass : float, optional
    Minimum relative range reduction for a further pass to be worth
    performing. Default is 0.0.

Returns
-------
IVar : list of Interval
    Tightened variable enclosures, in ``vVar`` order.
IDep : list of Interval
    Tightened dependent enclosures, in ``vDep`` order.

Notes
-----
New lists are returned; the input lists are not modified. For rigorous
results use a verified (outward-rounded) interval type; otherwise a
feasible constraint set may erroneously be contracted to empty.
)doc")
      .def(
          "reval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<I>& IDep, std::vector<mc::FFVar> const& vVar,
             std::vector<I>& IVar, I const& InfVal, unsigned const MaxPass,
             double const& ThresPass)
          {
            G.reval(vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass);
            return std::make_pair(IVar, IDep);
          },
          py::arg("vDep"), py::arg("IDep"), py::arg("vVar"), py::arg("IVar"),
          py::arg("InfVal"), py::arg("MaxPass") = 5, py::arg("ThresPass") = 0e0,
          // py::return_value_policy::take_ownership,
          R"doc(
Tighten variable and dependent bounds by forward/backward interval
propagation (constraint propagation).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.

Returns
-------
IVar : list of Interval
    Tightened variable enclosures, in ``vVar`` order.
IDep : list of Interval
    Tightened dependent enclosures, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<double>>& DVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = DVar1.size();
            std::vector<std::vector<double>> DDep(nSam,
                                                  std::vector<double>(nDep));
            G.veval(SgDep, vDep, DDep, vVar1, DVar1, vVar2, DVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return DDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("DVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("DVar2") = std::vector<double>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in double-precision arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

Performs one evaluation per entry of ``DVar1``, dispatching the
scenarios on up to ``options.MAXTHREAD`` parallel threads (0 means all
available hardware threads).

Parameters
----------
SgDep : FFSubgraph
    Subgraph of the dependents, as returned by `subgraph`; passing it
    avoids re-extracting it on every call.
vDep : list of FFVar
    Dependent nodes to evaluate.
vVar1 : list of FFVar
    Scenario-dependent variable nodes.
DVar1 : list of list of float
    One inner list of values of ``vVar1`` per scenario.
vVar2 : list of FFVar, optional
    Variable nodes shared by all scenarios. Default is [].
DVar2 : list of float, optional
    Values of the shared variables, in ``vVar2`` order. Default is [].
walltime : bool, optional
    If True, print the wall-clock evaluation time to standard error.
    Default is False.

Returns
-------
DDep : list of list of float
    One inner list per scenario holding the values of the dependents,
    in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<double>>& DVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = DVar1.size();
            std::vector<std::vector<double>> DDep(nSam,
                                                  std::vector<double>(nDep));
            G.veval(vDep, DDep, vVar1, DVar1, vVar2, DVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return DDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("DVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("DVar2") = std::vector<double>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in double-precision arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = [pymcpp.FFVar(DAG, "X" + str(i)) for i in range(2)]
>>> F = [X[0] * X[1] - 1, pymcpp.exp(X[0]) + X[1]]
>>> DAG.options.MAXTHREAD = 0  # use all available threads
>>> samX = [[1.0, 2.0], [3.0, 2.0]]  # one scenario per row
>>> DAG.veval(F, X, samX)
[[1.0, 4.718281828459045], [5.0, 22.085536923187668]]
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<I>>& IVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = IVar1.size();
            std::vector<std::vector<I>> IDep(nSam, std::vector<I>(nDep));
            G.veval(SgDep, vDep, IDep, vVar1, IVar1, vVar2, IVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return IDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("IVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("IVar2") = std::vector<I>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in interval arithmetic for multiple scenarios
(vectorized), using a precomputed subgraph.

``IVar1`` holds one list of `Interval` bounds on ``vVar1`` per
scenario; ``vVar2``/``IVar2`` are shared by all scenarios. Returns one
list of `Interval` enclosures of the dependents per scenario, in
``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<I>>& IVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = IVar1.size();
            std::vector<std::vector<I>> IDep(nSam, std::vector<I>(nDep));
            G.veval(vDep, IDep, vVar1, IVar1, vVar2, IVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return IDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("IVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("IVar2") = std::vector<I>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in interval arithmetic for multiple scenarios
(vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<MC>>& MCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = MCVar1.size();
            std::vector<std::vector<MC>> MCDep(nSam, std::vector<MC>(nDep));
            G.veval(SgDep, vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return MCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("MCVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("MCVar2") = std::vector<MC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in McCormick relaxation arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

``MCVar1`` holds one list of `McCormick` variables per scenario;
``vVar2``/``MCVar2`` are shared by all scenarios. Returns one list of
`McCormick` relaxations of the dependents per scenario, in ``vDep``
order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<MC>>& MCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = MCVar1.size();
            std::vector<std::vector<MC>> MCDep(nSam, std::vector<MC>(nDep));
            G.veval(vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return MCDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("MCVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("MCVar2") = std::vector<MC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in McCormick relaxation arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SB>>& SBVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SB>& SBVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SBVar1.size();
            std::vector<std::vector<SB>> SBDep(nSam, std::vector<SB>(nDep));
            G.veval(SgDep, vDep, SBDep, vVar1, SBVar1, vVar2, SBVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SBDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("SBVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("SBVar2") = std::vector<SB>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in spectral bound arithmetic for multiple scenarios
(vectorized), using a precomputed subgraph.

``SBVar1`` holds one list of `Specbnd` variables per scenario;
``vVar2``/``SBVar2`` are shared by all scenarios. Returns one list of
`Specbnd` results per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SB>>& SBVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SB>& SBVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SBVar1.size();
            std::vector<std::vector<SB>> SBDep(nSam, std::vector<SB>(nDep));
            G.veval(vDep, SBDep, vVar1, SBVar1, vVar2, SBVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SBDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("SBVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("SBVar2") = std::vector<SB>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in spectral bound arithmetic for multiple scenarios
(vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<T>>& TVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<T>& TVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = TVar1.size();
            std::vector<std::vector<T>> TDep(nSam, std::vector<T>(nDep));
            G.veval(SgDep, vDep, TDep, vVar1, TVar1, vVar2, TVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return TDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("TVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("TVar2") = std::vector<T>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in Taylor model arithmetic for multiple scenarios
(vectorized), using a precomputed subgraph.

``TVar1`` holds one list of `TVar` variables per scenario;
``vVar2``/``TVar2`` are shared by all scenarios. Returns one list of
`TVar` models of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<T>>& TVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<T>& TVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = TVar1.size();
            std::vector<std::vector<T>> TDep(nSam, std::vector<T>(nDep));
            G.veval(vDep, TDep, vVar1, TVar1, vVar2, TVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return TDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("TVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("TVar2") = std::vector<T>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in Taylor model arithmetic for multiple scenarios
(vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<C>>& CVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<C>& CVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = CVar1.size();
            std::vector<std::vector<C>> CDep(nSam, std::vector<C>(nDep));
            G.veval(SgDep, vDep, CDep, vVar1, CVar1, vVar2, CVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return CDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("CVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("CVar2") = std::vector<C>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dense Chebyshev model arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

``CVar1`` holds one list of `CVar` variables per scenario;
``vVar2``/``CVar2`` are shared by all scenarios. Returns one list of
`CVar` models of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<C>>& CVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<C>& CVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = CVar1.size();
            std::vector<std::vector<C>> CDep(nSam, std::vector<C>(nDep));
            G.veval(vDep, CDep, vVar1, CVar1, vVar2, CVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return CDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("CVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("CVar2") = std::vector<C>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dense Chebyshev model arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SC>>& SCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SC>& SCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SCVar1.size();
            std::vector<std::vector<SC>> SCDep(nSam, std::vector<SC>(nDep));
            G.veval(SgDep, vDep, SCDep, vVar1, SCVar1, vVar2, SCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("SCVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("SCVar2") = std::vector<SC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in sparse Chebyshev model arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

``SCVar1`` holds one list of `SCVar` variables per scenario;
``vVar2``/``SCVar2`` are shared by all scenarios. Returns one list of
`SCVar` models of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SC>>& SCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SC>& SCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SCVar1.size();
            std::vector<std::vector<SC>> SCDep(nSam, std::vector<SC>(nDep));
            G.veval(vDep, SCDep, vVar1, SCVar1, vVar2, SCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SCDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("SCVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("SCVar2") = std::vector<SC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in sparse Chebyshev model arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWLS>>& PWLSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWLS>& PWLSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWLSVar1.size();
            std::vector<std::vector<PWLS>> PWLSDep(nSam,
                                                   std::vector<PWLS>(nDep));
            G.veval(SgDep, vDep, PWLSDep, vVar1, PWLSVar1, vVar2, PWLSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWLSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PWLSVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PWLSVar2") = std::vector<PWLS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-linear univariate estimators for multiple scenarios
(vectorized), using a precomputed subgraph.

``PWLSVar1`` holds one list of `PWLSVar` variables per scenario;
``vVar2``/``PWLSVar2`` are shared by all scenarios. Returns one list of
`PWLSVar` relaxations of the dependents per scenario, in ``vDep``
order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWLS>>& PWLSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWLS>& PWLSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWLSVar1.size();
            std::vector<std::vector<PWLS>> PWLSDep(nSam,
                                                   std::vector<PWLS>(nDep));
            G.veval(vDep, PWLSDep, vVar1, PWLSVar1, vVar2, PWLSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWLSDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PWLSVar1"),
          py::arg("vVar2")    = std::vector<mc::FFVar>(),
          py::arg("PWLSVar2") = std::vector<PWLS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-linear univariate estimators for multiple scenarios
(vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWCS>>& PWCSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWCS>& PWCSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWCSVar1.size();
            std::vector<std::vector<PWCS>> PWCSDep(nSam,
                                                   std::vector<PWCS>(nDep));
            G.veval(SgDep, vDep, PWCSDep, vVar1, PWCSVar1, vVar2, PWCSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWCSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PWCSVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PWCSVar2") = std::vector<PWCS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-constant univariate estimators for multiple scenarios
(vectorized), using a precomputed subgraph.

``PWCSVar1`` holds one list of `PWCSVar` variables per scenario;
``vVar2``/``PWCSVar2`` are shared by all scenarios. Returns one list of
`PWCSVar` relaxations of the dependents per scenario, in ``vDep``
order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWCS>>& PWCSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWCS>& PWCSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWCSVar1.size();
            std::vector<std::vector<PWCS>> PWCSDep(nSam,
                                                   std::vector<PWCS>(nDep));
            G.veval(vDep, PWCSDep, vVar1, PWCSVar1, vVar2, PWCSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWCSDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PWCSVar1"),
          py::arg("vVar2")    = std::vector<mc::FFVar>(),
          py::arg("PWCSVar2") = std::vector<PWCS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in superposition model arithmetic with
piecewise-constant univariate estimators for multiple scenarios
(vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PV>>& PVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PV>& PVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PVVar1.size();
            std::vector<std::vector<PV>> PVDep(nSam, std::vector<PV>(nDep));
            G.veval(SgDep, vDep, PVDep, vVar1, PVVar1, vVar2, PVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PVVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PVVar2") = std::vector<PV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in polyhedral image arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

``PVVar1`` holds one list of `PolVar` variables per scenario;
``vVar2``/``PVVar2`` are shared by all scenarios. Returns one list of
`PolVar` results of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PV>>& PVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PV>& PVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PVVar1.size();
            std::vector<std::vector<PV>> PVDep(nSam, std::vector<PV>(nDep));
            G.veval(vDep, PVDep, vVar1, PVVar1, vVar2, PVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PVDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PVVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("PVVar2") = std::vector<PV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in polyhedral image arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<EV>>& EVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<EV>& EVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = EVVar1.size();
            std::vector<std::vector<EV>> EVDep(nSam, std::vector<EV>(nDep));
            G.veval(SgDep, vDep, EVDep, vVar1, EVVar1, vVar2, EVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return EVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("EVVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("EVVar2") = std::vector<EV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in ellipsoidal image arithmetic for multiple
scenarios (vectorized), using a precomputed subgraph.

``EVVar1`` holds one list of `EllVar` variables per scenario;
``vVar2``/``EVVar2`` are shared by all scenarios. Returns one list of
`EllVar` results of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<EV>>& EVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<EV>& EVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = EVVar1.size();
            std::vector<std::vector<EV>> EVDep(nSam, std::vector<EV>(nDep));
            G.veval(vDep, EVDep, vVar1, EVVar1, vVar2, EVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return EVDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("EVVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("EVVar2") = std::vector<EV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in ellipsoidal image arithmetic for multiple
scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FD>>& FDVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FD>& FDVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FDVar1.size();
            std::vector<std::vector<FD>> FDDep(nSam, std::vector<FD>(nDep));
            G.veval(SgDep, vDep, FDDep, vVar1, FDVar1, vVar2, FDVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FDDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("FDVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("FDVar2") = std::vector<FD>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dependency (sparsity pattern) arithmetic for
multiple scenarios (vectorized), using a precomputed subgraph.

``FDVar1`` holds one list of `FFDep` variables per scenario;
``vVar2``/``FDVar2`` are shared by all scenarios. Returns one list of
`FFDep` results of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FD>>& FDVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FD>& FDVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FDVar1.size();
            std::vector<std::vector<FD>> FDDep(nSam, std::vector<FD>(nDep));
            G.veval(vDep, FDDep, vVar1, FDVar1, vVar2, FDVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FDDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("FDVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("FDVar2") = std::vector<FD>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in dependency (sparsity pattern) arithmetic for
multiple scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FI>>& FIVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FI>& FIVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FIVar1.size();
            std::vector<std::vector<FI>> FIDep(nSam, std::vector<FI>(nDep));
            G.veval(SgDep, vDep, FIDep, vVar1, FIVar1, vVar2, FIVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FIDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("FIVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("FIVar2") = std::vector<FI>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in invariant (invertible structure detection)
arithmetic for multiple scenarios (vectorized), using a precomputed
subgraph.

``FIVar1`` holds one list of `FFInv` variables per scenario;
``vVar2``/``FIVar2`` are shared by all scenarios. Returns one list of
`FFInv` results of the dependents per scenario, in ``vDep`` order.
)doc")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FI>>& FIVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FI>& FIVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FIVar1.size();
            std::vector<std::vector<FI>> FIDep(nSam, std::vector<FI>(nDep));
            G.veval(vDep, FIDep, vVar1, FIVar1, vVar2, FIVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FIDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("FIVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("FIVar2") = std::vector<FI>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          R"doc(
Evaluate dependents in invariant (invertible structure detection)
arithmetic for multiple scenarios (vectorized).

Same as above but extracts the subgraph of ``vDep`` internally on each
call.
)doc");

  pyFFGraphOptions.def(py::init<>(), "Construct an option set with defaults.")
      .def(py::init<mc::FFGraph::Options const&>(), "Copy constructor.")
      .def_readwrite("DETECTSIGNOM", &mc::FFGraph::Options::DETECTSIGNOM,
                     "bool: Whether to detect signomial terms of the form "
                     "exp(d*log(x)) and handle them as powers x**d. Default "
                     "is True.")
      .def_readwrite("CHEBRECURS", &mc::FFGraph::Options::CHEBRECURS,
                     "bool: Whether to intersect Chebyshev variables with "
                     "their recursive expressions, building redundancy for "
                     "tighter relaxations. Default is False.")
      .def_readwrite("USEMOVE", &mc::FFGraph::Options::USEMOVE,
                     "bool: Whether to enable move semantics during DAG "
                     "evaluation. Default is False.")
      .def_readwrite("MAXTHREAD", &mc::FFGraph::Options::MAXTHREAD,
                     "int: Maximum number of threads for vectorized DAG "
                     "evaluation (`FFGraph.veval`); 0 means all concurrent "
                     "threads supported by the hardware. Default is 1.");
}
