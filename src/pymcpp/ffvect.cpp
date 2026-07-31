#include "ffvect.hpp"

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

namespace py = pybind11;

void
mc_ffvect(py::module_& m)
{
  py::class_<mc::Vect> pyVect(m, "Vect", R"doc(
Blocks of DAG expressions packaged for (parallel) vectorized evaluation.

A ``Vect`` object partitions a collection of dependent expressions of a
host :class:`FFGraph` into blocks (a list of lists of dependents), and
copies each block together with the participating variables into its own
private worker DAG. Embedded into a DAG as a single node via
:class:`FFVect`, the blocks are then evaluated independently -- on
separate threads when pymcpp is compiled with thread support -- in float,
interval, McCormick, Taylor/Chebyshev model and superposition model
arithmetics, which can substantially speed up the evaluation of large
block-structured functions.

Besides the participating variables, a list of constant variables may be
supplied; constants are passed to every block alongside the variables.

The host DAG is only referenced for bookkeeping; each block operates on
its own internal copy, so later changes to the host DAG do not affect a
``Vect`` object.

See Also
--------
FFVect : external DAG operation embedding a Vect as a single node.
)doc");

  py::class_<mc::Vect::Options> pyVectOptions(pyVect, "Options", R"doc(
Options controlling the differentiation of vectorized expressions.

Set the fields on the :attr:`Vect.options` attribute.
)doc");

  pyVect
      .def(py::init<>(), R"doc(
Construct an empty Vect; populate it with :meth:`set`.
)doc")
      .def(py::init<mc::FFGraph*, std::vector<mc::FFVar> const&,
                    std::vector<std::vector<mc::FFVar>> const&>(),
           py::arg("dag"), py::arg("var"), py::arg("fun"), R"doc(
Construct a Vect from expression blocks of a host DAG.

Each block in `fun` is copied, together with the variables in `var`, into
its own internal worker DAG.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
var : list of FFVar
    Participating (independent) variables.
fun : list of list of FFVar
    Blocks of dependent expressions; each inner list is evaluated as one
    unit, possibly on its own thread.
)doc")
      .def(py::init<mc::FFGraph*, std::vector<mc::FFVar> const&,
                    std::vector<mc::FFVar> const&,
                    std::vector<std::vector<mc::FFVar>> const&>(),
           py::arg("dag"), py::arg("var"), py::arg("cst"), py::arg("fun"),
           R"doc(
Construct a Vect with additional constant variables.

Same as the previous overload with an extra list of constant variables
that are passed to every block alongside the variables.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
var : list of FFVar
    Participating (independent) variables.
cst : list of FFVar
    Constant variables participating in the expressions.
fun : list of list of FFVar
    Blocks of dependent expressions.
)doc")
      .def(py::init<mc::Vect const&>(), R"doc(
Copy constructor; the internal worker DAGs are duplicated.
)doc")
      .def(
          "set",
          [](mc::Vect& self, mc::FFGraph* pDAG,
             std::vector<mc::FFVar> const& vVar,
             std::vector<std::vector<mc::FFVar>> const& vFun)
          { return self.set(pDAG, vVar, vFun); },
          py::arg("dag"), py::arg("var"), py::arg("fun"), R"doc(
Set the expression blocks, replacing any previous data.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
var : list of FFVar
    Participating (independent) variables.
fun : list of list of FFVar
    Blocks of dependent expressions.

Returns
-------
ok : bool
    True on success, False if copying an expression block failed.
)doc")
      .def(
          "set",
          [](mc::Vect& self, mc::FFGraph* pDAG,
             std::vector<mc::FFVar> const& vVar,
             std::vector<mc::FFVar> const& vCst,
             std::vector<std::vector<mc::FFVar>> const& vFun)
          { return self.set(pDAG, vVar, vCst, vFun); },
          py::arg("dag"), py::arg("var"), py::arg("cst"), py::arg("fun"),
          R"doc(
Set the expression blocks with additional constant variables.

Same as the first overload with an extra list of constant variables that
are passed to every block alongside the variables.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
var : list of FFVar
    Participating (independent) variables.
cst : list of FFVar
    Constant variables participating in the expressions.
fun : list of list of FFVar
    Blocks of dependent expressions.

Returns
-------
ok : bool
    True on success, False if copying an expression block failed.
)doc")
      .def_readwrite("options", &mc::Vect::options, R"doc(
Vect.Options: Options controlling differentiation of the vectorized
expressions.
)doc")
      .def_property_readonly(
          "DAG", [](mc::Vect& self) { return self.pDAG(); }, R"doc(
FFGraph: The host DAG the expression blocks were taken from.
)doc")
      .def_property_readonly(
          "Var", [](mc::Vect& self) { return self.vVar(); },
          R"doc(
list of FFVar: Participating (independent) variables in the host DAG.
)doc")
      .def_property_readonly(
          "Cst", [](mc::Vect& self) { return self.vCst(); },
          R"doc(
list of FFVar: Constant variables in the host DAG (empty if none were
supplied).
)doc")
      .def_property_readonly(
          "Dep", [](mc::Vect& self) { return self.vFun(); },
          R"doc(
list of list of FFVar: Blocks of dependent expressions in the host DAG.
)doc");


  pyVectOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::Vect::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def_readwrite(
          "AUTODIFF", &mc::Vect::Options::AUTODIFF,
          "Vect.Options.AD: Whether to apply forward (F) or backward (B) "
          "automatic differentiation when differentiating the vectorized "
          "expressions. Default is F.");

  py::enum_<mc::Vect::Options::AD>(pyVectOptions, "AD")
      .value("F", mc::Vect::Options::AD::F, "Forward differentiation")
      .value("B", mc::Vect::Options::AD::B, "Backward differentiation")
      .export_values();

  py::class_<mc::FFVect<I>, mc::FFOp> pyFFVect(m, "FFVect", R"doc(
Vectorized expression blocks as a single external DAG operation.

FFVect inserts a :class:`Vect` object into a host :class:`FFGraph` DAG as
one vector-valued node whose outputs are the concatenation of all
expression blocks. During DAG evaluation, the blocks are computed
independently -- on separate threads when pymcpp is compiled with thread
support -- in float, interval, McCormick, Taylor/Chebyshev model and
superposition model arithmetics. The node also supports symbolic
differentiation (with forward or backward AD per the Vect options) and
decomposition with :class:`SLift`.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(2)
>>> blocks = [[pymcpp.exp(X[0]) * X[1]], [X[0] - X[1], X[0] * X[1]]]
>>> fun = pymcpp.Vect(DAG, X, blocks)
>>> Op = pymcpp.FFVect()
>>> F = Op(fun)  # list of 3 dependents, one per blocked expression
>>> DAG.eval(F, X, [0., 2.])
[2.0, -2.0, 0.0]
)doc");

  pyFFVect
      .def(py::init<>(), R"doc(
Construct an unattached vectorized operation.

Call the resulting object with a :class:`Vect` to insert the operation
into the DAG referenced by that object.
)doc")
      .def(
          "__call__",
          [](mc::FFVect<I>& self, mc::Vect* pFun)
          {
            auto pDep = self(pFun);
            return std::vector<mc::FFVar*>(pDep, pDep + pFun->nFun());
          },
          py::return_value_policy::reference_internal, py::arg("fun"), R"doc(
Insert the vectorized operation into the DAG and return all outputs.

A single external operation node with the participating variables (and
constants, if any) of `fun` as operands is added to their host DAG.

Parameters
----------
fun : Vect
    Vectorized expression blocks to embed.

Returns
-------
dep : list of FFVar
    Dependent DAG variables of all blocked expressions, concatenated in
    block order.
)doc")
      .def_readwrite("type", &mc::FFVect<I>::type, R"doc(
int: Operation type identifier in the DAG.
)doc")
      .def_readwrite("info", &mc::FFVect<I>::info, R"doc(
int: Auxiliary operation identifier in the DAG.
)doc")
      .def_readwrite("varin", &mc::FFVect<I>::varin, R"doc(
list of FFVar: Operand DAG variables of this operation.
)doc")
      .def_readwrite("varout", &mc::FFVect<I>::varout, R"doc(
list of FFVar: Result DAG variables of this operation.
)doc")
      .def("name", &mc::FFVect<I>::name, R"doc(
Return the name of this operation as displayed in the DAG.
)doc")
      .def("__str__",
           [](mc::FFVect<I> const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFVect<I> const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           });
}
