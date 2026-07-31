// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

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

#include "ffdagext.hpp"

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>);

void
mc_ffdagext(py::module_& m)
{
  py::class_<mc::DAGEXT<I>> pyDAGEXT(m, "DAGEXT", R"doc(
Subgraph of a DAG packaged for embedding as a single external operation.

A ``DAGEXT`` object copies a set of expressions -- given by input
variables `varin` and dependent expressions `varout` of a host
:class:`FFGraph` -- into a self-contained internal DAG. The resulting
object defines a vector-valued operation mapping the ``nin`` inputs to
the ``nout`` outputs, which can then be embedded as a single node into
any DAG using :class:`FFDAGEXT`. This is useful for hiding a large or
repeatedly used subexpression behind one DAG node, while retaining
evaluation, differentiation, relaxation and constraint propagation of
the encapsulated expressions in all supported arithmetics.

Examples
--------
>>> import pymcpp
>>> SUB = pymcpp.FFGraph()
>>> X = SUB.add_vars(2)
>>> F = [X[0] * pymcpp.exp(X[1]), X[0] - X[1]]
>>> fun = pymcpp.DAGEXT(SUB, X, F)   # copies the subexpressions
>>> DAG = pymcpp.FFGraph()
>>> Y = DAG.add_vars(2)
>>> Op = pymcpp.FFDAGEXT()
>>> G = Op(Y, fun)  # embed as single node; list of 2 dependents
)doc");

  py::class_<mc::DAGEXT<I>::Options> pyDAGEXTOptions(
      pyDAGEXT, "Options", R"doc(
Options controlling evaluation of the encapsulated expression tree.

Set the fields on the :attr:`DAGEXT.options` attribute of an existing
object.
)doc");  //, py::module_local() );

  pyDAGEXT
      .def(py::init<>(), R"doc(
Construct an empty DAGEXT; populate it with :meth:`set`.
)doc")
      .def(py::init<mc::FFGraph*, std::vector<mc::FFVar> const&,
                    std::vector<mc::FFVar> const&>(),
           py::arg("dag"), py::arg("varin"), py::arg("varout"), R"doc(
Construct a DAGEXT from expressions of a host DAG.

The subexpressions defining `varout` in terms of `varin` are copied into
a new internal DAG; the host DAG is not modified and is no longer
referenced afterwards.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
varin : list of FFVar
    Input variables of the operation.
varout : list of FFVar
    Dependent expressions of the operation, in terms of `varin`.
)doc")
      .def(py::init<mc::DAGEXT<I> const&>(), R"doc(
Copy constructor.
)doc")
      .def_readwrite("options", &mc::DAGEXT<I>::options, R"doc(
DAGEXT.Options: Options controlling differentiation and constraint
propagation of the encapsulated expressions.
)doc")
      .def(
          "set",
          [](mc::DAGEXT<I>& self, mc::FFGraph* dag,
             std::vector<mc::FFVar> const& varin,
             std::vector<mc::FFVar> const& varout)
          { self.set(dag, varin, varout); },
          py::arg("dag"), py::arg("varin"), py::arg("varout"), R"doc(
Set the encapsulated expressions, replacing any previous data.

The subexpressions defining `varout` in terms of `varin` are copied from
the host DAG into a new internal DAG.

Parameters
----------
dag : FFGraph
    Host DAG holding the expressions.
varin : list of FFVar
    Input variables of the operation.
varout : list of FFVar
    Dependent expressions of the operation, in terms of `varin`.
)doc")
      .def_property_readonly(
          "dag", [](mc::DAGEXT<I> const& self) { return self.dag(); },
          py::return_value_policy::reference_internal, R"doc(
FFGraph: Internal DAG holding the copied expressions.
)doc")
      .def_property_readonly(
          "nin", [](mc::DAGEXT<I> const& self) { return self.nin(); },
          py::return_value_policy::reference_internal, R"doc(
int: Number of input variables of the operation.
)doc")
      .def_property_readonly(
          "varin", [](mc::DAGEXT<I> const& self) { return self.varin(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Input variables in the internal DAG.
)doc")
      .def_property_readonly(
          "nout", [](mc::DAGEXT<I> const& self) { return self.nout(); },
          py::return_value_policy::reference_internal, R"doc(
int: Number of output (dependent) variables of the operation.
)doc")
      .def_property_readonly(
          "varout", [](mc::DAGEXT<I> const& self) { return self.varout(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Output (dependent) variables in the internal DAG.
)doc");


  pyDAGEXTOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::DAGEXT<I>::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](mc::DAGEXT<I>::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("AUTODIFF", &mc::DAGEXT<I>::Options::AUTODIFF,
                     "DAGEXT.Options.AD_TYPE: Whether to apply forward (F) or "
                     "backward (B) automatic differentiation when computing "
                     "derivatives of the encapsulated expressions. Default is "
                     "F.")
      .def_readwrite("CPMAX", &mc::DAGEXT<I>::Options::CPMAX,
                     "int: Maximum rounds of interval constraint propagation "
                     "in reverse evaluation. Default is 1.")
      .def_readwrite(
          "CPTHRES", &mc::DAGEXT<I>::Options::CPTHRES,
          "float: Threshold for repeating constraint propagation (minimum "
          "relative reduction in any variable). Default is 1e4*DBL_EPSILON.")
      .def_readwrite("CPINF", &mc::DAGEXT<I>::Options::CPINF,
                     "float: Value used as infinity for unbounded variables "
                     "in constraint propagation. Default is 1e30.");

  py::enum_<mc::DAGEXT<I>::Options::AD_TYPE>(pyDAGEXTOptions, "AD_TYPE")
      .value("F", mc::DAGEXT<I>::Options::AD_TYPE::F, "Forward differentiation")
      .value("B", mc::DAGEXT<I>::Options::AD_TYPE::B,
             "Backward differentiation")
      .export_values();

  py::class_<mc::FFDAGEXT<I>, mc::FFOp> pyFFDAGEXT(m, "FFDAGEXT", R"doc(
Encapsulated expression tree as a single external DAG operation.

FFDAGEXT inserts a :class:`DAGEXT` object -- a self-contained subgraph
mapping ``nin`` inputs to ``nout`` outputs -- into a host
:class:`FFGraph` DAG as one vector-valued node. The resulting dependents
support the full DAG workflow: evaluation in float and validated
arithmetics (Interval, McCormick, Specbnd, sparse Chebyshev models,
superposition models), symbolic differentiation, polyhedral relaxation,
and reverse interval constraint propagation. The relaxation strategy for
the node is selected through :attr:`options` (field ``RELAX``).

The DAGEXT object passed to ``__call__`` is deep-copied into the DAG, so
later modifications of the original object do not affect the DAG node.

See Also
--------
DAGEXT : the encapsulated expression tree inserted by this operation.
)doc");

  py::class_<mc::FFDAGEXT<I>::Options> pyFFDAGEXTOptions(pyFFDAGEXT, "Options",
                                                         R"doc(
Options controlling the relaxation of an FFDAGEXT operation.

The field ``RELAX`` selects one or more relaxation strategies used when
building polyhedral cuts for the node; the remaining fields tune the
individual strategies. Set the fields on the :attr:`FFDAGEXT.options`
attribute.
)doc");

  pyFFDAGEXT
      .def(py::init<bool const>(), py::arg("sparse") = true, R"doc(
Construct an unattached DAGEXT operation.

Parameters
----------
sparse : bool, optional
    Whether to exploit derivative sparsity when differentiating the
    encapsulated expressions (only the structurally nonzero Jacobian
    entries are formed). Default is True.
)doc")
      .def(py::init<mc::FFDAGEXT<I> const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "__call__",
          [](mc::FFDAGEXT<I>& self, std::vector<mc::FFVar> const& vVar,
             mc::DAGEXT<I>* pDAGEXT)
          { return self(vVar, pDAGEXT, mc::FFDAGEXT<I>::COPY); },
          py::return_value_policy::reference_internal, py::arg("var"),
          py::arg("dag"), R"doc(
Insert the encapsulated operation into the DAG and return all outputs.

A single external operation node is added to the DAG that the variables
in `var` belong to; the DAGEXT object is deep-copied into the DAG.

Parameters
----------
var : list of FFVar
    Input variables of the operation; the length must equal ``dag.nin``,
    and all variables must belong to the same host DAG.
dag : DAGEXT
    Encapsulated expression tree to embed; copied into the DAG.

Returns
-------
dep : list of FFVar
    Dependent DAG variables of the ``dag.nout`` operation outputs.
)doc")
      .def(
          "__call__",
          [](mc::FFDAGEXT<I>& self, unsigned const idep,
             std::vector<mc::FFVar> const& vVar, mc::DAGEXT<I>* pDAGEXT)
          { return self(idep, vVar, pDAGEXT, mc::FFDAGEXT<I>::COPY); },
          py::return_value_policy::reference_internal, py::arg("dep"),
          py::arg("var"), py::arg("dag"), R"doc(
Insert the encapsulated operation into the DAG and return output `dep`.

Same as the first overload but returns the single output with index
`dep` instead of the full list.

Parameters
----------
dep : int
    Index of the output to return (0 <= dep < ``dag.nout``).
var : list of FFVar
    Input variables of the operation; the length must equal ``dag.nin``.
dag : DAGEXT
    Encapsulated expression tree to embed; copied into the DAG.

Returns
-------
dep : FFVar
    Dependent DAG variable of output `dep`.
)doc")
      .def_property_readonly(
          "name", [](mc::FFDAGEXT<I> const& self) { return self.name(); },
          py::return_value_policy::reference_internal, R"doc(
str: Name of this operation as displayed in the DAG.
)doc")
      .def("__str__",
           [](mc::FFDAGEXT<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFDAGEXT<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def_readwrite("options", &mc::FFDAGEXT<I>::options, R"doc(
FFDAGEXT.Options: Options controlling the relaxation of the node.
)doc");


  py::enum_<mc::FFDAGEXT<I>::Options::RELAX_TYPE>(pyFFDAGEXTOptions,
                                                  "RELAX_TYPE")
      .value("INT", mc::FFDAGEXT<I>::Options::RELAX_TYPE::INT,
             "Interval bounds")
      .value("AUX", mc::FFDAGEXT<I>::Options::RELAX_TYPE::AUX,
             "Auxiliary variable polyhedral relaxation")
      .value("MC", mc::FFDAGEXT<I>::Options::RELAX_TYPE::MC,
             "McCormick relaxation with interval bounds")
      .value("SB", mc::FFDAGEXT<I>::Options::RELAX_TYPE::SB, "Spectral bounds")
      .value("SCM", mc::FFDAGEXT<I>::Options::RELAX_TYPE::SCM,
             "Sparse Chebyshev model relaxation")
      .value("PWCS", mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWCS,
             "Piecewise-constant superposition bounds")
      .value("PWLS", mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWLS,
             "Piecewise-linear superposition bounds")
      .export_values();

  py::class_<std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>>(
      m, "RELAX_FFDAGEXT", py::module_local(), R"doc(
Set of relaxation strategies for FFDAGEXT operations.

A mutable, set-like container of :class:`FFDAGEXT.Options.RELAX_TYPE`
values backing the ``RELAX`` field of :class:`FFDAGEXT.Options`. Supports
``add``, ``remove``, ``discard``, membership tests, iteration and
``len``. The ``RELAX`` field can also be assigned directly from a Python
set, list or tuple of RELAX_TYPE values.
)doc")
      .def(py::init<>(), R"doc(
Construct an empty set of relaxation strategies.
)doc")
      .def(
          "add",
          [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s,
             mc::FFDAGEXT<I>::Options::RELAX_TYPE const& k) { s.insert(k); },
          py::arg("value"), R"doc(
Add a relaxation strategy to the set; no effect if already present.
)doc")
      .def(
          "remove",
          [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s,
             mc::FFDAGEXT<I>::Options::RELAX_TYPE const& k)
          {
            if (s.find(k) == s.end()) throw py::key_error();
            s.erase(k);
          },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set; raises KeyError if absent.
)doc")
      .def(
          "discard",
          [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s,
             mc::FFDAGEXT<I>::Options::RELAX_TYPE const& k) { s.erase(k); },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set if present.
)doc")
      .def("__contains__", [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s,
                              mc::FFDAGEXT<I>::Options::RELAX_TYPE const& k)
           { return s.find(k) != s.end(); })
      .def("__len__", [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s)
           { return s.size(); })
      .def(
          "__iter__", [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s)
          { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>())
      .def("__repr__",
           [](std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>& s)
           {
             std::ostringstream oss;
             oss << "{";
             bool first = true;
             for (auto const& v : s)
             {
               if (!first) oss << ", ";
               switch (v)
               {
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::INT:
                   oss << "INT";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::AUX:
                   oss << "AUX";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::MC:
                   oss << "MC";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::SB:
                   oss << "SB";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::SCM:
                   oss << "SCM";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWCS:
                   oss << "PWCS";
                   break;
                 case mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWLS:
                   oss << "PWLS";
                   break;
               }
               first = false;
             }
             oss << "}";
             return oss.str();
           });

  pyFFDAGEXTOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::FFDAGEXT<I>::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](mc::FFDAGEXT<I>::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_property(
          "RELAX",
          [](mc::FFDAGEXT<I>::Options& self)
              -> std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>&
          { return self.RELAX; },
          [](mc::FFDAGEXT<I>::Options& self, py::object const& val)
          {
            // 1. Optimized path: Assignment from another RELAX_FFDAGEXT (C++
            // copy)
            if (py::isinstance<std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>>(
                    val))
            {
              self.RELAX = val.cast<
                  std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> const&>();
              return;
            }
            // 2. Flexible path: Assignment from Python set/list/tuple
            else if (py::isinstance<py::iterable>(val))
            {
              self.RELAX.clear();
              for (auto item : val)
                // Cast items to the Enum type and insert
                self.RELAX.insert(
                    item.cast<mc::FFDAGEXT<I>::Options::RELAX_TYPE>());
              return;
            }
            throw py::type_error(
                "Cannot assign to RELAX: Expected RELAX_FFDAGEXT or iterable "
                "of RELAX_TYPE enums");
          },
          R"doc(
RELAX_FFDAGEXT: Relaxation strategies applied when building polyhedral
cuts for the node; any iterable of FFDAGEXT.Options.RELAX_TYPE values may
be assigned. Default is {MC}.
)doc")
      .def_readwrite("POLDEF", &mc::FFDAGEXT<I>::Options::POLDEF,
                     "bool: Whether the polyhedral relaxation options default "
                     "to those of the parent polyhedral relaxation "
                     "environment instead of the POLIMG field. Default is "
                     "True.")
      .def_readwrite("POLIMG", &mc::FFDAGEXT<I>::Options::POLIMG,
                     "PolImg.Options: Options for the internal polyhedral "
                     "relaxation environment (used when POLDEF is False).")
      .def_readwrite("SBLIN", &mc::FFDAGEXT<I>::Options::SBLIN,
                     "int: Number of linearization points of the spectral "
                     "bound relaxation. Default is 0, meaning 2*NIN+1 "
                     "points.")
      .def_readwrite(
          "SBSEED", &mc::FFDAGEXT<I>::Options::SBSEED,
          "int: Random number generator seed for the Latin hypercube sampler "
          "placing the linearization points of the spectral bound relaxation. "
          "Default is 42.")
      .def_readwrite("SCMODORD", &mc::FFDAGEXT<I>::Options::SCMODORD,
                     "int: Maximal degree of the sparse Chebyshev model "
                     "relaxation. Default is 3.")
      .def_readwrite(
          "SCBERNORD", &mc::FFDAGEXT<I>::Options::SCBERNORD,
          "int: Degree of the Bernstein basis conversion in the sparse "
          "Chebyshev model relaxation. Default is 0, meaning same as "
          "SCMODORD.")
      .def_readwrite("SCMODEL", &mc::FFDAGEXT<I>::Options::SCMODEL,
                     "SCModel.Options: Options for the sparse Chebyshev model "
                     "environment.")
      .def_readwrite("PWCDIV", &mc::FFDAGEXT<I>::Options::PWCDIV,
                     "int: Equipartition size in the piecewise-constant "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWCREL", &mc::FFDAGEXT<I>::Options::PWCREL,
                     "int: Representation of piecewise-constant univariates - "
                     "0: continuous relaxation; 1: binary encoding. Default "
                     "is 0.")
      .def_readwrite("PWCSLOPE", &mc::FFDAGEXT<I>::Options::PWCSLOPE,
                     "bool: Whether to append cuts from slopes in the "
                     "piecewise-constant superposition model relaxation. "
                     "Default is True.")
      .def_readwrite(
          "PWCSHADOW", &mc::FFDAGEXT<I>::Options::PWCSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-constant superposition model relaxation. Default is "
          "True.")
      .def_readwrite("PWCSUP", &mc::FFDAGEXT<I>::Options::PWCSUP,
                     "PWCSModel.Options: Options for the superposition model "
                     "with piecewise-constant univariate estimators.")
      .def_readwrite("PWLINI", &mc::FFDAGEXT<I>::Options::PWLINI,
                     "int: Initial partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWLMAX", &mc::FFDAGEXT<I>::Options::PWLMAX,
                     "int: Maximal partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite(
          "PWLREL", &mc::FFDAGEXT<I>::Options::PWLREL,
          "int: Representation of piecewise-linear univariates - 0: "
          "continuous relaxation; 1: binary encoding; 2: SOS2 encoding. "
          "Default is 0.")
      .def_readwrite(
          "PWLSHADOW", &mc::FFDAGEXT<I>::Options::PWLSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-linear superposition model relaxation. Default is True.")
      .def_readwrite("PWLSUP", &mc::FFDAGEXT<I>::Options::PWLSUP,
                     "PWLSModel.Options: Options for the superposition model "
                     "with adaptive piecewise-linear univariate estimators.");
}
