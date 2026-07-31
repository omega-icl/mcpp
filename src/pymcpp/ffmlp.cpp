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

#include "ffmlp.hpp"

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFMLP<I>::Options::RELAX_TYPE>);

void
mc_ffmlp(py::module_& m)
{
  py::class_<mc::MLP<I>> pyMLP(m, "MLP", R"doc(
Trained multilayer perceptron (MLP) for evaluation and relaxation.

An ``MLP`` object stores the weights, biases and activation functions of a
trained feedforward neural network and builds an internal :class:`FFGraph`
expression tree of the network, in which each neuron pre-activation is a
single :class:`FFLin` node. The network can then be evaluated in any MC++
arithmetic (float, Interval, McCormick, Specbnd, sparse Chebyshev models,
superposition models, polyhedral relaxation), differentiated, and
tightened by interval constraint propagation.

The network data can be populated in two ways:

- :meth:`read_data` imports a network exported from PyTorch as a
  TorchScript file (``torch.jit.script(model).save("model.pt")``);
- :meth:`set_data`/:meth:`append_data` supply the layers explicitly, one
  weight matrix per layer where each neuron row reads
  ``[bias, w_1, ..., w_n]`` with one weight per neuron of the previous
  layer.

By itself an ``MLP`` is not part of any user DAG; embed it into a host DAG
as a single node with :class:`FFMLP`.

Examples
--------
>>> import pymcpp
>>> NN = pymcpp.MLP()
>>> NN.nin, NN.nout, NN.nhid   # empty network so far
(0, 0, 0)
>>> NN.set_data([([[0., 1., 1.], [0., 1., -1.]], pymcpp.MLP.RELU),
...              ([[0., 1., 1.]], pymcpp.MLP.LINEAR)])
>>> NN.nin, NN.nout, NN.nhid   # 2 inputs, 1 output, 1 hidden layer
(2, 1, 1)
)doc");

  py::enum_<mc::MLP<I>::ACTIV_TYPE>(pyMLP, "ACTIV_TYPE")
      .value("LINEAR", mc::MLP<I>::ACTIV_TYPE::LINEAR,
             "Linear (identity) activation function")
      .value("RELU", mc::MLP<I>::ACTIV_TYPE::RELU,
             "ReLU activation function max(x, 0)")
      .value("TANH", mc::MLP<I>::ACTIV_TYPE::TANH, "Tanh activation function")
      .value("SIGMOID", mc::MLP<I>::ACTIV_TYPE::SIGMOID,
             "Sigmoid activation function 1/(1+exp(-x))")
      .export_values();

  py::class_<mc::MLP<I>::Options> pyMLPOptions(
      pyMLP, "Options", R"doc(
Options controlling how an MLP is expressed and evaluated.

Set the fields on the :attr:`MLP.options` attribute of an existing
network, e.g. ``NN.options.RELU2ABS = True``.
)doc");  //, py::module_local() );

  pyMLP
      .def(py::init<>(), R"doc(
Construct an empty MLP; populate it with read_data, set_data or
append_data.
)doc")
      .def_readwrite("options", &mc::MLP<I>::options, R"doc(
MLP.Options: Options controlling how the network is expressed and
evaluated.
)doc")
      .def(
          "reset_data", [](mc::MLP<I>& self) { self.reset_data(); },
          R"doc(
Discard all stored layers, resetting the network to an empty state.
)doc")
      .def(
          "set_data",
          [](mc::MLP<I>& self,
             std::vector<
                 std::pair<std::vector<std::vector<double>>, int>> const& mlp)
          { self.set_data(mlp); },
          py::arg("mlp"), R"doc(
Set the complete network data, replacing any previously stored layers.

Parameters
----------
mlp : list of tuple of (list of list of float, int)
    One tuple per layer (hidden layers first, output layer last). The
    first tuple element is the layer's weight matrix with one row per
    neuron, each row reading ``[bias, w_1, ..., w_n]`` where ``n`` is the
    number of neurons in the previous layer (or network inputs for the
    first layer). The second element is the activation applied to every
    neuron of the layer (a value of :class:`MLP.ACTIV_TYPE`). Weights
    with magnitude below ``options.ZEROTOL`` are zeroed out.
)doc")
      .def(
          "append_data",
          [](mc::MLP<I>& self, std::vector<std::vector<double>> const& layer,
             int const activ, bool const reset)
          { self.append_data(layer, activ, reset); },
          py::arg("layer"), py::arg("activ") = mc::MLP<I>::LINEAR,
          py::arg("reset") = false, R"doc(
Append a multi-neuron layer to the network data.

Parameters
----------
layer : list of list of float
    Weight matrix of the layer, one row per neuron; each row reads
    ``[bias, w_1, ..., w_n]`` with one weight per neuron of the previous
    layer (or network input for the first layer).
activ : MLP.ACTIV_TYPE, optional
    Activation function applied to every neuron of the layer. Default is
    LINEAR.
reset : bool, optional
    Whether to discard all previously stored layers first. Default is
    False.
)doc")
      .def(
          "append_data",
          [](mc::MLP<I>& self, std::vector<double> const& layer,
             int const activ, bool const reset)
          { self.append_data(layer, activ, reset); },
          py::arg("layer"), py::arg("activ") = mc::MLP<I>::LINEAR,
          py::arg("reset") = false, R"doc(
Append a single-neuron layer to the network data.

Same as the multi-neuron overload for a layer with exactly one neuron.

Parameters
----------
layer : list of float
    Weight row ``[bias, w_1, ..., w_n]`` of the single neuron.
activ : MLP.ACTIV_TYPE, optional
    Activation function of the neuron. Default is LINEAR.
reset : bool, optional
    Whether to discard all previously stored layers first. Default is
    False.
)doc")
      .def(
          "read_data",
          [](mc::MLP<I>& self, std::string const& filename, bool const disp)
          {
            if (!self.read_data(filename, disp))
              throw std::runtime_error(
                  "MLP.read_data: failed to load network data (was pymcpp "
                  "compiled with Torch support?)");
          },
          py::arg("filename"), py::arg("disp") = false, R"doc(
Load the network data from a TorchScript file.

Reads the weights, biases and activation functions of a feedforward
network saved with PyTorch (``torch.jit.script(model).save(...)``); the
supported activations are ReLU, Tanh, Sigmoid and Linear.

Parameters
----------
filename : str
    Path to the TorchScript (.pt) file.
disp : bool, optional
    Whether to print information about the loaded modules and tensors.
    Default is False.

Raises
------
RuntimeError
    If the file cannot be read, or if pymcpp was compiled without
    Torch support.
)doc")
      .def_property_readonly(
          "dag", [](mc::MLP<I>& self) { return self.dag(); },
          py::return_value_policy::reference_internal, R"doc(
FFGraph: Internal DAG expressing the network, with inputs named "X" and
outputs named "Y". Built (or refreshed) on access.
)doc")
      .def_property_readonly(
          "nin", [](mc::MLP<I> const& self) { return self.nin(); },
          py::return_value_policy::reference_internal, R"doc(
int: Number of network inputs.
)doc")
      .def_property_readonly(
          "varin", [](mc::MLP<I>& self) { return self.varin(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Input variables of the internal DAG.
)doc")
      .def_property_readonly(
          "nout", [](mc::MLP<I> const& self) { return self.nout(); },
          py::return_value_policy::reference_internal, R"doc(
int: Number of network outputs.
)doc")
      .def_property_readonly(
          "varout", [](mc::MLP<I>& self) { return self.varout(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Output variables of the internal DAG.
)doc")
      .def_property_readonly(
          "nhid", [](mc::MLP<I> const& self) { return self.nhid(); },
          py::return_value_policy::reference_internal, R"doc(
int: Number of hidden layers.
)doc");


  pyMLPOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::MLP<I>::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](mc::MLP<I>::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("EVALTORCH", &mc::MLP<I>::Options::EVALTORCH,
                     "bool: Whether to evaluate the network in floating-point "
                     "arithmetic using Torch instead of the internal "
                     "expression tree (only effective when the data was "
                     "loaded from a TorchScript file). Default is False.")
      .def_readwrite("ZEROTOL", &mc::MLP<I>::Options::ZEROTOL,
                     "float: Threshold below which network weights are "
                     "treated as zero. Default is DBL_EPSILON.")
      .def_readwrite("RELU2ABS", &mc::MLP<I>::Options::RELU2ABS,
                     "bool: Whether to express ReLU activations via abs, i.e. "
                     "(x+|x|)/2 (True), or via max(x,0) (False). Default is "
                     "False.")
      .def_readwrite("SIG2EXP", &mc::MLP<I>::Options::SIG2EXP,
                     "bool: Whether to express sigmoid activations via exp, "
                     "i.e. 1/(1+exp(-x)) (True), or via tanh (False). Default "
                     "is False.")
      .def_readwrite("AUTODIFF", &mc::MLP<I>::Options::AUTODIFF,
                     "MLP.Options.AD_TYPE: Whether to apply forward (F) or "
                     "backward (B) automatic differentiation when computing "
                     "network gradients. Default is F.")
      .def_readwrite("CPMAX", &mc::MLP<I>::Options::CPMAX,
                     "int: Maximum rounds of interval constraint propagation "
                     "in reverse evaluation. Default is 1.")
      .def_readwrite(
          "CPTHRES", &mc::MLP<I>::Options::CPTHRES,
          "float: Threshold for repeating constraint propagation (minimum "
          "relative reduction in any variable). Default is 1e4*DBL_EPSILON.")
      .def_readwrite("CPINF", &mc::MLP<I>::Options::CPINF,
                     "float: Value used as infinity for unbounded variables "
                     "in constraint propagation. Default is 1e30.");

  py::enum_<mc::MLP<I>::Options::AD_TYPE>(pyMLPOptions, "AD_TYPE")
      .value("F", mc::MLP<I>::Options::AD_TYPE::F, "Forward differentiation")
      .value("B", mc::MLP<I>::Options::AD_TYPE::B, "Backward differentiation")
      .export_values();

  py::class_<mc::FFMLP<I>, mc::FFOp> pyFFMLP(m, "FFMLP", R"doc(
Multilayer perceptron embedded as a single external DAG operation.

FFMLP inserts a trained :class:`MLP` into a host :class:`FFGraph` DAG as
one vector-valued node mapping the network inputs to its outputs. The
resulting dependents can be combined with any other DAG expressions and
support the full DAG workflow: evaluation in float and validated
arithmetics (Interval, McCormick, Specbnd, sparse Chebyshev models,
superposition models), symbolic differentiation (fdiff/bdiff, which
inserts a companion gradient operation), polyhedral relaxation, and
reverse interval constraint propagation.

The relaxation strategy used when constructing polyhedral cuts for the
network node is selected through :attr:`options` (field ``RELAX``).

The MLP object passed to ``__call__`` is deep-copied into the DAG, so
later modifications of the original network do not affect the DAG node.

Examples
--------
>>> import pymcpp
>>> NN = pymcpp.MLP()
>>> NN.set_data([([[0., 1., 1.], [0., 1., -1.]], pymcpp.MLP.RELU),
...              ([[0., 1., 1.]], pymcpp.MLP.LINEAR)])
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(NN.nin)
>>> OpNN = pymcpp.FFMLP()
>>> F = OpNN(X, NN)           # list with the network outputs
>>> DAG.eval(F, X, [1.0, 0.5])
[2.0]
>>> DAG.eval(F, X, [pymcpp.Interval(-1, 1)] * 2)  # bound
[[  0.0000000000000000e+00 :  4.0000000000000000e+00 ]]
)doc");

  py::class_<mc::FFMLP<I>::Options> pyFFMLPOptions(pyFFMLP, "Options", R"doc(
Options controlling the relaxation of an FFMLP operation.

The field ``RELAX`` selects one or more relaxation strategies used when
building polyhedral cuts for the network node; the remaining fields tune
the individual strategies. Set the fields on the :attr:`FFMLP.options`
attribute, e.g. ``OpNN.options.RELAX = {OpNN.Options.INT}``.
)doc");

  pyFFMLP
      .def(py::init<>(), R"doc(
Construct an unattached MLP operation.

Call the resulting object with DAG variables and an :class:`MLP` to
insert the operation into a DAG.
)doc")
      .def(py::init<mc::FFMLP<I> const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "__call__",
          [](mc::FFMLP<I>& self, std::vector<mc::FFVar> const& vVar,
             mc::MLP<I>* pMLP)
          {
            if (!pMLP || !pMLP->nin())
              throw std::runtime_error(
                  "FFMLP.__call__: empty network (populate the MLP with "
                  "read_data, set_data or append_data first)");
            if (vVar.size() != pMLP->nin())
              throw std::invalid_argument(
                  "FFMLP.__call__: var must have length mlp.nin");
            return self(vVar, pMLP, mc::FFMLP<I>::COPY);
          },
          py::return_value_policy::reference_internal, py::arg("var"),
          py::arg("mlp"), R"doc(
Insert the MLP operation into the DAG and return all network outputs.

A single external operation node is added to the DAG that the variables
in `var` belong to; the network `mlp` is deep-copied into the DAG.

Parameters
----------
var : list of FFVar
    Network input variables; the length must equal ``mlp.nin``, and all
    variables must belong to the same DAG.
mlp : MLP
    Trained network to embed; copied into the DAG.

Returns
-------
dep : list of FFVar
    Dependent DAG variables of the ``mlp.nout`` network outputs.

Raises
------
RuntimeError
    If `mlp` is empty (no layers loaded).
ValueError
    If the length of `var` differs from ``mlp.nin``.
)doc")
      .def(
          "__call__",
          [](mc::FFMLP<I>& self, unsigned const idep,
             std::vector<mc::FFVar> const& vVar, mc::MLP<I>* pMLP)
          {
            if (!pMLP || !pMLP->nin())
              throw std::runtime_error(
                  "FFMLP.__call__: empty network (populate the MLP with "
                  "read_data, set_data or append_data first)");
            if (vVar.size() != pMLP->nin())
              throw std::invalid_argument(
                  "FFMLP.__call__: var must have length mlp.nin");
            if (idep >= pMLP->nout())
              throw std::invalid_argument(
                  "FFMLP.__call__: dep must be smaller than mlp.nout");
            return self(idep, vVar, pMLP, mc::FFMLP<I>::COPY);
          },
          py::return_value_policy::reference_internal, py::arg("dep"),
          py::arg("var"), py::arg("mlp"), R"doc(
Insert the MLP operation into the DAG and return output `dep` only.

Same as the first overload but returns the single network output with
index `dep` instead of the full list.

Parameters
----------
dep : int
    Index of the network output to return (0 <= dep < ``mlp.nout``).
var : list of FFVar
    Network input variables; the length must equal ``mlp.nin``.
mlp : MLP
    Trained network to embed; copied into the DAG.

Returns
-------
dep : FFVar
    Dependent DAG variable of network output `dep`.

Raises
------
RuntimeError
    If `mlp` is empty (no layers loaded).
ValueError
    If the length of `var` differs from ``mlp.nin``, or if `dep` is
    not smaller than ``mlp.nout``.
)doc")
      .def_property_readonly(
          "name", [](mc::FFMLP<I> const& self) { return self.name(); },
          py::return_value_policy::reference_internal, R"doc(
str: Name of this operation as displayed in the DAG.
)doc")
      .def("__str__",
           [](mc::FFMLP<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFMLP<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def_readwrite("options", &mc::FFMLP<I>::options, R"doc(
FFMLP.Options: Options controlling the relaxation of the network node.
)doc");


  py::enum_<mc::FFMLP<I>::Options::RELAX_TYPE>(pyFFMLPOptions, "RELAX_TYPE")
      .value("INT", mc::FFMLP<I>::Options::RELAX_TYPE::INT, "Interval bounds")
      .value("AUX", mc::FFMLP<I>::Options::RELAX_TYPE::AUX,
             "Auxiliary variable polyhedral relaxation")
      .value("MC", mc::FFMLP<I>::Options::RELAX_TYPE::MC,
             "McCormick relaxation with interval bounds")
      .value("SB", mc::FFMLP<I>::Options::RELAX_TYPE::SB, "Spectral bounds")
      .value("SCM", mc::FFMLP<I>::Options::RELAX_TYPE::SCM,
             "Sparse Chebyshev model relaxation")
      .value("PWCS", mc::FFMLP<I>::Options::RELAX_TYPE::PWCS,
             "Piecewise-constant superposition bounds")
      .value("PWLS", mc::FFMLP<I>::Options::RELAX_TYPE::PWLS,
             "Piecewise-linear superposition bounds")
      .export_values();

  py::class_<std::set<mc::FFMLP<I>::Options::RELAX_TYPE>>(m, "RELAX_FFMLP",
                                                          py::module_local(),
                                                          R"doc(
Set of relaxation strategies for FFMLP operations.

A mutable, set-like container of :class:`FFMLP.Options.RELAX_TYPE` values
backing the ``RELAX`` field of :class:`FFMLP.Options`. Supports ``add``,
``remove``, ``discard``, membership tests, iteration and ``len``. The
``RELAX`` field can also be assigned directly from a Python set, list or
tuple of RELAX_TYPE values.
)doc")
      .def(py::init<>(), R"doc(
Construct an empty set of relaxation strategies.
)doc")
      .def(
          "add",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k) { s.insert(k); },
          py::arg("value"), R"doc(
Add a relaxation strategy to the set; no effect if already present.
)doc")
      .def(
          "remove",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k)
          {
            if (s.find(k) == s.end()) throw py::key_error();
            s.erase(k);
          },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set; raises KeyError if absent.
)doc")
      .def(
          "discard",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k) { s.erase(k); },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set if present.
)doc")
      .def("__contains__", [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
                              mc::FFMLP<I>::Options::RELAX_TYPE const& k)
           { return s.find(k) != s.end(); })
      .def("__len__", [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s)
           { return s.size(); })
      .def(
          "__iter__", [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s)
          { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>())
      .def("__repr__",
           [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s)
           {
             std::ostringstream oss;
             oss << "{";
             bool first = true;
             for (auto const& v : s)
             {
               if (!first) oss << ", ";
               switch (v)
               {
                 case mc::FFMLP<I>::Options::RELAX_TYPE::INT:
                   oss << "INT";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::AUX:
                   oss << "AUX";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::MC:
                   oss << "MC";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::SB:
                   oss << "SB";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::SCM:
                   oss << "SCM";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::PWCS:
                   oss << "PWCS";
                   break;
                 case mc::FFMLP<I>::Options::RELAX_TYPE::PWLS:
                   oss << "PWLS";
                   break;
               }
               first = false;
             }
             oss << "}";
             return oss.str();
           });

  pyFFMLPOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::FFMLP<I>::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](mc::FFMLP<I>::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_property(
          "RELAX",
          [](mc::FFMLP<I>::Options& self)
              -> std::set<mc::FFMLP<I>::Options::RELAX_TYPE>&
          { return self.RELAX; },
          [](mc::FFMLP<I>::Options& self, py::object const& val)
          {
            // 1. Optimized path: Assignment from another RELAX_FFMLP (C++ copy)
            if (py::isinstance<std::set<mc::FFMLP<I>::Options::RELAX_TYPE>>(
                    val))
            {
              self.RELAX = val.cast<
                  std::set<mc::FFMLP<I>::Options::RELAX_TYPE> const&>();
              return;
            }
            // 2. Flexible path: Assignment from Python set/list/tuple
            else if (py::isinstance<py::iterable>(val))
            {
              self.RELAX.clear();
              for (auto item : val)
                // Cast items to the Enum type and insert
                self.RELAX.insert(
                    item.cast<mc::FFMLP<I>::Options::RELAX_TYPE>());
              return;
            }
            throw py::type_error(
                "Cannot assign to RELAX: Expected RELAX_FFMLP or iterable of "
                "RELAX_TYPE enums");
          },
          R"doc(
RELAX_FFMLP: Relaxation strategies applied when building polyhedral cuts
for the network node; any iterable of FFMLP.Options.RELAX_TYPE values may
be assigned. Default is {MC}.
)doc")
      .def_readwrite("POLDEF", &mc::FFMLP<I>::Options::POLDEF,
                     "bool: Whether the polyhedral relaxation options default "
                     "to those of the parent polyhedral relaxation "
                     "environment instead of the POLIMG field. Default is "
                     "True.")
      .def_readwrite("POLIMG", &mc::FFMLP<I>::Options::POLIMG,
                     "PolImg.Options: Options for the internal polyhedral "
                     "relaxation environment (used when POLDEF is False).")
      .def_readwrite("SBLIN", &mc::FFMLP<I>::Options::SBLIN,
                     "int: Number of linearization points of the spectral "
                     "bound relaxation. Default is 0, meaning 2*NIN+1 "
                     "points.")
      .def_readwrite(
          "SBSEED", &mc::FFMLP<I>::Options::SBSEED,
          "int: Random number generator seed for the Latin hypercube sampler "
          "placing the linearization points of the spectral bound relaxation. "
          "Default is 42.")
      .def_readwrite("SCMODORD", &mc::FFMLP<I>::Options::SCMODORD,
                     "int: Maximal degree of the sparse Chebyshev model "
                     "relaxation. Default is 3.")
      .def_readwrite(
          "SCBERNORD", &mc::FFMLP<I>::Options::SCBERNORD,
          "int: Degree of the Bernstein basis conversion in the sparse "
          "Chebyshev model relaxation. Default is 0, meaning same as "
          "SCMODORD.")
      .def_readwrite("SCMODEL", &mc::FFMLP<I>::Options::SCMODEL,
                     "SCModel.Options: Options for the sparse Chebyshev model "
                     "environment.")
      .def_readwrite("PWCDIV", &mc::FFMLP<I>::Options::PWCDIV,
                     "int: Equipartition size in the piecewise-constant "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWCREL", &mc::FFMLP<I>::Options::PWCREL,
                     "int: Representation of piecewise-constant univariates - "
                     "0: continuous relaxation; 1: binary encoding. Default "
                     "is 0.")
      .def_readwrite("PWCSLOPE", &mc::FFMLP<I>::Options::PWCSLOPE,
                     "bool: Whether to append cuts from slopes in the "
                     "piecewise-constant superposition model relaxation. "
                     "Default is True.")
      .def_readwrite(
          "PWCSHADOW", &mc::FFMLP<I>::Options::PWCSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-constant superposition model relaxation. Default is "
          "True.")
      .def_readwrite("PWCSUP", &mc::FFMLP<I>::Options::PWCSUP,
                     "PWCSModel.Options: Options for the superposition model "
                     "with piecewise-constant univariate estimators.")
      .def_readwrite("PWLINI", &mc::FFMLP<I>::Options::PWLINI,
                     "int: Initial partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWLMAX", &mc::FFMLP<I>::Options::PWLMAX,
                     "int: Maximal partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite(
          "PWLREL", &mc::FFMLP<I>::Options::PWLREL,
          "int: Representation of piecewise-linear univariates - 0: "
          "continuous relaxation; 1: binary encoding; 2: SOS2 encoding. "
          "Default is 0.")
      .def_readwrite(
          "PWLSHADOW", &mc::FFMLP<I>::Options::PWLSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-linear superposition model relaxation. Default is True.")
      .def_readwrite("PWLSUP", &mc::FFMLP<I>::Options::PWLSUP,
                     "PWLSModel.Options: Options for the superposition model "
                     "with adaptive piecewise-linear univariate estimators.");
}
