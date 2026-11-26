// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

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

#include "ffkan.hpp"

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFKAN<I>::Options::RELAX_TYPE>);

void
mc_ffkan(py::module_& m)
{
  typedef mc::KAN<I> KAN;
  typedef mc::FFPWLU<I>::Options ENC;

  py::class_<KAN, std::shared_ptr<KAN>> pyKAN(m, "KAN", R"doc(
Kolmogorov-Arnold network (KAN) for evaluation and relaxation.

A ``KAN`` object stores the layers of a Kolmogorov-Arnold network: each
output of a layer is the sum over its inputs of an edge function of one
variable, ``y_j = sum_i f_ji(u_i)``. The edges are piecewise linear
(:class:`PWLU`), extended beyond their end breakpoints by their end
segments. The network builds an internal :class:`FFGraph` expression tree,
with one piecewise-linear operation per input of a layer, for the edges
leaving it, and one linear operation per output.

The network data can be populated in two ways:

- :meth:`read_data` imports a network saved as a dictionary by
  ``torch.save`` (model of type ``"KANmodel"``), with no other library;
- :meth:`set_data`/:meth:`append_data` supply the layers explicitly, each
  as its edges: ``layer[j][i]`` is the :class:`PWLU` from input ``i`` to
  output ``j``.

By itself a ``KAN`` is not part of any user DAG; embed it into a host DAG
as a single node with :class:`FFKAN`.

Examples
--------
>>> import pymcpp
>>> NN = pymcpp.KAN()
>>> NN.append_data([[pymcpp.PWLU([-1., 0., 1.], [1., 0., 1.], -1., 1.),
...                  pymcpp.PWLU.relu(0., 0., [1.], [1.], [0.])]])
>>> NN.nin, NN.nout, NN.nlay
(2, 1, 1)
)doc");

  py::class_<KAN::Options> pyKANOptions(pyKAN, "Options", R"doc(
Options controlling how a KAN is relaxed and differentiated.

Set the fields on the :attr:`KAN.options` attribute of an existing
network, e.g. ``NN.options.ENCODING = NN.Options.SOS2``.
)doc");

  pyKAN
      .def(py::init<>(), R"doc(
Construct an empty KAN; populate it with read_data, set_data or
append_data.
)doc")
      .def_readwrite("options", &KAN::options, R"doc(
KAN.Options: Options of the network.
)doc")
      .def(
          "reset_data", [](KAN& self) { self.reset_data(); }, R"doc(
Discard all stored layers and the internal DAG.
)doc")
      .def(
          "set_data",
          [](KAN& self, std::vector<std::vector<std::vector<mc::PWLU>>> const& layers)
          {
            std::vector<KAN::Layer> data(layers.size());
            for (size_t l = 0; l < layers.size(); ++l) data[l].pwl = layers[l];
            if (!self.set_data(data))
              throw std::invalid_argument(
                  "KAN.set_data: inconsistent layers (each layer is a list, "
                  "one entry per output, of the lists of its edges, one per "
                  "input; the number of inputs of a layer is the number of "
                  "outputs of the previous one)");
          },
          py::arg("layers"), R"doc(
Set all the layers of the network at once.

Parameters
----------
layers : list of list of list of PWLU
    ``layers[l][j][i]`` is the edge of layer ``l`` from its input ``i`` to
    its output ``j``.

Raises
------
ValueError
    If the layers are inconsistent. The network is then left as it was.
)doc")
      .def(
          "append_data",
          [](KAN& self, std::vector<std::vector<mc::PWLU>> const& layer, bool const reset)
          {
            if (!self.append_data(layer, reset))
              throw std::invalid_argument(
                  "KAN.append_data: inconsistent layer (a list, one entry "
                  "per output, of the lists of its edges, one per input; the "
                  "number of inputs is the number of outputs of the previous "
                  "layer)");
          },
          py::arg("layer"), py::arg("reset") = false, R"doc(
Append a layer to the network.

Parameters
----------
layer : list of list of PWLU
    ``layer[j][i]`` is the edge from input ``i`` to output ``j`` of the
    layer.
reset : bool, optional
    Whether to discard all previously stored layers first. Default is
    False.

Raises
------
ValueError
    If the layer is inconsistent. The network is then left as it was.
)doc")
      .def(
          "read_data",
          [](KAN& self, std::string const& filename, bool const disp)
          {
            if (!self.read_data(filename, disp))
              throw std::runtime_error(
                  "KAN.read_data: failed to load network data from " + filename +
                  " (call with disp=True for the reason)");
          },
          py::arg("filename"), py::arg("disp") = false, R"doc(
Load the network data from a dictionary saved by torch.save.

Reads a dictionary with entries ``model_type`` (``"KANmodel"``),
``config`` (``widths``, ``intervals``), ``state_dict`` -- for each layer
``l``, the tensors ``layers.l.edge_bias`` (c), ``linear_weight`` (v),
``relu_output_weight`` (a), ``relu_weight`` (w) and ``relu_bias`` (b) of
the edges ``c + v u + sum_r a_r max(w_r u + b_r, 0)`` -- and
``normalization``, with no other library. The min-max scalings of the
file are folded into the first and last layers: the network maps the raw
inputs to the raw outputs.

Parameters
----------
filename : str
    Path to the file.
disp : bool, optional
    Whether to print information about the network read, or the reason of
    a failure. Default is False.

Raises
------
RuntimeError
    If the file cannot be read or holds anything else than such a
    network. The KAN is then left as it was.
)doc")
      .def_property_readonly(
          "data",
          [](KAN const& self)
          {
            std::vector<std::vector<std::vector<mc::PWLU>>> layers;
            for (auto const& L : self.data()) layers.push_back(L.pwl);
            return layers;
          },
          R"doc(
list of list of list of PWLU: Copy of the layers, ``data[l][j][i]`` being
the edge of layer ``l`` from its input ``i`` to its output ``j``.
)doc")
      .def_property_readonly(
          "dag", [](KAN& self) { return self.dag(); },
          py::return_value_policy::reference_internal, R"doc(
FFGraph: Internal DAG expressing the network, with inputs named "X" and
outputs named "Y". Built (or refreshed) on access.
)doc")
      .def_property_readonly(
          "nin", [](KAN const& self) { return self.nin(); }, R"doc(
int: Number of inputs.
)doc")
      .def_property_readonly(
          "varin", [](KAN& self) { return self.varin(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Inputs of the internal DAG.
)doc")
      .def_property_readonly(
          "nout", [](KAN const& self) { return self.nout(); }, R"doc(
int: Number of outputs.
)doc")
      .def_property_readonly(
          "varout", [](KAN& self) { return self.varout(); },
          py::return_value_policy::reference_internal, R"doc(
list of FFVar: Outputs of the internal DAG.
)doc")
      .def_property_readonly(
          "nlay", [](KAN const& self) { return self.nlay(); }, R"doc(
int: Number of layers.
)doc")
      .def_property_readonly(
          "nwei", [](KAN const& self) { return self.nwei(); }, R"doc(
int: Total number of segments of the edges.
)doc");

  pyKANOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<KAN::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](KAN::Options& self) { self.reset(); }, R"doc(
Reset all options to their default values.
)doc")
      .def_property(
          "ENCODING", [](KAN::Options const& self) { return self.ENCODING; },
          [](KAN::Options& self, int const enc) { self.ENCODING = enc; },
          "KAN.Options.ENCODING_TYPE: Encoding of the edges in a polyhedral "
          "relaxation through the expression tree of the network: DEFAULT "
          "(as set by the BREAKPOINT_TYPE option of the polyhedral image), "
          "HULL (convex hull of each edge by cuts, no auxiliary variables), "
          "or a subdivision at the breakpoints with continuous (CONT), "
          "binary (BIN) or SOS2 (SOS2) auxiliary variables. Default is "
          "DEFAULT.")
      .def_readwrite("SEPARATE", &KAN::Options::SEPARATE,
                     "bool: Whether each edge leaving an input is encoded on "
                     "its own subdivision, rather than all on one, at the "
                     "union of their breakpoints. Default is False.")
      .def_property(
          "AUTODIFF", [](KAN::Options const& self) { return self.AUTODIFF; },
          [](KAN::Options& self, int const ad) { self.AUTODIFF = ad; },
          "KAN.Options.AD_TYPE: Whether the Jacobian is accumulated forward "
          "(F) or backward (B) through the layers. Default is F.")
      .def_readwrite("TANGENT", &KAN::Options::TANGENT,
                     "bool: Whether the forward symbolic differentiation of "
                     "an FFKAN operation propagates tangents through the "
                     "network, in place of its whole Jacobian. Default is "
                     "False.")
      .def_readwrite("ADJOINT", &KAN::Options::ADJOINT,
                     "bool: Whether the backward symbolic differentiation "
                     "of an FFKAN operation propagates adjoints through the "
                     "network, in place of its whole Jacobian. Default is "
                     "False.")
      .def_readwrite("CPMAX", &KAN::Options::CPMAX,
                     "int: Maximum number of passes in constraint "
                     "propagation. Default is 1.")
      .def_readwrite("CPTHRES", &KAN::Options::CPTHRES,
                     "float: Threshold of relative range reduction for "
                     "another pass in constraint propagation.")
      .def_readwrite("CPINF", &KAN::Options::CPINF,
                     "float: Bound standing for infinity in constraint "
                     "propagation. Default is 1e30.");

  py::enum_<ENC::ENCODING_TYPE>(pyKANOptions, "ENCODING_TYPE")
      .value("DEFAULT", ENC::DEFAULT, "As set by the polyhedral image")
      .value("HULL", ENC::HULL, "Convex hull of each edge by cuts")
      .value("CONT", ENC::CONT, "Subdivision with continuous variables")
      .value("BIN", ENC::BIN, "Subdivision with binary variables")
      .value("SOS2", ENC::SOS2, "Subdivision with SOS2 variables")
      .export_values();

  py::enum_<KAN::Options::AD_TYPE>(pyKANOptions, "AD_TYPE")
      .value("F", KAN::Options::AD_TYPE::F, "Forward differentiation")
      .value("B", KAN::Options::AD_TYPE::B, "Backward differentiation")
      .export_values();

  py::class_<mc::FFKAN<I>, mc::FFOp> pyFFKAN(m, "FFKAN", R"doc(
Kolmogorov-Arnold network embedded as a single external DAG operation.

FFKAN inserts a trained :class:`KAN` into a host :class:`FFGraph` DAG as
one vector-valued node mapping the network inputs to its outputs. The
resulting dependents can be combined with any other DAG expressions and
support the full DAG workflow: evaluation in float and validated
arithmetics (Interval, McCormick, Specbnd, sparse Chebyshev models,
superposition models), symbolic differentiation (fdiff/bdiff, which
inserts a companion gradient operation), polyhedral relaxation, and
reverse interval constraint propagation.

The relaxation strategy used when constructing polyhedral cuts for the
network node is selected through :attr:`options` (field ``RELAX``).

The KAN object passed to ``__call__`` is copied into the DAG, once per
DAG: every operation on that network in that DAG -- on whichever
variables, and the gradient operations that differentiation inserts --
shares the copy, which :attr:`KAN` returns. Later modifications of the
original network do not affect the DAG; a network modified with
``set_data`` or ``append_data`` and passed again is inserted as another
operation, with a copy of its own.

Options can be changed after the insertion, and the operations follow:
the relaxation options through :attr:`options`, the options of the
network in the DAG through ``KAN.options``.

Examples
--------
>>> import pymcpp
>>> NN = pymcpp.KAN()
>>> NN.append_data([[pymcpp.PWLU([-1., 0., 1.], [1., 0., 1.], -1., 1.),
...                  pymcpp.PWLU.relu(0., 0., [1.], [1.], [0.])]])
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(NN.nin)
>>> OpNN = pymcpp.FFKAN()
>>> F = OpNN(X, NN)           # list with the network outputs
>>> DAG.eval(F, X, [1.0, 0.5])
[1.5]
>>> DAG.eval(F, X, [pymcpp.Interval(-1, 1)] * 2)  # bound
[[  0.0000000000000000e+00 :  2.0000000000000000e+00 ]]
)doc");

  py::class_<mc::FFKAN<I>::Options, std::shared_ptr<mc::FFKAN<I>::Options>>
      pyFFKANOptions(pyFFKAN, "Options", R"doc(
Options controlling the relaxation of an FFKAN operation.

The field ``RELAX`` selects one or more relaxation strategies used when
building polyhedral cuts for the network node; the remaining fields tune
the individual strategies. Set the fields on the :attr:`FFKAN.options`
attribute, e.g. ``OpNN.options.RELAX = {OpNN.Options.INT}``, before or
after the operation is inserted in a DAG.
)doc");

  pyFFKAN
      .def(py::init<>(), R"doc(
Construct an unattached KAN operation.

Call the resulting object with DAG variables and an :class:`KAN` to
insert the operation into a DAG.
)doc")
      .def(py::init<mc::FFKAN<I> const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "__call__",
          [](mc::FFKAN<I>& self, std::vector<mc::FFVar> const& vVar,
             mc::KAN<I>* pKAN)
          {
            if (!pKAN || !pKAN->nin())
              throw std::runtime_error(
                  "FFKAN.__call__: empty network (populate the KAN with "
                  "read_data, set_data or append_data first)");
            if (vVar.size() != pKAN->nin())
              throw std::invalid_argument(
                  "FFKAN.__call__: var must have length kan.nin");
            return self(vVar, pKAN, mc::FFKAN<I>::COPY);
          },
          py::return_value_policy::reference_internal, py::arg("var"),
          py::arg("kan"), R"doc(
Insert the KAN operation into the DAG and return all network outputs.

A single external operation node is added to the DAG that the variables
in `var` belong to; the network `kan` is deep-copied into the DAG.

Parameters
----------
var : list of FFVar
    Network input variables; the length must equal ``kan.nin``, and all
    variables must belong to the same DAG.
kan : KAN
    Trained network to embed; copied into the DAG.

Returns
-------
dep : list of FFVar
    Dependent DAG variables of the ``kan.nout`` network outputs.

Raises
------
RuntimeError
    If `kan` is empty (no layers loaded).
ValueError
    If the length of `var` differs from ``kan.nin``.
)doc")
      .def(
          "__call__",
          [](mc::FFKAN<I>& self, unsigned const idep,
             std::vector<mc::FFVar> const& vVar, mc::KAN<I>* pKAN)
          {
            if (!pKAN || !pKAN->nin())
              throw std::runtime_error(
                  "FFKAN.__call__: empty network (populate the KAN with "
                  "read_data, set_data or append_data first)");
            if (vVar.size() != pKAN->nin())
              throw std::invalid_argument(
                  "FFKAN.__call__: var must have length kan.nin");
            if (idep >= pKAN->nout())
              throw std::invalid_argument(
                  "FFKAN.__call__: dep must be smaller than kan.nout");
            return self(idep, vVar, pKAN, mc::FFKAN<I>::COPY);
          },
          py::return_value_policy::reference_internal, py::arg("dep"),
          py::arg("var"), py::arg("kan"), R"doc(
Insert the KAN operation into the DAG and return output `dep` only.

Same as the first overload but returns the single network output with
index `dep` instead of the full list.

Parameters
----------
dep : int
    Index of the network output to return (0 <= dep < ``kan.nout``).
var : list of FFVar
    Network input variables; the length must equal ``kan.nin``.
kan : KAN
    Trained network to embed; copied into the DAG.

Returns
-------
dep : FFVar
    Dependent DAG variable of network output `dep`.

Raises
------
RuntimeError
    If `kan` is empty (no layers loaded).
ValueError
    If the length of `var` differs from ``kan.nin``, or if `dep` is
    not smaller than ``kan.nout``.
)doc")
      .def_property_readonly(
          "name", [](mc::FFKAN<I> const& self) { return self.name(); },
          py::return_value_policy::reference_internal, R"doc(
str: Name of this operation as displayed in the DAG.
)doc")
      .def("__str__",
           [](mc::FFKAN<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFKAN<I> const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def_property(
          "options",
          [](mc::FFKAN<I>& self) { return self.share_options(); },
          [](mc::FFKAN<I>& self, mc::FFKAN<I>::Options const& opt)
          { self.home_options() = opt; },
          R"doc(
FFKAN.Options: Options controlling the relaxation of the network node.

Before the first call, the options that the operations inserted through
this object will have. After a call, the options of the network inserted
last, in its DAG: they are shared by every operation on that network in
that DAG, the gradient operations included, and a change takes effect at
the next relaxation. (A network that another FFKAN object inserted first
in the DAG keeps the options of that object; they are the ones returned
here.)
)doc")
      .def_property_readonly(
          "map", [](mc::FFKAN<I>& self) { return self.share_object(); },
          R"doc(
KAN or None: The network of the operation in the DAG, as :attr:`KAN` --
the name that every operation on a network has for it.
)doc")
      .def_property_readonly(
          "KAN", [](mc::FFKAN<I>& self) { return self.share_object(); },
          R"doc(
KAN or None: The copy held by the DAG of the network inserted last, None
before the first call.

The copy is shared by every operation on that network in that DAG. Its
options can be changed, e.g. ``OpNN.KAN.options.TANGENT = True``: the
operations follow at their next evaluation or differentiation. Should its
layers be changed (``set_data``, ``append_data``), the operations on it
follow too, and the original network passed again is inserted as another
operation, with a copy of its own.
)doc")
      .def_property_readonly(
          "attached", [](mc::FFKAN<I> const& self) { return self.attached(); },
          R"doc(
bool: Whether the object has inserted an operation in a DAG.
)doc");


  py::enum_<mc::FFKAN<I>::Options::RELAX_TYPE>(pyFFKANOptions, "RELAX_TYPE")
      .value("INT", mc::FFKAN<I>::Options::RELAX_TYPE::INT, "Interval bounds")
      .value("AUX", mc::FFKAN<I>::Options::RELAX_TYPE::AUX,
             "Auxiliary variable polyhedral relaxation")
      .value("MC", mc::FFKAN<I>::Options::RELAX_TYPE::MC,
             "McCormick relaxation with interval bounds")
      .value("SB", mc::FFKAN<I>::Options::RELAX_TYPE::SB, "Spectral bounds")
      .value("SCM", mc::FFKAN<I>::Options::RELAX_TYPE::SCM,
             "Sparse Chebyshev model relaxation")
      .value("PWCS", mc::FFKAN<I>::Options::RELAX_TYPE::PWCS,
             "Piecewise-constant superposition bounds")
      .value("PWLS", mc::FFKAN<I>::Options::RELAX_TYPE::PWLS,
             "Piecewise-linear superposition bounds")
      .export_values();

  py::class_<std::set<mc::FFKAN<I>::Options::RELAX_TYPE>>(m, "RELAX_FFKAN",
                                                          py::module_local(),
                                                          R"doc(
Set of relaxation strategies for FFKAN operations.

A mutable, set-like container of :class:`FFKAN.Options.RELAX_TYPE` values
backing the ``RELAX`` field of :class:`FFKAN.Options`. Supports ``add``,
``remove``, ``discard``, membership tests, iteration and ``len``. The
``RELAX`` field can also be assigned directly from a Python set, list or
tuple of RELAX_TYPE values.
)doc")
      .def(py::init<>(), R"doc(
Construct an empty set of relaxation strategies.
)doc")
      .def(
          "add",
          [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s,
             mc::FFKAN<I>::Options::RELAX_TYPE const& k) { s.insert(k); },
          py::arg("value"), R"doc(
Add a relaxation strategy to the set; no effect if already present.
)doc")
      .def(
          "remove",
          [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s,
             mc::FFKAN<I>::Options::RELAX_TYPE const& k)
          {
            if (s.find(k) == s.end()) throw py::key_error();
            s.erase(k);
          },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set; raises KeyError if absent.
)doc")
      .def(
          "discard",
          [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s,
             mc::FFKAN<I>::Options::RELAX_TYPE const& k) { s.erase(k); },
          py::arg("value"), R"doc(
Remove a relaxation strategy from the set if present.
)doc")
      .def("__contains__", [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s,
                              mc::FFKAN<I>::Options::RELAX_TYPE const& k)
           { return s.find(k) != s.end(); })
      .def("__len__", [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s)
           { return s.size(); })
      .def(
          "__iter__", [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s)
          { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>())
      .def("__repr__",
           [](std::set<mc::FFKAN<I>::Options::RELAX_TYPE>& s)
           {
             std::ostringstream oss;
             oss << "{";
             bool first = true;
             for (auto const& v : s)
             {
               if (!first) oss << ", ";
               switch (v)
               {
                 case mc::FFKAN<I>::Options::RELAX_TYPE::INT:
                   oss << "INT";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::AUX:
                   oss << "AUX";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::MC:
                   oss << "MC";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::SB:
                   oss << "SB";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::SCM:
                   oss << "SCM";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::PWCS:
                   oss << "PWCS";
                   break;
                 case mc::FFKAN<I>::Options::RELAX_TYPE::PWLS:
                   oss << "PWLS";
                   break;
               }
               first = false;
             }
             oss << "}";
             return oss.str();
           });

  pyFFKANOptions
      .def(py::init<>(), R"doc(
Construct an option set initialized with the default values.
)doc")
      .def(py::init<mc::FFKAN<I>::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def(
          "reset", [](mc::FFKAN<I>::Options& self) { self.reset(); },
          R"doc(
Reset all options to their default values.
)doc")
      .def_property(
          "RELAX",
          [](mc::FFKAN<I>::Options& self)
              -> std::set<mc::FFKAN<I>::Options::RELAX_TYPE>&
          { return self.RELAX; },
          [](mc::FFKAN<I>::Options& self, py::object const& val)
          {
            // 1. Optimized path: Assignment from another RELAX_FFKAN (C++ copy)
            if (py::isinstance<std::set<mc::FFKAN<I>::Options::RELAX_TYPE>>(
                    val))
            {
              self.RELAX = val.cast<
                  std::set<mc::FFKAN<I>::Options::RELAX_TYPE> const&>();
              return;
            }
            // 2. Flexible path: Assignment from Python set/list/tuple
            else if (py::isinstance<py::iterable>(val))
            {
              self.RELAX.clear();
              for (auto item : val)
                // Cast items to the Enum type and insert
                self.RELAX.insert(
                    item.cast<mc::FFKAN<I>::Options::RELAX_TYPE>());
              return;
            }
            throw py::type_error(
                "Cannot assign to RELAX: Expected RELAX_FFKAN or iterable of "
                "RELAX_TYPE enums");
          },
          R"doc(
RELAX_FFKAN: Relaxation strategies applied when building polyhedral cuts
for the network node; any iterable of FFKAN.Options.RELAX_TYPE values may
be assigned. Default is {MC}.
)doc")
      .def_readwrite("POLDEF", &mc::FFKAN<I>::Options::POLDEF,
                     "bool: Whether the polyhedral relaxation options default "
                     "to those of the parent polyhedral relaxation "
                     "environment instead of the POLIMG field. Default is "
                     "True.")
      .def_readwrite("POLIMG", &mc::FFKAN<I>::Options::POLIMG,
                     "PolImg.Options: Options for the internal polyhedral "
                     "relaxation environment (used when POLDEF is False).")
      .def_readwrite("SBLIN", &mc::FFKAN<I>::Options::SBLIN,
                     "int: Number of linearization points of the spectral "
                     "bound relaxation. Default is 0, meaning 2*NIN+1 "
                     "points.")
      .def_readwrite(
          "SBSEED", &mc::FFKAN<I>::Options::SBSEED,
          "int: Random number generator seed for the Latin hypercube sampler "
          "placing the linearization points of the spectral bound relaxation. "
          "Default is 42.")
      .def_readwrite("SCMODORD", &mc::FFKAN<I>::Options::SCMODORD,
                     "int: Maximal degree of the sparse Chebyshev model "
                     "relaxation. Default is 3.")
      .def_readwrite(
          "SCBERNORD", &mc::FFKAN<I>::Options::SCBERNORD,
          "int: Degree of the Bernstein basis conversion in the sparse "
          "Chebyshev model relaxation. Default is 0, meaning same as "
          "SCMODORD.")
      .def_readwrite("SCMODEL", &mc::FFKAN<I>::Options::SCMODEL,
                     "SCModel.Options: Options for the sparse Chebyshev model "
                     "environment.")
      .def_readwrite("PWCDIV", &mc::FFKAN<I>::Options::PWCDIV,
                     "int: Equipartition size in the piecewise-constant "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWCREL", &mc::FFKAN<I>::Options::PWCREL,
                     "int: Representation of piecewise-constant univariates - "
                     "0: continuous relaxation; 1: binary encoding. Default "
                     "is 0.")
      .def_readwrite("PWCSLOPE", &mc::FFKAN<I>::Options::PWCSLOPE,
                     "bool: Whether to append cuts from slopes in the "
                     "piecewise-constant superposition model relaxation. "
                     "Default is True.")
      .def_readwrite(
          "PWCSHADOW", &mc::FFKAN<I>::Options::PWCSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-constant superposition model relaxation. Default is "
          "True.")
      .def_readwrite("PWCSUP", &mc::FFKAN<I>::Options::PWCSUP,
                     "PWCSModel.Options: Options for the superposition model "
                     "with piecewise-constant univariate estimators.")
      .def_readwrite("PWLINI", &mc::FFKAN<I>::Options::PWLINI,
                     "int: Initial partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite("PWLMAX", &mc::FFKAN<I>::Options::PWLMAX,
                     "int: Maximal partition size in the piecewise-linear "
                     "superposition model relaxation. Default is 16.")
      .def_readwrite(
          "PWLREL", &mc::FFKAN<I>::Options::PWLREL,
          "int: Representation of piecewise-linear univariates - 0: "
          "continuous relaxation; 1: binary encoding; 2: SOS2 encoding. "
          "Default is 0.")
      .def_readwrite(
          "PWLSHADOW", &mc::FFKAN<I>::Options::PWLSHADOW,
          "bool: Whether to append cuts from shadow estimators in the "
          "piecewise-linear superposition model relaxation. Default is True.")
      .def_readwrite("PWLSUP", &mc::FFKAN<I>::Options::PWLSUP,
                     "PWLSModel.Options: Options for the superposition model "
                     "with adaptive piecewise-linear univariate estimators.");
}
