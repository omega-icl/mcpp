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

#include "ffmlpreg.hpp"

namespace py = pybind11;

typedef mc::MLPREG<I> REG;
typedef mc::FFMLPREG<I> FFREG;
typedef FFREG::Options OPT;
typedef OPT::RELAX_TYPE RELAX;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<OPT::RELAX_TYPE>);

// operands of an MLPREG given as a list: checked against its size
static void
check_size(std::vector<double> const& v, size_t const n, char const* what)
{
  if (v.size() != n)
    throw std::invalid_argument(std::string("MLPREG.") + what + ": " +
                                std::to_string(v.size()) + " values where " +
                                std::to_string(n) + " are expected");
}

void
mc_ffmlpreg(py::module_& m)
{
  py::class_<REG, std::shared_ptr<REG>> pyREG(m, "MLPREG", R"doc(
Multilayer perceptron whose weights are variables, for training.

An ``MLPREG`` holds the structure of a feedforward network -- its number
of inputs, and the number of neurons and the activation function of each
layer -- and is a function of its *operands*: the network inputs followed
by its weights. The weights are ordered layer after layer and, in each
layer, neuron after neuron: the bias of the neuron followed by its weights,
one per neuron of the previous layer (the rows of :meth:`MLP.set_data`, end
to end). It mirrors :class:`MLP`, whose weights are data, and is embedded
in a DAG by :class:`FFMLPREG`.

The structure is set with :meth:`set_structure` or :meth:`append_layer`,
or read with the weights of a network from an ONNX file (or a TorchScript
file, where pymcpp was built with Torch) with :meth:`read_data`. The
weights read are in :attr:`weights`, and :meth:`data` gives the layers of
:meth:`MLP.set_data` for any weights -- a trained ``MLPREG`` becomes an
:class:`MLP`.

Examples
--------
>>> import pymcpp
>>> R = pymcpp.MLPREG(2, [(3, pymcpp.MLPREG.TANH), (1, pymcpp.MLPREG.LINEAR)])
>>> R.nin, R.nout, R.nwei, R.nvar
(2, 1, 13, 15)
)doc");

  py::enum_<REG::ACTIV_TYPE>(pyREG, "ACTIV_TYPE")
      .value("LINEAR", REG::ACTIV_TYPE::LINEAR, "Linear (identity) activation function")
      .value("RELU", REG::ACTIV_TYPE::RELU, "ReLU activation function max(x, 0)")
      .value("TANH", REG::ACTIV_TYPE::TANH, "Tanh activation function")
      .value("SIGMOID", REG::ACTIV_TYPE::SIGMOID, "Sigmoid activation function 1/(1+exp(-x))")
      .export_values();

  py::class_<REG::Options> pyREGOptions(pyREG, "Options", R"doc(
Options controlling how an MLPREG is expressed, evaluated and
differentiated. Set the fields on the :attr:`MLPREG.options` attribute,
e.g. ``R.options.ADJOINT = True``.
)doc");

  pyREG
      .def(py::init<>(), R"doc(
Construct an empty MLPREG; set its structure with set_structure,
append_layer or read_data.
)doc")
      .def(py::init(
               [](size_t const nin, std::vector<std::pair<size_t, int>> const& layers)
               {
                 auto reg = std::make_shared<REG>();
                 if (!reg->set_structure(nin, layers))
                   throw std::invalid_argument("MLPREG: invalid structure");
                 return reg;
               }),
           py::arg("nin"), py::arg("layers"), R"doc(
Construct an MLPREG of `nin` inputs and the layers `layers`: see
:meth:`set_structure`.
)doc")
      .def_readwrite("options", &REG::options, R"doc(
MLPREG.Options: Options controlling how the network is expressed,
evaluated and differentiated.
)doc")
      .def(
          "set_structure",
          [](REG& self, size_t const nin, std::vector<std::pair<size_t, int>> const& layers)
          {
            if (!self.set_structure(nin, layers))
              throw std::invalid_argument("MLPREG.set_structure: invalid structure");
          },
          py::arg("nin"), py::arg("layers"), R"doc(
Set the structure of the network.

Parameters
----------
nin : int
    Number of network inputs.
layers : list of tuple of (int, MLPREG.ACTIV_TYPE)
    Number of neurons and activation function of each layer, hidden layers
    first, output layer last.
)doc")
      .def(
          "append_layer",
          [](REG& self, size_t const nneu, int const activ)
          {
            if (!self.append_layer(nneu, activ))
              throw std::invalid_argument(
                  "MLPREG.append_layer: invalid layer, or no inputs (call "
                  "set_structure first)");
          },
          py::arg("nneu"), py::arg_v("activ", REG::LINEAR, "MLPREG.ACTIV_TYPE.LINEAR"),
          R"doc(
Append a layer of `nneu` neurons with the activation function `activ`,
after :meth:`set_structure` has given the number of inputs.
)doc")
      .def(
          "reset_structure", [](REG& self) { self.reset_structure(); },
          R"doc(
Discard the structure.
)doc")
      .def(
          "read_data",
          [](REG& self, std::string const& filename, bool const disp)
          {
            if (!self.read_data(filename, disp))
              throw std::runtime_error(
                  "MLPREG.read_data: failed to load network data from " + filename +
                  " (call with disp=True for the reason; a TorchScript file "
                  "needs pymcpp compiled with Torch support)");
          },
          py::arg("filename"), py::arg("disp") = false, R"doc(
Read the structure and the weights of a network from an ONNX file (name
ending with ``.onnx``) or from a TorchScript file (where pymcpp was built
with Torch), as :meth:`MLP.read_data` does. The weights read are then in
:attr:`weights`.

Raises
------
RuntimeError
    If the file cannot be read or holds anything else than a network of
    dense layers with ReLU, Tanh or Sigmoid activations. The MLPREG is
    then left as it was.
)doc")
      .def_property_readonly(
          "weights", [](REG const& self) { return self.weights(); },
          R"doc(
list of float: Weights read with the structure (read_data), in the ordering
of the operands; empty if none.
)doc")
      .def(
          "data",
          [](REG const& self, std::vector<double> const& wei)
          {
            check_size(wei, self.nwei(), "data");
            return self.data(wei.data());
          },
          py::arg("wei"), R"doc(
Layers of the network for the weights `wei`, in the form taken by
:meth:`MLP.set_data`: ``MLP().set_data(R.data(w))`` is the network with
these weights.
)doc")
      .def(
          "eval",
          [](REG& self, std::vector<double> const& valin)
          {
            check_size(valin, self.nvar(), "eval");
            std::vector<double> valout(self.nout());
            self.eval(valin.data(), valout.data());
            return valout;
          },
          py::arg("valin"), R"doc(
Outputs of the network at the operands `valin` (inputs, then weights).
)doc")
      .def(
          "grad",
          [](REG& self, std::vector<double> const& valin)
          {
            check_size(valin, self.nvar(), "grad");
            std::vector<double> gradout(self.nout() * self.nvar());
            self.grad(valin.data(), gradout.data());
            return gradout;
          },
          py::arg("valin"), R"doc(
Jacobian of the outputs in the operands at `valin`, row-major: the
derivative of output k in operand i at position k*nvar+i.
)doc")
      .def(
          "tang",
          [](REG& self, std::vector<double> const& valin, std::vector<double> const& dirin)
          {
            check_size(valin, self.nvar(), "tang");
            check_size(dirin, self.nvar(), "tang");
            std::vector<double> tanout(self.nout());
            self.tang(valin.data(), dirin.data(), tanout.data());
            return tanout;
          },
          py::arg("valin"), py::arg("dirin"), R"doc(
Derivatives of the outputs at `valin` along the direction `dirin` of the
operands (one forward pass).
)doc")
      .def(
          "adj",
          [](REG& self, std::vector<double> const& valin, std::vector<double> const& wout)
          {
            check_size(valin, self.nvar(), "adj");
            check_size(wout, self.nout(), "adj");
            std::vector<double> adjin(self.nvar());
            self.adj(valin.data(), wout.data(), adjin.data());
            return adjin;
          },
          py::arg("valin"), py::arg("wout"), R"doc(
Gradient in the operands of the outputs weighted by `wout`, at `valin`
(one backward pass).
)doc")
      .def_property_readonly("nin", &REG::nin, "int: Number of network inputs.")
      .def_property_readonly("nout", &REG::nout, "int: Number of network outputs.")
      .def_property_readonly("nhid", &REG::nhid, "int: Number of hidden layers.")
      .def_property_readonly("nwei", &REG::nwei, "int: Number of weights, biases included.")
      .def_property_readonly("nvar", &REG::nvar,
                             "int: Number of operands: the inputs, then the weights.")
      .def_property_readonly(
          "layers", [](REG const& self) { return self.layers(); },
          "list of tuple of (int, int): Number of neurons and activation "
          "function of each layer.")
      .def_property_readonly("uid", &REG::uid,
                             "int: Stamp of the structure, drawn again whenever it "
                             "changes and kept by the copies.")
      .def_property_readonly(
          "dag", [](REG& self) { return self.dag(); },
          py::return_value_policy::reference_internal, R"doc(
FFGraph: Internal DAG expressing the network, with variables named "X"
(inputs) and "W" (weights).
)doc")
      .def_property_readonly(
          "varin", [](REG& self) { return self.varin(); },
          "list of FFVar: Variables of the internal DAG: inputs, then weights.")
      .def_property_readonly(
          "varout", [](REG& self) { return self.varout(); },
          "list of FFVar: Outputs of the internal DAG.");

  pyREGOptions
      .def(py::init<>(), "Construct an option set initialized with the default values.")
      .def(py::init<REG::Options const&>(), "Copy constructor.")
      .def(
          "reset", [](REG::Options& self) { self.reset(); },
          "Reset all options to their default values.")
      .def_readwrite("EVALTORCH", &REG::Options::EVALTORCH,
                     "bool: Whether to evaluate the network and its derivatives with "
                     "Torch, where it was read from a TorchScript file whose layers "
                     "all have biases. Default is False.")
      .def_readwrite("RELU2ABS", &REG::Options::RELU2ABS,
                     "bool: Whether to express ReLU activations via abs (True) or "
                     "max (False). Default is False.")
      .def_readwrite("SIG2EXP", &REG::Options::SIG2EXP,
                     "bool: Whether to express sigmoid activations via exp (True) "
                     "or tanh (False). Default is False.")
      .def_readwrite("AUTODIFF", &REG::Options::AUTODIFF,
                     "MLPREG.Options.AD_TYPE: Mode of the Jacobian in arithmetics "
                     "other than float: a backward pass per output (B) or forward "
                     "AD (F). Default is B.")
      .def_readwrite("TANGENT", &REG::Options::TANGENT,
                     "bool: Whether the forward symbolic differentiation of the "
                     "network in a DAG inserts its directional derivative, one "
                     "forward pass per sweep (True), or its whole Jacobian (False). "
                     "Default is False.")
      .def_readwrite("ADJOINT", &REG::Options::ADJOINT,
                     "bool: Whether the backward symbolic differentiation of the "
                     "network in a DAG inserts its adjoint, one backward pass per "
                     "sweep (True), or its whole Jacobian (False). Default is False.")
      .def_readwrite("CPMAX", &REG::Options::CPMAX,
                     "int: Maximum rounds of interval constraint propagation. "
                     "Default is 1.")
      .def_readwrite("CPTHRES", &REG::Options::CPTHRES,
                     "float: Threshold for repeating constraint propagation.")
      .def_readwrite("CPINF", &REG::Options::CPINF,
                     "float: Value used as infinity in constraint propagation. "
                     "Default is 1e30.");

  py::enum_<REG::Options::AD_TYPE>(pyREGOptions, "AD_TYPE")
      .value("F", REG::Options::AD_TYPE::F, "Forward differentiation")
      .value("B", REG::Options::AD_TYPE::B, "Backward differentiation")
      .export_values();

  py::class_<FFREG, mc::FFOp> pyFFREG(m, "FFMLPREG", R"doc(
Multilayer perceptron whose weights are variables, embedded as a single
external DAG operation.

FFMLPREG inserts an :class:`MLPREG` into a host :class:`FFGraph` DAG as one
node mapping the network inputs and weights to its outputs. As
:class:`FFMLP`, it supports evaluation in float and validated arithmetics,
symbolic differentiation -- the Jacobian, or with ``MLPREG.options.TANGENT``
and ``ADJOINT`` one forward or backward pass per sweep -- polyhedral
relaxation and reverse constraint propagation. The MLPREG is copied into
the DAG, once per DAG, and returned by :attr:`MLPREG`; options may be
changed after the insertion.

Examples
--------
>>> import pymcpp
>>> R = pymcpp.MLPREG(1, [(2, pymcpp.MLPREG.TANH), (1, pymcpp.MLPREG.LINEAR)])
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(R.nin); W = DAG.add_vars(R.nwei)
>>> Y = pymcpp.FFMLPREG()(X, W, R)
)doc");

  py::class_<OPT, std::shared_ptr<OPT>> pyOPT(pyFFREG, "Options", R"doc(
Options controlling the relaxation of an FFMLPREG operation; see
:class:`FFMLP.Options`.
)doc");

  pyFFREG
      .def(py::init<>(), "Construct an unattached MLPREG operation.")
      .def(
          "__call__",
          [](FFREG& self, std::vector<mc::FFVar> const& vIn,
             std::vector<mc::FFVar> const& vWei, REG* pREG)
          {
            if (!pREG || !pREG->nout())
              throw std::runtime_error("FFMLPREG.__call__: empty network");
            if (vIn.size() != pREG->nin() || vWei.size() != pREG->nwei())
              throw std::invalid_argument(
                  "FFMLPREG.__call__: var and wei must have lengths mlpreg.nin "
                  "and mlpreg.nwei");
            return self(vIn, vWei, pREG, FFREG::COPY);
          },
          py::arg("var"), py::arg("wei"), py::arg("mlpreg"), R"doc(
Insert the MLPREG operation into the DAG and return all network outputs.

Parameters
----------
var : list of FFVar
    Network inputs (variables or constants); length ``mlpreg.nin``.
wei : list of FFVar
    Network weights, in the ordering of the operands; length
    ``mlpreg.nwei``.
mlpreg : MLPREG
    Network structure; copied into the DAG.

Returns
-------
dep : list of FFVar
    The ``mlpreg.nout`` network outputs.
)doc")
      .def(
          "__call__",
          [](FFREG& self, unsigned const idep, std::vector<mc::FFVar> const& vIn,
             std::vector<mc::FFVar> const& vWei, REG* pREG)
          {
            if (!pREG || !pREG->nout())
              throw std::runtime_error("FFMLPREG.__call__: empty network");
            if (vIn.size() != pREG->nin() || vWei.size() != pREG->nwei())
              throw std::invalid_argument(
                  "FFMLPREG.__call__: var and wei must have lengths mlpreg.nin "
                  "and mlpreg.nwei");
            if (idep >= pREG->nout())
              throw std::invalid_argument(
                  "FFMLPREG.__call__: dep must be smaller than mlpreg.nout");
            return self(idep, vIn, vWei, pREG, FFREG::COPY);
          },
          py::arg("dep"), py::arg("var"), py::arg("wei"), py::arg("mlpreg"),
          "Insert the MLPREG operation into the DAG and return output `dep` only.")
      .def_property_readonly(
          "name", [](FFREG const& self) { return self.name(); },
          "str: Name of this operation as displayed in the DAG.")
      .def("__str__",
           [](FFREG const& self)
           {
             std::ostringstream Oss;
             Oss << self;
             return Oss.str();
           })
      .def_property(
          "options", [](FFREG& self) { return self.share_options(); },
          [](FFREG& self, OPT const& opt) { self.home_options() = opt; },
          R"doc(
FFMLPREG.Options: Options controlling the relaxation of the operation, as
for :attr:`FFMLP.options`.
)doc")
      .def_property_readonly(
          "MLPREG", [](FFREG& self) { return self.share_object(); },
          R"doc(
MLPREG or None: The copy held by the DAG of the network inserted last,
None before the first call. Its options can be changed, e.g.
``Op.MLPREG.options.ADJOINT = True``: the operations follow.
)doc")
      .def_property_readonly(
          "attached", [](FFREG const& self) { return self.attached(); },
          "bool: Whether the object has inserted an operation in a DAG.");

  py::enum_<RELAX>(pyOPT, "RELAX_TYPE")
      .value("INT", RELAX::INT, "Interval bounds")
      .value("AUX", RELAX::AUX, "Auxiliary variable polyhedral relaxation")
      .value("MC", RELAX::MC, "McCormick relaxation with interval bounds")
      .value("SB", RELAX::SB, "Spectral bounds")
      .value("SCM", RELAX::SCM, "Sparse Chebyshev model relaxation")
      .value("PWCS", RELAX::PWCS, "Piecewise-constant superposition bounds")
      .value("PWLS", RELAX::PWLS, "Piecewise-linear superposition bounds")
      .export_values();

  py::class_<std::set<RELAX>>(m, "RELAX_FFMLPREG", py::module_local(), R"doc(
Set of relaxation strategies for FFMLPREG operations; see
:class:`RELAX_FFMLP`.
)doc")
      .def(py::init<>())
      .def("add", [](std::set<RELAX>& s, RELAX const& k) { s.insert(k); })
      .def("discard", [](std::set<RELAX>& s, RELAX const& k) { s.erase(k); })
      .def("__contains__", [](std::set<RELAX>& s, RELAX const& k) { return s.count(k) > 0; })
      .def("__len__", [](std::set<RELAX>& s) { return s.size(); })
      .def(
          "__iter__", [](std::set<RELAX>& s) { return py::make_iterator(s.begin(), s.end()); },
          py::keep_alive<0, 1>());

  pyOPT.def(py::init<>())
      .def(py::init<OPT const&>())
      .def("reset", [](OPT& self) { self.reset(); })
      .def_property(
          "RELAX", [](OPT& self) -> std::set<RELAX>& { return self.RELAX; },
          [](OPT& self, py::object const& val)
          {
            if (py::isinstance<std::set<RELAX>>(val))
            {
              self.RELAX = val.cast<std::set<RELAX> const&>();
              return;
            }
            if (py::isinstance<py::iterable>(val))
            {
              self.RELAX.clear();
              for (auto item : val) self.RELAX.insert(item.cast<RELAX>());
              return;
            }
            throw py::type_error("RELAX: expected an iterable of RELAX_TYPE values");
          },
          "RELAX_FFMLPREG: Relaxation strategies of the operation. Default is {MC}.")
      .def_readwrite("POLDEF", &OPT::POLDEF)
      .def_readwrite("POLIMG", &OPT::POLIMG)
      .def_readwrite("SBLIN", &OPT::SBLIN)
      .def_readwrite("SBSEED", &OPT::SBSEED)
      .def_readwrite("SCMODORD", &OPT::SCMODORD)
      .def_readwrite("SCBERNORD", &OPT::SCBERNORD)
      .def_readwrite("SCMODEL", &OPT::SCMODEL)
      .def_readwrite("PWCDIV", &OPT::PWCDIV)
      .def_readwrite("PWCREL", &OPT::PWCREL)
      .def_readwrite("PWCSLOPE", &OPT::PWCSLOPE)
      .def_readwrite("PWCSHADOW", &OPT::PWCSHADOW)
      .def_readwrite("PWCSUP", &OPT::PWCSUP)
      .def_readwrite("PWLINI", &OPT::PWLINI)
      .def_readwrite("PWLMAX", &OPT::PWLMAX)
      .def_readwrite("PWLREL", &OPT::PWLREL)
      .def_readwrite("PWLSHADOW", &OPT::PWLSHADOW)
      .def_readwrite("PWLSUP", &OPT::PWLSUP);
}
