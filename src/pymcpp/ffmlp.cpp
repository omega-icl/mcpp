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
  py::class_<mc::MLP<I>> pyMLP(m, "MLP");

  py::enum_<mc::MLP<I>::ACTIV_TYPE>(pyMLP, "ACTIV_TYPE")
      .value("LINEAR", mc::MLP<I>::ACTIV_TYPE::LINEAR,
             "Linear activation function")
      .value("RELU", mc::MLP<I>::ACTIV_TYPE::RELU, "ReLU activation function")
      .value("TANH", mc::MLP<I>::ACTIV_TYPE::TANH, "Tanh activation function")
      .value("SIGMOID", mc::MLP<I>::ACTIV_TYPE::SIGMOID,
             "Sigmoid activation function")
      .export_values();

  py::class_<mc::MLP<I>::Options> pyMLPOptions(
      pyMLP, "Options");  //, py::module_local() );

  pyMLP.def(py::init<>(), "default constructor")
      .def_readwrite("options", &mc::MLP<I>::options)
      .def(
          "reset_data", [](mc::MLP<I>& self) { self.reset_data(); },
          "reset MLP data")
      .def(
          "set_data",
          [](mc::MLP<I>& self,
             std::vector<
                 std::pair<std::vector<std::vector<double>>, int>> const& mlp)
          { self.set_data(mlp); },
          "set MLP data")
      .def(
          "append_data",
          [](mc::MLP<I>& self, std::vector<std::vector<double>> const& layer,
             int const activ, bool const reset)
          { self.append_data(layer, activ, reset); },
          py::arg("layer"), py::arg("activ") = mc::MLP<I>::LINEAR,
          py::arg("reset") = false, "Append multi-neuron layer to MLP data")
      .def(
          "append_data",
          [](mc::MLP<I>& self, std::vector<double> const& layer,
             int const activ, bool const reset)
          { self.append_data(layer, activ, reset); },
          py::arg("layer"), py::arg("activ") = mc::MLP<I>::LINEAR,
          py::arg("reset") = false, "Append single-layer layer to MLP data")
      .def(
          "read_data",
          [](mc::MLP<I>& self, std::string const& filename, bool const disp)
          { self.read_data(filename, disp); }, py::arg("filename"),
          py::arg("disp") = false, "read MLP data from Torch script")
      .def_property_readonly(
          "dag", [](mc::MLP<I>& self) { return self.dag(); },
          py::return_value_policy::reference_internal, "MLP internal DAG")
      .def_property_readonly(
          "nin", [](mc::MLP<I> const& self) { return self.nin(); },
          py::return_value_policy::reference_internal,
          "MLP number of input variables")
      .def_property_readonly(
          "varin", [](mc::MLP<I>& self) { return self.varin(); },
          py::return_value_policy::reference_internal, "MLP input variables")
      .def_property_readonly(
          "nout", [](mc::MLP<I> const& self) { return self.nout(); },
          py::return_value_policy::reference_internal,
          "MLP number of output variables")
      .def_property_readonly(
          "varout", [](mc::MLP<I>& self) { return self.varout(); },
          py::return_value_policy::reference_internal, "MLP output variables")
      .def_property_readonly(
          "nhid", [](mc::MLP<I> const& self) { return self.nhid(); },
          py::return_value_policy::reference_internal,
          "MLP number of hidden layers");


  pyMLPOptions.def(py::init<>())
      .def(py::init<mc::MLP<I>::Options const&>())
      .def(
          "reset", [](mc::MLP<I>::Options& self) { self.reset(); },
          "reset MLP options")
      .def_readwrite("EVALTORCH", &mc::MLP<I>::Options::EVALTORCH,
                     "Whether to evaluate ANN in floating-point arithemic "
                     "using Torch [Default: False]")
      .def_readwrite("ZEROTOL", &mc::MLP<I>::Options::ZEROTOL,
                     "Threshold for zero coefficient in neural network "
                     "[Default: DBL_EPSILON]")
      .def_readwrite("RELU2ABS", &mc::MLP<I>::Options::RELU2ABS,
                     "Whether to convert ReLU to abs (true) or max (false) "
                     "function [Default: False]")
      .def_readwrite("SIG2EXP", &mc::MLP<I>::Options::SIG2EXP,
                     "Whether to convert sigmoid to exp (true) or tanh (false) "
                     "[Default: False]")
      .def_readwrite("AUTODIFF", &mc::MLP<I>::Options::AUTODIFF,
                     "Automatic differentiation mode [Default: F]")
      .def_readwrite("CPMAX", &mc::MLP<I>::Options::CPMAX,
                     "Maximum rounds of constraint propagation [Default: 1]")
      .def_readwrite(
          "CPTHRES", &mc::MLP<I>::Options::CPTHRES,
          "Threshold for repeating constraint propagation (minimum relative "
          "reduction in any variable) [Default: 1e4*DBL_EPSILON]")
      .def_readwrite("CPINF", &mc::MLP<I>::Options::CPINF,
                     "Infinite value for unbounded variables in constraint "
                     "propagation [Default: 1e30]");

  py::enum_<mc::MLP<I>::Options::AD_TYPE>(pyMLPOptions, "AD_TYPE")
      .value("F", mc::MLP<I>::Options::AD_TYPE::F, "Forward differentiation")
      .value("B", mc::MLP<I>::Options::AD_TYPE::B, "Backward differentiation")
      .export_values();

  py::class_<mc::FFMLP<I>, mc::FFOp> pyFFMLP(m, "FFMLP");

  py::class_<mc::FFMLP<I>::Options> pyFFMLPOptions(pyFFMLP, "Options");

  pyFFMLP.def(py::init<>(), "default constructor")
      .def(py::init<mc::FFMLP<I> const&>(), "copy constructor")
      .def(
          "__call__",
          [](mc::FFMLP<I>& self, std::vector<mc::FFVar> const& vVar,
             mc::MLP<I>* pMLP) { return self(vVar, pMLP, mc::FFMLP<I>::COPY); },
          py::return_value_policy::reference_internal, py::arg("var"),
          py::arg("mlp"), "define MLP operation in DAG")
      .def(
          "__call__",
          [](mc::FFMLP<I>& self, unsigned const idep,
             std::vector<mc::FFVar> const& vVar, mc::MLP<I>* pMLP)
          { return self(idep, vVar, pMLP, mc::FFMLP<I>::COPY); },
          py::return_value_policy::reference_internal, py::arg("dep"),
          py::arg("var"), py::arg("mlp"), "define MLP operation in DAG")
      .def_property_readonly(
          "name", [](mc::FFMLP<I> const& self) { return self.name(); },
          py::return_value_policy::reference_internal, "MLP operation name")
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
      .def_readwrite("options", &mc::FFMLP<I>::options);


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
                                                          py::module_local())
      .def(py::init<>())
      .def(
          "add",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k) { s.insert(k); },
          "Add an element to the set")
      .def(
          "remove",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k)
          {
            if (s.find(k) == s.end()) throw py::key_error();
            s.erase(k);
          },
          "Remove an element from the set")
      .def(
          "discard",
          [](std::set<mc::FFMLP<I>::Options::RELAX_TYPE>& s,
             mc::FFMLP<I>::Options::RELAX_TYPE const& k) { s.erase(k); },
          "Remove an element if present")
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

  pyFFMLPOptions.def(py::init<>())
      .def(py::init<mc::FFMLP<I>::Options const&>())
      .def(
          "reset", [](mc::FFMLP<I>::Options& self) { self.reset(); },
          "reset options to defaults")
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
          "Type of relaxations [Default: {MC}]")
      .def_readwrite("POLDEF", &mc::FFMLP<I>::Options::POLDEF,
                     "Whether options default to the parent polyhedral "
                     "relaxation environment [Default: True]")
      .def_readwrite("POLIMG", &mc::FFMLP<I>::Options::POLIMG,
                     "Options for polyhedral relaxation environment")
      .def_readwrite("SBLIN", &mc::FFMLP<I>::Options::SBLIN,
                     "Number of linearization points of spectral bound "
                     "relaxation [Default: 0 (NVAR+1 points)]")
      .def_readwrite(
          "SBSEED", &mc::FFMLP<I>::Options::SBSEED,
          "Random number generator seed in LHS sampler for linearization "
          "points of spectral bound relaxation [Default: 42]")
      .def_readwrite("SCMODORD", &mc::FFMLP<I>::Options::SCMODORD,
                     "Maximal degree of sparse Chebyshev model [Default: 3]")
      .def_readwrite(
          "SCBERNORD", &mc::FFMLP<I>::Options::SCBERNORD,
          "Degree of Bernstein basis conversion [Default: 0 (same as SCORD)]")
      .def_readwrite("SCMODEL", &mc::FFMLP<I>::Options::SCMODEL,
                     "Options for sparse Chebyshev model")
      .def_readwrite("PWCDIV", &mc::FFMLP<I>::Options::PWCDIV,
                     "Equipartition size in piecewise-constant superposition "
                     "model [Default: 16]")
      .def_readwrite("PWCREL", &mc::FFMLP<I>::Options::PWCREL,
                     "Representation of piecewise-constant univariates - 0: "
                     "continuous relaxation; 1: binary encoding [Default: 0]")
      .def_readwrite("PWCSLOPE", &mc::FFMLP<I>::Options::PWCSLOPE,
                     "Whether to append cuts from slopes in piecewise-constant "
                     "superposition model relaxation [Default: True]")
      .def_readwrite(
          "PWCSHADOW", &mc::FFMLP<I>::Options::PWCSHADOW,
          "Whether to append cuts from shadow estimators in piecewise-constant "
          "superposition model relaxation [Default: True]")
      .def_readwrite("PWCSUP", &mc::FFMLP<I>::Options::PWCSUP,
                     "Options for superposition model with piecewise-constant "
                     "univariate estimators")
      .def_readwrite("PWLINI", &mc::FFMLP<I>::Options::PWLINI,
                     "Initial partition size in piecewise-linear superposition "
                     "model [Default: 16]")
      .def_readwrite("PWLMAX", &mc::FFMLP<I>::Options::PWLMAX,
                     "Maximal partition size in piecewise-linear superposition "
                     "model [Default: 16]")
      .def_readwrite(
          "PWLREL", &mc::FFMLP<I>::Options::PWLREL,
          "Representation of piecewise-linear univariates - 0: continuous "
          "relaxation; 1: binary encoding; 2: SOS2 encoding [Default: 0]")
      .def_readwrite(
          "PWLSHADOW", &mc::FFMLP<I>::Options::PWLSHADOW,
          "Whether to append cuts from shadow estimators in piecewise-linear "
          "superposition model relaxation [Default: True]")
      .def_readwrite("PWLSUP", &mc::FFMLP<I>::Options::PWLSUP,
                     "Options for superposition model with adaptive "
                     "piecewise-linear univariate estimators");
}
