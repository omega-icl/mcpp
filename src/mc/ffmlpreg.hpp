// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__FFMLPREG_HPP
#define MC__FFMLPREG_HPP

#include <algorithm>
#include <map>
#include <memory>
#include <mutex>
#include <thread>
#include <type_traits>
#if defined(MC__USE_TORCH)
#include <torch/autograd.h>
#include <torch/script.h>
#endif

#include "ffextern.hpp"
#include "interval.hpp"
#include "netio.hpp"

// #define MC__FFMLPREG_DEBUG
#define MC__FFMLPREG_CHECK

namespace mc
{

//! @brief C++ class for multilayer perceptrons whose weights are variables
////////////////////////////////////////////////////////////////////////
//! mc::MLPREG is a C++ class for the evaluation, differentiation and
//! relaxation of multilayer perceptrons (MLP) as functions of their
//! inputs and of their weights and biases -- for training, or for the
//! estimation of the weights -- that leverages expression trees and
//! arithmetics available through MC++. It mirrors mc::MLP, whose weights
//! are data: an MLPREG has the structure of a network -- its number of
//! inputs, and the number of neurons and the activation function of each
//! layer -- and is evaluated on a vector of operands made of its inputs
//! followed by its weights (nvar() = nin() + nwei()). The weights are
//! ordered layer after layer and, in each layer, neuron after neuron: the
//! bias of the neuron followed by its weights, one per neuron of the
//! previous layer -- the rows of mc::MLP::set_data, end to end.
//!
//! The structure is set (set_structure, append_layer) or read with the
//! weights of a trained network (read_data): from an ONNX file, or from a
//! TorchScript Module, which can then be evaluated by Torch
//! (Options::EVALTORCH). The weights read are kept (weights()), and
//! data() gives the layers of mc::MLP::set_data for any weights.
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class MLPREG
////////////////////////////////////////////////////////////////////////
{
 public:
  //! @brief Enumeration type for activation function
  enum ACTIV_TYPE
  {
    LINEAR = 0,  //!< Linear activation function
    RELU,        //!< ReLU activation function
    TANH,        //!< tanh activation function
    SIGMOID      //!< Sigmoid activation function
  };

  //! @brief Layers of a structure: number of neurons and activation function
  //! of each layer, the last one being the output layer
  typedef std::vector<std::pair<size_t, int>> t_Layers;

  //! @brief MLPREG options
  struct Options
  {
    //! @brief Default constructor
    Options() { reset(); }

    //! @brief Copy constructor
    Options(Options const& other) = default;

    //! @brief Assignment operator
    Options&
    operator=(Options const& other)
    {
      EVALTORCH = other.EVALTORCH;
      RELU2ABS  = other.RELU2ABS;
      SIG2EXP   = other.SIG2EXP;
      AUTODIFF  = other.AUTODIFF;
      TANGENT   = other.TANGENT;
      ADJOINT   = other.ADJOINT;
      CPMAX     = other.CPMAX;
      CPTHRES   = other.CPTHRES;
      CPINF     = other.CPINF;
      return *this;
    }

    //! @brief Reset options
    void
    reset()
    {
      EVALTORCH = false;
      RELU2ABS  = false;
      SIG2EXP   = false;
      AUTODIFF  = B;
      TANGENT   = false;
      ADJOINT   = false;
      CPMAX     = 1;
      CPTHRES   = 1e4 * DBL_EPSILON;
      CPINF     = 1e30;
    }

    //! @brief Enumeration type for AD strategy
    enum AD_TYPE
    {
      F = 0,  //!< Forward differentiation
      B       //!< Backward differentiation
    };

    //! @brief Whether to evaluate the network and its derivatives in
    //! floating-point arithmetic using Torch, in double precision, where it
    //! was read from a TorchScript Module whose layers all have biases. Each
    //! thread works on a copy of its own of the Module, into which the
    //! weights of each evaluation are written
    bool EVALTORCH;
    //! @brief Whether to convert ReLU to abs (true) or max (false) function. It
    //! may be changed at any time, also on the object of an operation in a DAG:
    //! the expression tree of the network is built again
    bool RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false). It may
    //! be changed at any time, as RELU2ABS
    bool SIG2EXP;
    //! @brief Mode of the Jacobian (grad) in arithmetics other than double:
    //! one backward pass through the layers per output (B, default), or the
    //! forward mode of FADType with one direction per operand (F). In double
    //! arithmetic, the reverse mode is applied where the outputs are no more
    //! than the operands -- that is, nearly always -- and the forward mode
    //! otherwise with F. (With the weights among the operands, there are far
    //! more operands than outputs: mc::MLP has F as its default.)
    int AUTODIFF;
    //! @brief Whether the forward symbolic differentiation of an MLPREG
    //! operation in a DAG (FAD, SFAD, DFAD) inserts its directional derivative
    //! along the tangents of its operands, mc::FFTanMLPREG -- one forward pass
    //! per sweep -- in place of its whole Jacobian, mc::FFGradMLPREG (false).
    //! With the weights among the operands, the Jacobian has nout()*nvar()
    //! entries, and the passes are always the cheaper
    bool TANGENT;
    //! @brief Whether the backward symbolic differentiation of an MLPREG
    //! operation in a DAG (BAD, SBAD, DBAD) inserts the product of its
    //! transposed Jacobian with the adjoints of its results,
    //! mc::FFAdjMLPREG -- one backward pass per sweep -- in place of its whole
    //! Jacobian, mc::FFGradMLPREG (false). The passes are always the cheaper
    bool ADJOINT;
    //! @brief Maximum rounds of constraint propagation
    size_t CPMAX;
    //! @brief Threshold for repeating constraint propagation (minimum relative
    //! reduction in any variable)
    double CPTHRES;
    //! @brief Infinite value for unbounded variables in constraint propagation
    double CPINF;
  } options;

 private:
  //! @brief Number of inputs
  size_t _nin;
  //! @brief Number of outputs
  size_t _nout;
  //! @brief Number of weights and biases
  size_t _nwei;
  //! @brief Layers
  t_Layers _layers;
  //! @brief Position of the first weight of each layer among the weights
  std::vector<size_t> _off;
  //! @brief Weights read with the structure (empty if none)
  std::vector<double> _wei;
  //! @brief Stamp of the structure: see uid()
  size_t _uid;

  //! @brief Expression tree
  FFGraph* _dag;
  //! @brief Expression tree needs updating
  bool _update;
  //! @brief Options RELU2ABS and SIG2EXP with which the expression tree was built: the tree is
  //! built again when they change
  std::pair<bool, bool> _dagopt{false, false};
  //! @brief Variables of the expression tree: inputs, then weights
  std::vector<FFVar> _varin;
  //! @brief Dependents of the expression tree
  std::vector<FFVar> _varout;
  //! @brief Codelist
  FFSubgraph _codelist;
  //! @brief Work arrays of the evaluations through the expression tree
  std::vector<PolVar<T>> _wkPOL;
  std::vector<T> _wkCPI;

  //! @brief Torch Module available for the evaluation
  bool _hasmodule;
#if defined(MC__USE_TORCH)
  //! @brief Copy of the Torch Module of a thread, and its Linear layers
  struct TorchSlot
  {
    torch::jit::script::Module module;
    std::vector<std::pair<at::Tensor, at::Tensor>> lin;
    std::vector<at::Tensor> param;
  };
  //! @brief Copies of the Torch Module, one per thread
  struct TorchPool
  {
    std::mutex mtx;
    std::map<std::thread::id, std::unique_ptr<TorchSlot>> slots;
  };
  //! @brief Torch Module read, never evaluated itself
  torch::jit::script::Module _module;
  //! @brief Copies of the Torch Module, one per thread -- not shared with
  //! the copies of the MLPREG
  std::shared_ptr<TorchPool> _pool;

  //! @brief Copy of the Torch Module of the calling thread
  TorchSlot& _slot() const;
  //! @brief Write the weights wei into the copy of the Torch Module of a
  //! thread, with gradients if grad
  void _torch_load(TorchSlot& slot, double const* wei, bool const grad) const;
  //! @brief Forward pass through the copy of the Torch Module of a thread,
  //! at the inputs and weights of valin -- with gradients if grad, the
  //! inputs then also returned in tensin
  torch::Tensor _torch_forward(TorchSlot& slot, double const* valin,
                               bool const grad, torch::Tensor* tensin) const;
  //! @brief Gradients returned by Torch, written in the ordering of the
  //! operands
  void _torch_row(std::vector<torch::Tensor> const& g, double* row) const;
  static torch::Tensor
  _tensor(double const* val, size_t const n)
  {
    return torch::from_blob(const_cast<double*>(val), {(int64_t)n},
                            torch::TensorOptions().dtype(torch::kFloat64))
        .clone();
  }
  static void
  _copy(torch::Tensor const& t, size_t const n, double* val)
  {
    if (!t.defined() || (size_t)t.numel() != n)
      throw std::runtime_error(
          "MLPREG: Torch Module returns " +
          std::to_string(t.defined() ? t.numel() : 0) + " values where " +
          std::to_string(n) + " are expected\n");
    torch::Tensor const td =
        t.detach().to(torch::kCPU, torch::kFloat64).contiguous();
    std::copy(td.data_ptr<double>(), td.data_ptr<double>() + n, val);
  }
  bool
  _torch() const
  {
    return _hasmodule && options.EVALTORCH;
  }
#else
  bool
  _torch() const
  {
    return false;
  }
#endif

  //! @brief ReLU activation
  template <typename U>
  U
  ReLU(U const& x) const
  {
    return options.RELU2ABS ? (x + Op<U>::fabs(x)) * 0.5 : Op<U>::max(x, 0.);
  }
  template <typename U>
  FADType<U>
  ReLU(FADType<U> const& x) const
  {
    FADType<U> z = ReLU(x.val());
    z.setDepend(x);
    for (unsigned j = 0; j < z.size(); ++j) z[j] = Op<U>::fstep(x.val()) * x[j];
    return z;
  }

  //! @brief Activation of val
  template <typename U>
  void
  _act(int const activ, U& val) const
  {
    switch (activ)
    {
      case LINEAR:
      default:
        break;
      case RELU:
        val = ReLU(val);
        break;
      case TANH:
        val = Op<U>::tanh(val);
        break;
      case SIGMOID:
        val = (options.SIG2EXP ? Op<U>::inv(Op<U>::exp(-val) + 1.)
                               : Op<U>::tanh(val * 0.5) * 0.5 + 0.5);
        break;
    }
  }

  //! @brief Activation of val and its derivative der, for a nonlinear
  //! activation function (false for LINEAR, der then left as it was)
  template <typename U>
  bool
  _actder(int const activ, U& val, U& der) const
  {
    switch (activ)
    {
      case LINEAR:
      default:
        return false;
      case RELU:
        der = Op<U>::fstep(val);
        val = ReLU(val);
        return true;
      case TANH:
        val = Op<U>::tanh(val);
        der = -Op<U>::sqr(val) + 1.;
        return true;
      case SIGMOID:
        _act(SIGMOID, val);
        // s (1-s), written with a single occurrence of s
        der = (-Op<U>::sqr(val * 2. - 1.) + 1.) * 0.25;
        return true;
    }
  }

  //! @brief Set the sizes and positions of the weights after a change of the
  //! layers
  bool _resize();

  //! @brief Set expression tree
  void _set_dag();

  //! @brief Evaluate MLPREG
  template <typename U>
  void _eval(U const* valin, U* valout) const;

  //! @brief Evaluate MLPREG directional derivative: forward pass with the
  //! values and the tangents of the neurons
  template <typename U>
  void _tang(U const* valin, U const* dirin, U* tanout, U* valout) const;

  //! @brief Evaluate MLPREG adjoint: forward pass keeping the derivative of
  //! each activation, then backward pass through the layers
  template <typename U>
  void _adj(U const* valin, U const* wout, U* adjin, U* valout) const;

  //! @brief Evaluate MLPREG Jacobian in reverse mode: forward pass keeping the
  //! derivative of each activation, then all the rows of the Jacobian backward
  //! through the layers together
  template <typename U>
  void _jacrev(U const* valin, U* gradout, U* valout) const;

  //! @brief Evaluate MLPREG Jacobian in forward mode, with the FADType of U
  template <typename U>
  void _jacfad(U const* valin, U* gradout, U* valout) const;

 public:
  //! @brief Default constructor
  MLPREG()
      : _nin(0),
        _nout(0),
        _nwei(0),
        _uid(ffhome_uid()),
        _dag(nullptr),
        _update(false),
        _hasmodule(false)
#if defined(MC__USE_TORCH)
        ,
        _pool(std::make_shared<TorchPool>())
#endif
  {
  }

  //! @brief Constructor of a structure: see set_structure
  MLPREG(size_t const nin, t_Layers const& layers) : MLPREG()
  {
    set_structure(nin, layers);
  }

  //! @brief Copy constructor
  MLPREG(MLPREG const& other)
      : options(other.options),
        _nin(other._nin),
        _nout(other._nout),
        _nwei(other._nwei),
        _layers(other._layers),
        _off(other._off),
        _wei(other._wei),
        _uid(other._uid),
        _dag(nullptr),
        _update(false),
        _hasmodule(other._hasmodule)
#if defined(MC__USE_TORCH)
        ,
        _module(other._module),
        _pool(std::make_shared<TorchPool>())
#endif
  {
  }

  //! @brief Assignment operator: takes the structure, the weights read, the
  //! options and the Torch Module of another MLPREG; the expression tree is
  //! built again when it is needed
  MLPREG&
  operator=(MLPREG const& other)
  {
    if (this == &other) return *this;
    options = other.options;
    _nin    = other._nin;
    _nout   = other._nout;
    _nwei   = other._nwei;
    _layers = other._layers;
    _off    = other._off;
    _wei    = other._wei;
    _uid    = other._uid;
    delete _dag;
    _dag    = nullptr;
    _update = false;
    _varin.clear();
    _varout.clear();
    _codelist.clear();
    _wkPOL.clear();
    _wkCPI.clear();
    _hasmodule = other._hasmodule;
#if defined(MC__USE_TORCH)
    _module = other._module;
    _pool   = std::make_shared<TorchPool>();
#endif
    return *this;
  }

  virtual ~MLPREG() { delete _dag; }

  //! @brief Reset the structure
  void
  reset_structure()
  {
    _nin = _nout = _nwei = 0;
    _layers.clear();
    _off.clear();
    _wei.clear();
    _hasmodule = false;
    _update    = true;
    _uid       = ffhome_uid();
  }

  //! @brief Set the structure: the number of inputs nin, and the number of
  //! neurons and the activation function of each layer -- {{number of
  //! neurons of the first hidden layer, activation}, ..., {number of outputs,
  //! activation of the output layer}}. False, with the structure unchanged,
  //! for a layer without neuron or with an unknown activation, or without
  //! input. No layer leaves a structure with inputs only, for append_layer
  bool
  set_structure(size_t const nin, t_Layers const& layers)
  {
    if (!nin) return false;
    for (auto const& [nneu, activ] : layers)
      if (!nneu || activ < LINEAR || activ > SIGMOID) return false;
    _nin    = nin;
    _layers = layers;
    _wei.clear();
    _hasmodule = false;
    return _resize();
  }

  //! @brief Append a layer of nneu neurons with the activation activ, after
  //! set_structure has given the number of inputs
  bool
  append_layer(size_t const nneu, int const activ = LINEAR)
  {
    if (!_nin || !nneu || activ < LINEAR || activ > SIGMOID) return false;
    _layers.push_back(std::make_pair(nneu, activ));
    _wei.clear();
    _hasmodule = false;
    return _resize();
  }

  //! @brief Read the structure and the weights of a network from an ONNX
  //! file (name ending with .onnx), or from a TorchScript Module (with
  //! Torch), as mc::MLP::read_data does: the weights read are in weights().
  //! A Module whose layers all have biases may then be evaluated by Torch
  //! (Options::EVALTORCH). False, with the MLPREG unchanged, where the file
  //! cannot be read or holds anything else than dense layers with ReLU, tanh
  //! or sigmoid activations (see mc::NETIO); the reason is displayed with
  //! disp
  bool read_data(std::string const& filename, bool const disp = false);

  //! @brief Weights read with the structure (read_data), in the ordering of
  //! the operands; empty if none
  std::vector<double> const&
  weights() const
  {
    return _wei;
  }

  //! @brief Layers of the network for the weights wei, in the form of
  //! mc::MLP::set_data
  NETIO::t_Data
  data(double const* wei) const
  {
    NETIO::t_Data dat;
    size_t nprv = _nin;
    for (size_t l = 0; l < _layers.size(); ++l)
    {
      std::vector<std::vector<double>> layer(_layers[l].first);
      double const* w = wei + _off[l];
      for (auto& neuron : layer)
      {
        neuron.assign(w, w + nprv + 1);
        w += nprv + 1;
      }
      dat.push_back(std::make_pair(layer, _layers[l].second));
      nprv = _layers[l].first;
    }
    return dat;
  }

  //! @brief Number of inputs
  size_t
  nin() const
  {
    return _nin;
  }
  //! @brief Number of outputs
  size_t
  nout() const
  {
    return _nout;
  }
  //! @brief Number of hidden layers
  size_t
  nhid() const
  {
    return _layers.empty() ? 0 : _layers.size() - 1;
  }
  //! @brief Number of weights, biases included
  size_t
  nwei() const
  {
    return _nwei;
  }
  //! @brief Number of operands: the inputs, then the weights
  size_t
  nvar() const
  {
    return _nin + _nwei;
  }
  //! @brief Layers: number of neurons and activation function of each
  t_Layers const&
  layers() const
  {
    return _layers;
  }
  //! @brief Stamp of the structure -- its layers, the weights read and its
  //! Torch Module if any: a number that no other structure has had, drawn
  //! again whenever the structure is set, appended to, read or reset, and
  //! kept by the copies and by an MLPREG it is assigned to. The options are
  //! not part of it. An MLPREG operation that holds a copy of the structure
  //! is identified in a DAG by this stamp: see mc::FFHOME
  size_t
  uid() const
  {
    return _uid;
  }

  //! @brief Variables of the expression tree: the inputs, then the weights
  std::vector<FFVar> const&
  varin()
  {
    _set_dag();
    return _varin;
  }
  //! @brief Dependents of the expression tree
  std::vector<FFVar> const&
  varout()
  {
    _set_dag();
    return _varout;
  }
  //! @brief Expression tree
  FFGraph*
  dag()
  {
    _set_dag();
    return _dag;
  }

  //! @brief MLPREG evaluation: valout (size nout) at the operands valin (size
  //! nvar: the inputs, then the weights)
  void
  eval(double const* valin, double* valout)
  {
#if defined(MC__USE_TORCH)
    if (_torch())
    {
      TorchSlot& slot = _slot();
      torch::NoGradGuard nograd;
      _copy(_torch_forward(slot, valin, false, nullptr), _nout, valout);
      return;
    }
#endif
    _eval(valin, valout);
  }
  void
  eval(PolVar<T> const* valin, PolVar<T>* valout)
  {
    _set_dag();
    _dag->eval(_codelist, _wkPOL, _varout.size(), _varout.data(), valout,
               _varin.size(), _varin.data(), valin);
  }
  template <typename U>
  void
  eval(U const* valin, U* valout)
  {
    _eval(valin, valout);
  }

  //! @brief MLPREG Jacobian: gradout (size nout*nvar, row-major) at the
  //! operands valin (size nvar), and the outputs valout if not null
  void grad(double const* valin, double* gradout, double* const valout = nullptr);
  template <typename U>
  void
  grad(U const* valin, U* gradout, U* const valout = nullptr)
  {
    if (options.AUTODIFF == Options::AD_TYPE::F)
      _jacfad(valin, gradout, valout);
    else
      _jacrev(valin, gradout, valout);
  }

  //! @brief MLPREG directional derivative: tanout (size nout) is the
  //! derivative of the outputs at valin along the direction dirin (size nvar)
  void tang(double const* valin, double const* dirin, double* tanout,
            double* const valout = nullptr);
  template <typename U>
  void
  tang(U const* valin, U const* dirin, U* tanout, U* const valout = nullptr)
  {
    _tang(valin, dirin, tanout, valout);
  }

  //! @brief MLPREG adjoint: adjin (size nvar) is the product of the
  //! transposed Jacobian of the outputs at valin with the weights wout (size
  //! nout), i.e. the gradient of the weighted sum of the outputs in the
  //! inputs and the weights
  void adj(double const* valin, double const* wout, double* adjin,
           double* const valout = nullptr);
  template <typename U>
  void
  adj(U const* valin, U const* wout, U* adjin, U* const valout = nullptr)
  {
    _adj(valin, wout, adjin, valout);
  }

  //! @brief MLPREG reverse evaluation (constraint propagation)
  bool
  reval(T* valin, T* valout)
  {
    _set_dag();
    int const flag = _dag->reval(_codelist, _wkCPI, _varout.size(), _varout.data(),
                                 valout, _varin.size(), _varin.data(), valin,
                                 options.CPINF * T(-1, 1), options.CPMAX,
                                 options.CPTHRES);
    return flag >= 0;
  }
};

template <typename T>
inline bool
MLPREG<T>::_resize()
{
  _off.resize(_layers.size());
  _nwei       = 0;
  size_t nprv = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    _off[l] = _nwei;
    _nwei += (nprv + 1) * _layers[l].first;  // +1 for the bias
    nprv = _layers[l].first;
  }
  _nout   = (_layers.empty() ? 0 : _layers.back().first);
  _update = true;
  _uid    = ffhome_uid();
  return true;
}

template <typename T>
inline bool
MLPREG<T>::read_data(std::string const& filename, bool const disp)
{
  NETIO::t_Data dat;
  std::string msg;
  bool hasmodule = false;
#if defined(MC__USE_TORCH)
  torch::jit::script::Module module;
#endif

  if (NETIO::is_onnx(filename))
  {
    if (!NETIO::read_onnx(filename, dat, msg))
    {
      if (disp) std::cerr << "Error reading ONNX model: " << msg << "\n";
      return false;
    }
  }
  else
  {
#if defined(MC__USE_TORCH)
    try
    {
      module = torch::jit::load(filename);
      module.to(torch::kFloat64);
      module.eval();
      if (!NETIO::read_torch_net(module, dat, msg, disp))
      {
        if (disp) std::cerr << msg << "\n";
        return false;
      }
      // Torch evaluates the network with the weights of the operands where
      // every layer has a bias to write them into
      hasmodule = true;
      for (auto const& [weight, bias] : NETIO::torch_linear(module))
        if (!bias.defined()) hasmodule = false;
    }
    catch (std::exception const& e)  // c10::Error among them
    {
      if (disp)
        std::cerr << "Error reading Torch Module from " << filename << ":\n"
                  << e.what() << std::endl;
      return false;
    }
#else
    if (disp)
      std::cout << "Program compiled without Torch library: " << filename
                << " cannot be read.\n";
    return false;
#endif
  }

  // Structure and weights of the network read
  t_Layers layers;
  std::vector<double> wei;
  for (auto const& [layer, activ] : dat)
  {
    layers.push_back(std::make_pair(layer.size(), activ));
    for (auto const& neuron : layer) wei.insert(wei.end(), neuron.begin(), neuron.end());
  }
  MLPREG<T> reg;
  reg.options           = options;
  reg.options.EVALTORCH = false;
  if (!NETIO::check(dat, msg) ||
      !reg.set_structure(dat[0].first[0].size() - 1, layers))
  {
    if (disp) std::cerr << "Inconsistent network: " << msg << "\n";
    return false;
  }

#if defined(MC__USE_TORCH)
  // Check that the Module computes what was read from it, at two points: its
  // forward method may do more than go through 'net'
  if (!NETIO::is_onnx(filename))
  {
    try
    {
      size_t const nin = reg.nin(), nout = reg.nout();
      std::vector<double> valin(nin + wei.size()), valout(nout), refout(nout);
      std::copy(wei.begin(), wei.end(), valin.begin() + nin);
      torch::NoGradGuard nograd;
      for (size_t k = 0; k < 2; ++k)
      {
        for (size_t i = 0; i < nin; ++i) valin[i] = (k ? 0.5 * std::cos(i + 1.) : 0.);
        reg.eval(valin.data(), refout.data());
        _copy(module.forward({_tensor(valin.data(), nin)}).toTensor(), nout,
              valout.data());
        for (size_t j = 0; j < nout; ++j)
          if (!(std::fabs(valout[j] - refout[j]) <= 1e-8 * (1. + std::fabs(refout[j]))))
          {
            if (disp)
              std::cerr << "Module does not evaluate as its 'net' submodule: "
                        << valout[j] << " against " << refout[j] << "\n";
            return false;
          }
      }
    }
    catch (std::exception const& e)
    {
      if (disp) std::cerr << "Error evaluating Torch Module: " << e.what() << "\n";
      return false;
    }
  }
#endif

  // Commit
  set_structure(reg.nin(), layers);
  _wei       = wei;
  _hasmodule = hasmodule;
#if defined(MC__USE_TORCH)
  _module = module;
  _pool   = std::make_shared<TorchPool>();
#endif
  if (disp)
    std::cout << "Network read from " << filename << ": " << _nin << " inputs, "
              << _layers.size() << " layers, " << _nout << " outputs, " << _nwei
              << " weights" << (_hasmodule ? ", Torch Module" : "") << "\n";
  return true;
}

template <typename T>
inline void
MLPREG<T>::_set_dag()
{
  // Nothing to be updated: the same network, and the same forms of ReLU and
  // sigmoid activations
  std::pair<bool, bool> const dagopt(options.RELU2ABS, options.SIG2EXP);
  if (_dag && !_update && _dagopt == dagopt) return;
  _dagopt = dagopt;

  delete _dag;
  _dag = nullptr;
  _varin.clear();
  _varout.clear();
  _codelist.clear();
  _update = false;
  if (!_nin || !_nout) return;

  _dag = new FFGraph;
  _varin = _dag->add_vars(_nin, "X");
  std::vector<FFVar> const W = _dag->add_vars(_nwei, "W");
  _varin.insert(_varin.end(), W.begin(), W.end());

  std::vector<FFVar> cur, prv;
  size_t nprv = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    size_t const nneu = _layers[l].first;
    prv.swap(cur);
    cur.assign(nneu, FFVar());
    FFVar const* h = (l ? prv.data() : _varin.data());
    FFVar const* w = W.data() + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1)
    {
      FFVar val = w[0];  // bias term
      for (size_t j = 0; j < nprv; ++j) val += w[1 + j] * h[j];
      switch (_layers[l].second)
      {
        case LINEAR:
        default:
          break;
        case RELU:
          val = ReLU(val);
          break;
        case TANH:
          val = tanh(val);
          break;
        case SIGMOID:
          val = (options.SIG2EXP ? 1. / (exp(-val) + 1.) : tanh(val * 0.5) * 0.5 + 0.5);
          break;
      }
      cur[i] = val;
    }
    nprv = nneu;
  }
  _varout = cur;
}

template <typename T>
template <typename U>
inline void
MLPREG<T>::_eval(U const* valin, U* valout) const
{
  // Work arrays: values of the neurons of two consecutive layers
  static thread_local std::vector<U> wkcur, wkprv;
  U const* wei = valin + _nin;
  size_t nprv  = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    size_t const nneu = _layers[l].first;
    wkcur.swap(wkprv);
    wkcur.resize(nneu);
    U const* h = (l ? wkprv.data() : valin);
    U const* w = wei + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1)
    {
      U val = w[0];  // bias term
      for (size_t j = 0; j < nprv; ++j) val += h[j] * w[1 + j];
      _act(_layers[l].second, val);
      wkcur[i] = val;
    }
    nprv = nneu;
  }
  for (size_t k = 0; k < _nout; ++k) valout[k] = wkcur[k];
  // the values held by the work arrays may refer to others -- the nodes of a
  // backward AD tape, say: released
  if constexpr (!std::is_arithmetic<U>::value)
  {
    wkcur.clear();
    wkprv.clear();
  }
}

template <typename T>
template <typename U>
inline void
MLPREG<T>::_tang(U const* valin, U const* dirin, U* tanout, U* valout) const
{
  // Work arrays: values and tangents of the neurons of two consecutive layers
  static thread_local std::vector<U> wkval, wktan, wkvalprv, wktanprv;
  U const* wei  = valin + _nin;
  U const* dwei = dirin + _nin;
  size_t nprv   = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    size_t const nneu = _layers[l].first;
    wkval.swap(wkvalprv);
    wktan.swap(wktanprv);
    wkval.resize(nneu);
    wktan.resize(nneu);
    U const* h  = (l ? wkvalprv.data() : valin);
    U const* dh = (l ? wktanprv.data() : dirin);
    U const* w  = wei + _off[l];
    U const* dw = dwei + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1, dw += nprv + 1)
    {
      U val = w[0];   // bias term
      U tan = dw[0];  // and its tangent
      for (size_t j = 0; j < nprv; ++j)
      {
        val += h[j] * w[1 + j];
        tan += dh[j] * w[1 + j] + h[j] * dw[1 + j];
      }
      U der;
      if (_actder(_layers[l].second, val, der)) tan *= der;
      wkval[i] = val;
      wktan[i] = tan;
    }
    nprv = nneu;
  }
  for (size_t k = 0; k < _nout; ++k)
  {
    if (valout) valout[k] = wkval[k];
    tanout[k] = wktan[k];
  }
  if constexpr (!std::is_arithmetic<U>::value)
  {
    wkval.clear();
    wktan.clear();
    wkvalprv.clear();
    wktanprv.clear();
  }
}

template <typename T>
template <typename U>
inline void
MLPREG<T>::_adj(U const* valin, U const* wout, U* adjin, U* valout) const
{
  // Work arrays: values and activation derivatives of the neurons, layer by
  // layer, and the adjoints of two consecutive layers
  static thread_local std::vector<std::vector<U>> wkval, wkder;
  static thread_local std::vector<U> wkadj, wkprv;
  size_t const nlay = _layers.size();
  wkval.resize(nlay);
  wkder.resize(nlay);
  U const* wei = valin + _nin;

  // Forward pass, as in _eval
  size_t nprv = _nin;
  for (size_t l = 0; l < nlay; ++l)
  {
    size_t const nneu = _layers[l].first;
    wkval[l].resize(nneu);
    wkder[l].resize(nneu);
    U const* h = (l ? wkval[l - 1].data() : valin);
    U const* w = wei + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1)
    {
      U val = w[0];  // bias term
      for (size_t j = 0; j < nprv; ++j) val += h[j] * w[1 + j];
      _actder(_layers[l].second, val, wkder[l][i]);
      wkval[l][i] = val;
    }
    nprv = nneu;
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[nlay - 1][k];

  // Backward pass: the adjoints of the outputs of a layer give those of its
  // weights and of its inputs
  wkadj.assign(wout, wout + _nout);
  for (size_t l = nlay; l-- > 0;)
  {
    size_t const nneu = _layers[l].first;
    nprv              = (l ? _layers[l - 1].first : _nin);
    int const activ   = _layers[l].second;
    if (activ == RELU || activ == TANH || activ == SIGMOID)
      for (size_t i = 0; i < nneu; ++i) wkadj[i] *= wkder[l][i];
    U const* h = (l ? wkval[l - 1].data() : valin);
    U const* w = wei + _off[l];
    U* g       = adjin + _nin + _off[l];
    wkprv.assign(nprv, U(0.));
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1, g += nprv + 1)
    {
      g[0] = wkadj[i];
      for (size_t j = 0; j < nprv; ++j)
      {
        g[1 + j] = wkadj[i] * h[j];
        wkprv[j] += wkadj[i] * w[1 + j];
      }
    }
    wkadj.swap(wkprv);
  }
  for (size_t i = 0; i < _nin; ++i) adjin[i] = wkadj[i];
  if constexpr (!std::is_arithmetic<U>::value)
  {
    wkval.clear();
    wkder.clear();
    wkadj.clear();
    wkprv.clear();
  }
}

template <typename T>
template <typename U>
inline void
MLPREG<T>::_jacrev(U const* valin, U* gradout, U* valout) const
{
  // Work arrays: values and activation derivatives of the neurons, layer by
  // layer, and the adjoints of two consecutive layers for every output -- a
  // row per output
  static thread_local std::vector<std::vector<U>> wkval, wkder;
  static thread_local std::vector<U> wkadj, wkprv;
  size_t const nlay = _layers.size(), nvar = _nin + _nwei;
  wkval.resize(nlay);
  wkder.resize(nlay);
  U const* wei = valin + _nin;

  // Forward pass, as in _eval
  size_t nprv = _nin;
  for (size_t l = 0; l < nlay; ++l)
  {
    size_t const nneu = _layers[l].first;
    wkval[l].resize(nneu);
    wkder[l].resize(nneu);
    U const* h = (l ? wkval[l - 1].data() : valin);
    U const* w = wei + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1)
    {
      U val = w[0];  // bias term
      for (size_t j = 0; j < nprv; ++j) val += h[j] * w[1 + j];
      _actder(_layers[l].second, val, wkder[l][i]);
      wkval[l][i] = val;
    }
    nprv = nneu;
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[nlay - 1][k];

  // Backward pass: row k holds the adjoints of the outputs of a layer for
  // output k -- the identity to begin with -- and gives those of its weights
  // and of its inputs
  size_t ncur = _nout;
  wkadj.assign(_nout * _nout, U(0.));
  for (size_t k = 0; k < _nout; ++k) wkadj[k * _nout + k] = 1.;
  for (size_t l = nlay; l-- > 0;)
  {
    size_t const nneu = _layers[l].first;
    nprv              = (l ? _layers[l - 1].first : _nin);
    int const activ   = _layers[l].second;
    bool const nonlin = (activ == RELU || activ == TANH || activ == SIGMOID);
    U const* h        = (l ? wkval[l - 1].data() : valin);
    wkprv.assign(_nout * nprv, U(0.));
    for (size_t k = 0; k < _nout; ++k)
    {
      U const* adj = wkadj.data() + k * ncur;
      U* prv       = wkprv.data() + k * nprv;
      U* g         = gradout + k * nvar + _nin + _off[l];
      U const* w   = wei + _off[l];
      for (size_t i = 0; i < nneu; ++i, w += nprv + 1, g += nprv + 1)
      {
        U a = adj[i];
        if (nonlin) a *= wkder[l][i];
        g[0] = a;
        if constexpr (std::is_arithmetic<U>::value)
          if (a == 0.)  // a neuron that does not reach output k, or flat
          {
            for (size_t j = 0; j < nprv; ++j) g[1 + j] = 0.;
            continue;
          }
        for (size_t j = 0; j < nprv; ++j)
        {
          g[1 + j] = a * h[j];
          prv[j] += a * w[1 + j];
        }
      }
    }
    wkadj.swap(wkprv);
    ncur = nprv;
  }
  for (size_t k = 0; k < _nout; ++k)
    for (size_t i = 0; i < _nin; ++i) gradout[k * nvar + i] = wkadj[k * _nin + i];
  if constexpr (!std::is_arithmetic<U>::value)
  {
    wkval.clear();
    wkder.clear();
    wkadj.clear();
    wkprv.clear();
  }
}

template <typename T>
template <typename U>
inline void
MLPREG<T>::_jacfad(U const* valin, U* gradout, U* valout) const
{
  size_t const nvar = _nin + _nwei;
  std::vector<FADType<U>> vFin(nvar), vFout(_nout);
  for (size_t i = 0; i < nvar; ++i)
  {
    vFin[i] = valin[i];
    vFin[i].diff(i, nvar);
  }
  _eval(vFin.data(), vFout.data());
  for (size_t k = 0; k < _nout; ++k)
  {
    if (valout) valout[k] = vFout[k].x();
    for (size_t i = 0; i < nvar; ++i) gradout[k * nvar + i] = vFout[k].d(i);
  }
}

template <typename T>
inline void
MLPREG<T>::grad(double const* valin, double* gradout, double* const valout)
{
#if defined(MC__USE_TORCH)
  if (_torch())
  {
    // One backward pass per output, into the inputs and the parameters
    TorchSlot& slot = _slot();
    torch::Tensor tensin;
    torch::Tensor tensout = _torch_forward(slot, valin, true, &tensin);
    if (valout) _copy(tensout, _nout, valout);
    std::vector<torch::Tensor> inputs(1, tensin);
    inputs.insert(inputs.end(), slot.param.begin(), slot.param.end());
    size_t const nvar = _nin + _nwei;
    for (size_t k = 0; k < _nout; ++k)
      _torch_row(torch::autograd::grad({tensout[k]}, inputs, {},
                                       /*retain_graph=*/(k + 1 < _nout),
                                       /*create_graph=*/false,
                                       /*allow_unused=*/true),
                 gradout + k * nvar);
    return;
  }
#endif
  // The Jacobian by passes through the layers: the reverse mode costs a
  // pass per output, the forward mode a pass per operand
  size_t const nvar = _nin + _nwei;
  if (_nout <= nvar || options.AUTODIFF == Options::AD_TYPE::B)
  {
    _jacrev(valin, gradout, valout);
    return;
  }
  std::vector<double> dir(nvar, 0.), tan(_nout);
  for (size_t i = 0; i < nvar; ++i)
  {
    dir[i] = 1.;
    _tang(valin, dir.data(), tan.data(), i ? nullptr : valout);
    for (size_t k = 0; k < _nout; ++k) gradout[k * nvar + i] = tan[k];
    dir[i] = 0.;
  }
}

template <typename T>
inline void
MLPREG<T>::tang(double const* valin, double const* dirin, double* tanout,
                double* const valout)
{
#if defined(MC__USE_TORCH)
  if (_torch())
  {
    // The product of the Jacobian with the direction
    size_t const nvar = _nin + _nwei;
    std::vector<double> J(_nout * nvar);
    grad(valin, J.data(), valout);
    for (size_t k = 0; k < _nout; ++k)
    {
      tanout[k] = 0.;
      for (size_t i = 0; i < nvar; ++i) tanout[k] += J[k * nvar + i] * dirin[i];
    }
    return;
  }
#endif
  _tang(valin, dirin, tanout, valout);
}

template <typename T>
inline void
MLPREG<T>::adj(double const* valin, double const* wout, double* adjin,
               double* const valout)
{
#if defined(MC__USE_TORCH)
  if (_torch())
  {
    // One backward pass from the outputs, with the weights as their gradient
    TorchSlot& slot = _slot();
    torch::Tensor tensin;
    torch::Tensor tensout = _torch_forward(slot, valin, true, &tensin);
    if (valout) _copy(tensout, _nout, valout);
    std::vector<torch::Tensor> inputs(1, tensin);
    inputs.insert(inputs.end(), slot.param.begin(), slot.param.end());
    _torch_row(torch::autograd::grad({tensout}, inputs, {_tensor(wout, _nout)},
                                     /*retain_graph=*/false,
                                     /*create_graph=*/false,
                                     /*allow_unused=*/true),
               adjin);
    return;
  }
#endif
  _adj(valin, wout, adjin, valout);
}

#if defined(MC__USE_TORCH)
template <typename T>
inline typename MLPREG<T>::TorchSlot&
MLPREG<T>::_slot() const
{
  std::lock_guard<std::mutex> lock(_pool->mtx);
  std::unique_ptr<TorchSlot>& slot = _pool->slots[std::this_thread::get_id()];
  if (!slot)
  {
    slot         = std::make_unique<TorchSlot>();
    slot->module = _module.clone();
    slot->module.eval();
    slot->lin = NETIO::torch_linear(slot->module);
    for (auto const& [weight, bias] : slot->lin)
    {
      slot->param.push_back(weight);
      slot->param.push_back(bias);
    }
    if (slot->lin.size() != _layers.size())
      throw std::runtime_error("MLPREG: Torch Module does not match the structure\n");
  }
  return *slot;
}

template <typename T>
inline void
MLPREG<T>::_torch_load(TorchSlot& slot, double const* wei, bool const grad) const
{
  torch::NoGradGuard nograd;
  static thread_local std::vector<double> wkW, wkB;
  size_t nprv = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    size_t const nneu = _layers[l].first;
    wkW.resize(nneu * nprv);
    wkB.resize(nneu);
    double const* w = wei + _off[l];
    for (size_t i = 0; i < nneu; ++i, w += nprv + 1)
    {
      wkB[i] = w[0];
      std::copy(w + 1, w + 1 + nprv, wkW.data() + i * nprv);
    }
    auto& [weight, bias] = slot.lin[l];
    weight.copy_(torch::from_blob(wkW.data(), {(int64_t)nneu, (int64_t)nprv},
                                  torch::TensorOptions().dtype(torch::kFloat64)));
    bias.copy_(torch::from_blob(wkB.data(), {(int64_t)nneu},
                                torch::TensorOptions().dtype(torch::kFloat64)));
    nprv = nneu;
  }
  for (auto& p : slot.param) p.set_requires_grad(grad);
}

template <typename T>
inline torch::Tensor
MLPREG<T>::_torch_forward(TorchSlot& slot, double const* valin, bool const grad,
                          torch::Tensor* tensin) const
{
  _torch_load(slot, valin + _nin, grad);
  torch::Tensor tin = _tensor(valin, _nin);
  if (grad) tin.set_requires_grad(true);
  torch::Tensor tout = slot.module.forward({tin}).toTensor();
  if (!tout.defined() || (size_t)tout.numel() != _nout)
    throw std::runtime_error("MLPREG: Torch Module returns " +
                             std::to_string(tout.defined() ? tout.numel() : 0) +
                             " values where " + std::to_string(_nout) +
                             " are expected\n");
  if (tensin) *tensin = tin;
  return tout.view({-1});
}

template <typename T>
inline void
MLPREG<T>::_torch_row(std::vector<torch::Tensor> const& g, double* row) const
{
  auto get = [](torch::Tensor const& t, size_t const n, std::vector<double>& v)
  {
    v.assign(n, 0.);
    if (!t.defined()) return;  // unused: zero
    torch::Tensor const td = t.detach().to(torch::kCPU, torch::kFloat64).contiguous();
    std::copy(td.data_ptr<double>(), td.data_ptr<double>() + n, v.data());
  };
  static thread_local std::vector<double> gW, gB;
  get(g[0], _nin, gB);
  std::copy(gB.begin(), gB.end(), row);
  size_t nprv = _nin;
  for (size_t l = 0; l < _layers.size(); ++l)
  {
    size_t const nneu = _layers[l].first;
    get(g[1 + 2 * l], nneu * nprv, gW);
    get(g[2 + 2 * l], nneu, gB);
    double* r = row + _nin + _off[l];
    for (size_t i = 0; i < nneu; ++i, r += nprv + 1)
    {
      r[0] = gB[i];
      std::copy(gW.data() + i * nprv, gW.data() + (i + 1) * nprv, r + 1);
    }
    nprv = nneu;
  }
}
#endif

//! @brief C++ base class of the external DAG operations on an mc::MLPREG
////////////////////////////////////////////////////////////////////////
//! mc::FFMLPREGOp is the base class of mc::FFMLPREG, mc::FFGradMLPREG,
//! mc::FFTanMLPREG and mc::FFAdjMLPREG: the home of the MLPREG in a DAG
//! and the options of the relaxations (mc::FFEXTERN), the definition of
//! the operations, and their evaluations other than numerical ones. The
//! class OP of the operation gives its numbers of results and operands
//! (nres, nopd), its numerical evaluation (evalnum) and its name (tag).
////////////////////////////////////////////////////////////////////////
template <typename T, typename OP>
class FFMLPREGOp
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, mc::MLPREG<T>>
{
 protected:
  typedef FFEXTERN<T, mc::MLPREG<T>> t_Base;
  using t_Base::_ptrObj;

  OP const&
  _op() const
  {
    return static_cast<OP const&>(*this);
  }

  void
  _check(size_t const nRes, size_t const nVar) const
  {
#ifdef MC__FFMLPREG_CHECK
    if (!_ptrObj || nRes != OP::nres(*_ptrObj) || nVar != OP::nopd(*_ptrObj))
      throw std::runtime_error(
          "FFMLPREG: " + std::to_string(nVar) + " operands and " +
          std::to_string(nRes) + " results for " +
          std::to_string(_ptrObj ? OP::nopd(*_ptrObj) : 0) + " operands and " +
          std::to_string(_ptrObj ? OP::nres(*_ptrObj) : 0) + " results\n");
#endif
  }

  // set the object and related operation in DAG: home of the MLPREG in the
  // DAG of the operands, made on its first insertion there, with a copy of
  // the MLPREG (policy > 0), found again afterwards -- see mc::FFHOME::_attach
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, mc::MLPREG<T>* pMLP, int const policy)
  {
    if (!pMLP || nVar != OP::nopd(*pMLP))
      throw std::runtime_error("FFMLPREG: " + std::to_string(nVar) +
                               " operands where the MLPREG takes " +
                               std::to_string(pMLP ? OP::nopd(*pMLP) : 0) + "\n");
    this->owndata = false;
    this->_attach(this->_dag_of(nVar, pVar), pMLP, policy);
    return this->insert_external_operation(_op(), OP::nres(*_ptrObj), nVar, pVar);
  }

  // set the operation in DAG on the MLPREG of another operation, whose home
  // it shares
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, t_Base const& parent)
  {
    this->owndata = false;
    this->_join(parent);
    return this->insert_external_operation(_op(), OP::nres(*_ptrObj), nVar, pVar);
  }

 public:
  //! @brief Enumeration type for copy policy of MLPREG object
  enum POLICY_TYPE
  {
    SHALLOW  = 0,  //!< The DAG refers to the MLPREG object, which it does not own
    COPY     = 1,  //!< The DAG holds a copy of the MLPREG object: one per DAG
    TRANSFER = -1  //!< The DAG refers to the MLPREG object, and deletes it in the end
  };

  FFMLPREGOp() : t_Base() {}
  FFMLPREGOp(FFMLPREGOp const& other) : t_Base(other) {}
  virtual ~FFMLPREGOp() {}

  //! @brief MLPREG object of the DAG: after an insertion, the copy that the
  //! DAG holds (policy COPY) or the object given, shared by every operation
  //! on that structure in that DAG -- see mc::FFEXTERN::object. Throws before
  //! any insertion
  mc::MLPREG<T>&
  MLPREG() const
  {
    return this->object();
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type(idU, typeid(FFVar)))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFDep)))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFInv)))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(FFExpr)))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(SLiftVar)))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(double)))
      return evalnum(nRes, static_cast<double*>(vRes), nVar,
                     static_cast<double const*>(vVar));
    else if (mc::same_type(idU, typeid(FADType<double>)))
      return evalnum(nRes, static_cast<FADType<double>*>(vRes), nVar,
                     static_cast<FADType<double> const*>(vVar));
    else if (mc::same_type(idU, typeid(T)))
      return evalnum(nRes, static_cast<T*>(vRes), nVar, static_cast<T const*>(vVar));
    else if (mc::same_type(idU, typeid(McCormick<T>)))
      return evalnum(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                     static_cast<McCormick<T> const*>(vVar));
    else if (mc::same_type(idU, typeid(Specbnd<T>)))
      return evalnum(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                     static_cast<Specbnd<T> const*>(vVar));
    else if (mc::same_type(idU, typeid(SCVar<T>)))
      return evalnum(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                     static_cast<SCVar<T> const*>(vVar));
    else if (mc::same_type(idU, typeid(SupVar<PWCU>)))
      return evalnum(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                     static_cast<SupVar<PWCU> const*>(vVar));
    else if (mc::same_type(idU, typeid(SupVar<PWLU>)))
      return evalnum(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                     static_cast<SupVar<PWLU> const*>(vVar));

    throw std::runtime_error(_op().name() +
                             "::feval: **ERROR** No evaluation method with type " +
                             std::string(idU.name()) + "\n");
  }

  //! @brief Numerical evaluation, in any arithmetic of the MLPREG
  template <typename U>
  void
  evalnum(size_t const nRes, U* vRes, size_t const nVar, U const* vVar) const
  {
    _check(nRes, nVar);
    _op().evalobj(*_ptrObj, vRes, vVar);
  }

  void
  eval(size_t const nRes, FFVar* vRes, size_t const nVar, FFVar const* vVar,
       unsigned const* mVar) const
  {
    _check(nRes, nVar);
    // in the DAG of the operands, with the home of the MLPREG there
    FFVar** ppRes = this->_reinsert(_op(), nRes, nVar, vVar);
    for (size_t j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
  }

  void
  eval(size_t const nRes, FFDep* vRes, size_t const nVar, FFDep const* vVar,
       unsigned const* mVar) const
  {
    _check(nRes, nVar);
    vRes[0] = 0;
    for (size_t i = 0; i < nVar; ++i) vRes[0] += vVar[i];
    vRes[0].update(FFDep::TYPE::N);
    for (size_t j = 1; j < nRes; ++j) vRes[j] = vRes[0];
  }

  void
  eval(size_t const nRes, FFInv* vRes, size_t const nVar, FFInv const* vVar,
       unsigned const* mVar) const
  {
    _check(nRes, nVar);
    vRes[0] = 0;
    for (size_t i = 0; i < nVar; ++i) vRes[0] += vVar[i];
    vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
    for (size_t j = 1; j < nRes; ++j) vRes[j] = vRes[0];
  }

  void
  eval(size_t const nRes, FFExpr* vRes, size_t const nVar, FFExpr const* vVar,
       unsigned const* mVar) const
  {
    switch (FFExpr::options.LANG)
    {
      case FFExpr::Options::DAG:
        for (size_t j = 0; j < nRes; ++j)
        {
          std::ostringstream os;
          os << _op().name() << "[" << j << "]";
          vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
        }
        break;
      case FFExpr::Options::GAMS:
      default:
        throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
    }
  }

  void
  eval(size_t const nRes, SLiftVar* vRes, size_t const nVar, SLiftVar const* vVar,
       unsigned const* mVar) const
  {
    _check(nRes, nVar);
    vVar->env()->lift(nRes, vRes, nVar, vVar);
  }

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(_op().name() +
                             "::reval: **ERROR** No evaluation method with type " +
                             std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return std::string(OP::tag()) + "[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T>
class FFGradMLPREG;
template <typename T>
class FFTanMLPREG;
template <typename T>
class FFAdjMLPREG;

//! @brief C++ class defining multilayer perceptrons whose weights are
//! variables as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLPREG is a C++ class for defining an mc::MLPREG as an external
//! DAG operation in MC++: its operands are the inputs of the network
//! followed by its weights, and its results are the outputs of the
//! network. As mc::FFMLP, it derives from mc::FFEXTERN: the MLPREG has
//! one home per DAG, and the relaxation of the operation follows the
//! options of the handle (options) -- interval bounds, McCormick or
//! polyhedral relaxations, and so on; its differentiation inserts the
//! Jacobian (mc::FFGradMLPREG), or with mc::MLPREG::Options::TANGENT and
//! ADJOINT the directional derivative (mc::FFTanMLPREG) and the adjoint
//! (mc::FFAdjMLPREG).
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class FFMLPREG
    ////////////////////////////////////////////////////////////////////////
    : public FFMLPREGOp<T, FFMLPREG<T>>
{
  friend class FFMLPREGOp<T, FFMLPREG<T>>;
  typedef FFMLPREGOp<T, FFMLPREG<T>> t_Op;
  using t_Op::_ptrObj;
  using t_Op::_check;

 public:
  using typename t_Op::POLICY_TYPE;
  using t_Op::COPY;
  using t_Op::SHALLOW;
  using t_Op::TRANSFER;
  using t_Op::eval;
  using t_Op::evalnum;
  using t_Op::reval;
  using t_Op::name;

  static size_t
  nres(mc::MLPREG<T> const& reg)
  {
    return reg.nout();
  }
  static size_t
  nopd(mc::MLPREG<T> const& reg)
  {
    return reg.nvar();
  }
  static char const*
  tag()
  {
    return "MLPREG";
  }
  template <typename U>
  static void
  evalobj(mc::MLPREG<T>& reg, U* vRes, U const* vVar)
  {
    reg.eval(vVar, vRes);
  }

  FFMLPREG() : t_Op() {}
  FFMLPREG(FFMLPREG const& other) : t_Op(other) {}
  virtual ~FFMLPREG() {}

  //! @brief Define the operation on the inputs vIn and the weights vWei:
  //! the outputs of the network
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vIn, std::vector<FFVar> const& vWei,
             mc::MLPREG<T>* pMLP, int const policy = COPY)
  {
    std::vector<FFVar> vOpd(vIn);
    vOpd.insert(vOpd.end(), vWei.begin(), vWei.end());
    FFVar** ppRes = this->_set(vOpd.size(), vOpd.data(), pMLP, policy);
    std::vector<FFVar> vRes(_ptrObj->nout());
    for (size_t k = 0; k < vRes.size(); ++k) vRes[k] = *ppRes[k];
    return vRes;
  }

  //! @brief Define the operation on the inputs vIn and the weights vWei:
  //! output idep of the network
  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vIn,
             std::vector<FFVar> const& vWei, mc::MLPREG<T>* pMLP,
             int const policy = COPY)
  {
    std::vector<FFVar> vOpd(vIn);
    vOpd.insert(vOpd.end(), vWei.begin(), vWei.end());
    if (!pMLP || idep >= pMLP->nout())
      throw std::runtime_error("FFMLPREG: no output " + std::to_string(idep) + "\n");
    return *(this->_set(vOpd.size(), vOpd.data(), pMLP, policy)[idep]);
  }

  //! @brief Define the operation on the nIn inputs pIn and the nWei weights
  //! pWei: the outputs of the network
  FFVar**
  operator()(size_t const nIn, FFVar const* pIn, size_t const nWei,
             FFVar const* pWei, mc::MLPREG<T>* pMLP, int const policy = COPY)
  {
    std::vector<FFVar> vOpd(pIn, pIn + nIn);
    vOpd.insert(vOpd.end(), pWei, pWei + nWei);
    return this->_set(vOpd.size(), vOpd.data(), pMLP, policy);
  }

  //! @brief Define the operation on the nIn inputs pIn and the nWei weights
  //! pWei: output idep of the network
  FFVar&
  operator()(size_t const idep, size_t const nIn, FFVar const* pIn,
             size_t const nWei, FFVar const* pWei, mc::MLPREG<T>* pMLP,
             int const policy = COPY)
  {
    if (!pMLP || idep >= pMLP->nout())
      throw std::runtime_error("FFMLPREG: no output " + std::to_string(idep) + "\n");
    return *((*this)(nIn, pIn, nWei, pWei, pMLP, policy)[idep]);
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type(idU, typeid(FADType<FFVar>)))
      return eval(nRes, static_cast<FADType<FFVar>*>(vRes), nVar,
                  static_cast<FADType<FFVar> const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(PolVar<T>)))
      return eval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                  static_cast<PolVar<T> const*>(vVar), mVar);
    else if (mc::same_type(idU, typeid(BADType<double>)))
      return evalbad(nRes, static_cast<BADType<double>*>(vRes), nVar,
                     static_cast<BADType<double> const*>(vVar));
    return t_Op::feval(idU, nRes, vRes, nVar, vVar, mVar);
  }

  //! @brief Evaluation in the reverse mode of AD: each output is a single
  //! node whose partial derivatives with respect to the inputs and the
  //! weights are the Jacobian of the network, computed in double precision by
  //! passes through the layers (see MLPREG::grad), in place of a tape of every
  //! neuron. The gradient of a training loss in the weights then costs one
  //! evaluation and one reverse pass per sample
  void
  evalbad(size_t const nRes, BADType<double>* vRes, size_t const nVar,
          BADType<double> const* vVar) const
  {
    this->_check(nRes, nVar);
    std::vector<double> x(nVar), y(nRes), J(nRes * nVar);
    for (size_t i = 0; i < nVar; ++i) x[i] = vVar[i].val();
    this->_ptrObj->grad(x.data(), J.data(), y.data());
    for (size_t k = 0; k < nRes; ++k)
      vRes[k] = bad_node(y[k], nVar, vVar, J.data() + k * nVar);
  }

  void eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
            FADType<FFVar> const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
            PolVar<T> const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    if (mc::same_type(idU, typeid(T)))
      return reval(nRes, static_cast<T*>(vRes), nVar, static_cast<T*>(vVar));
    else if (mc::same_type(idU, typeid(PolVar<T>)))
      return reval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                   static_cast<PolVar<T>*>(vVar));
    return t_Op::reval(idU, nRes, vRes, nVar, vVar);
  }

  bool reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const;

  bool reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
             PolVar<T>* vVar) const;

  // Derivatives: the Jacobian operation mc::FFGradMLPREG
  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
             FFVar const* vVar, FFVar** vDer) const;

  // Tangents: with mc::MLPREG::Options::TANGENT, the forward symbolic
  // differentiation inserts the directional derivative operation
  // FFTanMLPREG along the tangents of the operands, in place of the Jacobian
  // FFGradMLPREG. With the weights among the operands, the Jacobian has
  // nout*nvar entries -- more than the weights for every output -- and the
  // products that the differentiation forms with them in each sweep cost
  // more than a tangent pass, which costs about two evaluations of the
  // network: the passes are the cheaper for any number of sweeps
  bool
  has_tangent(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.TANGENT;
  }
  void tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vTan, FFVar* vResTan) const;

  // Adjoints: with mc::MLPREG::Options::ADJOINT, the backward symbolic
  // differentiation inserts the adjoint operation FFAdjMLPREG for the
  // adjoints of the results, in place of the Jacobian FFGradMLPREG -- the
  // cheaper for any number of sweeps, as for the tangents
  bool
  has_adjoint(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.ADJOINT;
  }
  void adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vResAdj, FFVar* vVarAdj) const;
};

//! @brief C++ class defining the Jacobian of an mc::MLPREG as an external
//! DAG operation in MC++: its operands are those of mc::FFMLPREG, and its
//! results are the derivatives of the outputs in the operands, that of
//! output k in operand i at position k*nvar+i
template <typename T = mc::Interval>
class FFGradMLPREG
    ////////////////////////////////////////////////////////////////////////
    : public FFMLPREGOp<T, FFGradMLPREG<T>>
{
  friend class FFMLPREG<T>;
  friend class FFMLPREGOp<T, FFGradMLPREG<T>>;
  typedef FFMLPREGOp<T, FFGradMLPREG<T>> t_Op;

 public:
  static size_t
  nres(mc::MLPREG<T> const& reg)
  {
    return reg.nout() * reg.nvar();
  }
  static size_t
  nopd(mc::MLPREG<T> const& reg)
  {
    return reg.nvar();
  }
  static char const*
  tag()
  {
    return "GradMLPREG";
  }
  template <typename U>
  static void
  evalobj(mc::MLPREG<T>& reg, U* vRes, U const* vVar)
  {
    reg.grad(vVar, vRes);
  }

  FFGradMLPREG() : t_Op() {}
  FFGradMLPREG(FFGradMLPREG const& other) : t_Op(other) {}
  virtual ~FFGradMLPREG() {}

  //! @brief Define the operation on the operands vVar (inputs, then
  //! weights): the Jacobian, row-major
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, mc::MLPREG<T>* pMLP,
             int const policy = t_Op::COPY)
  {
    FFVar** ppRes = this->_set(vVar.size(), vVar.data(), pMLP, policy);
    std::vector<FFVar> vRes(nres(*pMLP));
    for (size_t k = 0; k < vRes.size(); ++k) vRes[k] = *ppRes[k];
    return vRes;
  }
};

//! @brief C++ class defining the directional derivative of an mc::MLPREG as
//! an external DAG operation in MC++: its operands are those of
//! mc::FFMLPREG followed by a direction of the same size, and its results
//! are the derivatives of the outputs along that direction, evaluated in a
//! single forward pass. It is inserted by the forward symbolic
//! differentiation of mc::FFMLPREG when mc::MLPREG::Options::TANGENT is
//! set, and can also be defined directly. It has no derivative of its own.
template <typename T = mc::Interval>
class FFTanMLPREG
    ////////////////////////////////////////////////////////////////////////
    : public FFMLPREGOp<T, FFTanMLPREG<T>>
{
  friend class FFMLPREG<T>;
  friend class FFMLPREGOp<T, FFTanMLPREG<T>>;
  typedef FFMLPREGOp<T, FFTanMLPREG<T>> t_Op;

 public:
  static size_t
  nres(mc::MLPREG<T> const& reg)
  {
    return reg.nout();
  }
  static size_t
  nopd(mc::MLPREG<T> const& reg)
  {
    return 2 * reg.nvar();
  }
  static char const*
  tag()
  {
    return "TanMLPREG";
  }
  template <typename U>
  static void
  evalobj(mc::MLPREG<T>& reg, U* vRes, U const* vVar)
  {
    reg.tang(vVar, vVar + reg.nvar(), vRes);
  }

  FFTanMLPREG() : t_Op() {}
  FFTanMLPREG(FFTanMLPREG const& other) : t_Op(other) {}
  virtual ~FFTanMLPREG() {}

  //! @brief Define the operation: derivatives of the outputs at the operands
  //! vVar (inputs, then weights) along vDir
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vDir,
             mc::MLPREG<T>* pMLP, int const policy = t_Op::COPY)
  {
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vDir.begin(), vDir.end());
    FFVar** ppRes = this->_set(vOpd.size(), vOpd.data(), pMLP, policy);
    std::vector<FFVar> vRes(nres(*pMLP));
    for (size_t k = 0; k < vRes.size(); ++k) vRes[k] = *ppRes[k];
    return vRes;
  }
};

//! @brief C++ class defining the adjoint of an mc::MLPREG as an external DAG
//! operation in MC++: its operands are those of mc::FFMLPREG followed by
//! one weight per output, and its results are the derivatives of the
//! weighted sum of the outputs in the operands -- the product of the
//! transposed Jacobian with the weights -- evaluated in a single backward
//! pass. It is inserted by the backward symbolic differentiation of
//! mc::FFMLPREG when mc::MLPREG::Options::ADJOINT is set, and can also be
//! defined directly. It has no derivative of its own.
template <typename T = mc::Interval>
class FFAdjMLPREG
    ////////////////////////////////////////////////////////////////////////
    : public FFMLPREGOp<T, FFAdjMLPREG<T>>
{
  friend class FFMLPREG<T>;
  friend class FFMLPREGOp<T, FFAdjMLPREG<T>>;
  typedef FFMLPREGOp<T, FFAdjMLPREG<T>> t_Op;

 public:
  static size_t
  nres(mc::MLPREG<T> const& reg)
  {
    return reg.nvar();
  }
  static size_t
  nopd(mc::MLPREG<T> const& reg)
  {
    return reg.nvar() + reg.nout();
  }
  static char const*
  tag()
  {
    return "AdjMLPREG";
  }
  template <typename U>
  static void
  evalobj(mc::MLPREG<T>& reg, U* vRes, U const* vVar)
  {
    reg.adj(vVar, vVar + reg.nvar(), vRes);
  }

  FFAdjMLPREG() : t_Op() {}
  FFAdjMLPREG(FFAdjMLPREG const& other) : t_Op(other) {}
  virtual ~FFAdjMLPREG() {}

  //! @brief Define the operation: derivatives in the operands vVar (inputs,
  //! then weights) of the sum of the outputs weighted by vWei
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vWei,
             mc::MLPREG<T>* pMLP, int const policy = t_Op::COPY)
  {
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vWei.begin(), vWei.end());
    FFVar** ppRes = this->_set(vOpd.size(), vOpd.data(), pMLP, policy);
    std::vector<FFVar> vRes(nres(*pMLP));
    for (size_t k = 0; k < vRes.size(); ++k) vRes[k] = *ppRes[k];
    return vRes;
  }
};

template <typename T>
inline void
FFMLPREG<T>::eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
                  FADType<FFVar> const* vVar, unsigned const* mVar) const
{
  _check(nRes, nVar);
  std::vector<FFVar> vVarVal(nVar);
  for (size_t i = 0; i < nVar; ++i) vVarVal[i] = vVar[i].val();
  FFVar const* const* vResVal = this->_reinsert(*this, nRes, nVar, vVarVal.data());
  FFGradMLPREG<T> ResDer;
  FFVar const* const* vResDer = ResDer._set(nVar, vVarVal.data(), *this);

  for (size_t k = 0; k < nRes; ++k)
  {
    vRes[k] = *vResVal[k];
    for (size_t i = 0; i < nVar; ++i) vRes[k].setDepend(vVar[i]);
    for (size_t j = 0; j < vRes[k].size(); ++j)
    {
      vRes[k][j] = 0.;
      for (size_t i = 0; i < nVar; ++i)
      {
        if (vVar[i][j].cst() && vVar[i][j].num().val() == 0.) continue;
        vRes[k][j] += *vResDer[k * nVar + i] * vVar[i][j];
      }
    }
  }
}

template <typename T>
inline void
FFMLPREG<T>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                   FFVar const* vVar, FFVar** vDer) const
{
  _check(nRes, nVar);
  FFGradMLPREG<T> ResDer;
  FFVar const* const* vResDer = ResDer._set(nVar, vVar, *this);
  // the Jacobian is row-major: the derivative of result k in operand i at
  // position k*nVar+i
  for (unsigned k = 0; k < nRes; ++k)
    for (unsigned i = 0; i < nVar; ++i) vDer[k][i] = *vResDer[k * nVar + i];
}

template <typename T>
inline void
FFMLPREG<T>::tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                     FFVar const* vVar, FFVar const* const* vTan,
                     FFVar* vResTan) const
{
  _check(nRes, nVar);
  // Operands of the directional derivative: the operands, then their
  // tangents -- zero where the tangent of an operand is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(2 * nVar);
  for (unsigned i = 0; i < nVar; ++i) vOpd.push_back(vTan[i] ? *vTan[i] : FFVar(0.));
  FFTanMLPREG<T> ResTan;
  FFVar const* const* vResDer = ResTan._set(vOpd.size(), vOpd.data(), *this);
  for (unsigned k = 0; k < nRes; ++k) vResTan[k] = *vResDer[k];
}

template <typename T>
inline void
FFMLPREG<T>::adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                     FFVar const* vVar, FFVar const* const* vResAdj,
                     FFVar* vVarAdj) const
{
  _check(nRes, nVar);
  // Operands of the adjoint operation: the operands, then the adjoints of the
  // results -- zero where the adjoint of a result is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(nVar + nRes);
  for (unsigned k = 0; k < nRes; ++k)
    vOpd.push_back(vResAdj[k] ? *vResAdj[k] : FFVar(0.));
  FFAdjMLPREG<T> ResAdj;
  FFVar const* const* vResDer = ResAdj._set(vOpd.size(), vOpd.data(), *this);
  for (unsigned i = 0; i < nVar; ++i) vVarAdj[i] = *vResDer[i];
}

template <typename T>
inline void
FFMLPREG<T>::eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                  PolVar<T> const* vVar, unsigned const* mVar) const
{
  _check(nRes, nVar);
  PolImg<T>* img = vVar[0].image();
  FFBase* dag    = vVar[0].var().dag();
#ifdef MC__FFMLPREG_CHECK
  assert(img && dag);
#endif
  FFVar** ppRes = dag->curOp()->varout.data();
#ifdef MC__FFMLPREG_CHECK
  assert(nRes == dag->curOp()->varout.size());
#endif
  this->_resize_relax(img);
  this->_propagate_relax(img, ppRes, vRes, vVar);
}

template <typename T>
inline bool
FFMLPREG<T>::reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                   PolVar<T>* vVar) const
{
  _check(nRes, nVar);
  PolImg<T>* img = vVar[0].image();
  FFOp* pop      = vVar[0].var().opdef().first;
#ifdef MC__FFMLPREG_CHECK
  assert(img && pop);
#endif
  this->_backpropagate_relax(img, pop, vRes, vVar);
  return true;
}

template <typename T>
inline bool
FFMLPREG<T>::reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const
{
  _check(nRes, nVar);
  return _ptrObj->reval(vVar, vRes);
}

}  // end namespace mc

#endif
