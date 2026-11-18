// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__FFMLP_HPP
#define MC__FFMLP_HPP

#include <algorithm>
#if defined(MC__USE_TORCH)
#include <torch/script.h>  // One-stop header for TorchScript
                           // #include <ATen/cuda/CUDAContext.h>
#include <torch/autograd.h>  // torch::autograd::grad, forward-mode AD
#endif

#include "ffextern.hpp"
#include "fflin.hpp"
#include "interval.hpp"

// #define MC__FFMLP_DEBUG
#define MC__FFMLP_CHECK

namespace mc
{

//! @brief C++ class for evaluation and relaxation of multilayer perceptrons
////////////////////////////////////////////////////////////////////////
//! mc::MLP is a C++ class for evaluation and relaxation of (trained)
//! multilayer perceptrons (MLP) that leverages expression trees and
//! arithmetics available through MC++
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class MLP
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

  //! @brief Reset MLP data
  void reset_data();

  //! @brief Set MLP data
  bool set_data(
      std::vector<std::pair<std::vector<std::vector<double>>, int>> const& mlp);

  //! @brief Read MLP data from TorchScript
  bool read_data(std::string const& filename, bool const disp = false);

  //! @brief Append multi-neuron layer to MLP data: outer vector size is
  //! number of neurons in hidden layer; inner vector has scalar weights multiplying each
  //! neuron in previous (input or hidden) layer and a bias term and size
  //! 1 + the number of neurons in previous layer
  bool append_data(std::vector<std::vector<double>> const& layer,
                   int const activ = LINEAR, bool const reset = false);

  //! @brief Append single-neuron layer to MLP data: vector has scalar weights
  //! multiplying each neuron in previous (input or hidden) layer and a bias
  //! term and size 1 + the number of neurons in previous layer
  bool append_data(std::vector<double> const& layer, int const activ = LINEAR,
                   bool const reset = false);

  //! @brief MLP options
  struct Options
  {
    //! @brief Default constructor
    Options() { reset(); }

    //! @brief Assignment operator
    Options&
    operator=(Options const& other)
    {
      EVALTORCH = other.EVALTORCH;
      ZEROTOL   = other.ZEROTOL;
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
      ZEROTOL   = DBL_EPSILON;
      RELU2ABS  = false;
      SIG2EXP   = false;
      AUTODIFF  = F;
      TANGENT   = false;
      ADJOINT   = false;
      CPMAX =
          1;  // Leave it to the outer DAG to iterate as necessary by default
      CPTHRES = 1e4 * DBL_EPSILON;
      CPINF   = 1e30;
    }

    //! @brief Enumeration type for AD strategy
    enum AD_TYPE
    {
      F = 0,  //!< Forward differentiation
      B       //!< Backward differentiation
    };

    //! @brief Whether to evaluate ANN in floating-point arithemic using Torch,
    //! in double precision whatever the precision of the Module that was read
    bool EVALTORCH;
    //! @brief Threshold for zero coefficient in neural network
    double ZEROTOL;
    //! @brief Whether to convert ReLU to abs (true) or max (false) function
    bool RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false)
    bool SIG2EXP;
    //! @brief Whether to apply forward or reverse automatic differentiation
    //! for the Jacobian (grad). In double arithmetic the reverse mode is also
    //! applied with the default (F) where the outputs are no more than the
    //! inputs, both modes giving the same Jacobian to rounding errors
    int AUTODIFF;
    //! @brief Whether the forward symbolic differentiation of an MLP operation
    //! in a DAG (FAD, SFAD, DFAD) inserts its directional derivative along the
    //! tangents of its operands, mc::FFTanMLP -- one operation with as many
    //! results as the MLP has outputs per sweep, evaluated in one forward pass
    //! -- in place of its whole Jacobian, mc::FFGradMLP (false). The Jacobian
    //! is still inserted where it is estimated to be the cheaper, from the
    //! sizes of the MLP and the number of sweeps through the operation -- as
    //! many as the independents it depends on, one for a directional
    //! derivative: see mc::FFMLP::has_tangent
    bool TANGENT;
    //! @brief Whether the backward symbolic differentiation of an MLP operation
    //! in a DAG (BAD, SBAD, DBAD) inserts the product of its transposed
    //! Jacobian with the adjoints of its results, mc::FFAdjMLP -- one operation
    //! with as many results as the MLP has inputs per sweep, evaluated in one
    //! backward pass -- in place of its whole Jacobian, mc::FFGradMLP (false).
    //! The Jacobian is still inserted where it is estimated to be the cheaper,
    //! from the sizes of the MLP and the number of sweeps through the
    //! operation -- as many as the dependents computed from it, one for a
    //! directional derivative: see mc::FFMLP::has_adjoint
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
  //! @brief Expression tree
  FFGraph* _dag;
  //! @brief Expression tree needs updating
  bool _update;
  //! @brief Number of variables
  size_t _nin;
  //! @brief Number of dependents
  size_t _nout;
  //! @brief Number of hidden layers
  size_t _nhid;
  //! @brief Variables
  std::vector<FFVar> _varin;
  //! @brief Dependents
  std::vector<FFVar> _varout;
  //! @brief Codelist
  FFSubgraph _codelist;

  //! @brief Torch Module avaialble
  bool _hasmodule;
#if defined(MC__USE_TORCH)
  //! @brief Torch Device
  torch::Device _device;
  //! @brief Torch Module
  torch::jit::script::Module _module;
  //! @brief Torch data type
  torch::Dtype _dtype;
  //! @brief Whether the Torch Module takes a batch of inputs, one per row
  bool _batched;
#endif

  //! @brief MLP data
  std::vector<std::pair<std::vector<std::vector<double>>, int>> _data;

  //! @brief Intermediate storage for DAG evaluation
  std::vector<std::vector<FFVar>> _wkFF;
  std::vector<std::vector<double>> _wkD;
  std::vector<std::vector<FADType<double>>> _wkFD;
  std::vector<std::vector<BADType<double>>> _wkBD;
  std::vector<std::vector<T>> _wkI;
  std::vector<std::vector<McCormick<T>>> _wkMC;
  std::vector<std::vector<Specbnd<T>>> _wkSB;
  std::vector<std::vector<SupVar<PWCU>>> _wkPWCS;
  std::vector<std::vector<SupVar<PWLU>>> _wkPWLS;
  std::vector<std::vector<SCVar<T>>> _wkSC;
  std::vector<std::vector<OCVar<double>>> _wkOCD;
  std::vector<std::vector<OCVar<FADType<double>>>> _wkOCFD;
  std::vector<PolVar<T>> _wkPOL;
  std::vector<T> _wkCPI;

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

  //! @brief Set expression tree
  void _set_dag(std::vector<std::vector<FFVar>>& wkhid);

  //! @brief Evaluate MLP
  template <typename U>
  void _eval(U const* valin, U* valout,
             std::vector<std::vector<U>>& wkhid) const;

  //! @brief Evaluate MLP through expression tree
  template <typename U>
  void _eval(U const* valin, U* valout, std::vector<U>& wk);

  //! @brief Evaluate MLP directional derivative: forward pass with the values
  //! and the tangents of the neurons
  template <typename U>
  void _tang(U const* valin, U const* dirin, U* tanout, U* valout) const;

  //! @brief Evaluate MLP adjoint: forward pass keeping the derivative of each
  //! activation, then backward pass through the layers
  template <typename U>
  void _adj(U const* valin, U const* wout, U* adjin, U* valout) const;

  //! @brief Evaluate MLP Jacobian in reverse mode: forward pass keeping the
  //! derivative of each activation, then all the rows of the Jacobian
  //! backward through the layers together
  void _jacrev(double const* valin, double* gradout, double* valout) const;

  //! @brief Evaluate MLP Jacobian in forward mode: one pass through the
  //! layers with the values of the neurons and their derivatives in all the
  //! inputs together
  void _jacfwd(double const* valin, double* gradout, double* valout) const;

  //! @brief Value and derivative of an activation function in double: val
  //! holds the argument on entry and the value on exit, as _eval computes it
  double
  _activate(int const activ, double& val) const
  {
    switch (activ)
    {
      case LINEAR:
      default:
        return 1.;
      case RELU:
      {
        double const der = Op<double>::fstep(val);
        val              = ReLU(val);
        return der;
      }
      case TANH:
        val = Op<double>::tanh(val);
        return 1. - val * val;
      case SIGMOID:
        val = (options.SIG2EXP ? Op<double>::inv(Op<double>::exp(-val) + 1.)
                               : Op<double>::tanh(val * 0.5) * 0.5 + 0.5);
        return val * (1. - val);
    }
  }
  //    const;

  //! @brief Reverse evaluate MLP through expression tree
  template <typename U>
  bool _reval(U* valin, U* valout, std::vector<U>& wk, U const& inf);
  //    const;

#if defined(MC__USE_TORCH)
  //! @brief Convert Torch tensor to std::vector
  bool _tensor_to_vector(std::vector<std::vector<double>>& v,
                         at::Tensor const& t, size_t offset, bool disp = false);

  //! @brief Torch tensor of n elements in a given data type, from an array of
  //! doubles -- which it shares where that type is double precision
  static torch::Tensor
  _tensor(double const* val, size_t const n, torch::Dtype const dtype)
  {
    return torch::from_blob(const_cast<double*>(val), {(int64_t)n},
                            torch::TensorOptions().dtype(torch::kFloat64))
        .to(dtype);
  }

  //! @brief Torch tensor of n elements in the data type of the Module
  torch::Tensor
  _tensor(double const* val, size_t const n) const
  {
    return _tensor(val, n, _dtype);
  }

  //! @brief Check that a Torch tensor has n elements
  static void
  _check(torch::Tensor const& t, size_t const n)
  {
    if (!t.defined() || (size_t)t.numel() != n)
      throw std::runtime_error(
          "MLP: Torch Module returns " +
          std::to_string(t.defined() ? t.numel() : 0) + " values where " +
          std::to_string(n) + " are expected\n");
  }

  //! @brief Copy the n elements of a Torch tensor into an array of doubles
  static void
  _copy(torch::Tensor const& t, size_t const n, double* val)
  {
    _check(t, n);
    torch::Tensor const td =
        t.detach().to(torch::kCPU, torch::kFloat64).contiguous();
    std::copy(td.data_ptr<double>(), td.data_ptr<double>() + n, val);
  }
#endif

 public:
  //! @brief Default constructor
  MLP()
      : _dag(nullptr),
        _update(false),
        _nin(0),
        _nout(0),
        _nhid(0),
        _hasmodule(false)
#if defined(MC__USE_TORCH)
        ,
        _device(torch::kCPU),
        _dtype(torch::kFloat32),
        _batched(false)
#endif
  {
    // Decide device once
    // if( torch::cuda::is_available() ){
    //  _device = torch::kCUDA;
    //}
    // else{
    //  _device = torch::kCPU;
    //}
  }

  //! @brief Copy constructor
  MLP(MLP const& other)
      : options(other.options),
        _dag(nullptr),
        _update(false),
        _nin(other._nin),
        _nout(other._nout),
        _nhid(other._nhid),
        _hasmodule(other._hasmodule),
#if defined(MC__USE_TORCH)
        _device(other._device),
        _module(other._module),
        _dtype(other._dtype),
        _batched(other._batched),
#endif
        _data(other._data)
  {
  }

  virtual ~MLP() { delete _dag; }

  //! @brief Assignment operator: takes the network, the options and the Torch
  //! Module of another MLP. The expression tree is not copied but built again
  //! when it is needed, as in the copy constructor. (The implicit operator is
  //! deleted, FFSubgraph not being copy-assignable.)
  MLP&
  operator=(MLP const& other)
  {
    if (this == &other) return *this;
    options = other.options;
    delete _dag;
    _dag    = nullptr;
    _update = false;
    _nin    = other._nin;
    _nout   = other._nout;
    _nhid   = other._nhid;
    _varin.clear();
    _varout.clear();
    _codelist.clear();
    _hasmodule = other._hasmodule;
#if defined(MC__USE_TORCH)
    _device  = other._device;
    _module  = other._module;
    _dtype   = other._dtype;
    _batched = other._batched;
#endif
    _data = other._data;
    // The work arrays hold variables of the previous expression tree, and
    // have the sizes of the previous network
    _wkFF.clear();
    _wkD.clear();
    _wkFD.clear();
    _wkBD.clear();
    _wkI.clear();
    _wkMC.clear();
    _wkSB.clear();
    _wkPWCS.clear();
    _wkPWLS.clear();
    _wkSC.clear();
    _wkOCD.clear();
    _wkOCFD.clear();
    _wkPOL.clear();
    _wkCPI.clear();
    return *this;
  }

  //! @brief Set MLP DAG
  // void set_dag
  //   ()
  //   {
  //     _set_dag( _wkFF );
  //   }

  //! @brief MLP automatic differentiation
  void
  grad(double const* valin, double* gradout, double* const valout = nullptr)
  {
#if defined(MC__USE_TORCH)
    if (_hasmodule && options.EVALTORCH)
    {
      // The parameters of the Module do not require gradients (read_data), so
      // that the backward passes stop at the inputs
      if (_batched && _nout > 1)
      {
        // One forward pass on as many copies of the inputs as there are
        // outputs, and one backward pass from the sum of output i of copy i:
        // the gradient in copy i is the gradient of output i
        torch::Tensor tensin =
            _tensor(valin, _nin).repeat({(int64_t)_nout, 1});
        tensin.set_requires_grad(true);
        torch::Tensor tensout = _module.forward({tensin}).toTensor();
        _check(tensout, _nout * _nout);
        tensout = tensout.view({(int64_t)_nout, (int64_t)_nout});
        if (valout) _copy(tensout[0], _nout, valout);
        _copy(torch::autograd::grad({tensout.diagonal().sum()}, {tensin})[0],
              _nout * _nin, gradout);  // row-major: [i * _nin + j]
      }
      else
      {
        // One backward pass per output
        torch::Tensor tensin = _tensor(valin, _nin);
        tensin.set_requires_grad(true);
        torch::Tensor tensout = _module.forward({tensin}).toTensor();
        _check(tensout, _nout);
        tensout = tensout.view({-1});
        if (valout) _copy(tensout, _nout, valout);
        for (size_t i = 0; i < _nout; ++i)
          _copy(torch::autograd::grad({tensout[i]}, {tensin}, {},
                                      /*retain_graph=*/(i + 1 < _nout))[0],
                _nin, gradout + i * _nin);
      }

#ifdef MC__FFMLP_DEBUG_EVAL
      for (size_t k = 0; k < _nout; ++k)
      {
        if (valout)
          std::cout << "valout[" << k << "] = " << valout[k] << std::endl;
        for (size_t i = 0; i < _nin; ++i)
          std::cout << "gradout[" << k << "][" << i
                    << "] = " << gradout[k * _nin + i] << std::endl;
      }
#endif
      return;
    }
#endif

#ifdef MC__FFMLP_DEBUG_EVAL
    std::cout << "Gradient evaluation separate from Torch\n";
#endif
    // In double arithmetic the Jacobian is computed by passes through the
    // layers, and the mode follows the sizes: the forward mode costs about a
    // pass per input, the reverse mode a pass per output. AUTODIFF = B asks
    // for the reverse mode whatever the sizes -- it used to select the tape
    // of BADType here, a hundred times or so slower. The two modes agree to
    // rounding errors
    if (_nout <= _nin || options.AUTODIFF == Options::AD_TYPE::B)
      _jacrev(valin, gradout, valout);
    else
      _jacfwd(valin, gradout, valout);

#ifdef MC__FFMLP_DEBUG_EVAL
    for (size_t k = 0; k < _nout; ++k)
    {
      for (size_t i = 0; i < _nin; ++i)
        std::cout << "gradout[" << k << "][" << i
                  << "] = " << gradout[k * _nin + i] << std::endl;
    }
#endif
  }
  template <typename U>
  void
  grad(U const* valin, U* gradout, U* const valout = nullptr)
  {
#ifdef MC__FFMLP_DEBUG_EVAL
    std::cout << "Gradient evaluation separate from Torch\n";
#endif
    switch (options.AUTODIFF)
    {
      default:
      case Options::AD_TYPE::F:
      {
        static thread_local std::vector<std::vector<FADType<U>>> wkFU;
        std::vector<FADType<U>> vFvalin(_nin);
        for (size_t i = 0; i < _nin; ++i)
        {
          vFvalin[i] = valin[i];
          vFvalin[i].diff(i, _nin);
        }
        std::vector<FADType<U>> vFvalout(_nout);
        _eval(vFvalin.data(), vFvalout.data(), wkFU);
        for (size_t k = 0; k < _nout; ++k)
        {
          if (valout) valout[k] = vFvalout[k].x();
          for (size_t i = 0; i < _nin; ++i)
            gradout[k * _nin + i] = vFvalout[k].d(i);
          // gradout[i*_nout+k] = vFvalout[k].d(i);
        }
        break;
      }

      case Options::AD_TYPE::B:
      {
        static thread_local std::vector<std::vector<BADType<U>>> wkBU;
        std::vector<BADType<U>> vBvalin(_nin);
        for (size_t i = 0; i < _nin; ++i) vBvalin[i] = valin[i];
        std::vector<BADType<U>> vBvalout(_nout);
        _eval(vBvalin.data(), vBvalout.data(), wkBU);
        for (size_t k = 0; k < _nout; ++k) vBvalout[k].diff(k, _nout);
        // Release hidden-layer reverse-AD nodes before reading input adjoints.
        wkBU.clear();
        for (size_t k = 0; k < _nout; ++k)
        {
          if (valout) valout[k] = vBvalout[k].x();
          for (size_t i = 0; i < _nin; ++i)
            gradout[k * _nin + i] = vBvalin[i].d(k);
          // gradout[i*_nout+k] = vBvalin[i].d(k);
        }
        break;
      }
    }

#ifdef MC__FFMLP_DEBUG_EVAL
    for (size_t k = 0; k < _nout; ++k)
    {
      for (size_t i = 0; i < _nin; ++i)
        std::cout << "gradout[" << k << "][" << i
                  << "] = " << gradout[k * _nin + i] << std::endl;
    }
#endif
  }

  //! @brief MLP directional derivative: tanout (size nout) is the derivative
  //! of the outputs at valin along the direction dirin (size nin)
  void
  tang(double const* valin, double const* dirin, double* tanout,
       double* const valout = nullptr)
  {
#if defined(MC__USE_TORCH)
    if (_hasmodule && options.EVALTORCH)
    {
      // By backward passes only, which are safe across threads -- the dual
      // level of the forward-mode AD of Torch is global to the process. With
      // few outputs, the product with the Jacobian
      if (_nout == 1 || (_batched && _nout <= 8))
      {
        std::vector<double> J(_nout * _nin);
        grad(valin, J.data(), valout);
        for (size_t k = 0; k < _nout; ++k)
        {
          tanout[k] = 0.;
          for (size_t i = 0; i < _nin; ++i)
            tanout[k] += J[k * _nin + i] * dirin[i];
        }
        return;
      }
      // Otherwise two backward passes: the adjoint J^T v as a function of the
      // weights v, then its derivative in v along the direction, which is J d
      torch::Tensor tensin = _tensor(valin, _nin);
      tensin.set_requires_grad(true);
      torch::Tensor tensout = _module.forward({tensin}).toTensor();
      _check(tensout, _nout);
      tensout = tensout.view({-1});
      if (valout) _copy(tensout, _nout, valout);
      torch::Tensor v = torch::zeros({(int64_t)_nout},
                                     torch::TensorOptions().dtype(_dtype));
      v.set_requires_grad(true);
      torch::Tensor const adj = torch::autograd::grad(
          {tensout}, {tensin}, {v}, /*retain_graph=*/true,
          /*create_graph=*/true)[0];
      _copy(torch::autograd::grad({adj}, {v}, {_tensor(dirin, _nin)})[0], _nout,
            tanout);
      return;
    }
#endif
    _tang(valin, dirin, tanout, valout);
  }
  template <typename U>
  void
  tang(U const* valin, U const* dirin, U* tanout, U* const valout = nullptr)
  {
    _tang(valin, dirin, tanout, valout);
  }

  //! @brief MLP adjoint: adjin (size nin) is the product of the transposed
  //! Jacobian of the outputs at valin with the weights wout (size nout), i.e.
  //! the gradient of the weighted sum of the outputs
  void
  adj(double const* valin, double const* wout, double* adjin,
      double* const valout = nullptr)
  {
#if defined(MC__USE_TORCH)
    if (_hasmodule && options.EVALTORCH)
    {
      // One backward pass from the outputs, with the weights as their
      // gradient
      torch::Tensor tensin = _tensor(valin, _nin);
      tensin.set_requires_grad(true);
      torch::Tensor tensout = _module.forward({tensin}).toTensor();
      _check(tensout, _nout);
      tensout = tensout.view({-1});
      if (valout) _copy(tensout, _nout, valout);
      _copy(torch::autograd::grad({tensout}, {tensin},
                                  {_tensor(wout, _nout)})[0],
            _nin, adjin);
      return;
    }
#endif
    _adj(valin, wout, adjin, valout);
  }
  template <typename U>
  void
  adj(U const* valin, U const* wout, U* adjin, U* const valout = nullptr)
  {
    _adj(valin, wout, adjin, valout);
  }

  //! @brief MLP evaluation
  void
  eval(double const* valin, double* valout)
  {
#if defined(MC__USE_TORCH)
    if (_hasmodule && options.EVALTORCH)
    {
      // No graph is needed for the value alone
      torch::NoGradGuard nograd;
      _copy(_module.forward({_tensor(valin, _nin)}).toTensor(), _nout, valout);

#ifdef MC__FFMLP_DEBUG_EVAL
      for (size_t i = 0; i < _nout; ++i)
        std::cout << "valout[" << i << "] = " << valout[i] << std::endl;
#endif
      return;
    }
#endif

    _eval(valin, valout, _wkD);
#ifdef MC__FFMLP_DEBUG_EVAL
    for (size_t i = 0; i < _nout; ++i)
      std::cout << "valout[" << i << "] = " << valout[i] << std::endl;
#endif
  }
  void
  eval(FADType<double> const* valin, FADType<double>* valout)
  {
    _eval(valin, valout, _wkFD);
  }
  void
  eval(BADType<double> const* valin, BADType<double>* valout)
  {
    _eval(valin, valout, _wkBD);
  }
  void
  eval(T const* valin, T* valout)
  {
    _eval(valin, valout, _wkI);
  }
  void
  eval(McCormick<T> const* valin, McCormick<T>* valout)
  {
    _eval(valin, valout, _wkMC);
  }
  void
  eval(Specbnd<T> const* valin, Specbnd<T>* valout)
  {
    _eval(valin, valout, _wkSB);
  }
  void
  eval(SCVar<T> const* valin, SCVar<T>* valout)
  {
    _eval(valin, valout, _wkSC);
  }
  void
  eval(SupVar<PWCU> const* valin, SupVar<PWCU>* valout)
  {
    _eval(valin, valout, _wkPWCS);
  }
  void
  eval(SupVar<PWLU> const* valin, SupVar<PWLU>* valout)
  {
    _eval(valin, valout, _wkPWLS);
  }
  void
  eval(PolVar<T> const* valin, PolVar<T>* valout)
  {
    // Update MLP DAG first
    _set_dag(_wkFF);
    _eval(valin, valout, _wkPOL);
  }
  void
  eval(OCVar<double> const* valin, OCVar<double>* valout)
  {
    _eval(valin, valout, _wkOCD);
  }
  void
  eval(OCVar<FADType<double>> const* valin, OCVar<FADType<double>>* valout)
  {
    _eval(valin, valout, _wkOCFD);
  }
  template <typename U>
  void
  eval(U const* valin, U* valout)
  {
    static thread_local std::vector<std::vector<U>> wkU;
    _eval(valin, valout, wkU);
  }

  //! @brief MLP reverse evaluation
  bool
  reval(T* valin, T* valout)
  {
    // Update MLP DAG first
    _set_dag(_wkFF);
    return _reval(valin, valout, _wkCPI, options.CPINF * T(-1, 1));
  }

  //! @brief Query members
  size_t
  nin() const
  {
    return _nin;
  }
  size_t
  nout() const
  {
    return _nout;
  }
  //! @brief Number of weights, biases aside
  size_t
  nwei() const
  {
    size_t n = 0;
    for (auto const& [layer, activ] : _data)
      n += layer.size() * (layer[0].size() - 1);
    return n;
  }
  size_t
  nhid() const
  {
    return _nhid;
  }
  std::vector<FFVar> const&
  varin()
  // const
  {
    _set_dag(_wkFF);
    return _varin;
  }
  std::vector<FFVar> const&
  varout()
  // const
  {
    _set_dag(_wkFF);
    return _varout;
  }
  FFGraph*
  dag()
  // const
  {
    _set_dag(_wkFF);
    return _dag;
  }
};

#if defined(MC__USE_TORCH)
template <typename T>
inline bool
MLP<T>::read_data(std::string const& filename, bool const disp)
{
  // The Module is read into local objects, and the MLP is only changed once
  // it has been read in full and checked: a failure leaves the MLP as it was
  torch::jit::script::Module module;
  torch::Dtype const dtype = torch::kFloat64;
  std::vector<std::pair<std::vector<std::vector<double>>, int>> data;
  bool batched = false;

  try
  {
    if (disp) std::cout << "Loading Torch Module from " << filename << "...\n";
    module = torch::jit::load(filename);
    module.to(_device);
    // Double precision, whatever the precision of the Module: the conversion
    // of its parameters is exact, and Torch then evaluates the very network
    // that MC++ bounds and relaxes, to rounding errors -- in single
    // precision its values lie some 1e-7 off, outside the bounds of MC++ at
    // the same point
    if (disp)
      for (auto const& p : module.parameters())
      {
        std::cout << "Model parameter dtype: " << p.scalar_type()
                  << ", converted to " << dtype << "\n";
        break;
      }
    module.to(dtype);
    // Evaluation mode, and no gradients with respect to the parameters: the
    // backward passes of grad and adj stop at the inputs, and a forward pass
    // records no graph unless the inputs require one
    module.eval();
    for (auto const& p : module.parameters()) p.set_requires_grad(false);
    if (disp) std::cout << "Model loaded successfully.\n\n";

    if (!module.hasattr("net"))
    {
      if (disp)
        std::cerr << "Module has no 'net' submodule; cannot inspect layers "
                     "and activation functions.\n";
      return false;
    }

    // 'net' is the Sequential defined in Python: Linear submodules, each
    // followed by at most one activation function
    torch::jit::script::Module net = module.attr("net").toModule();
    if (disp) std::cout << "ANN Parameters (weights and biases)\n";
    bool pending = false;  // a Linear layer without activation function yet
    for (auto const& child : net.named_children())
    {
      // Class name of the submodule, without its qualifiers: a name that
      // merely contains "ReLU" or "Tanh" (LeakyReLU, Tanhshrink, ...) is
      // another function
      std::string const type_str = child.value.type()->str();
      std::string const type = type_str.substr(type_str.rfind('.') + 1);
      if (disp)
        std::cout << "Submodule '" << child.name << "'  type: " << type_str
                  << "\n";

      if (type == "Linear")
      {
        if (!child.value.hasattr("weight"))
        {
          if (disp) std::cerr << "Linear submodule without weight\n";
          return false;
        }
        at::Tensor const weight = child.value.attr("weight").toTensor();
        if (weight.dim() != 2 || !weight.size(0) || !weight.size(1) ||
            (!data.empty() &&
             data.back().first.size() != (size_t)weight.size(1)))
        {
          if (disp) std::cerr << "Linear submodule with inconsistent size\n";
          return false;
        }

        // Each neuron holds its bias followed by its weights; a layer without
        // bias has zero biases
        std::vector<std::vector<double>> layer;
        if (!_tensor_to_vector(layer, weight, 1, disp)) return false;
        if (child.value.hasattr("bias"))
        {
          c10::IValue const bias = child.value.attr("bias");
          if (bias.isTensor() &&
              !_tensor_to_vector(layer, bias.toTensor(), 0, disp))
            return false;
        }
        data.push_back(std::make_pair(layer, (int)LINEAR));
        pending = true;
      }

      else if (type == "Tanh" || type == "ReLU" || type == "Sigmoid")
      {
        if (!pending)
        {
          if (disp)
            std::cerr << "Activation function without a Linear layer before "
                         "it\n";
          return false;
        }
        data.back().second =
            (type == "Tanh" ? TANH : type == "ReLU" ? RELU : SIGMOID);
        pending = false;
      }

      // No effect in evaluation mode
      else if (type == "Identity" || type == "Dropout")
        continue;

      else
      {
        if (disp) std::cerr << "Non-supported submodule: " << type_str << "\n";
        return false;
      }
    }
    if (data.empty())
    {
      if (disp) std::cerr << "Module has no Linear submodule\n";
      return false;
    }

    // Check that the Module computes what was read from it, at two points:
    // its forward method may do more than go through 'net'. The second check
    // is whether it takes a batch of inputs, one per row
    MLP<T> mlp;
    mlp.options           = options;
    mlp.options.EVALTORCH = false;
    if (!mlp.set_data(data)) return false;
    size_t const nin = mlp.nin(), nout = mlp.nout();
    std::vector<double> valin(2 * nin, 0.), valout(2 * nout), refout(2 * nout);
    for (size_t i = 0; i < nin; ++i) valin[nin + i] = 0.5 * std::cos(i + 1.);
    torch::NoGradGuard nograd;
    for (size_t k = 0; k < 2; ++k)
    {
      mlp.eval(valin.data() + k * nin, refout.data() + k * nout);
      _copy(module.forward({_tensor(valin.data() + k * nin, nin, dtype)})
                .toTensor(),
            nout, valout.data() + k * nout);
    }
    for (size_t k = 0; k < 2 * nout; ++k)
      if (!(std::fabs(valout[k] - refout[k]) <=
            1e-8 * (1. + std::fabs(refout[k]))))
      {
        if (disp)
          std::cerr << "Module does not evaluate as its 'net' submodule: "
                    << valout[k] << " against " << refout[k] << "\n";
        return false;
      }
    try
    {
      torch::Tensor const tensin =
          _tensor(valin.data(), 2 * nin, dtype).view({2, (int64_t)nin});
      std::vector<double> batout(2 * nout);
      _copy(module.forward({tensin}).toTensor(), 2 * nout, batout.data());
      batched = true;
      for (size_t k = 0; k < 2 * nout; ++k)
        if (!(std::fabs(batout[k] - refout[k]) <=
              1e-8 * (1. + std::fabs(refout[k]))))
          batched = false;
    }
    catch (std::exception const&)
    {
      batched = false;
    }
    if (disp)
      std::cout << "Model " << (batched ? "takes" : "does not take")
                << " batches of inputs\n";
  }

  catch (std::exception const& e)  // c10::Error among them
  {
    if (disp)
      std::cerr << "Error reading Torch Module from " << filename << ":\n"
                << e.what() << std::endl;
    return false;
  }

  // Commit
  if (!set_data(data)) return false;
  _module    = module;
  _dtype     = dtype;
  _batched   = batched;
  _hasmodule = true;
  return true;
}

template <typename T>
inline bool
MLP<T>::_tensor_to_vector(std::vector<std::vector<double>>& v,
                          at::Tensor const& t, size_t offset, bool disp)
{
  if (disp)
  {
    std::cout << "  Size: [";
    for (size_t i = 0; i < t.sizes().size(); ++i)
    {
      std::cout << t.sizes()[i];
      if (i + 1 < t.sizes().size()) std::cout << ", ";
    }
    std::cout << "]\n";
    std::cout << t << "\n\n";
  }
  if (t.dim() < 1 || t.dim() > 2 || !t.is_floating_point())
  {
    if (disp)
      std::cerr << "Unsupported tensor: " << t.dim() << " dimensions, type "
                << t.scalar_type() << "\n";
    return false;
  }

  // In double precision on the CPU, whatever the type and device of the tensor
  at::Tensor const tc =
      t.detach().to(torch::kCPU, torch::kFloat64).contiguous();
  double const* data_ptr = tc.data_ptr<double>();

  // One row of v per element along dimension 0; rows that exist already --
  // from the other tensor of the layer -- must be as many
  size_t const n0 = t.size(0), n1 = (t.dim() == 2 ? t.size(1) : 1);
  if (v.empty())
    v.resize(n0);
  else if (v.size() != n0)
    return false;

  for (size_t i0 = 0; i0 < n0; ++i0)
  {
    if (v[i0].size() < offset + n1) v[i0].resize(offset + n1, 0.);
    for (size_t i1 = 0; i1 < n1; ++i1)
      v[i0][offset + i1] = data_ptr[i0 * n1 + i1];
  }
  return true;
}

#else
template <typename T>
inline bool
MLP<T>::read_data(std::string const& filename, bool const disp)
{
  if (disp) std::cout << "Program compiled without Torch library.\n";
  return false;
}
#endif

template <typename T>
inline bool
MLP<T>::set_data(
    std::vector<std::pair<std::vector<std::vector<double>>, int>> const& mlp)
{
  if (!mlp.size() || !mlp[0].first.size() || (mlp[0].first)[0].size() <= 1)
    return false;

  // Set data
  _data   = mlp;
  _nin    = _data.front().first.front().size() - 1;
  _nout   = _data.back().first.size();
  _nhid   = _data.size() - 1;
  _update = true;

  // Cleanse data
  for (auto& [layer, activ] : _data)
    for (auto& neuron : layer)
      for (auto& weight : neuron)
        if (std::fabs(weight) < options.ZEROTOL) weight = 0.;

  return true;
}

template <typename T>
inline void
MLP<T>::reset_data()
{
  // Reset data
  _nin = _nout = _nhid = 0;
  _data.clear();
  _update = true;
}

template <typename T>
inline bool
MLP<T>::append_data(std::vector<std::vector<double>> const& layer,
                    int const activ, bool const reset)
{
  if (reset)
  {
    reset_data();
    _update = true;
  }
  if (!layer.size() || layer[0].size() <= 1 ||
      (!_data.empty() && _data.back().first.size() != layer[0].size() - 1))
    return false;

  // Set data
  _data.push_back(std::make_pair(layer, activ));
  _nin    = _data.front().first.front().size() - 1;
  _nout   = _data.back().first.size();
  _nhid   = _data.size() - 1;
  _update = true;

  // Cleanse data
  for (auto& neuron : _data.back().first)
    for (auto& weight : neuron)
      if (std::fabs(weight) < options.ZEROTOL) weight = 0.;

  return true;
}

template <typename T>
inline bool
MLP<T>::append_data(std::vector<double> const& layer, int const activ,
                    bool const reset)
{
  return append_data(std::vector<std::vector<double>>({layer}), activ, reset);
}

template <typename T>
inline void
MLP<T>::_set_dag(std::vector<std::vector<FFVar>>& wkhid)
{
  // Nothing to be updated
  if (_dag && !_update) return;

  delete _dag;
  _dag = nullptr;
  if (!_nin || !_nout) return;

  // Created new DAG for MLP
  _dag    = new FFGraph;
  _varin  = _dag->add_vars(_nin, "X");
  _varout = _dag->add_vars(_nout, "Y");

  // Propagate DAG through hidden layers
  wkhid.resize(_nhid);
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << _nhid << std::endl;
#endif
  for (unsigned l = 0; l < _nhid; ++l)
  {
    size_t const nneu = _data[l].first.size();
#ifdef MC__FFMLP_CHECK
    assert(nneu);  // number of neurons in hidden layer l+1
#endif
    wkhid[l].resize(nneu);
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for (unsigned i = 0; i < nneu; ++i)
    {
      FFLin<T> sum;
      // Need to clean data beforehand
      wkhid[l][i] = sum((_data[l].first)[i].size() - 1,
                        l ? wkhid[l - 1].data() : _varin.data(),
                        (_data[l].first)[i].data() + 1, (_data[l].first)[i][0],
                        FFLin<T>::SHALLOW);
      // #ifdef MC__FFMLP_DEBUG
      //       std::cerr << "No inputs to neuron " << i << " in layer " << l <<
      //       ": " << (_data[l].first)[i].size()-1 << std::endl;
      // #endif
      switch (_data[l].second)
      {
        case LINEAR:
        default:
          break;
        case RELU:
          wkhid[l][i] = ReLU(wkhid[l][i]);
          break;
        case TANH:
          wkhid[l][i] = tanh(wkhid[l][i]);
          break;
        case SIGMOID:
          wkhid[l][i] = (options.SIG2EXP ? 1. / (exp(-wkhid[l][i]) + 1.)
                                         : tanh(wkhid[l][i] * 0.5) * 0.5 + 0.5);
          break;
      }
    }
  }

  // Propagate DAG through output layer
#ifdef MC__FFMLP_CHECK
  assert(_data.back().first.size() ==
         _nout);  // number of neurons in output layer
#endif
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << _nhid << ": " << _nout << std::endl;
#endif
  for (unsigned i = 0; i < _nout; ++i)
  {
    FFLin<T> sum;
    // Need to clean data beforehand
    _varout[i] = sum((_data.back().first)[i].size() - 1,
                     _nhid ? wkhid[_nhid - 1].data() : _varin.data(),
                     (_data.back().first)[i].data() + 1,
                     (_data.back().first)[i][0], FFLin<T>::SHALLOW);
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << _nhid << ": "
              << (_data.back().first)[i].size() - 1 << std::endl;
#endif
    switch (_data.back().second)
    {
      case LINEAR:
      default:
        break;
      case RELU:
        _varout[i] = ReLU(_varout[i]);
        break;
      case TANH:
        _varout[i] = tanh(_varout[i]);
        break;
      case SIGMOID:
        _varout[i] = (options.SIG2EXP ? 1. / (exp(-_varout[i]) + 1.)
                                      : tanh(_varout[i] * 0.5) * 0.5 + 0.5);
        break;
    }
  }

#ifdef MC__FFMLP_DEBUG
  _codelist = _dag->subgraph(_varout);
  //  _dag->output( _codelist );
  std::vector<FFExpr> strout = FFExpr::subgraph(_dag, _codelist);
  for (size_t i = 0; i < _nout; ++i)
    std::cout << "F" << i << ": " << strout[i] << std::endl;
//  std::ofstream oFile( "MLP.dot", std::ios_base::out );
//  _dag->dot_script( _varout, oFile );
//  oFile.close();
#else
  _codelist.clear();
#endif
  _update = false;
}

template <typename T>
template <typename U>
inline void
MLP<T>::_eval(U const* valin, U* valout,
              std::vector<std::vector<U>>& wkhid) const
{
  // Propagate through hidden layers
  wkhid.resize(_nhid);
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No hidden layers: " << _nhid << std::endl;
#endif
  for (unsigned l = 0; l < _nhid; ++l)
  {
#ifdef MC__FFMLP_CHECK
    assert(_data[l].first.size());  // number of neurons in layer l+1
#endif
    size_t const nneu = _data[l].first.size();
    wkhid[l].resize(nneu);
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No neurons in layer " << l << ": " << nneu << std::endl;
#endif
    for (unsigned i = 0; i < nneu; ++i)
    {
      wkhid[l][i] = (_data[l].first)[i][0];  // bias term
#ifdef MC__FFMLP_DEBUG
      std::cerr << "No inputs to neuron " << i << " in layer " << l << ": "
                << (_data[l].first)[i].size() - 1 << std::endl;
#endif
      for (unsigned j = 0; j < (_data[l].first)[i].size() - 1; ++j)
      {
        // #ifdef MC__FFMLP_DEBUG
        //         std::cout << "layer:" << l << " neuron:" << i << " input:" <<
        //         j << std::endl;
        // #endif
        if (std::fabs((_data[l].first)[i][1 + j]) < options.ZEROTOL) continue;
        wkhid[l][i] +=
            (l ? wkhid[l - 1][j] : valin[j]) * (_data[l].first)[i][1 + j];
      }
      switch (_data[l].second)
      {
        case LINEAR:
        default:
          break;
        case RELU:
          wkhid[l][i] = ReLU(wkhid[l][i]);
          break;
        case TANH:
          wkhid[l][i] = Op<U>::tanh(wkhid[l][i]);
          break;
        case SIGMOID:
          wkhid[l][i] =
              (options.SIG2EXP ? 1. / (Op<U>::exp(-wkhid[l][i]) + 1.)
                               : Op<U>::tanh(wkhid[l][i] * 0.5) * 0.5 + 0.5);
          break;
      }
    }
  }

  // Propagate through output layers
#ifdef MC__FFMLP_CHECK
  assert(_data.back().first.size());  // number of neurons in layer l+1
#endif
  size_t const nneu = _data.back().first.size();
#ifdef MC__FFMLP_DEBUG
  std::cerr << "No neurons in layer " << _nhid << ": " << nneu << std::endl;
#endif
  for (unsigned i = 0; i < nneu; ++i)
  {
    valout[i] = (_data.back().first)[i][0];  // bias term
#ifdef MC__FFMLP_DEBUG
    std::cerr << "No inputs to neuron " << i << " in layer " << _nhid << ": "
              << (_data.back().first)[i].size() - 1 << std::endl;
#endif
    for (unsigned j = 0; j < (_data.back().first)[i].size() - 1; ++j)
    {
      // #ifdef MC__FFMLP_DEBUG
      //       std::cout << "layer:" << _nhid << " neuron:" << i << " input:" <<
      //       j << std::endl;
      // #endif
      if (std::fabs((_data.back().first)[i][1 + j]) < options.ZEROTOL) continue;
      valout[i] += (_nhid ? wkhid[_nhid - 1][j] : valin[j]) *
                   (_data.back().first)[i][1 + j];
    }
    switch (_data.back().second)
    {
      case LINEAR:
      default:
        break;
      case RELU:
        valout[i] = ReLU(valout[i]);
        break;
      case TANH:
        valout[i] = Op<U>::tanh(valout[i]);
        break;
      case SIGMOID:
        valout[i] =
            (options.SIG2EXP ? Op<U>::inv(Op<U>::exp(-valout[i]) + 1.)
                             : Op<U>::tanh(valout[i] * 0.5) * 0.5 + 0.5);
        break;
    }
  }
}

template <typename T>
template <typename U>
inline void
MLP<T>::_tang(U const* valin, U const* dirin, U* tanout, U* valout) const
{
  // Work arrays: values and tangents of the neurons of two consecutive layers
  static thread_local std::vector<U> wkval, wktan, wkvalprv, wktanprv;
  size_t const nlay = _data.size();

  for (size_t l = 0; l < nlay; ++l)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size();
    wkval.swap(wkvalprv);
    wktan.swap(wktanprv);
    wkval.resize(nneu);
    wktan.resize(nneu);
    for (size_t i = 0; i < nneu; ++i)
    {
      U& val = wkval[i];
      U& tan = wktan[i];
      val    = layer[i][0];  // bias term
      tan    = 0.;
      for (size_t j = 0; j + 1 < layer[i].size(); ++j)
      {
        if (std::fabs(layer[i][1 + j]) < options.ZEROTOL) continue;
        val += (l ? wkvalprv[j] : valin[j]) * layer[i][1 + j];
        tan += (l ? wktanprv[j] : dirin[j]) * layer[i][1 + j];
      }
      switch (activ)
      {
        case LINEAR:
        default:
          break;
        case RELU:
          tan *= Op<U>::fstep(val);
          val = ReLU(val);
          break;
        case TANH:
          val = Op<U>::tanh(val);
          tan *= -Op<U>::sqr(val) + 1.;
          break;
        case SIGMOID:
          val = (options.SIG2EXP ? Op<U>::inv(Op<U>::exp(-val) + 1.)
                                 : Op<U>::tanh(val * 0.5) * 0.5 + 0.5);
          // s (1-s), written with a single occurrence of s
          tan *= (-Op<U>::sqr(val * 2. - 1.) + 1.) * 0.25;
          break;
      }
    }
  }
  for (size_t k = 0; k < _nout; ++k)
  {
    if (valout) valout[k] = wkval[k];
    tanout[k] = wktan[k];
  }
}

template <typename T>
template <typename U>
inline void
MLP<T>::_adj(U const* valin, U const* wout, U* adjin, U* valout) const
{
  // Work arrays: values and activation derivatives of the neurons, layer by
  // layer, and the adjoints of two consecutive layers
  static thread_local std::vector<std::vector<U>> wkval, wkder;
  static thread_local std::vector<U> wkadj, wkprv;
  size_t const nlay = _data.size();
  wkval.resize(nlay);
  wkder.resize(nlay);

  // Forward pass, as in _eval
  for (size_t l = 0; l < nlay; ++l)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size();
    wkval[l].resize(nneu);
    wkder[l].resize(nneu);
    for (size_t i = 0; i < nneu; ++i)
    {
      U& val = wkval[l][i];
      val    = layer[i][0];  // bias term
      for (size_t j = 0; j + 1 < layer[i].size(); ++j)
      {
        if (std::fabs(layer[i][1 + j]) < options.ZEROTOL) continue;
        val += (l ? wkval[l - 1][j] : valin[j]) * layer[i][1 + j];
      }
      switch (activ)
      {
        case LINEAR:
        default:
          break;
        case RELU:
          wkder[l][i] = Op<U>::fstep(val);
          val         = ReLU(val);
          break;
        case TANH:
          val         = Op<U>::tanh(val);
          wkder[l][i] = -Op<U>::sqr(val) + 1.;
          break;
        case SIGMOID:
          val = (options.SIG2EXP ? Op<U>::inv(Op<U>::exp(-val) + 1.)
                                 : Op<U>::tanh(val * 0.5) * 0.5 + 0.5);
          // s (1-s), written with a single occurrence of s
          wkder[l][i] = (-Op<U>::sqr(val * 2. - 1.) + 1.) * 0.25;
          break;
      }
    }
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[nlay - 1][k];

  // Backward pass: the adjoints of the outputs of a layer give those of its
  // inputs
  wkadj.assign(wout, wout + _nout);
  for (size_t l = nlay; l-- > 0;)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size(), nprv = layer[0].size() - 1;
    if (activ == RELU || activ == TANH || activ == SIGMOID)
      for (size_t i = 0; i < nneu; ++i) wkadj[i] *= wkder[l][i];
    wkprv.resize(nprv);
    for (size_t j = 0; j < nprv; ++j) wkprv[j] = 0.;
    // one neuron after the other, along its weights: each adjoint is still
    // summed over the neurons in their order
    for (size_t i = 0; i < nneu; ++i)
      for (size_t j = 0; j < nprv; ++j)
      {
        if (std::fabs(layer[i][1 + j]) < options.ZEROTOL) continue;
        wkprv[j] += wkadj[i] * layer[i][1 + j];
      }
    wkadj.swap(wkprv);
  }
  for (size_t i = 0; i < _nin; ++i) adjin[i] = wkadj[i];
}

template <typename T>
inline void
MLP<T>::_jacrev(double const* valin, double* gradout, double* valout) const
{
  // Work arrays: values and activation derivatives of the neurons, layer by
  // layer, and the adjoints of two consecutive layers for every output -- a
  // row per output
  static thread_local std::vector<std::vector<double>> wkval, wkder;
  static thread_local std::vector<double> wkadj, wkprv;
  size_t const nlay = _data.size();
  wkval.resize(nlay);
  wkder.resize(nlay);

  // Forward pass, with the operations of _eval so that the values are the same
  for (size_t l = 0; l < nlay; ++l)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size();
    wkval[l].resize(nneu);
    wkder[l].resize(nneu);
    double const* prv = (l ? wkval[l - 1].data() : valin);
    for (size_t i = 0; i < nneu; ++i)
    {
      double val = layer[i][0];  // bias term
      for (size_t j = 0; j + 1 < layer[i].size(); ++j)
      {
        if (std::fabs(layer[i][1 + j]) < options.ZEROTOL) continue;
        val += prv[j] * layer[i][1 + j];
      }
      wkder[l][i] = _activate(activ, val);
      wkval[l][i] = val;
    }
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[nlay - 1][k];

  // Backward pass: row k holds the adjoints of the outputs of a layer for
  // output k -- the identity to begin with -- and gives those of its inputs
  size_t ncur = _nout;
  wkadj.assign(_nout * _nout, 0.);
  for (size_t k = 0; k < _nout; ++k) wkadj[k * _nout + k] = 1.;
  for (size_t l = nlay; l-- > 0;)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size(), nprv = layer[0].size() - 1;
    wkprv.assign(_nout * nprv, 0.);
    for (size_t k = 0; k < _nout; ++k)
    {
      double const* adj = wkadj.data() + k * ncur;
      double* prv       = wkprv.data() + k * nprv;
      for (size_t i = 0; i < nneu; ++i)
      {
        // a neuron that does not reach output k, or whose activation is flat
        double const a = adj[i] * wkder[l][i];
        if (a == 0.) continue;
        double const* weight = layer[i].data() + 1;
        for (size_t j = 0; j < nprv; ++j) prv[j] += a * weight[j];
      }
    }
    wkadj.swap(wkprv);
    ncur = nprv;
  }
  std::copy(wkadj.begin(), wkadj.begin() + _nout * _nin, gradout);
}

template <typename T>
inline void
MLP<T>::_jacfwd(double const* valin, double* gradout, double* valout) const
{
  // Work arrays: values of the neurons of two consecutive layers, and their
  // derivatives in the inputs -- a row per neuron
  static thread_local std::vector<double> wkval, wkprv, wktan, wktanprv;
  size_t const nlay = _data.size();

  for (size_t l = 0; l < nlay; ++l)
  {
    auto const& [layer, activ] = _data[l];
    size_t const nneu = layer.size(), nprv = layer[0].size() - 1;
    wkval.swap(wkprv);
    wktan.swap(wktanprv);
    wkval.resize(nneu);
    wktan.assign(nneu * _nin, 0.);
    double const* prv = (l ? wkprv.data() : valin);
    for (size_t i = 0; i < nneu; ++i)
    {
      double val  = layer[i][0];  // bias term
      double* tan = wktan.data() + i * _nin;
      double const* weight = layer[i].data() + 1;
      for (size_t j = 0; j < nprv; ++j)
      {
        if (std::fabs(weight[j]) < options.ZEROTOL) continue;
        val += prv[j] * weight[j];
        // the derivatives of the inputs themselves are the identity
        if (!l)
        {
          tan[j] = weight[j];
          continue;
        }
        double const* tanprv = wktanprv.data() + j * _nin;
        for (size_t c = 0; c < _nin; ++c) tan[c] += weight[j] * tanprv[c];
      }
      double const der = _activate(activ, val);
      if (der != 1.)
        for (size_t c = 0; c < _nin; ++c) tan[c] *= der;
      wkval[i] = val;
    }
  }
  if (valout) std::copy(wkval.begin(), wkval.begin() + _nout, valout);
  std::copy(wktan.begin(), wktan.begin() + _nout * _nin, gradout);
}

template <typename T>
template <typename U>
inline void
MLP<T>::_eval(U const* valin, U* valout, std::vector<U>& wk)
// const
{
  // Run eval on MLP DAG
  _dag->eval(_codelist, wk, _varout.size(), _varout.data(), valout,
             _varin.size(), _varin.data(), valin);
}

template <typename T>
template <typename U>
inline bool
MLP<T>::_reval(U* valin, U* valout, std::vector<U>& wk, U const& inf)
// const
{
  // Run reval on MLP DAG
  int flag = _dag->reval(_codelist, wk, _varout.size(), _varout.data(), valout,
                         _varin.size(), _varin.data(), valin, inf,
                         options.CPMAX, options.CPTHRES);

#ifdef MC__FFMLP_DEBUG
  std::cout << "MLP:: Work array: " << flag << " passes\n";
  for (unsigned i = 0; i < _codelist.len_tap - _codelist.len_wrk; i++)
    std::cout << "wk[" << i << "] = " << wk[i] << std::endl;
  for (unsigned i = 0; i < _varin.size(); i++)
    std::cout << "valin[" << i << "] = " << valin[i] << std::endl;
#endif

  return (flag < 0 ? false : true);
}

//! @brief C++ class defining neural networks as external DAG operations in
//! MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFMLP is a C++ class for defining neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class FFMLP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MLP<T>>
{
 protected:
  using FFEXTERN<T, MLP<T>>::_ptrObj;
  using FFEXTERN<T, MLP<T>>::_ownObj;

  // set the object and related operation in DAG
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin());
#endif
    if (_ownObj && _ptrObj) delete _ptrObj;
    _ownObj = (policy > 0 ? true : false);  // copy;
    //_ownObj = true;
    this->owndata = false;
    this->data = _ptrObj = pMLP;

    FFVar** ppRes =
        this->insert_external_operation(*this, pMLP->nout(), nVar, pVar);

    _ownObj   = false;
    FFOp* pOp = (*ppRes)->opdef().first;
    if (policy > 0)
      _ptrObj =
          static_cast<FFMLP<T>*>(pOp)->_ptrObj;  // set MLP pointer to DAG copy
    else if (policy < 0)
      static_cast<FFMLP<T>*>(pOp)->_ownObj = true;  // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFMLP_DEBUG
    std::cerr << "MLP operation address: " << this << std::endl;
    std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
    return ppRes;
  }

 public:
  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY    = 1,  //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER =
        -1  //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFMLP() : FFEXTERN<T, MLP<T>>() {}

  // Destructor
  virtual ~FFMLP() {}

  // Copy constructor
  FFMLP(FFMLP<T> const& Other) : FFEXTERN<T, MLP<T>>(Other) {}

  // Define operation
  // FFVar** operator()
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin());
#endif
    // return _set( vVar.size(), vVar.data(), pMLP, policy );
    FFVar** ppDer = _set(vVar.size(), vVar.data(), pMLP, policy);
    std::vector<FFVar> vDer(pMLP->nout());
    for (size_t i = 0; i < vDer.size(); ++i) vDer[i] = *ppDer[i];
    return vDer;  // std::move( vDer );
  }

  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vVar, MLP<T>* pMLP,
             int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin() && idep < pMLP->nout());
#endif
    return *(_set(vVar.size(), vVar.data(), pMLP, policy)[idep]);
  }

  FFVar**
  operator()(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP,
             int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin());
#endif
    return _set(nVar, pVar, pMLP, policy);
  }

  FFVar&
  operator()(size_t const idep, size_t const nVar, FFVar const* pVar,
             MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin() && idep < pMLP->nout());
#endif
    return *(_set(nVar, pVar, pMLP, policy)[idep]);
  }

  // MLP pointer
  MLP<T>*
  pMLP() const
  {
    // std::cerr << "MLP address retreived: " << _ptrObj << std::endl;
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FADType<FFVar>) ))
      return eval(nRes, static_cast<FADType<FFVar>*>(vRes), nVar,
                  static_cast<FADType<FFVar> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FADType<double>) ))
      return eval(nRes, static_cast<FADType<double>*>(vRes), nVar,
                  static_cast<FADType<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(PolVar<T>) ))
      return eval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                  static_cast<PolVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFMLP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
            FADType<FFVar> const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
            PolVar<T> const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    if (mc::same_type( idU, typeid(T) ))
      return reval(nRes, static_cast<T*>(vRes), nVar, static_cast<T*>(vVar));
    else if (mc::same_type( idU, typeid(PolVar<T>) ))
      return reval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                   static_cast<PolVar<T>*>(vVar));

    throw std::runtime_error(
        "FFMLP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  bool reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const;

  bool reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
             PolVar<T>* vVar) const;

  // Derivatives
  void deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
             FFVar const* vVar, FFVar** vDer) const;

  // Whether nSweep passes through the network, each costing npass
  // evaluations, are estimated to be cheaper in a DAG than the Jacobian
  // operation and the products that the differentiation forms with its
  // entries in each sweep. In units of one weight in an evaluation of the
  // network, with W weights, nin inputs and nout outputs, timings of the
  // derivative DAG in double arithmetic give about
  //  - a pass:       npass W + 100, the latter for the operation itself
  //  - the Jacobian: W (1/2 + min(nin,nout)/4)
  //  - the products: 35 nin nout, when the Jacobian is evaluated and again
  //                  in each sweep
  // The passes are always the cheaper where the Jacobian has many entries
  // for the size of the network, and never, even for a single sweep, where
  // it has few -- a network with a single output above all. Over networks
  // of 2 to 200 inputs and 1 to 50 outputs and 1 to 32 sweeps, the choice
  // came within 4% of the better of the two in total; where it was the
  // wrong one it lost tens of microseconds, and 0.2 to 0.7 ms (16 to 34%)
  // at worst, with 200 inputs for 5 outputs and 12 to 16 sweeps
  bool
  passes_pay(unsigned const nSweep, double const npass) const
  {
    double const W = _ptrObj->nwei(), E = _ptrObj->nin() * _ptrObj->nout(),
                 R = std::min(_ptrObj->nin(), _ptrObj->nout());
    return nSweep * (npass * W + 100. - 35. * E) <=
           W * (0.5 + 0.25 * R) + 35. * E;
  }

  // Tangents: with MLP<T>::Options::TANGENT, the forward symbolic
  // differentiation inserts the directional derivative operation FFTanMLP
  // along the tangents of the operands, in place of the Jacobian FFGradMLP
  // -- where that is estimated to be the cheaper for nSweep sweeps, a
  // tangent pass costing about 1.1 evaluations of the network
  bool
  has_tangent(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.TANGENT && passes_pay(nSweep, 1.1);
  }
  void tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vTan,
               FFVar* vResTan) const;

  // Adjoints: with MLP<T>::Options::ADJOINT, the backward symbolic
  // differentiation inserts the adjoint operation FFAdjMLP for the adjoints
  // of the results, in place of the Jacobian FFGradMLP -- where that is
  // estimated to be the cheaper for nSweep sweeps, an adjoint pass costing
  // about 1.3 evaluations of the network
  bool
  has_adjoint(unsigned const nSweep) const
  {
    return _ptrObj && _ptrObj->options.ADJOINT && passes_pay(nSweep, 1.3);
  }
  void adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
               FFVar const* vVar, FFVar const* const* vResAdj,
               FFVar* vVarAdj) const;

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return "MLP[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

//! @brief C++ class defining gradient of neural networks as external DAG
//! operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFGradMLP is a C++ class for defining gradient of neural
//! networks as external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FFGradMLP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MLP<T>>
{
  friend class FFMLP<T>;

 protected:
  using FFEXTERN<T, MLP<T>>::_ptrObj;
  using FFEXTERN<T, MLP<T>>::_ownObj;

  // set the object and related operation in DAG
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin());
#endif
    if (_ownObj && _ptrObj) delete _ptrObj;
    _ownObj = (policy > 0 ? true : false);  // copy;
    //_ownObj = true;
    this->owndata = false;
    this->data = _ptrObj = pMLP;

    FFVar** ppRes =
        this->insert_external_operation(*this, nVar * pMLP->nout(), nVar, pVar);

    _ownObj   = false;
    FFOp* pOp = (*ppRes)->opdef().first;
    if (policy > 0)
      _ptrObj = static_cast<FFGradMLP<T>*>(pOp)
                    ->_ptrObj;  // set MLP pointer to DAG copy
    else if (policy < 0)
      static_cast<FFGradMLP<T>*>(pOp)->_ownObj = true;  // transfer ownership
      // nothing to do if policy = 0
#ifdef MC__FFMLP_DEBUG
    std::cerr << "GradMLP operation address: " << this << std::endl;
    std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
    return ppRes;
  }

 public:
  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY    = 1,  //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER =
        -1  //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFGradMLP() : FFEXTERN<T, MLP<T>>() {}

  // Destructor
  virtual ~FFGradMLP() {}

  // Copy constructor
  FFGradMLP(FFGradMLP<T> const& Other) : FFEXTERN<T, MLP<T>>(Other) {}

  // Define operation
  FFVar**
  operator()(std::vector<FFVar> const& vVar, MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin());
#endif
    return _set(vVar.size(), vVar.data(), pMLP, policy);
  }

  FFVar&
  operator()(size_t const idep, std::vector<FFVar> const& vVar, MLP<T>* pMLP,
             int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin() && idep < pMLP->nout() * pMLP->nin());
#endif
    return *(_set(vVar.size(), vVar.data(), pMLP, policy)[idep]);
  }

  FFVar**
  operator()(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP,
             int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin());
#endif
    return _set(nVar, pVar, pMLP, policy);
  }

  FFVar&
  operator()(size_t const idep, size_t const nVar, FFVar const* pVar,
             MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin() && idep < pMLP->nout() * pMLP->nin());
#endif
    return *(_set(nVar, pVar, pMLP, policy)[idep]);
  }

  // MLP pointer
  MLP<T>*
  pMLP() const
  {
    // std::cerr << "MLP address retreived: " << _ptrObj << std::endl;
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFGradMLP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFGradMLP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return "GradMLP[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

//! @brief C++ class defining directional derivatives of neural networks as
//! external DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFTanMLP is a C++ class for defining the directional derivative
//! of a neural network as an external DAG operation in MC++: its
//! operands are the inputs of the network followed by a direction of
//! the same size, and its results are the derivatives of the outputs
//! of the network along that direction, evaluated in a single forward
//! pass. It is inserted by the forward symbolic differentiation of
//! mc::FFMLP when mc::MLP::Options::TANGENT is set (FFOp::tangent),
//! and can also be defined directly. It has no derivative of its own.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FFTanMLP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MLP<T>>
{
  friend class FFMLP<T>;

 protected:
  using FFEXTERN<T, MLP<T>>::_ptrObj;
  using FFEXTERN<T, MLP<T>>::_ownObj;

  // set the object and related operation in DAG; the operands are the inputs
  // followed by the direction
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == 2 * pMLP->nin());
#endif
    if (_ownObj && _ptrObj) delete _ptrObj;
    _ownObj = (policy > 0 ? true : false);  // copy;
    this->owndata = false;
    this->data = _ptrObj = pMLP;

    FFVar** ppRes =
        this->insert_external_operation(*this, pMLP->nout(), nVar, pVar);

    _ownObj   = false;
    FFOp* pOp = (*ppRes)->opdef().first;
    if (policy > 0)
      _ptrObj = static_cast<FFTanMLP<T>*>(pOp)
                    ->_ptrObj;  // set MLP pointer to DAG copy
    else if (policy < 0)
      static_cast<FFTanMLP<T>*>(pOp)->_ownObj = true;  // transfer ownership
      // nothing to do if policy = 0
    return ppRes;
  }

 public:
  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY    = 1,  //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER =
        -1  //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFTanMLP() : FFEXTERN<T, MLP<T>>() {}

  // Destructor
  virtual ~FFTanMLP() {}

  // Copy constructor
  FFTanMLP(FFTanMLP<T> const& Other) : FFEXTERN<T, MLP<T>>(Other) {}

  // Define operation: derivatives of the outputs at vVar along vDir
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vDir,
             MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin() && vDir.size() == pMLP->nin());
#endif
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vDir.begin(), vDir.end());
    FFVar** ppTan = _set(vOpd.size(), vOpd.data(), pMLP, policy);
    std::vector<FFVar> vTan(pMLP->nout());
    for (size_t k = 0; k < vTan.size(); ++k) vTan[k] = *ppTan[k];
    return vTan;
  }

  // MLP pointer
  MLP<T>*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFTanMLP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFTanMLP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return "TanMLP[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T>
template <typename U>
inline void
FFTanMLP<T>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                  unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFTanMLP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  _ptrObj->tang(vVar, vVar + _ptrObj->nin(), vRes);
}

template <typename T>
inline void
FFTanMLP<T>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                  FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFTanMLP<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                  FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFTanMLP<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                  FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  FFVar** ppRes = this->insert_external_operation(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T>
inline void
FFTanMLP<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                  SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == 2 * _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFTanMLP<T>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                  FFExpr const* vVar, unsigned const* mVar) const
{
  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
      {
        std::ostringstream os;
        os << name() << "[" << j << "]";
        vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
      }
      break;
    case FFExpr::Options::GAMS:
    default:
      throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
  }
}

//! @brief C++ class defining adjoints of neural networks as external
//! DAG operations in MC++.
////////////////////////////////////////////////////////////////////////
//! mc::FFAdjMLP is a C++ class for defining the adjoint of a neural
//! network as an external DAG operation in MC++: its operands are the
//! inputs of the network followed by one weight per output, and its
//! results are the derivatives of the weighted sum of the outputs with
//! respect to the inputs -- the product of the transposed Jacobian with
//! the weights -- evaluated in a single backward pass. It is inserted
//! by the backward symbolic differentiation of mc::FFMLP when
//! mc::MLP::Options::ADJOINT is set (FFOp::adjoint), and can also be
//! defined directly. It has no derivative of its own.
////////////////////////////////////////////////////////////////////////
template <typename T>
class FFAdjMLP
    ////////////////////////////////////////////////////////////////////////
    : public FFEXTERN<T, MLP<T>>
{
  friend class FFMLP<T>;

 protected:
  using FFEXTERN<T, MLP<T>>::_ptrObj;
  using FFEXTERN<T, MLP<T>>::_ownObj;

  // set the object and related operation in DAG; the operands are the inputs
  // followed by the weights of the outputs
  FFVar**
  _set(size_t const nVar, FFVar const* pVar, MLP<T>* pMLP, int policy)
  {
#ifdef MC__FFMLP_CHECK
    assert(nVar == pMLP->nin() + pMLP->nout());
#endif
    if (_ownObj && _ptrObj) delete _ptrObj;
    _ownObj = (policy > 0 ? true : false);  // copy;
    this->owndata = false;
    this->data = _ptrObj = pMLP;

    FFVar** ppRes =
        this->insert_external_operation(*this, pMLP->nin(), nVar, pVar);

    _ownObj   = false;
    FFOp* pOp = (*ppRes)->opdef().first;
    if (policy > 0)
      _ptrObj = static_cast<FFAdjMLP<T>*>(pOp)
                    ->_ptrObj;  // set MLP pointer to DAG copy
    else if (policy < 0)
      static_cast<FFAdjMLP<T>*>(pOp)->_ownObj = true;  // transfer ownership
      // nothing to do if policy = 0
    return ppRes;
  }

 public:
  //! @brief Enumeration type for copy policy of MLP object
  enum POLICY_TYPE
  {
    SHALLOW = 0,  //!< Shallow copy of MLP object in FFGraph (without ownership)
    COPY    = 1,  //!< Deep copy of MLP object in FFGraph (with ownership)
    TRANSFER =
        -1  //!< Shallow copy of MLP object in FFGraph (with ownership transfer)
  };

  //! @brief Default constructor
  FFAdjMLP() : FFEXTERN<T, MLP<T>>() {}

  // Destructor
  virtual ~FFAdjMLP() {}

  // Copy constructor
  FFAdjMLP(FFAdjMLP<T> const& Other) : FFEXTERN<T, MLP<T>>(Other) {}

  // Define operation: derivatives in vVar of the sum of the outputs weighted
  // by vWei
  std::vector<FFVar>
  operator()(std::vector<FFVar> const& vVar, std::vector<FFVar> const& vWei,
             MLP<T>* pMLP, int policy = COPY)
  {
#ifdef MC__FFMLP_CHECK
    assert(vVar.size() == pMLP->nin() && vWei.size() == pMLP->nout());
#endif
    std::vector<FFVar> vOpd(vVar);
    vOpd.insert(vOpd.end(), vWei.begin(), vWei.end());
    FFVar** ppAdj = _set(vOpd.size(), vOpd.data(), pMLP, policy);
    std::vector<FFVar> vAdj(pMLP->nin());
    for (size_t i = 0; i < vAdj.size(); ++i) vAdj[i] = *ppAdj[i];
    return vAdj;
  }

  // MLP pointer
  MLP<T>*
  pMLP() const
  {
    return _ptrObj;
  }

  // Forward evaluation overloads
  virtual void
  feval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void const* vVar, unsigned const* mVar) const
  {
    if (mc::same_type( idU, typeid(FFVar) ))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFDep) ))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFInv) ))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(double) ))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(T) ))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(Specbnd<T>) ))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SCVar<T>) ))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(McCormick<T>) ))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWCU>) ))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SupVar<PWLU>) ))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<double>) ))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FADType<double>>) ))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(OCVar<FFDep>) ))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(SLiftVar) ))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (mc::same_type( idU, typeid(FFExpr) ))
      return eval(nRes, static_cast<FFExpr*>(vRes), nVar,
                  static_cast<FFExpr const*>(vVar), mVar);

    throw std::runtime_error(
        "FFAdjMLP::feval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  template <typename U>
  void eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
            unsigned const* mVar) const;

  void eval(size_t const nRes, FFDep* vRes, size_t const nVar,
            FFDep const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFInv* vRes, size_t const nVar,
            FFInv const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
            FFExpr const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, FFVar* vRes, size_t const nVar,
            FFVar const* vVar, unsigned const* mVar) const;

  void eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
            SLiftVar const* vVar, unsigned const* mVar) const;

  // Backward evaluation overloads
  virtual bool
  reval(std::type_info const& idU, unsigned const nRes, void* vRes,
        unsigned const nVar, void* vVar) const
  {
    throw std::runtime_error(
        "FFAdjMLP::reval: **ERROR** No evaluation method with type" +
        std::string(idU.name()) + "\n");
  }

  // Properties
  std::string
  name() const
  {
    std::ostringstream oss;
    oss << this->data;
    return "AdjMLP[" + oss.str() + "]";
  }

  //! @brief Return whether or not operation is commutative
  bool
  commutative() const
  {
    return false;
  }
};

template <typename T>
template <typename U>
inline void
FFAdjMLP<T>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                  unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFAdjMLP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  _ptrObj->adj(vVar, vVar + _ptrObj->nin(), vRes);
}

template <typename T>
inline void
FFAdjMLP<T>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                  FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFAdjMLP<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                  FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFAdjMLP<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                  FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  FFVar** ppRes = this->insert_external_operation(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T>
inline void
FFAdjMLP<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                  SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nin() &&
         nVar == _ptrObj->nin() + _ptrObj->nout());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFAdjMLP<T>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                  FFExpr const* vVar, unsigned const* mVar) const
{
  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
      {
        std::ostringstream os;
        os << name() << "[" << j << "]";
        vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
      }
      break;
    case FFExpr::Options::GAMS:
    default:
      throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
  }
}

template <typename T>
inline void
FFMLP<T>::tangent(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                  FFVar const* vVar, FFVar const* const* vTan,
                  FFVar* vResTan) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::tangent\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  // Operands of the directional derivative: the inputs, then their tangents
  // -- zero where the tangent of an input is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(2 * nVar);
  for (unsigned i = 0; i < nVar; ++i)
    vOpd.push_back(vTan[i] ? *vTan[i] : FFVar(0.));

  FFTanMLP<T> ResTan;
  FFVar const* const* vResDer = ResTan._set(vOpd.size(), vOpd.data(), _ptrObj, COPY);
  for (unsigned k = 0; k < nRes; ++k) vResTan[k] = *vResDer[k];
}

template <typename T>
inline void
FFMLP<T>::adjoint(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                  FFVar const* vVar, FFVar const* const* vResAdj,
                  FFVar* vVarAdj) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::adjoint\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  // Operands of the adjoint operation: the inputs, then the adjoints of the
  // outputs -- zero where the adjoint of an output is a structural zero
  std::vector<FFVar> vOpd(vVar, vVar + nVar);
  vOpd.reserve(nVar + nRes);
  for (unsigned k = 0; k < nRes; ++k)
    vOpd.push_back(vResAdj[k] ? *vResAdj[k] : FFVar(0.));

  FFAdjMLP<T> ResAdj;
  FFVar const* const* vResDer = ResAdj._set(vOpd.size(), vOpd.data(), _ptrObj, COPY);
  for (unsigned i = 0; i < nVar; ++i) vVarAdj[i] = *vResDer[i];
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
               FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFDep\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
               FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFInv\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFGradMLP<T>::eval(size_t const nRes, FFDep* vRes, size_t const nVar,
                   FFDep const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFDep\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFDep::TYPE::N);
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFGradMLP<T>::eval(size_t const nRes, FFInv* vRes, size_t const nVar,
                   FFInv const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFInv\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vRes[0] = 0;
  for (unsigned i = 0; i < nVar; ++i) vRes[0] += vVar[i];
  vRes[0].update(FFInv::TYPE::U);  // Not a candidate for inversion
  for (unsigned j = 1; j < nRes; ++j) vRes[j] = vRes[0];
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
               FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFVar\n";
  std::cerr << "FFMLP operation address: " << this << std::endl;
  std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  FFVar** ppRes = this->insert_external_operation(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T>
inline void
FFGradMLP<T>::eval(size_t const nRes, FFVar* vRes, size_t const nVar,
                   FFVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFVar\n";
  std::cerr << "FFGradMLP operation address: " << this << std::endl;
  std::cerr << "MLP address in DAG: " << _ptrObj << std::endl;
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  FFVar** ppRes = this->insert_external_operation(*this, nRes, nVar, vVar);
  for (unsigned j = 0; j < nRes; ++j) vRes[j] = *(ppRes[j]);
}

template <typename T>
template <typename U>
inline void
FFMLP<T>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
               unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  _ptrObj->eval(vVar, vRes);
}

template <typename T>
template <typename U>
inline void
FFGradMLP<T>::eval(size_t const nRes, U* vRes, size_t const nVar, U const* vVar,
                   unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: " << typeid(vRes[0]).name() << " (generic)\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  _ptrObj->grad(vVar, vRes);
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, FADType<FFVar>* vRes, size_t const nVar,
               FADType<FFVar> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FADType<FFVar>\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  std::vector<FFVar> vVarVal(nVar);
  for (unsigned i = 0; i < nVar; ++i) vVarVal[i] = vVar[i].val();
  FFVar const* const* vResVal =
      this->insert_external_operation(*this, nRes, nVar, vVarVal.data());

  FFGradMLP<T> ResDer;
  FFVar const* const* vResDer =
      ResDer._set(nVar, vVarVal.data(), _ptrObj, COPY);

  for (unsigned k = 0; k < nRes; ++k)
  {
    vRes[k] = *vResVal[k];
    for (unsigned i = 0; i < nVar; ++i) vRes[k].setDepend(vVar[i]);
    for (unsigned j = 0; j < vRes[k].size(); ++j)
    {
      vRes[k][j] = 0.;
      for (unsigned i = 0; i < nVar; ++i)
      {
        if (vVar[i][j].cst() && vVar[i][j].num().val() == 0.) continue;
        // vRes[k][j] += *vResDer[k+nRes*i] * vVar[i][j];
        vRes[k][j] += *vResDer[k * nVar + i] * vVar[i][j];
      }
    }
  }
}

template <typename T>
inline void
FFMLP<T>::deriv(unsigned const nRes, FFVar const* vRes, unsigned const nVar,
                FFVar const* vVar, FFVar** vDer) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::deriv\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  FFGradMLP<T> ResDer;
  FFVar const* const* vResDer = ResDer._set(nVar, vVar, _ptrObj, COPY);
  for (unsigned k = 0; k < nRes; ++k)
    for (unsigned i = 0; i < nVar; ++i)
      // vDer[k][i] = *vResDer[k+nRes*i];
      vDer[k][i] = *vResDer[k * nVar + i];
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
               SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFGradMLP<T>::eval(size_t const nRes, SLiftVar* vRes, size_t const nVar,
                   SLiftVar const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: SLiftVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() * _ptrObj->nin() &&
         nVar == _ptrObj->nin());
#endif

  vVar->env()->lift(nRes, vRes, nVar, vVar);
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
               FFExpr const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: FFExpr\n";
#endif

  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
      {
        std::ostringstream os;
        os << name() << "[" << j << "]";
        vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
      }
      break;
    case FFExpr::Options::GAMS:
    default:
      throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
  }
}

template <typename T>
inline void
FFGradMLP<T>::eval(size_t const nRes, FFExpr* vRes, size_t const nVar,
                   FFExpr const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFGradMLP::eval: FFExpr\n";
#endif

  switch (FFExpr::options.LANG)
  {
    case FFExpr::Options::DAG:
      for (unsigned j = 0; j < nRes; ++j)
      {
        std::ostringstream os;
        os << name() << "[" << j << "]";
        vRes[j] = FFExpr::compose(os.str(), nVar, vVar);
      }
      break;
    case FFExpr::Options::GAMS:
    default:
      throw typename FFExpr::Exceptions(FFExpr::Exceptions::UNDEF);
  }
}

template <typename T>
inline void
FFMLP<T>::eval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
               PolVar<T> const* vVar, unsigned const* mVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::eval: PolVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFBase* dag    = vVar[0].var().dag();
#ifdef MC__FFMLP_CHECK
  assert(img && dag);
#endif
  FFVar** ppRes = dag->curOp()->varout.data();
#ifdef MC__FFMLP_CHECK
  assert(nRes == dag->curOp()->varout.size());
#endif

  this->_resize_relax(img);
  this->_propagate_relax(img, ppRes, vRes, vVar);
}

template <typename T>
inline bool
FFMLP<T>::reval(size_t const nRes, PolVar<T>* vRes, size_t const nVar,
                PolVar<T>* vVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: PolVar\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  PolImg<T>* img = vVar[0].image();
  FFOp* pop      = vVar[0].var().opdef().first;
#ifdef MC__FFMLP_CHECK
  assert(img && pop);
#endif

  this->_backpropagate_relax(img, pop, vRes, vVar);
  return true;
}

template <typename T>
inline bool
FFMLP<T>::reval(size_t const nRes, T* vRes, size_t const nVar, T* vVar) const
{
#ifdef MC__FFMLP_TRACE
  std::cout << "FFMLP::reval: T\n";
#endif
#ifdef MC__FFMLP_CHECK
  assert(_ptrObj && nRes == _ptrObj->nout() && nVar == _ptrObj->nin());
#endif

  return _ptrObj->reval(vVar, vRes);
}

}  // end namespace mc

#endif
