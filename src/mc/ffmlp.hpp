// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#ifndef MC__FFMLP_HPP
#define MC__FFMLP_HPP

#include <algorithm>
#if defined(MC__USE_TORCH)
#include <torch/script.h>  // One-stop header for TorchScript
                           // #include <ATen/cuda/CUDAContext.h>
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
  //! #neurons in hidden layer; inner vector has scalar weights multiplying each
  //! neuron in previous (input or hidden) layer and a bias term and size
  //! 1+#neuron in previous layer
  bool append_data(std::vector<std::vector<double>> const& layer,
                   int const activ = LINEAR, bool const reset = false);

  //! @brief Append single-neuron layer to MLP data: vector has scalar weights
  //! multiplying each neuron in previous (input or hidden) layer and a bias
  //! term and size 1+#neuron in previous layer
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

    //! @brief Whether to evaluate ANN in floating-point arithemic using Torch
    bool EVALTORCH;
    //! @brief Threshold for zero coefficient in neural network
    double ZEROTOL;
    //! @brief Whether to convert ReLU to abs (true) or max (false) function
    bool RELU2ABS;
    //! @brief Whether to convert sigmoid to exp (true) or tanh (false)
    bool SIG2EXP;
    //! @brief Whether to apply forward or reverse automatic differentiation
    int AUTODIFF;
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
  //    const;

  //! @brief Reverse evaluate MLP through expression tree
  template <typename U>
  bool _reval(U* valin, U* valout, std::vector<U>& wk, U const& inf);
  //    const;

#if defined(MC__USE_TORCH)
  //! @brief Convert Torch tensor to std::vector
  bool _tensor_to_vector(std::vector<std::vector<double>>& v,
                         at::Tensor const& t, size_t offset, bool disp = false);
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
        _dtype(torch::kFloat32)
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
#endif
        _data(other._data)
  {
  }

  virtual ~MLP() { delete _dag; }

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
      // 1) Wrap input and make it a leaf with grad
      auto base =
          torch::from_blob(const_cast<double*>(valin), {(int64_t)_nin},
                           torch::TensorOptions().dtype(torch::kFloat64));
      torch::Tensor tensin = base.to(_dtype).detach();
      tensin.set_requires_grad(true);  // now a leaf

      // 2) Forward
      auto tensout = _module.forward({tensin}).toTensor().to(
          torch::kFloat64);          // shape (_nout) or (_nout,1)
      tensout = tensout.view({-1});  // ( _nout )

      if (valout)
      {
        double* ptensout = tensout.data_ptr<double>();
        std::copy(ptensout, ptensout + _nout, valout);
      }

      // Storage for full Jacobian d y_i / d x_j
      // Row-major: [i * _nin + j]
      std::vector<double> J(_nout * _nin);

      for (size_t i = 0; i < _nout; ++i)
      {
        // Clear previous grad on tensin
        if (tensin.grad().defined()) tensin.grad().zero_();

        // 3) Take scalar output y_i
        torch::Tensor yi = tensout[i];  // scalar

        // 4) Backward from this scalar
        //    gradient (first arg) left empty, keep_graph set explicitly
        yi.backward(
            /*gradient=*/torch::Tensor(),    // equivalent to default {}
            /*keep_graph=*/(i + 1 < _nout),  // keep graph until last iteration
            /*create_graph=*/false);

        // 5) Gradient wrt input: d y_i / d x
        torch::Tensor grad = tensin.grad().to(torch::kFloat64);  // shape (_nin)

        double* gptr = grad.data_ptr<double>();
        std::copy(gptr, gptr + _nin, gradout + i * _nin);
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
    switch (options.AUTODIFF)
    {
      default:
      case Options::AD_TYPE::F:
      {
        std::vector<FADType<double>> vFvalin(_nin);
        for (size_t i = 0; i < _nin; ++i)
        {
          vFvalin[i] = valin[i];
          vFvalin[i].diff(i, _nin);
        }
        std::vector<FADType<double>> vFvalout(_nout);
        _eval(vFvalin.data(), vFvalout.data(), _wkFD);
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
        std::vector<BADType<double>> vBvalin(_nin);
        for (size_t i = 0; i < _nin; ++i) vBvalin[i] = valin[i];
        std::vector<BADType<double>> vBvalout(_nout);
        _eval(vBvalin.data(), vBvalout.data(), _wkBD);
        for (size_t k = 0; k < _nout; ++k) vBvalout[k].diff(k, _nout);
        // FADBAD++ reverse mode propagates adjoints when expression nodes
        // are released. The persistent hidden-layer work array holds
        // references to those nodes, so release it before reading input
        // derivatives.
        _wkBD.clear();
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

  //! @brief MLP evaluation
  void
  eval(double const* valin, double* valout)
  {
#if defined(MC__USE_TORCH)
    if (_hasmodule && options.EVALTORCH)
    {
      // 1) Wrap input and make it a leaf
      auto base =
          torch::from_blob(const_cast<double*>(valin), {(int64_t)_nin},
                           torch::TensorOptions().dtype(torch::kFloat64));
      torch::Tensor tensin = base.to(_dtype).detach();

      // 2) Forward
      auto tensout = _module.forward({tensin}).toTensor().to(
          torch::kFloat64);          // shape (_nout) or (_nout,1)
      tensout = tensout.view({-1});  // ( _nout )

      double* ptensout = tensout.data_ptr<double>();
      std::copy(ptensout, ptensout + _nout, valout);

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
  _hasmodule = false;
  try
  {
    if (disp) std::cout << "Loading Torch Module from " << filename << "...\n";
    _module = torch::jit::load(filename);
    _module.to(_device);
    _hasmodule = true;
    if (disp) std::cout << "Model loaded successfully.\n\n";
  }

  catch (const c10::Error& e)
  {
    if (disp)
      std::cerr << "Error loading Torch Module from " << filename << ":\n"
                << e.msg() << std::endl;
    return false;
  }

  if (!_module.hasattr("net"))
  {
    if (disp)
      std::cout << "Module has no 'net' submodule; cannot inspect activation "
                   "functions.\n";
    return false;
  }

  // Get data type
  for (const auto& p : _module.named_parameters())
  {
    _dtype = p.value.scalar_type();
    if (disp) std::cout << "Model parameter dtype: " << _dtype << "\n";
    break;
  }

  // 'net' is the Sequential defined in Python
  torch::jit::script::Module net = _module.attr("net").toModule();

  reset_data();
  if (disp) std::cout << "ANN Parameters (weights and biases)\n";

  // Populate _data
  std::vector<std::vector<double>> layer;
  size_t offset;
  bool bias = false, weight = false;
  for (const auto& p : _module.named_parameters())
  {
    if (disp) std::cout << "Parameter name: " << p.name << "\n";
    if (p.name.find("bias") != std::string::npos)
    {
      if (bias)
      {
        if (disp) std::cerr << "Module has multiple bias as parameter name\n";
        return false;
      }
      offset = 0;
      bias   = true;
    }
    else if (p.name.find("weight") != std::string::npos)
    {
      if (weight)
      {
        if (disp) std::cerr << "Module has multiple weight as parameter name\n";
        return false;
      }
      offset = 1;
      weight = true;
    }
    else
    {
      if (disp)
        std::cerr << "Module has unrecognized parameter name: " << p.name
                  << "\n";
      return false;
    }

    // p.value is an at::Tensor (PyTorch tensor in the C++/ATen API).
    if (!_tensor_to_vector(layer, p.value, offset, disp)) return false;
    if (bias && weight)
    {
      append_data(layer);
      bias = weight = false;
      layer.clear();
    }
  }

  // Retrieve activation functions
  size_t i0   = 0;
  bool linear = false;
  for (const auto& child : net.named_children())
  {
    if (i0 == _data.size())
    {
      if (disp)
        std::cerr << "Module has too many submodules for " << _data.size()
                  << " layers\n";
      return false;
    }

    const auto& submod   = child.value;
    std::string type_str = submod.type()->str();
    if (disp)
      std::cout << "Submodule '" << child.name << "'  type: " << type_str;
    if (type_str.find("Linear") != std::string::npos)
    {
      if (disp) std::cout << "  (Linear layer)\n";
      if (linear || i0 + 1 == _data.size()) _data[i0++].second = LINEAR;
    }
    else if (type_str.find("Tanh") != std::string::npos)
    {
      if (disp) std::cout << "  (Activation: Tanh)\n";
      _data[i0++].second = TANH;
      linear             = false;
    }
    else if (type_str.find("ReLU") != std::string::npos)
    {
      if (disp) std::cout << "  (Activation: ReLU)\n";
      _data[i0++].second = RELU;
      linear             = false;
    }
    else if (type_str.find("Sigmoid") != std::string::npos)
    {
      if (disp) std::cout << "  (Activation: Sigmoid)\n";
      _data[i0++].second = SIGMOID;
      linear             = false;
    }
    else
    {
      if (disp)
        std::cerr << "\nNon-supported activation function: " << type_str
                  << "\n";
      return false;
    }
  }
  std::cout << "\n";

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
  if (t.sizes().size() > 2) return false;
  at::Tensor tc = t.contiguous();

  // 1D tensor
  if (t.sizes().size() == 1)
  {
    size_t n0 = t.sizes()[0];  // size of dim 0
    if (v.size() < n0) v.resize(n0);

    switch (tc.scalar_type())
    {
      case at::kFloat:
      {
        float* data_ptr = tc.data_ptr<float>();
        for (size_t i0 = 0; i0 < n0; ++i0)
        {
          if (v[i0].size() < offset) v[i0].resize(offset);
          v[i0][offset] = static_cast<double>(data_ptr[i0]);
        }
        break;
      }
      case at::kDouble:
      {
        double* data_ptr = tc.data_ptr<double>();
        for (size_t i0 = 0; i0 < n0; ++i0)
        {
          if (v[i0].size() < offset) v[i0].resize(offset);
          v[i0][offset] = data_ptr[i0];
        }
        break;
      }
      // add kInt, kLong, ... as needed
      default:
        if (disp)
          std::cerr << "Unsupported tensor type " << tc.scalar_type() << "\n";
        return false;
    }

    return true;
  }

  // 2D tensor
  size_t n0 = t.sizes()[0];  // size of dim 0
  size_t n1 = t.sizes()[1];  // size of dim 1
  if (v.size() < n0) v.resize(n0);

  switch (tc.scalar_type())
  {
    case at::kFloat:
    {
      float* data_ptr = tc.data_ptr<float>();
      for (size_t i0 = 0, i01 = 0; i0 < n0; ++i0, i01 += n1)
      {
        if (v[i0].size() < offset + n1) v[i0].resize(offset + n1);
        for (size_t i1 = 0; i1 < n1; ++i1)
          v[i0][offset + i1] = static_cast<double>(data_ptr[i01 + i1]);
      }
      break;
    }
    case at::kDouble:
    {
      double* data_ptr = tc.data_ptr<double>();
      for (size_t i0 = 0, i01 = 0; i0 < n0; ++i0, i01 += n1)
      {
        if (v[i0].size() < offset + n1) v[i0].resize(offset + n1);
        for (size_t i1 = 0; i1 < n1; ++i1)
          v[i0][offset + i1] = data_ptr[i01 + i1];
      }
      break;
    }
    // add kInt, kLong, ... as needed
    default:
      if (disp)
        std::cerr << "Unsupported tensor type " << tc.scalar_type() << "\n";
      return false;
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
    if (idU == typeid(FFVar))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (idU == typeid(FADType<FFVar>))
      return eval(nRes, static_cast<FADType<FFVar>*>(vRes), nVar,
                  static_cast<FADType<FFVar> const*>(vVar), mVar);
    else if (idU == typeid(FFDep))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (idU == typeid(FFInv))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (idU == typeid(double))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (idU == typeid(FADType<double>))
      return eval(nRes, static_cast<FADType<double>*>(vRes), nVar,
                  static_cast<FADType<double> const*>(vVar), mVar);
    else if (idU == typeid(T))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (idU == typeid(Specbnd<T>))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (idU == typeid(SCVar<T>))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (idU == typeid(McCormick<T>))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (idU == typeid(SupVar<PWCU>))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (idU == typeid(SupVar<PWLU>))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (idU == typeid(PolVar<T>))
      return eval(nRes, static_cast<PolVar<T>*>(vRes), nVar,
                  static_cast<PolVar<T> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<double>))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<FADType<double>>))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<FFDep>))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (idU == typeid(SLiftVar))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (idU == typeid(FFExpr))
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
    if (idU == typeid(T))
      return reval(nRes, static_cast<T*>(vRes), nVar, static_cast<T*>(vVar));
    else if (idU == typeid(PolVar<T>))
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
    if (idU == typeid(FFVar))
      return eval(nRes, static_cast<FFVar*>(vRes), nVar,
                  static_cast<FFVar const*>(vVar), mVar);
    else if (idU == typeid(FFDep))
      return eval(nRes, static_cast<FFDep*>(vRes), nVar,
                  static_cast<FFDep const*>(vVar), mVar);
    else if (idU == typeid(FFInv))
      return eval(nRes, static_cast<FFInv*>(vRes), nVar,
                  static_cast<FFInv const*>(vVar), mVar);
    else if (idU == typeid(double))
      return eval(nRes, static_cast<double*>(vRes), nVar,
                  static_cast<double const*>(vVar), mVar);
    else if (idU == typeid(T))
      return eval(nRes, static_cast<T*>(vRes), nVar,
                  static_cast<T const*>(vVar), mVar);
    else if (idU == typeid(Specbnd<T>))
      return eval(nRes, static_cast<Specbnd<T>*>(vRes), nVar,
                  static_cast<Specbnd<T> const*>(vVar), mVar);
    else if (idU == typeid(SCVar<T>))
      return eval(nRes, static_cast<SCVar<T>*>(vRes), nVar,
                  static_cast<SCVar<T> const*>(vVar), mVar);
    else if (idU == typeid(McCormick<T>))
      return eval(nRes, static_cast<McCormick<T>*>(vRes), nVar,
                  static_cast<McCormick<T> const*>(vVar), mVar);
    else if (idU == typeid(SupVar<PWCU>))
      return eval(nRes, static_cast<SupVar<PWCU>*>(vRes), nVar,
                  static_cast<SupVar<PWCU> const*>(vVar), mVar);
    else if (idU == typeid(SupVar<PWLU>))
      return eval(nRes, static_cast<SupVar<PWLU>*>(vRes), nVar,
                  static_cast<SupVar<PWLU> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<double>))
      return eval(nRes, static_cast<OCVar<double>*>(vRes), nVar,
                  static_cast<OCVar<double> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<FADType<double>>))
      return eval(nRes, static_cast<OCVar<FADType<double>>*>(vRes), nVar,
                  static_cast<OCVar<FADType<double>> const*>(vVar), mVar);
    else if (idU == typeid(OCVar<FFDep>))
      return eval(nRes, static_cast<OCVar<FFDep>*>(vRes), nVar,
                  static_cast<OCVar<FFDep> const*>(vVar), mVar);
    else if (idU == typeid(SLiftVar))
      return eval(nRes, static_cast<SLiftVar*>(vRes), nVar,
                  static_cast<SLiftVar const*>(vVar), mVar);
    else if (idU == typeid(FFExpr))
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
