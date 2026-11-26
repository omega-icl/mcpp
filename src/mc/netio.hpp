// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__NETIO_HPP
#define MC__NETIO_HPP

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#if defined(MC__USE_TORCH)
#include <torch/script.h>
#endif

namespace mc
{

//! @brief C++ class reading the layers of a network -- a multilayer perceptron
//! or a Kolmogorov-Arnold network -- from a file
////////////////////////////////////////////////////////////////////////
//! mc::NETIO reads the dense layers of a multilayer perceptron, and their
//! activation functions, for mc::MLP and mc::MLPREG:
//! - from an ONNX file (read_onnx), with a decoder of the protobuf format
//!   of its own and no other library. The graph must be a chain of nodes
//!   from its single input to its single output: Gemm (transA = 0), or
//!   MatMul followed by Add, for each dense layer, with the weights and
//!   biases as initializers or Constant nodes; BatchNormalization after a
//!   dense layer, folded into it; Relu, Tanh or Sigmoid after a dense
//!   layer; and Identity, Dropout, Flatten, Reshape, Squeeze, Unsqueeze
//!   and Cast, which leave a vector of features as it is. Any other node
//!   is refused, as are weights stored outside the file and weights in
//!   other types than float, double and float16. This covers the networks
//!   that PyTorch (torch.onnx.export), scikit-learn (skl2onnx), Keras
//!   (tf2onnx) and others export;
//! - from a TorchScript Module (read_torch_net, with Torch), whose
//!   submodule 'net' is a Sequential of Linear layers, each followed by at
//!   most one activation function (Tanh, ReLU, Sigmoid), Identity and
//!   Dropout aside.
//! - from a dictionary saved by torch.save (read_dict), with a reader of
//!   the archive and of the pickle format of its own and without Torch: a
//!   dense network (ANNmodel) or a Kolmogorov-Arnold network with
//!   piecewise-linear edges (KANmodel, for mc::KAN), with the scalings of
//!   its inputs and outputs.
//! A network is returned as the data of mc::MLP::set_data: one entry per
//! layer, holding one row per neuron -- its bias followed by its weights --
//! and its activation function (LINEAR, RELU, TANH, SIGMOID, with the
//! values of mc::MLP::ACTIV_TYPE and mc::MLPREG::ACTIV_TYPE).
////////////////////////////////////////////////////////////////////////
class NETIO
////////////////////////////////////////////////////////////////////////
{
 public:
  //! @brief Activation functions, with the values of mc::MLP::ACTIV_TYPE
  enum ACTIV_TYPE
  {
    LINEAR = 0,  //!< Linear activation function
    RELU,        //!< ReLU activation function
    TANH,        //!< tanh activation function
    SIGMOID      //!< Sigmoid activation function
  };

  //! @brief Layers of a network: neurons -- bias, then weights -- and
  //! activation function of each layer
  typedef std::vector<std::pair<std::vector<std::vector<double>>, int>> t_Data;

  //! @brief Read the layers of the network in the ONNX file filename into
  //! data: false, with the reason in msg, where the file cannot be read or
  //! holds anything else than a network of dense layers (see mc::NETIO)
  static bool
  read_onnx(std::string const& filename, t_Data& data, std::string& msg);

  //! @brief Read the layers of the network in an ONNX model held in memory,
  //! the bytes of an ONNX file
  static bool
  read_onnx_bytes(std::string const& bytes, t_Data& data, std::string& msg);

#if defined(MC__USE_TORCH)
  //! @brief Read the layers of the submodule 'net' of a TorchScript Module
  //! into data: false, with the reason in msg, where it is anything else than
  //! a Sequential of Linear layers each followed by at most one activation
  //! function (Tanh, ReLU, Sigmoid), Identity and Dropout aside. With disp,
  //! the submodules and parameters are displayed
  static bool
  read_torch_net(torch::jit::script::Module const& module, t_Data& data,
                 std::string& msg, bool const disp = false);

  //! @brief Linear submodules of the submodule 'net' of a TorchScript Module,
  //! in their order: their weight and bias (undefined where there is none)
  static std::vector<std::pair<at::Tensor, at::Tensor>>
  torch_linear(torch::jit::script::Module const& module);

  //! @brief Rows of the tensor t, a vector or a matrix, written in the
  //! columns of v from column offset -- the rows of v being made, or as many
  //! already
  static bool
  tensor_to_rows(std::vector<std::vector<double>>& v, at::Tensor const& t,
                 size_t const offset, bool const disp = false);
#endif

  //! @brief Whether the name of the file filename ends with .onnx (whatever
  //! the case)
  static bool
  is_onnx(std::string const& filename)
  {
    if (filename.size() < 5) return false;
    std::string ext = filename.substr(filename.size() - 5);
    for (auto& c : ext) c = (char)std::tolower((unsigned char)c);
    return ext == ".onnx";
  }

  //! @brief Check the layers of a network: false, with the reason in msg,
  //! where a layer is empty or its neurons do not have one bias and one
  //! weight per neuron of the previous layer
  static bool
  check(t_Data const& data, std::string& msg);

  //! @brief Value held in a file saved by torch.save: none, a boolean, an
  //! integer, a real, a string, a list or tuple, a dictionary, a tensor, or
  //! any other object (kept by its name only)
  struct PValue
  {
    //! @brief Kinds of values
    enum KIND
    {
      NONE = 0,  //!< None
      BOOL,      //!< Boolean, in b
      INT,       //!< Integer, in n
      REAL,      //!< Real, in x
      STR,       //!< String, in s
      LIST,      //!< List or tuple, in list
      DICT,      //!< Dictionary, in dict, in the order of the file
      TENSOR,    //!< Tensor: dimensions in dims, entries in val in row-major order
      OTHER      //!< Any other object: its name in s
    };
    int kind = NONE;
    bool b = false;
    long long n = 0;
    double x = 0.;
    std::string s;
    std::vector<std::shared_ptr<PValue>> list;
    std::vector<std::pair<std::shared_ptr<PValue>, std::shared_ptr<PValue>>> dict;
    std::vector<size_t> dims;
    std::vector<double> val;

    //! @brief Entry of a dictionary with the string key, nullptr if this is
    //! not a dictionary or has no such entry
    PValue const*
    at(std::string const& key) const
    {
      if (kind != DICT) return nullptr;
      for (auto const& [k, v] : dict)
        if (k && k->kind == STR && k->s == key) return v.get();
      return nullptr;
    }
    //! @brief Whether this is an integer or a real, with its value in v
    bool
    number(double& v) const
    {
      if (kind != INT && kind != REAL) return false;
      v = (kind == INT ? (double)n : x);
      return true;
    }
    //! @brief Whether this is a list of numbers, with their values in v
    bool
    numbers(std::vector<double>& v) const
    {
      if (kind != LIST) return false;
      v.resize(list.size());
      for (size_t i = 0; i < list.size(); ++i)
        if (!list[i] || !list[i]->number(v[i])) return false;
      return true;
    }
  };

  //! @brief Network saved as a dictionary by torch.save: see read_dict
  struct Dict
  {
    //! @brief Layer of a Kolmogorov-Arnold network with edges \f$c + v u +
    //! \sum_{r<R} a_r \max(w_r u + b_r, 0)\f$: c and v indexed [j*nin+i], a, w
    //! and b indexed [(j*nin+i)*R+r], for the edge from input i to output j
    struct KANLayer
    {
      size_t nin = 0, nout = 0, R = 0;
      std::vector<double> c, v, a, w, b;
    };
    //! @brief Type of the model: "KANmodel" or "ANNmodel"
    std::string type;
    //! @brief Numbers of inputs, of nodes in each hidden layer, and of outputs
    std::vector<size_t> widths;
    //! @brief Layers of a KANmodel
    std::vector<KANLayer> kan;
    //! @brief Layers of an ANNmodel, as for mc::MLP::set_data
    t_Data ann;
    //! @brief Min-max scaling of the inputs, u = (x - xlower) / (xupper -
    //! xlower): empty where the inputs are not scaled
    std::vector<double> xlower, xupper;
    //! @brief Min-max scaling of the outputs, z = zlower + (zupper - zlower)
    //! y: empty where the outputs are not scaled
    std::vector<double> zlower, zupper;
    //! @brief Whole content of the file, e.g. root->at("bounds")
    std::shared_ptr<PValue> root;
  };

  //! @brief Whether the file filename is an archive of torch.save holding a
  //! pickled object (data.pkl) and no TorchScript Module (constants.pkl)
  static bool is_torch_dict(std::string const& filename);

  //! @brief Read the object saved by torch.save in the file filename into
  //! root, without Torch: false, with the reason in msg, where the file is not
  //! such an archive, is compressed, or holds values other than those of
  //! mc::NETIO::PValue -- tensors of floats, doubles, half floats, integers
  //! and booleans, as contiguous or strided views of their storage
  static bool read_pickle(std::string const& filename, std::shared_ptr<PValue>& root,
                          std::string& msg);

  //! @brief Read the network saved as a dictionary by torch.save in the file
  //! filename, without Torch. The dictionary holds
  //! - model_type: "KANmodel" or "ANNmodel";
  //! - config: widths, the numbers of inputs, of nodes of each hidden layer
  //!   and of outputs; for a KANmodel, intervals, the number R+1 of intervals
  //!   of the edges of each layer; for an ANNmodel, optionally activation,
  //!   the activation function of the hidden layers ("relu", "tanh" or
  //!   "sigmoid");
  //! - state_dict: for each layer l of a KANmodel, the tensors
  //!   layers.l.edge_bias (c, nout x nin), layers.l.linear_weight (v, nout x
  //!   nin x 1), layers.l.relu_output_weight (a), layers.l.relu_weight (w)
  //!   and layers.l.relu_bias (b, nout x nin x R each); for each layer of an
  //!   ANNmodel, layers.l.weight (nout x nin) and layers.l.bias (nout), the
  //!   output layer being linear;
  //! - normalization, optionally: x and z, each None or a dictionary with
  //!   the lists lower and upper of a min-max scaling.
  //! Where an ANNmodel does not name its activation function, it is the one
  //! given by activ (RELU, TANH or SIGMOID), or, with activ negative, the one
  //! that matches the bounds of the first hidden layer in the file (entry
  //! bounds, tensors pre and post of the first layer of a bound set). False,
  //! with the reason in msg, where the file cannot be read or holds anything
  //! else
  static bool read_dict(std::string const& filename, Dict& dict, std::string& msg,
                        int const activ = -1);

  //! @brief Fold the min-max scalings of dict into the first and last layers
  //! of the network data, which then maps raw inputs to raw outputs: false,
  //! with the reason in msg, where the last layer is not linear
  static bool fold_scaling(t_Data& data, Dict const& dict, std::string& msg);

 private:
  //! @brief Entries of a zip archive held in memory: offset and size of the
  //! data of each, the offset being npos for a compressed entry
  static bool _zip_entries(std::string const& bytes,
                           std::map<std::string, std::pair<size_t, size_t>>& entries,
                           std::string& msg);
  //! @brief Bytes of a file
  static bool _file_bytes(std::string const& filename, std::string& bytes);

  //! @brief Reader of the protobuf wire format over a range of bytes
  struct PB
  {
    uint8_t const* p;
    uint8_t const* e;

    static void
    fail(std::string const& what)
    {
      throw std::runtime_error("corrupted protobuf message (" + what + ")");
    }
    bool
    done() const
    {
      return p >= e;
    }
    void
    need(size_t const n) const
    {
      if ((size_t)(e - p) < n) fail("truncated field");
    }
    uint64_t
    varint()
    {
      uint64_t v = 0;
      for (int s = 0;; s += 7)
      {
        if (p >= e || s > 63) fail("varint");
        uint8_t const b = *p++;
        v |= uint64_t(b & 0x7f) << s;
        if (!(b & 0x80)) break;
      }
      return v;
    }
    uint32_t
    fixed32()
    {
      need(4);
      uint32_t v = 0;
      for (int k = 3; k >= 0; --k) v = (v << 8) | p[k];
      p += 4;
      return v;
    }
    uint64_t
    fixed64()
    {
      need(8);
      uint64_t v = 0;
      for (int k = 7; k >= 0; --k) v = (v << 8) | p[k];
      p += 8;
      return v;
    }
    // next field: its number and wire type; false at the end
    bool
    field(uint32_t& num, uint32_t& wt)
    {
      if (p >= e) return false;
      uint64_t const key = varint();
      num                = (uint32_t)(key >> 3);
      wt                 = (uint32_t)(key & 7);
      return true;
    }
    // length-delimited field, as a reader of its own
    PB
    sub()
    {
      uint64_t const n = varint();
      need(n);
      PB s{p, p + n};
      p += n;
      return s;
    }
    std::string
    str()
    {
      PB s = sub();
      return std::string((char const*)s.p, s.e - s.p);
    }
    void
    skip(uint32_t const wt)
    {
      switch (wt)
      {
        case 0:
          varint();
          break;
        case 1:
          need(8);
          p += 8;
          break;
        case 2:
          sub();
          break;
        case 5:
          need(4);
          p += 4;
          break;
        default:
          fail("wire type " + std::to_string(wt));
      }
    }
  };

  static double
  _f32(uint32_t const u)
  {
    float f;
    std::memcpy(&f, &u, 4);
    return f;
  }
  static double
  _f64(uint64_t const u)
  {
    double d;
    std::memcpy(&d, &u, 8);
    return d;
  }
  static double
  _f16(uint16_t const h)
  {
    int const s = (h >> 15) & 1, x = (h >> 10) & 0x1f, m = h & 0x3ff;
    double v = (x == 0    ? std::ldexp((double)m, -24)
                : x == 31 ? (m ? NAN : INFINITY)
                          : std::ldexp((double)(m | 0x400), x - 25));
    return s ? -v : v;
  }

  // ONNX TensorProto, values converted to double
  struct Tensor
  {
    std::string name;
    std::vector<int64_t> dims;
    int dtype     = 0;
    bool external = false;
    std::vector<double> val;
    size_t
    numel() const
    {
      size_t n = 1;
      for (auto d : dims) n *= (size_t)d;
      return n;
    }
  };

  // ONNX AttributeProto
  struct Attr
  {
    std::string name;
    double f  = 0.;
    int64_t i = 0;
    std::string s;
    std::vector<int64_t> ints;
    std::vector<double> floats;
    Tensor t;
    bool hast = false;
  };

  // ONNX NodeProto
  struct Node
  {
    std::string op, domain, name;
    std::vector<std::string> in, out;
    std::vector<Attr> attr;
    Attr const*
    get(std::string const& n) const
    {
      for (auto const& a : attr)
        if (a.name == n) return &a;
      return nullptr;
    }
    int64_t
    geti(std::string const& n, int64_t const def) const
    {
      Attr const* a = get(n);
      return a ? a->i : def;
    }
    double
    getf(std::string const& n, double const def) const
    {
      Attr const* a = get(n);
      return a ? a->f : def;
    }
    std::string
    label() const
    {
      return op + (name.empty() ? std::string() : " '" + name + "'");
    }
  };

  // ONNX GraphProto
  struct Graph
  {
    std::vector<Node> nodes;
    std::map<std::string, Tensor> init;
    std::vector<std::string> inputs, outputs;
    std::map<std::string, std::vector<int64_t>> shapes;  // -1: unknown
    bool sparse = false;
  };

  static void _read_tensor(PB pb, Tensor& t);
  static void _read_attr(PB pb, Attr& a);
  static void _read_node(PB pb, Node& n);
  static void _read_valueinfo(PB pb, std::string& name,
                              std::vector<int64_t>& shape);
  static void _read_graph(PB pb, Graph& g);
  static bool _chain(Graph const& g, t_Data& data, std::string& msg);
};

inline bool
NETIO::check(t_Data const& data, std::string& msg)
{
  if (data.empty())
  {
    msg = "no dense layer";
    return false;
  }
  size_t nprv = 0;
  for (size_t l = 0; l < data.size(); ++l)
  {
    auto const& [layer, activ] = data[l];
    if (layer.empty() || layer[0].size() < 2 || (l && layer[0].size() != nprv + 1))
    {
      msg = "layer " + std::to_string(l) + " is empty or does not match the previous one";
      return false;
    }
    for (auto const& neuron : layer)
      if (neuron.size() != layer[0].size())
      {
        msg = "layer " + std::to_string(l) + " has neurons of different sizes";
        return false;
      }
    if (activ < LINEAR || activ > SIGMOID)
    {
      msg = "layer " + std::to_string(l) + " has an unknown activation function";
      return false;
    }
    nprv = layer.size();
  }
  return true;
}

inline void
NETIO::_read_tensor(PB pb, Tensor& t)
{
  std::string raw;
  std::vector<double> fdata, ddata, idata;
  uint32_t num, wt;
  while (pb.field(num, wt))
  {
    switch (num)
    {
      case 1:  // dims
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) t.dims.push_back((int64_t)s.varint());
        }
        else
          t.dims.push_back((int64_t)pb.varint());
        break;
      case 2:  // data_type
        t.dtype = (int)pb.varint();
        break;
      case 4:  // float_data
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) fdata.push_back(_f32(s.fixed32()));
        }
        else
          fdata.push_back(_f32(pb.fixed32()));
        break;
      case 5:  // int32_data -- float16 values among them
      case 7:  // int64_data
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) idata.push_back((double)(int64_t)s.varint());
        }
        else
          idata.push_back((double)(int64_t)pb.varint());
        break;
      case 8:  // name
        t.name = pb.str();
        break;
      case 9:  // raw_data
        raw = pb.str();
        break;
      case 10:  // double_data
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) ddata.push_back(_f64(s.fixed64()));
        }
        else
          ddata.push_back(_f64(pb.fixed64()));
        break;
      case 13:  // external_data
        t.external = true;
        pb.skip(wt);
        break;
      case 14:  // data_location
        if (pb.varint() == 1) t.external = true;
        break;
      default:
        pb.skip(wt);
    }
  }

  size_t const n = t.numel();
  // little-endian bytes, whatever the host
  auto rd = [&](size_t const k, size_t const w) -> uint64_t
  {
    uint64_t v = 0;
    for (size_t b = w; b-- > 0;) v = (v << 8) | (uint8_t)raw[k * w + b];
    return v;
  };
  t.val.clear();
  switch (t.dtype)
  {
    case 1:  // FLOAT
      if (!raw.empty())
      {
        if (raw.size() != 4 * n) PB::fail("size of tensor '" + t.name + "'");
        for (size_t k = 0; k < n; ++k) t.val.push_back(_f32((uint32_t)rd(k, 4)));
      }
      else
        t.val = fdata;
      break;
    case 11:  // DOUBLE
      if (!raw.empty())
      {
        if (raw.size() != 8 * n) PB::fail("size of tensor '" + t.name + "'");
        for (size_t k = 0; k < n; ++k) t.val.push_back(_f64(rd(k, 8)));
      }
      else
        t.val = ddata;
      break;
    case 10:  // FLOAT16
      if (!raw.empty())
      {
        if (raw.size() != 2 * n) PB::fail("size of tensor '" + t.name + "'");
        for (size_t k = 0; k < n; ++k) t.val.push_back(_f16((uint16_t)rd(k, 2)));
      }
      else
        for (double const v : idata) t.val.push_back(_f16((uint16_t)v));
      break;
    case 6:  // INT32
    case 7:  // INT64
      if (!raw.empty())
      {
        size_t const w = (t.dtype == 6 ? 4 : 8);
        if (raw.size() != w * n) PB::fail("size of tensor '" + t.name + "'");
        for (size_t k = 0; k < n; ++k)
          t.val.push_back(w == 4 ? (double)(int32_t)rd(k, 4)
                                 : (double)(int64_t)rd(k, 8));
      }
      else
        t.val = idata;
      break;
    default:  // other types: no values, refused if used as weights
      break;
  }
}

inline void
NETIO::_read_attr(PB pb, Attr& a)
{
  uint32_t num, wt;
  while (pb.field(num, wt))
  {
    switch (num)
    {
      case 1:
        a.name = pb.str();
        break;
      case 2:
        a.f = _f32(pb.fixed32());
        break;
      case 3:
        a.i = (int64_t)pb.varint();
        break;
      case 4:
        a.s = pb.str();
        break;
      case 5:
        _read_tensor(pb.sub(), a.t);
        a.hast = true;
        break;
      case 7:
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) a.floats.push_back(_f32(s.fixed32()));
        }
        else
          a.floats.push_back(_f32(pb.fixed32()));
        break;
      case 8:
        if (wt == 2)
        {
          PB s = pb.sub();
          while (!s.done()) a.ints.push_back((int64_t)s.varint());
        }
        else
          a.ints.push_back((int64_t)pb.varint());
        break;
      default:
        pb.skip(wt);
    }
  }
}

inline void
NETIO::_read_node(PB pb, Node& n)
{
  uint32_t num, wt;
  while (pb.field(num, wt))
  {
    switch (num)
    {
      case 1:
        n.in.push_back(pb.str());
        break;
      case 2:
        n.out.push_back(pb.str());
        break;
      case 3:
        n.name = pb.str();
        break;
      case 4:
        n.op = pb.str();
        break;
      case 5:
        n.attr.emplace_back();
        _read_attr(pb.sub(), n.attr.back());
        break;
      case 7:
        n.domain = pb.str();
        break;
      default:
        pb.skip(wt);
    }
  }
}

inline void
NETIO::_read_valueinfo(PB pb, std::string& name, std::vector<int64_t>& shape)
{
  uint32_t num, wt;
  while (pb.field(num, wt))
  {
    if (num == 1)
      name = pb.str();
    else if (num == 2)
    {
      // TypeProto: tensor_type (1) -> shape (2) -> dim (1) -> dim_value (1)
      PB tp = pb.sub();
      uint32_t n1, w1;
      while (tp.field(n1, w1))
      {
        if (n1 != 1)
        {
          tp.skip(w1);
          continue;
        }
        PB tt = tp.sub();
        uint32_t n2, w2;
        while (tt.field(n2, w2))
        {
          if (n2 != 2)
          {
            tt.skip(w2);
            continue;
          }
          PB sh = tt.sub();
          uint32_t n3, w3;
          while (sh.field(n3, w3))
          {
            if (n3 != 1)
            {
              sh.skip(w3);
              continue;
            }
            PB dim = sh.sub();
            int64_t d = -1;
            uint32_t n4, w4;
            while (dim.field(n4, w4))
              if (n4 == 1)
                d = (int64_t)dim.varint();
              else
                dim.skip(w4);
            shape.push_back(d);
          }
        }
      }
    }
    else
      pb.skip(wt);
  }
}

inline void
NETIO::_read_graph(PB pb, Graph& g)
{
  uint32_t num, wt;
  while (pb.field(num, wt))
  {
    switch (num)
    {
      case 1:
        g.nodes.emplace_back();
        _read_node(pb.sub(), g.nodes.back());
        break;
      case 5:
      {
        Tensor t;
        _read_tensor(pb.sub(), t);
        std::string const name = t.name;
        g.init[name] = std::move(t);
        break;
      }
      case 11:
      case 12:
      {
        std::string name;
        std::vector<int64_t> shape;
        _read_valueinfo(pb.sub(), name, shape);
        (num == 11 ? g.inputs : g.outputs).push_back(name);
        g.shapes[name] = shape;
        break;
      }
      case 15:  // sparse_initializer
        g.sparse = true;
        pb.skip(wt);
        break;
      default:
        pb.skip(wt);
    }
  }
}

inline bool
NETIO::_chain(Graph const& g, t_Data& data, std::string& msg)
{
  // Constant tensors: the initializers, and the outputs of Constant nodes
  std::map<std::string, Tensor> cst = g.init;
  for (auto const& n : g.nodes)
  {
    if (n.op != "Constant" || n.out.size() != 1) continue;
    Attr const* a = n.get("value");
    Tensor t;
    if (a && a->hast)
      t = a->t;
    else if ((a = n.get("value_float")))
    {
      t.dtype = 1;
      t.val   = {a->f};
    }
    else if ((a = n.get("value_floats")))
    {
      t.dtype = 1;
      t.dims  = {(int64_t)a->floats.size()};
      t.val   = a->floats;
    }
    else if ((a = n.get("value_int")))
    {
      t.dtype = 7;
      t.val   = {(double)a->i};
    }
    else if ((a = n.get("value_ints")))
    {
      t.dtype = 7;
      t.dims  = {(int64_t)a->ints.size()};
      for (auto v : a->ints) t.val.push_back((double)v);
    }
    else
    {
      msg = "Constant node without a value that can be read";
      return false;
    }
    cst[n.out[0]] = t;
  }

  // The input of the network: the input of the graph that is no constant
  std::string cur;
  for (auto const& name : g.inputs)
  {
    if (cst.count(name)) continue;
    if (!cur.empty())
    {
      msg = "the graph has several inputs ('" + cur + "', '" + name + "')";
      return false;
    }
    cur = name;
  }
  if (cur.empty() || g.outputs.size() != 1)
  {
    msg = "the graph must have one input and one output";
    return false;
  }
  int64_t nin = -1;
  {
    auto it = g.shapes.find(cur);
    if (it != g.shapes.end() && !it->second.empty()) nin = it->second.back();
  }

  // Values of a constant tensor that are weights
  auto weights = [&](std::string const& name, Tensor const*& t) -> bool
  {
    auto it = cst.find(name);
    if (it == cst.end())
    {
      msg = "'" + name + "' is not a constant tensor";
      return false;
    }
    t = &it->second;
    if (t->external)
    {
      msg = "tensor '" + name + "' is stored outside the file";
      return false;
    }
    if (t->dtype != 1 && t->dtype != 10 && t->dtype != 11)
    {
      msg = "tensor '" + name + "' is not of a floating-point type";
      return false;
    }
    if (t->val.size() != t->numel())
    {
      msg = "tensor '" + name + "' has " + std::to_string(t->val.size()) +
            " values for " + std::to_string(t->numel()) + " elements";
      return false;
    }
    return true;
  };

  // Number of features after the latest layer: nin before any
  auto nfeat = [&]() -> int64_t
  { return data.empty() ? nin : (int64_t)data.back().first.size(); };

  // Along the chain: whether the latest layer still takes a bias, a batch
  // normalization or an activation
  bool open = false;
  data.clear();
  for (auto const& n : g.nodes)
  {
    if (n.op == "Constant") continue;
    if (!n.domain.empty() && n.domain != "ai.onnx")
    {
      msg = "node " + n.label() + " of domain '" + n.domain + "' is not supported";
      return false;
    }
    if (n.in.empty() || n.out.empty())
    {
      msg = "node " + n.label() + " without input or output";
      return false;
    }

    // The data path enters each node through its first input, except for an
    // Add whose first input is the bias
    size_t const ipath = (n.op == "Add" && n.in.size() == 2 && n.in[1] == cur ? 1 : 0);
    if (n.in[ipath] != cur)
    {
      msg = "node " + n.label() + " is not on the path from the input to the output: "
            "only a chain of layers is supported";
      return false;
    }

    if (n.op == "Gemm" || n.op == "MatMul")
    {
      Tensor const* B = nullptr;
      if (n.in.size() < 2 || !weights(n.in[1], B)) return false;
      if (B->dims.size() != 2)
      {
        msg = "node " + n.label() + ": weights of " + std::to_string(B->dims.size()) + " dimensions";
        return false;
      }
      bool const transB  = (n.op == "Gemm" && n.geti("transB", 0));
      double const alpha = (n.op == "Gemm" ? n.getf("alpha", 1.) : 1.);
      double const beta  = (n.op == "Gemm" ? n.getf("beta", 1.) : 1.);
      if (n.op == "Gemm" && n.geti("transA", 0))
      {
        msg = "node " + n.label() + " with transA = 1";
        return false;
      }
      size_t const nneu = (size_t)(transB ? B->dims[0] : B->dims[1]);
      size_t const nprv = (size_t)(transB ? B->dims[1] : B->dims[0]);
      if (!nneu || !nprv || (nfeat() >= 0 && (size_t)nfeat() != nprv))
      {
        msg = "node " + n.label() + " takes " + std::to_string(nprv) +
              " features where there are " + std::to_string(nfeat());
        return false;
      }
      std::vector<std::vector<double>> layer(nneu, std::vector<double>(nprv + 1, 0.));
      for (size_t i = 0; i < nneu; ++i)
        for (size_t j = 0; j < nprv; ++j)
          layer[i][1 + j] = alpha * B->val[transB ? i * nprv + j : j * nneu + i];
      if (n.op == "Gemm" && n.in.size() > 2 && !n.in[2].empty())
      {
        Tensor const* C = nullptr;
        if (!weights(n.in[2], C)) return false;
        if (C->val.size() != nneu && C->val.size() != 1)
        {
          msg = "node " + n.label() + ": bias of " + std::to_string(C->val.size()) +
                " values for " + std::to_string(nneu) + " neurons";
          return false;
        }
        for (size_t i = 0; i < nneu; ++i)
          layer[i][0] = beta * C->val[C->val.size() == 1 ? 0 : i];
      }
      data.push_back(std::make_pair(layer, (int)LINEAR));
      open = true;
    }

    else if (n.op == "Add")
    {
      Tensor const* C = nullptr;
      if (n.in.size() == 2 && !cst.count(n.in[1 - ipath]))
      {
        msg = "node " + n.label() + " adds two values computed by the network: "
              "only a chain of layers is supported";
        return false;
      }
      if (n.in.size() != 2 || !weights(n.in[1 - ipath], C)) return false;
      size_t const nneu = (open ? data.back().first.size() : 0);
      if (!open || (C->val.size() != nneu && C->val.size() != 1))
      {
        msg = "node " + n.label() + " is not the bias of a dense layer";
        return false;
      }
      for (size_t i = 0; i < nneu; ++i)
        data.back().first[i][0] += C->val[C->val.size() == 1 ? 0 : i];
    }

    else if (n.op == "BatchNormalization")
    {
      Tensor const *S = nullptr, *Bb = nullptr, *M = nullptr, *V = nullptr;
      if (n.in.size() < 5 || !weights(n.in[1], S) || !weights(n.in[2], Bb) ||
          !weights(n.in[3], M) || !weights(n.in[4], V))
        return false;
      size_t const nneu = (open ? data.back().first.size() : 0);
      if (!open || S->val.size() != nneu || Bb->val.size() != nneu ||
          M->val.size() != nneu || V->val.size() != nneu)
      {
        msg = "node " + n.label() + " does not follow a dense layer of its size";
        return false;
      }
      double const eps = n.getf("epsilon", 1e-5);
      for (size_t i = 0; i < nneu; ++i)
      {
        double const s = S->val[i] / std::sqrt(V->val[i] + eps);
        auto& neuron   = data.back().first[i];
        for (auto& w : neuron) w *= s;
        neuron[0] += Bb->val[i] - M->val[i] * s;
      }
    }

    else if (n.op == "Relu" || n.op == "Tanh" || n.op == "Sigmoid")
    {
      if (!open)
      {
        msg = "activation " + n.label() + " without a dense layer before it";
        return false;
      }
      data.back().second = (n.op == "Relu" ? RELU : n.op == "Tanh" ? TANH : SIGMOID);
      open = false;
    }

    else if (n.op == "Identity" || n.op == "Dropout" || n.op == "Flatten" ||
             n.op == "Reshape" || n.op == "Squeeze" || n.op == "Unsqueeze" ||
             n.op == "Cast")
    {
      // a vector of features stays as it is; the other inputs are constants
      for (size_t k = 1; k < n.in.size(); ++k)
        if (!n.in[k].empty() && !cst.count(n.in[k]))
        {
          msg = "node " + n.label() + " with an input that is no constant";
          return false;
        }
      if (n.op == "Cast")
      {
        int64_t const to = n.geti("to", 1);
        if (to != 1 && to != 10 && to != 11)
        {
          msg = "node " + n.label() + " casts to a type that is not floating-point";
          return false;
        }
      }
    }

    else
    {
      msg = "node " + n.label() + " is not supported";
      return false;
    }

    cur = n.out[0];
  }

  if (cur != g.outputs[0])
  {
    msg = "the output of the graph, '" + g.outputs[0] + "', is not the end of the chain";
    return false;
  }
  return check(data, msg);
}

inline bool
NETIO::read_onnx_bytes(std::string const& bytes, t_Data& data, std::string& msg)
{
  t_Data res;
  try
  {
    PB pb{(uint8_t const*)bytes.data(), (uint8_t const*)bytes.data() + bytes.size()};
    Graph g;
    bool hasgraph = false;
    uint32_t num, wt;
    while (pb.field(num, wt))
    {
      if (num == 7 && wt == 2)
      {
        _read_graph(pb.sub(), g);
        hasgraph = true;
      }
      else
        pb.skip(wt);
    }
    if (!hasgraph)
    {
      msg = "no graph in the ONNX model";
      return false;
    }
    if (g.sparse)
    {
      msg = "sparse initializers are not supported";
      return false;
    }
    if (!_chain(g, res, msg)) return false;
  }
  catch (std::exception const& e)
  {
    msg = e.what();
    return false;
  }
  data.swap(res);
  return true;
}

inline bool
NETIO::read_onnx(std::string const& filename, t_Data& data, std::string& msg)
{
  std::ifstream is(filename, std::ios::binary);
  if (!is)
  {
    msg = "cannot open " + filename;
    return false;
  }
  std::string const bytes((std::istreambuf_iterator<char>(is)),
                          std::istreambuf_iterator<char>());
  if (!read_onnx_bytes(bytes, data, msg))
  {
    msg = filename + ": " + msg;
    return false;
  }
  return true;
}

#if defined(MC__USE_TORCH)
inline bool
NETIO::tensor_to_rows(std::vector<std::vector<double>>& v, at::Tensor const& t,
                      size_t const offset, bool const disp)
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
  at::Tensor const tc = t.detach().to(torch::kCPU, torch::kFloat64).contiguous();
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
    for (size_t i1 = 0; i1 < n1; ++i1) v[i0][offset + i1] = data_ptr[i0 * n1 + i1];
  }
  return true;
}

inline bool
NETIO::read_torch_net(torch::jit::script::Module const& module, t_Data& data,
                      std::string& msg, bool const disp)
{
  data.clear();
  if (!module.hasattr("net"))
  {
    msg = "Module has no 'net' submodule; cannot inspect layers and activation functions";
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
    std::string const type     = type_str.substr(type_str.rfind('.') + 1);
    if (disp)
      std::cout << "Submodule '" << child.name << "'  type: " << type_str << "\n";

    if (type == "Linear")
    {
      if (!child.value.hasattr("weight"))
      {
        msg = "Linear submodule without weight";
        return false;
      }
      at::Tensor const weight = child.value.attr("weight").toTensor();
      if (weight.dim() != 2 || !weight.size(0) || !weight.size(1) ||
          (!data.empty() && data.back().first.size() != (size_t)weight.size(1)))
      {
        msg = "Linear submodule with inconsistent size";
        return false;
      }

      // Each neuron holds its bias followed by its weights; a layer without
      // bias has zero biases
      std::vector<std::vector<double>> layer;
      if (!tensor_to_rows(layer, weight, 1, disp))
      {
        msg = "Linear submodule with an unsupported weight tensor";
        return false;
      }
      if (child.value.hasattr("bias"))
      {
        c10::IValue const bias = child.value.attr("bias");
        if (bias.isTensor() && !tensor_to_rows(layer, bias.toTensor(), 0, disp))
        {
          msg = "Linear submodule with an unsupported bias tensor";
          return false;
        }
      }
      data.push_back(std::make_pair(layer, (int)LINEAR));
      pending = true;
    }

    else if (type == "Tanh" || type == "ReLU" || type == "Sigmoid")
    {
      if (!pending)
      {
        msg = "Activation function without a Linear layer before it";
        return false;
      }
      data.back().second = (type == "Tanh" ? TANH : type == "ReLU" ? RELU : SIGMOID);
      pending = false;
    }

    // No effect in evaluation mode
    else if (type == "Identity" || type == "Dropout")
      continue;

    else
    {
      msg = "Non-supported submodule: " + type_str;
      return false;
    }
  }
  if (data.empty())
  {
    msg = "Module has no Linear submodule";
    return false;
  }
  return true;
}

inline std::vector<std::pair<at::Tensor, at::Tensor>>
NETIO::torch_linear(torch::jit::script::Module const& module)
{
  std::vector<std::pair<at::Tensor, at::Tensor>> lin;
  torch::jit::script::Module net = module.attr("net").toModule();
  for (auto const& child : net.named_children())
  {
    std::string const type_str = child.value.type()->str();
    if (type_str.substr(type_str.rfind('.') + 1) != "Linear") continue;
    at::Tensor bias;
    if (child.value.hasattr("bias") && child.value.attr("bias").isTensor())
      bias = child.value.attr("bias").toTensor();
    lin.push_back(std::make_pair(child.value.attr("weight").toTensor(), bias));
  }
  return lin;
}
#endif

////////////////////////////////////////////////////////////////////////
// Dictionaries saved by torch.save

inline bool
NETIO::_file_bytes(std::string const& filename, std::string& bytes)
{
  std::ifstream is(filename, std::ios::binary);
  if (!is) return false;
  bytes.assign(std::istreambuf_iterator<char>(is), std::istreambuf_iterator<char>());
  return true;
}

inline bool
NETIO::_zip_entries(std::string const& bytes,
                    std::map<std::string, std::pair<size_t, size_t>>& entries, std::string& msg)
{
  auto const* p   = reinterpret_cast<uint8_t const*>(bytes.data());
  size_t const sz = bytes.size();
  auto u16 = [&](size_t const o) { return (uint64_t)p[o] | ((uint64_t)p[o + 1] << 8); };
  auto u32 = [&](size_t const o) { return u16(o) | (u16(o + 2) << 16); };
  auto u64 = [&](size_t const o) { return u32(o) | (u32(o + 4) << 32); };

  // End of central directory record, from the end of the file
  if (sz < 22)
  {
    msg = "not a zip archive";
    return false;
  }
  size_t eocd = sz;
  for (size_t o = sz - 22;; --o)
  {
    if (u32(o) == 0x06054b50)
    {
      eocd = o;
      break;
    }
    if (o == 0 || sz - o > 22 + 65535) break;
  }
  if (eocd == sz)
  {
    msg = "not a zip archive";
    return false;
  }
  uint64_t nent = u16(eocd + 10), cdoff = u32(eocd + 16);
  // Zip64 record, where the counts do not fit
  if (nent == 0xFFFF || cdoff == 0xFFFFFFFF)
  {
    if (eocd < 20 || u32(eocd - 20) != 0x07064b50)
    {
      msg = "zip archive without zip64 locator";
      return false;
    }
    uint64_t const o64 = u64(eocd - 20 + 8);
    if (o64 + 56 > sz || u32(o64) != 0x06064b50)
    {
      msg = "corrupted zip64 record";
      return false;
    }
    nent  = u64(o64 + 32);
    cdoff = u64(o64 + 48);
  }

  entries.clear();
  uint64_t o = cdoff;
  for (uint64_t k = 0; k < nent; ++k)
  {
    if (o + 46 > sz || u32(o) != 0x02014b50)
    {
      msg = "corrupted zip directory";
      return false;
    }
    uint64_t const method = u16(o + 10), nlen = u16(o + 28), xlen = u16(o + 30),
                   clen = u16(o + 32);
    uint64_t csize = u32(o + 20), usize = u32(o + 24), loc = u32(o + 42);
    if (o + 46 + nlen + xlen + clen > sz)
    {
      msg = "corrupted zip directory";
      return false;
    }
    std::string const name(bytes, o + 46, nlen);
    // Zip64 extra field: the values that do not fit, in this order
    for (uint64_t x = o + 46 + nlen, xe = x + xlen; x + 4 <= xe;)
    {
      uint64_t const id = u16(x), len = u16(x + 2);
      if (id == 0x0001)
      {
        uint64_t q = x + 4;
        for (uint64_t* v : {&usize, &csize, &loc})
          if (*v == 0xFFFFFFFF && q + 8 <= x + 4 + len)
          {
            *v = u64(q);
            q += 8;
          }
      }
      x += 4 + len;
    }
    if (method != 0 || csize != usize)
    {
      entries[name] = {std::string::npos, (size_t)usize};
      o += 46 + nlen + xlen + clen;
      continue;
    }
    if (loc + 30 > sz || u32(loc) != 0x04034b50)
    {
      msg = "corrupted zip entry '" + name + "'";
      return false;
    }
    uint64_t const data = loc + 30 + u16(loc + 26) + u16(loc + 28);
    if (data + usize > sz)
    {
      msg = "truncated zip entry '" + name + "'";
      return false;
    }
    entries[name] = {(size_t)data, (size_t)usize};
    o += 46 + nlen + xlen + clen;
  }
  return true;
}

inline bool
NETIO::is_torch_dict(std::string const& filename)
{
  std::string bytes, msg;
  std::map<std::string, std::pair<size_t, size_t>> entries;
  if (!_file_bytes(filename, bytes) || !_zip_entries(bytes, entries, msg)) return false;
  auto ends = [](std::string const& s, std::string const& e) {
    return s == e || (s.size() > e.size() && s.compare(s.size() - e.size() - 1, e.size() + 1, "/" + e) == 0);
  };
  bool pkl = false;
  for (auto const& [name, pos] : entries)
  {
    if (ends(name, "constants.pkl")) return false;
    if (ends(name, "data.pkl")) pkl = true;
  }
  return pkl;
}

inline bool
NETIO::read_pickle(std::string const& filename, std::shared_ptr<PValue>& root, std::string& msg)
{
  typedef std::shared_ptr<PValue> t_P;
  std::string bytes;
  if (!_file_bytes(filename, bytes))
  {
    msg = "cannot open file " + filename;
    return false;
  }
  std::map<std::string, std::pair<size_t, size_t>> entries;
  if (!_zip_entries(bytes, entries, msg)) return false;

  // The pickle and the directory of the archive that holds it
  std::string prefix;
  std::pair<size_t, size_t> pkl{0, 0};
  bool found = false;
  for (auto const& [name, pos] : entries)
    if (name == "data.pkl" || (name.size() > 9 && name.compare(name.size() - 9, 9, "/data.pkl") == 0))
    {
      prefix = name.substr(0, name.size() - 8);
      pkl    = pos;
      found  = true;
      break;
    }
  if (!found)
  {
    msg = "no data.pkl in the archive: not a file saved by torch.save";
    return false;
  }
  if (pkl.first == std::string::npos)
  {
    msg = "compressed data.pkl in the archive";
    return false;
  }
  auto itbo = entries.find(prefix + "byteorder");
  if (itbo != entries.end() && (itbo->second.first == std::string::npos ||
                                bytes.compare(itbo->second.first, itbo->second.second, "little") != 0))
  {
    msg = "tensors not stored in little-endian byte order";
    return false;
  }

  try
  {
    auto const* p = reinterpret_cast<uint8_t const*>(bytes.data()) + pkl.first;
    auto const* e = p + pkl.second;
    auto need = [&](size_t const n) {
      if ((size_t)(e - p) < n) throw std::runtime_error("truncated pickle");
    };
    auto le = [&](size_t const n) {
      need(n);
      uint64_t v = 0;
      for (size_t k = 0; k < n; ++k) v |= (uint64_t)p[k] << (8 * k);
      p += n;
      return v;
    };
    auto str = [&](size_t const n) {
      need(n);
      std::string s(reinterpret_cast<char const*>(p), n);
      p += n;
      return s;
    };
    auto line = [&]() {
      std::string s;
      for (;; ++p)
      {
        need(1);
        if (*p == '\n') break;
        s.push_back((char)*p);
      }
      ++p;
      return s;
    };
    auto make = [](int const kind) {
      t_P v   = std::make_shared<PValue>();
      v->kind = kind;
      return v;
    };

    std::vector<t_P> stack;
    std::vector<size_t> marks;
    std::map<uint64_t, t_P> memo;
    auto pop = [&]() {
      if (stack.empty()) throw std::runtime_error("corrupted pickle (empty stack)");
      t_P v = stack.back();
      stack.pop_back();
      return v;
    };
    auto popmark = [&]() {
      if (marks.empty() || marks.back() > stack.size())
        throw std::runtime_error("corrupted pickle (no mark)");
      std::vector<t_P> v(stack.begin() + marks.back(), stack.end());
      stack.resize(marks.back());
      marks.pop_back();
      return v;
    };
    auto top = [&]() -> PValue& {
      if (stack.empty()) throw std::runtime_error("corrupted pickle (empty stack)");
      return *stack.back();
    };

    // Tensor from its storage, offset, sizes and strides
    auto tensor = [&](std::vector<t_P> const& args) {
      if (args.size() < 4 || !args[0] || args[0]->kind != PValue::OTHER || args[0]->list.size() != 2 ||
          args[1]->kind != PValue::INT || args[2]->kind != PValue::LIST || args[3]->kind != PValue::LIST ||
          args[2]->list.size() != args[3]->list.size())
        throw std::runtime_error("unsupported tensor");
      std::string const& type = args[0]->s;
      size_t const elsize =
          (type == "torch.DoubleStorage" || type == "torch.LongStorage" ? 8
           : type == "torch.FloatStorage" || type == "torch.IntStorage" ? 4
           : type == "torch.HalfStorage" || type == "torch.ShortStorage" ? 2
           : type == "torch.CharStorage" || type == "torch.ByteStorage" || type == "torch.BoolStorage" ? 1
                                                                                                         : 0);
      if (!elsize) throw std::runtime_error("unsupported tensor storage " + type);
      auto its = entries.find(prefix + "data/" + args[0]->list[0]->s);
      if (its == entries.end() || its->second.first == std::string::npos)
        throw std::runtime_error("tensor storage not found in the archive, or compressed");
      auto const* d      = reinterpret_cast<uint8_t const*>(bytes.data()) + its->second.first;
      size_t const numel = its->second.second / elsize;

      t_P t          = make(PValue::TENSOR);
      size_t const nd = args[2]->list.size();
      std::vector<size_t> stride(nd);
      size_t total = 1;
      for (size_t k = 0; k < nd; ++k)
      {
        if (args[2]->list[k]->kind != PValue::INT || args[3]->list[k]->kind != PValue::INT ||
            args[2]->list[k]->n < 0 || args[3]->list[k]->n < 0)
          throw std::runtime_error("unsupported tensor");
        t->dims.push_back((size_t)args[2]->list[k]->n);
        stride[k] = (size_t)args[3]->list[k]->n;
        total *= t->dims.back();
        if (total > numel) throw std::runtime_error("tensor larger than its storage");
      }
      t->val.resize(total);
      std::vector<size_t> idx(nd, 0);
      for (size_t c = 0; c < total; ++c)
      {
        size_t o = (size_t)args[1]->n;
        for (size_t k = 0; k < nd; ++k) o += idx[k] * stride[k];
        if (o >= numel) throw std::runtime_error("tensor larger than its storage");
        uint8_t const* q = d + o * elsize;
        double v         = 0.;
        if (type == "torch.DoubleStorage")
          std::memcpy(&v, q, 8);
        else if (type == "torch.FloatStorage")
        {
          float f;
          std::memcpy(&f, q, 4);
          v = f;
        }
        else if (type == "torch.HalfStorage")
        {
          uint16_t h;
          std::memcpy(&h, q, 2);
          int const ex   = (h >> 10) & 0x1f, fr = h & 0x3ff;
          double const m = (ex == 0 ? std::ldexp((double)fr, -24)
                            : ex == 31 ? (fr ? std::nan("") : INFINITY)
                                       : std::ldexp((double)(fr | 0x400), ex - 25));
          v = (h & 0x8000 ? -m : m);
        }
        else if (type == "torch.LongStorage")
        {
          int64_t i;
          std::memcpy(&i, q, 8);
          v = (double)i;
        }
        else if (type == "torch.IntStorage")
        {
          int32_t i;
          std::memcpy(&i, q, 4);
          v = i;
        }
        else if (type == "torch.ShortStorage")
        {
          int16_t i;
          std::memcpy(&i, q, 2);
          v = i;
        }
        else if (type == "torch.CharStorage")
          v = (int8_t)*q;
        else
          v = *q;
        t->val[c] = v;
        for (size_t k = nd; k-- > 0;)
        {
          if (++idx[k] < t->dims[k]) break;
          idx[k] = 0;
        }
      }
      return t;
    };

    for (bool stop = false; !stop;)
    {
      need(1);
      uint8_t const op = *p++;
      switch (op)
      {
        case 0x80: le(1); break;                      // PROTO
        case 0x95: le(8); break;                      // FRAME
        case '.': stop = true; break;                   // STOP
        case 'N': stack.push_back(make(PValue::NONE)); break;
        case 0x88:
        case 0x89:
          stack.push_back(make(PValue::BOOL));
          top().b = (op == 0x88);
          break;
        case 'J':                                       // BININT (signed)
          stack.push_back(make(PValue::INT));
          top().n = (int32_t)le(4);
          break;
        case 'K':
        case 'M':                                       // BININT1, BININT2
          stack.push_back(make(PValue::INT));
          top().n = (long long)le(op == 'K' ? 1 : 2);
          break;
        case 0x8a:                                      // LONG1
        {
          size_t const n = le(1);
          if (n > 8) throw std::runtime_error("integer too large");
          uint64_t v = le(n);
          if (n && n < 8 && (v >> (8 * n - 1)) & 1) v |= ~uint64_t(0) << (8 * n);
          stack.push_back(make(PValue::INT));
          top().n = (long long)v;
          break;
        }
        case 'G':                                       // BINFLOAT, big-endian
        {
          need(8);
          uint64_t v = 0;
          for (size_t k = 0; k < 8; ++k) v = (v << 8) | p[k];
          p += 8;
          stack.push_back(make(PValue::REAL));
          std::memcpy(&top().x, &v, 8);
          break;
        }
        case 'X':
        case 'T':
        case 'B':                                       // BINUNICODE, BINSTRING, BINBYTES
          stack.push_back(make(PValue::STR));
          top().s = str(le(4));
          break;
        case 0x8c:
        case 'U':
        case 'C':                                       // SHORT_BINUNICODE, SHORT_BINSTRING, SHORT_BINBYTES
          stack.push_back(make(PValue::STR));
          top().s = str(le(1));
          break;
        case 0x8d:                                      // BINUNICODE8
          stack.push_back(make(PValue::STR));
          top().s = str(le(8));
          break;
        case ']':
        case ')': stack.push_back(make(PValue::LIST)); break;
        case '}': stack.push_back(make(PValue::DICT)); break;
        case '(': marks.push_back(stack.size()); break;
        case 't':
        case 'l':                                       // TUPLE, LIST
        {
          std::vector<t_P> items = popmark();
          stack.push_back(make(PValue::LIST));
          top().list = items;
          break;
        }
        case 0x85:
        case 0x86:
        case 0x87:                                      // TUPLE1, TUPLE2, TUPLE3
        {
          size_t const n = op - 0x84;
          if (stack.size() < n) throw std::runtime_error("corrupted pickle (empty stack)");
          std::vector<t_P> items(stack.end() - n, stack.end());
          stack.resize(stack.size() - n);
          stack.push_back(make(PValue::LIST));
          top().list = items;
          break;
        }
        case 'd':                                       // DICT
        {
          std::vector<t_P> items = popmark();
          stack.push_back(make(PValue::DICT));
          for (size_t k = 0; k + 1 < items.size(); k += 2) top().dict.push_back({items[k], items[k + 1]});
          break;
        }
        case 'a':                                       // APPEND
        {
          t_P v = pop();
          if (top().kind != PValue::LIST) throw std::runtime_error("corrupted pickle (append)");
          top().list.push_back(v);
          break;
        }
        case 'e':                                       // APPENDS
        {
          std::vector<t_P> items = popmark();
          if (top().kind != PValue::LIST) throw std::runtime_error("corrupted pickle (append)");
          top().list.insert(top().list.end(), items.begin(), items.end());
          break;
        }
        case 's':                                       // SETITEM
        {
          t_P v = pop(), k = pop();
          if (top().kind == PValue::DICT) top().dict.push_back({k, v});
          break;
        }
        case 'u':                                       // SETITEMS
        {
          std::vector<t_P> items = popmark();
          if (top().kind == PValue::DICT)
            for (size_t k = 0; k + 1 < items.size(); k += 2) top().dict.push_back({items[k], items[k + 1]});
          break;
        }
        case 'q':
        case 'r':                                       // BINPUT, LONG_BINPUT
        {
          uint64_t const k = le(op == 'q' ? 1 : 4);
          if (stack.empty()) throw std::runtime_error("corrupted pickle (empty stack)");
          memo[k] = stack.back();
          break;
        }
        case 0x94:                                      // MEMOIZE
          if (stack.empty()) throw std::runtime_error("corrupted pickle (empty stack)");
          memo[memo.size()] = stack.back();
          break;
        case 'h':
        case 'j':                                       // BINGET, LONG_BINGET
        {
          auto it = memo.find(le(op == 'h' ? 1 : 4));
          if (it == memo.end()) throw std::runtime_error("corrupted pickle (memo)");
          stack.push_back(it->second);
          break;
        }
        case 'c':                                       // GLOBAL
        {
          std::string const mod = line(), name = line();
          stack.push_back(make(PValue::OTHER));
          top().s = mod + "." + name;
          break;
        }
        case 0x93:                                      // STACK_GLOBAL
        {
          t_P name = pop(), mod = pop();
          stack.push_back(make(PValue::OTHER));
          top().s = mod->s + "." + name->s;
          break;
        }
        case 'Q':                                       // BINPERSID: storage of a tensor
        {
          t_P pid = pop();
          if (pid->kind != PValue::LIST || pid->list.size() < 5 || pid->list[0]->s != "storage" ||
              pid->list[1]->kind != PValue::OTHER || pid->list[2]->kind != PValue::STR)
            throw std::runtime_error("unsupported persistent object");
          stack.push_back(make(PValue::OTHER));
          top().s    = pid->list[1]->s;
          top().list = {pid->list[2], pid->list[4]};
          break;
        }
        case 'R':                                       // REDUCE
        {
          t_P args = pop(), fct = pop();
          if (fct->kind != PValue::OTHER || args->kind != PValue::LIST)
            throw std::runtime_error("corrupted pickle (reduce)");
          if (fct->s == "collections.OrderedDict")
            stack.push_back(make(PValue::DICT));
          else if (fct->s == "torch._utils._rebuild_tensor_v2")
            stack.push_back(tensor(args->list));
          else if ((fct->s == "torch._utils._rebuild_parameter" ||
                    fct->s == "torch._utils._rebuild_parameter_with_state") &&
                   !args->list.empty())
            stack.push_back(args->list[0]);
          else
          {
            stack.push_back(make(PValue::OTHER));
            top().s = fct->s;
          }
          break;
        }
        case 0x81:                                      // NEWOBJ
        {
          t_P args = pop(), cls = pop();
          stack.push_back(make(PValue::OTHER));
          top().s = cls->s;
          break;
        }
        case 'b': pop(); break;                         // BUILD: state ignored
        case '0': pop(); break;                         // POP
        default:
          throw std::runtime_error("unsupported pickle opcode " + std::to_string((int)op));
      }
    }
    if (stack.size() != 1) throw std::runtime_error("corrupted pickle (stack at the end)");
    root = stack.back();
  }
  catch (std::exception const& ex)
  {
    msg = ex.what();
    return false;
  }
  return true;
}

inline bool
NETIO::read_dict(std::string const& filename, Dict& dict, std::string& msg, int const activ)
{
  Dict D;
  if (!read_pickle(filename, D.root, msg)) return false;
  PValue const& root = *D.root;
  auto fail = [&](std::string const& what) {
    msg = what;
    return false;
  };

  PValue const *type = root.at("model_type"), *config = root.at("config"),
               *state = root.at("state_dict");
  if (!type || type->kind != PValue::STR || !config || !state || state->kind != PValue::DICT)
    return fail("not a dictionary with entries model_type, config and state_dict");
  D.type = type->s;
  if (D.type != "KANmodel" && D.type != "ANNmodel") return fail("unknown model type '" + D.type + "'");
  std::vector<double> widths, intervals;
  if (!config->at("widths") || !config->at("widths")->numbers(widths) || widths.size() < 2)
    return fail("no list widths of at least two numbers in config");
  for (double const w : widths)
  {
    if (!(w >= 1.)) return fail("width below 1 in config");
    D.widths.push_back((size_t)w);
  }
  size_t const nlay = D.widths.size() - 1;

  // Tensor of the state dictionary with the given dimensions
  auto tensor = [&](std::string const& name, std::vector<size_t> const& dims,
                    std::vector<double>& val) {
    PValue const* t = state->at(name);
    if (!t || t->kind != PValue::TENSOR) return fail("no tensor " + name + " in state_dict");
    if (t->dims != dims) return fail("tensor " + name + " of unexpected dimensions");
    for (double const v : t->val)
      if (!std::isfinite(v)) return fail("tensor " + name + " with entries that are not finite");
    val = t->val;
    return true;
  };

  if (D.type == "KANmodel")
  {
    if (!config->at("intervals") || !config->at("intervals")->numbers(intervals) ||
        intervals.size() != nlay)
      return fail("no list intervals with one number per layer in config");
    D.kan.resize(nlay);
    for (size_t l = 0; l < nlay; ++l)
    {
      if (!(intervals[l] >= 1.)) return fail("number of intervals below 1 in config");
      Dict::KANLayer& K = D.kan[l];
      K.nin             = D.widths[l];
      K.nout            = D.widths[l + 1];
      K.R               = (size_t)intervals[l] - 1;
      std::string const k = "layers." + std::to_string(l) + ".";
      if (!tensor(k + "edge_bias", {K.nout, K.nin}, K.c) ||
          !tensor(k + "linear_weight", {K.nout, K.nin, 1}, K.v) ||
          !tensor(k + "relu_output_weight", {K.nout, K.nin, K.R}, K.a) ||
          !tensor(k + "relu_weight", {K.nout, K.nin, K.R}, K.w) ||
          !tensor(k + "relu_bias", {K.nout, K.nin, K.R}, K.b))
        return false;
    }
  }
  else
  {
    // Activation of the hidden layers: named in config, given, or matching the bounds of the
    // first hidden layer
    int act = activ;
    if (PValue const* a = config->at("activation"); a && a->kind == PValue::STR)
    {
      std::string s = a->s;
      for (auto& c : s) c = (char)std::tolower((unsigned char)c);
      act = (s == "relu" ? (int)RELU : s == "tanh" ? (int)TANH : s == "sigmoid" ? (int)SIGMOID : -2);
      if (act == -2) return fail("activation '" + a->s + "' is not supported");
    }
    else if (act != RELU && act != TANH && act != SIGMOID && nlay > 1)
    {
      act = -1;
      PValue const* bounds = root.at("bounds");
      for (size_t s = 0; act < 0 && bounds && bounds->kind == PValue::DICT && s < bounds->dict.size(); ++s)
      {
        PValue const* b = bounds->dict[s].second ? bounds->dict[s].second->at("bounds") : nullptr;
        PValue const* L = b ? b->at("layers") : nullptr;
        if (!L || L->kind != PValue::LIST || L->list.empty() || !L->list[0]) continue;
        PValue const *pre = L->list[0]->at("pre"), *post = L->list[0]->at("post");
        if (!pre || !post || pre->kind != PValue::TENSOR || post->kind != PValue::TENSOR ||
            pre->val.size() != post->val.size() || pre->val.empty())
          continue;
        for (int cand : {(int)RELU, (int)TANH, (int)SIGMOID})
        {
          bool match = true;
          for (size_t k = 0; match && k < pre->val.size(); ++k)
          {
            double const z = pre->val[k],
                         g = (cand == RELU ? std::max(z, 0.) : cand == TANH ? std::tanh(z)
                                                                            : 1. / (1. + std::exp(-z)));
            match = std::fabs(g - post->val[k]) <= 1e-6 * (1. + std::fabs(g));
          }
          if (match)
          {
            act = cand;
            break;
          }
        }
      }
      if (act < 0)
        return fail("activation of the hidden layers neither named in config nor recognised from the "
                    "bounds in the file");
    }
    for (size_t l = 0; l < nlay; ++l)
    {
      size_t const nin = D.widths[l], nout = D.widths[l + 1];
      std::string const k = "layers." + std::to_string(l) + ".";
      std::vector<double> W, b;
      if (!tensor(k + "weight", {nout, nin}, W) || !tensor(k + "bias", {nout}, b)) return false;
      std::vector<std::vector<double>> layer(nout);
      for (size_t j = 0; j < nout; ++j)
      {
        layer[j].push_back(b[j]);
        layer[j].insert(layer[j].end(), W.begin() + j * nin, W.begin() + (j + 1) * nin);
      }
      D.ann.push_back({layer, l + 1 < nlay ? act : (int)LINEAR});
    }
  }

  // Scalings
  if (PValue const* norm = root.at("normalization"); norm && norm->kind == PValue::DICT)
  {
    auto scaling = [&](std::string const& side, size_t const n, std::vector<double>& lower,
                       std::vector<double>& upper) {
      PValue const* sc = norm->at(side);
      if (!sc || sc->kind == PValue::NONE) return true;
      if (!sc->at("lower") || !sc->at("upper") || !sc->at("lower")->numbers(lower) ||
          !sc->at("upper")->numbers(upper) || lower.size() != n || upper.size() != n)
        return fail("normalization of " + side + " without lists lower and upper of " +
                    std::to_string(n) + " numbers");
      for (size_t i = 0; i < n; ++i)
        if (!(upper[i] > lower[i]) || !std::isfinite(lower[i]) || !std::isfinite(upper[i]))
          return fail("normalization of " + side + " with upper not above lower");
      return true;
    };
    if (!scaling("x", D.widths.front(), D.xlower, D.xupper) ||
        !scaling("z", D.widths.back(), D.zlower, D.zupper))
      return false;
  }

  dict = D;
  return true;
}

inline bool
NETIO::fold_scaling(t_Data& data, Dict const& dict, std::string& msg)
{
  if (data.empty()) return true;
  // Inputs: u = (x - lower) / (upper - lower) in the neurons of the first layer
  if (!dict.xlower.empty())
    for (auto& neuron : data.front().first)
    {
      if (neuron.size() != dict.xlower.size() + 1)
      {
        msg = "scaling of the inputs of another size than the first layer";
        return false;
      }
      for (size_t i = 0; i < dict.xlower.size(); ++i)
      {
        neuron[1 + i] /= dict.xupper[i] - dict.xlower[i];
        neuron[0] -= neuron[1 + i] * dict.xlower[i];
      }
    }
  // Outputs: z = lower + (upper - lower) y in the neurons of the last layer
  if (!dict.zlower.empty())
  {
    if (data.back().second != LINEAR)
    {
      msg = "scaling of the outputs cannot be folded into a layer that is not linear";
      return false;
    }
    if (data.back().first.size() != dict.zlower.size())
    {
      msg = "scaling of the outputs of another size than the last layer";
      return false;
    }
    for (size_t j = 0; j < dict.zlower.size(); ++j)
    {
      for (auto& w : data.back().first[j]) w *= dict.zupper[j] - dict.zlower[j];
      data.back().first[j][0] += dict.zlower[j];
    }
  }
  return true;
}

}  // end namespace mc

#endif
