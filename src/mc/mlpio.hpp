// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#ifndef MC__MLPIO_HPP
#define MC__MLPIO_HPP

#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
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

//! @brief C++ class reading the layers of a multilayer perceptron from a file
////////////////////////////////////////////////////////////////////////
//! mc::MLPIO reads the dense layers of a multilayer perceptron, and their
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
//! A network is returned as the data of mc::MLP::set_data: one entry per
//! layer, holding one row per neuron -- its bias followed by its weights --
//! and its activation function (LINEAR, RELU, TANH, SIGMOID, with the
//! values of mc::MLP::ACTIV_TYPE and mc::MLPREG::ACTIV_TYPE).
////////////////////////////////////////////////////////////////////////
class MLPIO
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
  //! holds anything else than a network of dense layers (see mc::MLPIO)
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

 private:
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
MLPIO::check(t_Data const& data, std::string& msg)
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
MLPIO::_read_tensor(PB pb, Tensor& t)
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
MLPIO::_read_attr(PB pb, Attr& a)
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
MLPIO::_read_node(PB pb, Node& n)
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
MLPIO::_read_valueinfo(PB pb, std::string& name, std::vector<int64_t>& shape)
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
MLPIO::_read_graph(PB pb, Graph& g)
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
MLPIO::_chain(Graph const& g, t_Data& data, std::string& msg)
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
MLPIO::read_onnx_bytes(std::string const& bytes, t_Data& data, std::string& msg)
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
MLPIO::read_onnx(std::string const& filename, t_Data& data, std::string& msg)
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
MLPIO::tensor_to_rows(std::vector<std::vector<double>>& v, at::Tensor const& t,
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
MLPIO::read_torch_net(torch::jit::script::Module const& module, t_Data& data,
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
MLPIO::torch_linear(torch::jit::script::Module const& module)
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

}  // end namespace mc

#endif
