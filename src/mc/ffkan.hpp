// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

/*!
\page page_KAN Kolmogorov-Arnold networks
\author Benoit Chachuat

The class mc::KAN holds a Kolmogorov-Arnold network (KAN): a sequence of layers, each mapping its
inputs \f$u_1,\ldots,u_n\f$ to the outputs
\f[ y_j = \sum_{i=1}^{n} f_{ji}(u_i), \quad j=1,\ldots,m, \f]
where every edge function \f$f_{ji}\f$ is a function of one variable with finitely many parameters.
It mirrors mc::MLP: the network is set layer by layer (mc::KAN::append_data), all at once
(mc::KAN::set_data) or read from a file (mc::KAN::read_data), and is then evaluated in the
arithmetics of MC++, differentiated, relaxed and propagated backward.

The edge functions of a layer are of one type (mc::KAN::EDGE_TYPE). The type available is PWL,
continuous piecewise-linear edges held as mc::PWLU objects extended beyond their domain (see
\ref page_PWLU), which mc::PWLU::relu builds from the form \f$c + v u + \sum_r a_r \max(w_r u + b_r,
0)\f$. Other types -- splines, polynomials -- have their place in mc::KAN::Layer.

\code
  #include "ffkan.hpp"
  typedef mc::Interval I;

  mc::KAN<I> NN;
  // layer of 2 inputs and 1 output: layer[j][i] is the edge from input i to output j
  NN.append_data( { { mc::PWLU( { 0., 1. }, { 0., 2. }, -1., 0.5 ),
                      mc::PWLU::relu( 0.1, -0.5, { 1.2 }, { 1. }, { 0.3 } ) } } );
  double x[2] = { 0.3, -0.2 }, y;
  NN.eval( x, &y );
  I IX[2] = { I(0.,1.), I(-1.,1.) }, IY;
  NN.eval( IX, &IY );                      // exact range of each edge, summed
\endcode

mc::KAN::read_data reads a network saved as a dictionary by torch.save (mc::NETIO::read_dict),
without Torch.

The network is evaluated edge by edge with the evaluations of mc::FFPWLU (mc::FFPWLU::value): exact
ranges in interval arithmetic, convex and concave envelopes in McCormick arithmetic, exact summands
of a variable in superposition models. Its expression tree (mc::KAN::dag), which the polyhedral
relaxations and the backward propagation go through, holds one mc::FFPWLU operation per input of
a layer, for the edges leaving it, and one mc::FFLin operation per output; the encoding of the
edges in a polyhedral relaxation is set by mc::KAN::Options::ENCODING and SEPARATE.

mc::FFKAN defines the network as an external operation of a DAG -- it is mc::FFMAP on an mc::KAN,
as mc::FFMLP is mc::FFMAP on an mc::MLP:

\code
  mc::FFGraph DAG;
  std::vector<mc::FFVar> X(2);
  for (auto& Xi : X) Xi.set(&DAG);
  mc::FFKAN<I> Op;
  std::vector<mc::FFVar> Y = Op(X, &NN);   // the DAG holds a copy of NN: Op.map()
  Op.options.RELAX = {Op.options.AUX};     // polyhedral relaxations through the expression tree
                                           // of the network, with the encoding of the edges
\endcode

The operation is evaluated by the network in the arithmetics of MC++; its polyhedral relaxation
is by McCormick relaxations by default, and through the expression tree of the network -- the
encodings of the edges -- with mc::FFEXTERN::Options::RELAX set to AUX. The differentiation of
the DAG inserts the Jacobian of the network (mc::FFGradKAN), or, with mc::KAN::Options::TANGENT
and ADJOINT, its directional derivative (mc::FFTanKAN) and its adjoint (mc::FFAdjKAN).
*/

#ifndef MC__FFKAN_HPP
#define MC__FFKAN_HPP

#include <algorithm>
#include <type_traits>
#include <cfloat>
#include <string>
#include <vector>

#include "ffhome.hpp"
#include "fflin.hpp"
#include "ffmap.hpp"
#include "ffpwlu.hpp"
#include "interval.hpp"
#include "netio.hpp"
#include "specbnd.hpp"

namespace mc
{

//! @brief C++ class for the evaluation, differentiation and relaxation of
//! Kolmogorov-Arnold networks
////////////////////////////////////////////////////////////////////////
//! mc::KAN is a C++ class for Kolmogorov-Arnold networks (KAN), mirroring
//! mc::MLP: evaluation in the arithmetics of MC++, Jacobian, directional
//! derivatives and adjoints, expression tree, backward propagation. See
//! \ref page_KAN.
////////////////////////////////////////////////////////////////////////
template <typename T = mc::Interval>
class KAN
////////////////////////////////////////////////////////////////////////
{
 public:
  //! @brief Enumeration type for the edge functions of a layer
  enum EDGE_TYPE
  {
    PWL = 0  //!< Continuous piecewise-linear functions (mc::PWLU)
  };

  //! @brief Layer of a network: its outputs are the sums over its inputs of
  //! the edge functions, of one type per layer. A further type of edge adds
  //! its value to EDGE_TYPE, the container of its descriptors here, and its
  //! case to the evaluations of an edge (_edge, _dedge) and to its insertion
  //! in the expression tree (_build)
  struct Layer
  {
    //! @brief Type of the edge functions
    int type = PWL;
    //! @brief Piecewise-linear edges: pwl[j][i] from input i to output j
    std::vector<std::vector<PWLU>> pwl;

    //! @brief Number of inputs
    size_t
    nin() const
    {
      return pwl.empty() ? 0 : pwl.front().size();
    }
    //! @brief Number of outputs
    size_t
    nout() const
    {
      return pwl.size();
    }
  };

  //! @brief Reset KAN data
  void reset_data();

  //! @brief Set KAN data: false, with the network unchanged, where a layer is
  //! empty, has outputs with different numbers of edges, an edge without
  //! segment, or fewer or more inputs than the previous layer has outputs
  bool set_data(std::vector<Layer> const& kan);

  //! @brief Set KAN data from a dictionary saved by torch.save, as read by
  //! mc::NETIO::read_dict (model of type "KANmodel"): the file is then read
  //! once, where its type, sizes or scalings are also wanted. The min-max
  //! scalings of the dictionary are folded into the first and last layers, as
  //! by read_data. False, with the network unchanged and the reason in
  //! <a>msg</a> if given, where the dictionary holds anything else
  //! (a template only so that set_data({}) remains a call with layers)
  template <typename DICT, typename = typename std::enable_if<std::is_same<DICT, NETIO::Dict>::value>::type>
  bool
  set_data(DICT const& dict, std::string* msg = nullptr)
  {
    return _set_dict(dict, msg);
  }

  //! @brief Append a layer of piecewise-linear edges to KAN data: layer[j][i]
  //! is the edge from input i to output j, extended beyond its domain (see
  //! mc::PWLU); false, with the network unchanged, as for set_data
  bool append_data(std::vector<std::vector<PWLU>> const& layer, bool const reset = false);

  //! @brief Read KAN data from a file: a dictionary saved by torch.save, with
  //! entries model_type ("KANmodel"), config (widths, intervals), state_dict
  //! (for each layer l, the tensors layers.l.edge_bias c, linear_weight v,
  //! relu_output_weight a, relu_weight w and relu_bias b of the edges \f$c + v
  //! u + \sum_r a_r \max(w_r u + b_r, 0)\f$) and normalization (see
  //! mc::NETIO::read_dict), read without Torch. The min-max scalings of the
  //! file are folded into the first and last layers: the network maps the raw
  //! inputs to the raw outputs. False, with the network unchanged, where the
  //! file cannot be read or holds anything else; the reason is displayed with
  //! <a>disp</a>
  bool read_data(std::string const& filename, bool const disp = false);

  //! @brief Layer of piecewise-linear edges from the parameters of the form
  //! \f$c + v u + \sum_{r<R} a_r \max(w_r u + b_r, 0)\f$: c and v indexed
  //! [j*nin+i], a, w and b indexed [(j*nin+i)*R+r]
  static std::vector<std::vector<PWLU>> relu_layer(size_t const nout, size_t const nin,
                                                   size_t const R, double const* c,
                                                   double const* v, double const* a,
                                                   double const* w, double const* b);

  //! @brief KAN options
  struct Options
  {
    //! @brief Default constructor
    Options() { reset(); }

    //! @brief Reset options
    void
    reset()
    {
      ENCODING = FFPWLU<T>::Options::DEFAULT;
      SEPARATE = false;
      AUTODIFF = F;
      TANGENT  = false;
      ADJOINT  = false;
      CPMAX    = 1;  // Leave it to the outer DAG to iterate as necessary by default
      CPTHRES  = 1e4 * DBL_EPSILON;
      CPINF    = 1e30;
    }

    //! @brief Enumeration type for AD strategy
    enum AD_TYPE
    {
      F = 0,  //!< Forward differentiation
      B       //!< Backward differentiation
    };

    //! @brief Encoding of the piecewise-linear edges in a polyhedral
    //! relaxation of the expression tree, see mc::FFPWLU::Options::ENCODING
    //! - Default: DEFAULT. It may be changed at any time: the expression tree
    //! of the network is built again
    int ENCODING;
    //! @brief Whether each piecewise-linear edge is encoded on its own
    //! subdivision, rather than the edges leaving an input on one, see
    //! mc::FFPWLU::Options::SEPARATE - Default: false. It may be changed at
    //! any time, as ENCODING
    bool SEPARATE;
    //! @brief Whether the Jacobian (grad) is accumulated forward or backward
    //! through the layers. With the default (F) it is still accumulated
    //! backward where the outputs are no more than the inputs
    int AUTODIFF;
    //! @brief Whether the forward symbolic differentiation of an operation on
    //! the network in a DAG inserts its directional derivative, as
    //! mc::MLP::Options::TANGENT
    bool TANGENT;
    //! @brief Whether the backward symbolic differentiation of an operation
    //! on the network in a DAG inserts its adjoint, as
    //! mc::MLP::Options::ADJOINT
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
  //! @brief Options ENCODING and SEPARATE with which the expression tree was
  //! built: the tree is built again when they change
  std::pair<int, bool> _dagopt{0, false};
  //! @brief Number of variables
  size_t _nin;
  //! @brief Number of dependents
  size_t _nout;
  //! @brief Variables
  std::vector<FFVar> _varin;
  //! @brief Dependents
  std::vector<FFVar> _varout;
  //! @brief Codelist
  FFSubgraph _codelist;
  //! @brief KAN data
  std::vector<Layer> _data;
  //! @brief Stamp of the network: see uid()
  size_t _uid;
  //! @brief Intermediate storage for the expression tree
  std::vector<PolVar<T>> _wkPOL;
  std::vector<T> _wkCPI;

  //! @brief Value of the edge from input i to output j of a layer at x, in
  //! the arithmetic of x
  template <typename U>
  static U
  _edge(Layer const& L, size_t const j, size_t const i, U const& x)
  {
    switch (L.type)
    {
      case PWL:
      default:
        return FFPWLU<T>::value(L.pwl[j][i], x);
    }
  }

  //! @brief Derivative of the edge from input i to output j of a layer at x,
  //! in the arithmetic of x
  template <typename U>
  static U
  _dedge(Layer const& L, size_t const j, size_t const i, U const& x)
  {
    switch (L.type)
    {
      case PWL:
      default:
        return FFPWLU<T>::slope(L.pwl[j][i], x);
    }
  }

  //! @brief Spectral bounds are not available for piecewise-linear edges: the
  //! values and derivatives of the network in that arithmetic -- which an
  //! external operation on the network (mc::FFMAP) may be asked for -- throw
  static Specbnd<T>
  _edge(Layer const& L, size_t const j, size_t const i, Specbnd<T> const& x)
  {
    throw std::runtime_error("mc::KAN ** No spectral bounds for piecewise-linear edges\n");
  }
  static Specbnd<T>
  _dedge(Layer const& L, size_t const j, size_t const i, Specbnd<T> const& x)
  {
    throw std::runtime_error("mc::KAN ** No spectral bounds for piecewise-linear edges\n");
  }

  //! @brief Check a layer
  static bool _check(Layer const& L);

  //! @brief Insert the network in the DAG of its inputs
  void _build(FFVar const* varin, FFVar* varout) const;

  //! @brief Set KAN data from a dictionary read by mc::NETIO::read_dict
  bool _set_dict(NETIO::Dict const& dict, std::string* msg);

  //! @brief Set expression tree
  void _set_dag();

  //! @brief Evaluate KAN, layer by layer
  template <typename U>
  void _eval(U const* valin, U* valout) const;

  //! @brief Evaluate KAN Jacobian, accumulated forward or backward
  template <typename U>
  void _jac(U const* valin, U* gradout, U* valout, bool const backward) const;

  //! @brief Evaluate KAN directional derivative: forward pass with the values
  //! and the tangents of the nodes
  template <typename U>
  void _tang(U const* valin, U const* dirin, U* tanout, U* valout) const;

  //! @brief Evaluate KAN adjoint: forward pass keeping the derivative of each
  //! edge, then backward pass through the layers
  template <typename U>
  void _adj(U const* valin, U const* wout, U* adjin, U* valout) const;

 public:
  //! @brief Default constructor
  KAN() : _dag(nullptr), _nin(0), _nout(0), _uid(ffhome_uid()) {}

  //! @brief Copy constructor: the copy keeps the stamp of the network and
  //! builds an expression tree of its own
  KAN(KAN const& other)
      : options(other.options),
        _dag(nullptr),
        _nin(other._nin),
        _nout(other._nout),
        _data(other._data),
        _uid(other._uid)
  {
  }

  virtual ~KAN() { delete _dag; }

  //! @brief Assignment operator
  KAN&
  operator=(KAN const& other)
  {
    if (this == &other) return *this;
    options = other.options;
    delete _dag;
    _dag  = nullptr;
    _nin  = other._nin;
    _nout = other._nout;
    _varin.clear();
    _varout.clear();
    _codelist.clear();
    _data = other._data;
    _uid  = other._uid;
    _wkPOL.clear();
    _wkCPI.clear();
    return *this;
  }

  //! @brief Evaluate the outputs at valin, in the arithmetic of valin: layer
  //! by layer, each edge with mc::FFPWLU::value
  template <typename U>
  void
  eval(U const* valin, U* valout) const
  {
    _eval(valin, valout);
  }

  //! @brief Polyhedral image of the outputs, through the expression tree
  void
  eval(PolVar<T> const* valin, PolVar<T>* valout)
  {
    _set_dag();
    _dag->eval(_codelist, _wkPOL, _varout.size(), _varout.data(), valout, _varin.size(),
               _varin.data(), valin);
  }

  //! @brief Insert the network in the DAG of the inputs valin: one mc::FFPWLU
  //! operation per input of a layer and one mc::FFLin operation per output
  void
  eval(FFVar const* valin, FFVar* valout) const
  {
    _build(valin, valout);
  }

  //! @brief Spectral bounds are not available for piecewise-linear edges
  void
  eval(Specbnd<T> const* valin, Specbnd<T>* valout) const
  {
    throw std::runtime_error("mc::KAN::eval ** No spectral bounds for piecewise-linear edges\n");
  }

  //! @brief Evaluate the Jacobian at valin, gradout[k*nin()+i] being the
  //! derivative of output k in input i, and the outputs in valout if given;
  //! right derivatives at the breakpoints in double precision
  template <typename U>
  void
  grad(U const* valin, U* gradout, U* const valout = nullptr) const
  {
    _jac(valin, gradout, valout, _nout <= _nin || options.AUTODIFF == Options::AD_TYPE::B);
  }

  //! @brief Evaluate the directional derivative of the outputs at valin along
  //! dirin, and the outputs in valout if given
  template <typename U>
  void
  tang(U const* valin, U const* dirin, U* tanout, U* const valout = nullptr) const
  {
    _tang(valin, dirin, tanout, valout);
  }

  //! @brief Evaluate the product of the transposed Jacobian at valin with
  //! wout, and the outputs in valout if given
  template <typename U>
  void
  adj(U const* valin, U const* wout, U* adjin, U* const valout = nullptr) const
  {
    _adj(valin, wout, adjin, valout);
  }

  //! @brief Contract the inputs valin and outputs valout by forward and
  //! backward propagation through the expression tree: false if the network
  //! has no point with inputs in valin and outputs in valout
  bool
  reval(T* valin, T* valout)
  {
    _set_dag();
    int const flag = _dag->reval(_codelist, _wkCPI, _varout.size(), _varout.data(), valout,
                                 _varin.size(), _varin.data(), valin, options.CPINF * T(-1, 1),
                                 options.CPMAX, options.CPTHRES);
    return flag >= 0;
  }

  //! @brief Name of the external operations on the network (mc::FFMAP)
  static std::string
  tag()
  {
    return "KAN";
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
  //! @brief Number of layers
  size_t
  nlay() const
  {
    return _data.size();
  }
  //! @brief Stamp of the network: a number that no other network has had,
  //! drawn again whenever the network is set, appended to, read or reset,
  //! and kept by the copies and by a KAN the network is assigned to. The
  //! options are not part of it. See mc::MLP::uid
  size_t
  uid() const
  {
    return _uid;
  }
  //! @brief Number of parameters: segments of the edges
  size_t
  nwei() const
  {
    size_t n = 0;
    for (auto const& L : _data)
      for (auto const& row : L.pwl)
        for (auto const& f : row) n += f.size();
    return n;
  }
  //! @brief Layers of the network
  std::vector<Layer> const&
  data() const
  {
    return _data;
  }
  //! @brief Variables of the expression tree
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
};

////////////////////////////////////////////////////////////////////////

template <typename T>
inline bool
KAN<T>::_check(Layer const& L)
{
  if (L.type != PWL || !L.nout() || !L.nin()) return false;
  for (auto const& row : L.pwl)
  {
    if (row.size() != L.nin()) return false;
    for (auto const& f : row)
      if (!f.size() || f.size() != f.dy().size()) return false;
  }
  return true;
}

template <typename T>
inline void
KAN<T>::reset_data()
{
  _nin = _nout = 0;
  _data.clear();
  delete _dag;
  _dag = nullptr;
  _uid = ffhome_uid();
}

template <typename T>
inline bool
KAN<T>::set_data(std::vector<Layer> const& kan)
{
  if (kan.empty()) return false;
  for (size_t l = 0; l < kan.size(); ++l)
    if (!_check(kan[l]) || (l && kan[l].nin() != kan[l - 1].nout())) return false;

  _data = kan;
  _nin  = _data.front().nin();
  _nout = _data.back().nout();
  delete _dag;
  _dag = nullptr;
  _uid = ffhome_uid();
  return true;
}

template <typename T>
inline bool
KAN<T>::append_data(std::vector<std::vector<PWLU>> const& layer, bool const reset)
{
  Layer L;
  L.type = PWL;
  L.pwl  = layer;
  if (!_check(L) || (!reset && !_data.empty() && L.nin() != _data.back().nout())) return false;

  if (reset) _data.clear();
  _data.push_back(L);
  _nin  = _data.front().nin();
  _nout = _data.back().nout();
  delete _dag;
  _dag = nullptr;
  _uid = ffhome_uid();
  return true;
}

template <typename T>
inline std::vector<std::vector<PWLU>>
KAN<T>::relu_layer(size_t const nout, size_t const nin, size_t const R, double const* c,
                   double const* v, double const* a, double const* w, double const* b)
{
  std::vector<std::vector<PWLU>> layer(nout);
  for (size_t j = 0; j < nout; ++j)
    for (size_t i = 0; i < nin; ++i)
    {
      size_t const e = j * nin + i;
      layer[j].push_back(PWLU::relu(c[e], v[e], R, a + e * R, w + e * R, b + e * R));
    }
  return layer;
}

template <typename T>
inline bool
KAN<T>::_set_dict(NETIO::Dict const& dict, std::string* msg)
{
  if (dict.type != "KANmodel" || dict.kan.empty())
  {
    if (msg) *msg = "model of type '" + dict.type + "' where KANmodel is expected";
    return false;
  }

  std::vector<Layer> data(dict.kan.size());
  for (size_t l = 0; l < data.size(); ++l)
  {
    auto const& K = dict.kan[l];
    data[l].pwl   = relu_layer(K.nout, K.nin, K.R, K.c.data(), K.v.data(), K.a.data(), K.w.data(),
                               K.b.data());
  }
  // Scaling of the inputs, u = (x - lower) / (upper - lower), folded into the edges of the first
  // layer
  if (!dict.xlower.empty())
    for (auto& row : data.front().pwl)
      for (size_t i = 0; i < row.size(); ++i)
      {
        double const s = dict.xupper[i] - dict.xlower[i];
        row[i].xL()    = dict.xlower[i] + s * row[i].xL();
        for (auto& dx : row[i].dx()) dx *= s;
        for (auto& dy : row[i].dy()) dy /= s;
      }
  // Scaling of the outputs, z = lower + (upper - lower) y, folded into the edges of the last
  // layer, the constant in the first edge
  if (!dict.zlower.empty())
    for (size_t j = 0; j < data.back().pwl.size(); ++j)
    {
      for (auto& f : data.back().pwl[j]) f *= dict.zupper[j] - dict.zlower[j];
      data.back().pwl[j].front() += dict.zlower[j];
    }
  if (!set_data(data))
  {
    if (msg) *msg = "inconsistent layers";
    return false;
  }
  return true;
}

template <typename T>
inline bool
KAN<T>::read_data(std::string const& filename, bool const disp)
{
  NETIO::Dict dict;
  std::string msg;
  if (!NETIO::read_dict(filename, dict, msg) || !set_data(dict, &msg))
  {
    if (disp) std::cerr << "Error reading " << filename << ": " << msg << "\n";
    return false;
  }
  if (disp)
  {
    std::cout << "KAN read from " << filename << ": widths";
    for (size_t const n : dict.widths) std::cout << " " << n;
    std::cout << ", " << nwei() << " segments\n";
  }
  return true;
}

template <typename T>
inline void
KAN<T>::_build(FFVar const* varin, FFVar* varout) const
{
  std::vector<FFVar> in(varin, varin + _nin), out;
  for (auto const& L : _data)
  {
    size_t const nin = L.nin(), nout = L.nout();
    // Edges leaving each input, then their sum for each output
    std::vector<std::vector<FFVar>> edge(nout, std::vector<FFVar>(nin));
    for (size_t i = 0; i < nin; ++i)
      switch (L.type)
      {
        case PWL:
        default:
        {
          std::vector<PWLU> F(nout);
          for (size_t j = 0; j < nout; ++j) F[j] = L.pwl[j][i];
          FFPWLU<T> op;
          op.options.ENCODING     = options.ENCODING;
          op.options.SEPARATE     = options.SEPARATE;
          std::vector<FFVar> const y = op(in[i], F);
          for (size_t j = 0; j < nout; ++j) edge[j][i] = y[j];
          break;
        }
      }
    out.resize(nout);
    for (size_t j = 0; j < nout; ++j)
    {
      if (nin == 1)
      {
        out[j] = edge[j][0];
        continue;
      }
      FFLin<T> sum;
      out[j] = sum(nin, edge[j].data(), 1., 0.);
    }
    in.swap(out);
  }
  for (size_t k = 0; k < _nout; ++k) varout[k] = in[k];
}

template <typename T>
inline void
KAN<T>::_set_dag()
{
  std::pair<int, bool> const dagopt(options.ENCODING, options.SEPARATE);
  if (_dag && _dagopt == dagopt) return;
  _dagopt = dagopt;

  delete _dag;
  _dag = nullptr;
  _varin.clear();
  _varout.clear();
  _codelist.clear();
  if (!_nin || !_nout) return;

  _dag   = new FFGraph;
  _varin = _dag->add_vars(_nin, "X");
  _varout.resize(_nout);
  _build(_varin.data(), _varout.data());
}

template <typename T>
template <typename U>
inline void
KAN<T>::_eval(U const* valin, U* valout) const
{
  static thread_local std::vector<U> wk, wkprv;
  for (size_t l = 0; l < _data.size(); ++l)
  {
    Layer const& L = _data[l];
    size_t const nin = L.nin(), nout = L.nout();
    wk.swap(wkprv);
    wk.resize(nout);
    U const* in = (l ? wkprv.data() : valin);
    for (size_t j = 0; j < nout; ++j)
    {
      wk[j] = _edge(L, j, 0, in[0]);
      for (size_t i = 1; i < nin; ++i) wk[j] += _edge(L, j, i, in[i]);
    }
  }
  for (size_t k = 0; k < _nout; ++k) valout[k] = wk[k];
  // the work arrays keep their capacity but none of their entries: a reverse-mode AD type that
  // propagates when the last handle of a value is released (FADBAD++) needs them gone
  wk.clear();
  wkprv.clear();
}

template <typename T>
template <typename U>
inline void
KAN<T>::_jac(U const* valin, U* gradout, U* valout, bool const backward) const
{
  static thread_local std::vector<U> wkval, wkprv, wkjac, wkjacprv;
  static thread_local std::vector<std::vector<U>> wkder;
  size_t const nlay = _data.size();
  if (backward) wkder.resize(nlay);

  // Forward pass: values, and either the Jacobian in the inputs of the network (forward) or the
  // derivatives of the edges (backward)
  for (size_t l = 0; l < nlay; ++l)
  {
    Layer const& L = _data[l];
    size_t const nin = L.nin(), nout = L.nout();
    wkval.swap(wkprv);
    wkval.resize(nout);
    U const* in = (l ? wkprv.data() : valin);
    for (size_t j = 0; j < nout; ++j)
    {
      wkval[j] = _edge(L, j, 0, in[0]);
      for (size_t i = 1; i < nin; ++i) wkval[j] += _edge(L, j, i, in[i]);
    }
    if (backward)
    {
      wkder[l].resize(nout * nin);
      for (size_t j = 0; j < nout; ++j)
        for (size_t i = 0; i < nin; ++i) wkder[l][j * nin + i] = _dedge(L, j, i, in[i]);
      continue;
    }
    wkjac.swap(wkjacprv);
    wkjac.resize(nout * _nin);
    for (size_t j = 0; j < nout; ++j)
      for (size_t c = 0; c < _nin; ++c)
      {
        U& jac = wkjac[j * _nin + c];
        if (!l)
        {
          jac = _dedge(L, j, c, in[c]);
          continue;
        }
        jac = _dedge(L, j, 0, in[0]) * wkjacprv[c];
        for (size_t i = 1; i < nin; ++i) jac += _dedge(L, j, i, in[i]) * wkjacprv[i * _nin + c];
      }
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[k];
  wkval.clear();
  wkprv.clear();
  if (!backward)
  {
    for (size_t e = 0; e < _nout * _nin; ++e) gradout[e] = wkjac[e];
    wkjac.clear();
    wkjacprv.clear();
    return;
  }

  // Backward pass: the rows of the Jacobian through the layers together
  size_t ncur = _nout;
  wkjac.resize(_nout * _nout);
  for (size_t k = 0; k < _nout; ++k)
    for (size_t j = 0; j < _nout; ++j) wkjac[k * _nout + j] = (j == k ? 1. : 0.);
  for (size_t l = nlay; l-- > 0;)
  {
    size_t const nin = _data[l].nin(), nout = _data[l].nout();
    wkjacprv.resize(_nout * nin);
    for (size_t k = 0; k < _nout; ++k)
      for (size_t i = 0; i < nin; ++i)
      {
        U& jac = wkjacprv[k * nin + i];
        jac    = wkjac[k * ncur] * wkder[l][i];
        for (size_t j = 1; j < nout; ++j) jac += wkjac[k * ncur + j] * wkder[l][j * nin + i];
      }
    wkjac.swap(wkjacprv);
    ncur = nin;
  }
  for (size_t e = 0; e < _nout * _nin; ++e) gradout[e] = wkjac[e];
  wkjac.clear();
  wkjacprv.clear();
  for (auto& der : wkder) der.clear();
}

template <typename T>
template <typename U>
inline void
KAN<T>::_tang(U const* valin, U const* dirin, U* tanout, U* valout) const
{
  static thread_local std::vector<U> wkval, wktan, wkvalprv, wktanprv;
  for (size_t l = 0; l < _data.size(); ++l)
  {
    Layer const& L = _data[l];
    size_t const nin = L.nin(), nout = L.nout();
    wkval.swap(wkvalprv);
    wktan.swap(wktanprv);
    wkval.resize(nout);
    wktan.resize(nout);
    U const* in  = (l ? wkvalprv.data() : valin);
    U const* dir = (l ? wktanprv.data() : dirin);
    for (size_t j = 0; j < nout; ++j)
    {
      wkval[j] = _edge(L, j, 0, in[0]);
      wktan[j] = _dedge(L, j, 0, in[0]) * dir[0];
      for (size_t i = 1; i < nin; ++i)
      {
        wkval[j] += _edge(L, j, i, in[i]);
        wktan[j] += _dedge(L, j, i, in[i]) * dir[i];
      }
    }
  }
  for (size_t k = 0; k < _nout; ++k)
  {
    if (valout) valout[k] = wkval[k];
    tanout[k] = wktan[k];
  }
  wkval.clear();
  wktan.clear();
  wkvalprv.clear();
  wktanprv.clear();
}

template <typename T>
template <typename U>
inline void
KAN<T>::_adj(U const* valin, U const* wout, U* adjin, U* valout) const
{
  static thread_local std::vector<U> wkval, wkprv, wkadj, wkadjprv;
  static thread_local std::vector<std::vector<U>> wkder;
  size_t const nlay = _data.size();
  wkder.resize(nlay);
  for (size_t l = 0; l < nlay; ++l)
  {
    Layer const& L = _data[l];
    size_t const nin = L.nin(), nout = L.nout();
    wkval.swap(wkprv);
    wkval.resize(nout);
    wkder[l].resize(nout * nin);
    U const* in = (l ? wkprv.data() : valin);
    for (size_t j = 0; j < nout; ++j)
    {
      wkval[j] = _edge(L, j, 0, in[0]);
      for (size_t i = 1; i < nin; ++i) wkval[j] += _edge(L, j, i, in[i]);
      for (size_t i = 0; i < nin; ++i) wkder[l][j * nin + i] = _dedge(L, j, i, in[i]);
    }
  }
  if (valout)
    for (size_t k = 0; k < _nout; ++k) valout[k] = wkval[k];

  wkadj.assign(wout, wout + _nout);
  for (size_t l = nlay; l-- > 0;)
  {
    size_t const nin = _data[l].nin(), nout = _data[l].nout();
    wkadjprv.resize(nin);
    for (size_t i = 0; i < nin; ++i)
    {
      wkadjprv[i] = wkadj[0] * wkder[l][i];
      for (size_t j = 1; j < nout; ++j) wkadjprv[i] += wkadj[j] * wkder[l][j * nin + i];
    }
    wkadj.swap(wkadjprv);
  }
  for (size_t i = 0; i < _nin; ++i) adjin[i] = wkadj[i];
  wkval.clear();
  wkprv.clear();
  wkadj.clear();
  wkadjprv.clear();
  for (auto& der : wkder) der.clear();
}

//! @brief Kolmogorov-Arnold network as an external DAG operation: mc::FFMAP
//! on an mc::KAN
template <typename T = mc::Interval>
using FFKAN = FFMAP<T, KAN<T>>;

//! @brief Jacobian of a Kolmogorov-Arnold network as an external DAG
//! operation
template <typename T = mc::Interval>
using FFGradKAN = FFGradMAP<T, KAN<T>>;

//! @brief Directional derivative of a Kolmogorov-Arnold network as an
//! external DAG operation
template <typename T = mc::Interval>
using FFTanKAN = FFTanMAP<T, KAN<T>>;

//! @brief Adjoint of a Kolmogorov-Arnold network as an external DAG operation
template <typename T = mc::Interval>
using FFAdjKAN = FFAdjMAP<T, KAN<T>>;

}  // namespace mc

#endif
