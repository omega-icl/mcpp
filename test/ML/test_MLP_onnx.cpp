// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

// Networks read from ONNX files by mc::MLP::read_data (mc::NETIO::read_onnx), against their values computed
// by ONNX Runtime (onnx_values.txt, from test_MLP_onnx.py), and the models that must be refused.

#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "ffmlp.hpp"
#include "interval.hpp"

typedef mc::Interval I;

static int npass = 0, nfail = 0;
static void
check(bool const ok, std::string const& what)
{
  std::cout << (ok ? "  ok    " : "  FAIL  ") << what << std::endl;
  (ok ? npass : nfail)++;
}

int
main()
{
  std::cout << "ONNX networks against ONNX Runtime\n";
  std::ifstream is("onnx_values.txt");
  if (!is)
  {
    std::cout << "onnx_values.txt not found: run test_MLP_onnx.py\n";
    return 1;
  }
  std::string line, last;
  mc::MLP<I> mlp;
  while (std::getline(is, line))
  {
    std::istringstream ls(line);
    std::string name;
    size_t nin, nout;
    ls >> name >> nin >> nout;
    std::vector<double> x(nin), y(nout), v(nout);
    for (auto& xi : x) ls >> xi;
    for (auto& yi : y) ls >> yi;
    if (name != last)
    {
      bool const ok = mlp.read_data(name);
      check(ok && mlp.nin() == nin && mlp.nout() == nout,
            name + ": read, " + std::to_string(mlp.nin()) + " inputs, " +
                std::to_string(mlp.nhid()) + " hidden layers, " +
                std::to_string(mlp.nout()) + " outputs");
      last = name;
    }
    if (mlp.nin() != nin || mlp.nout() != nout) continue;
    mlp.eval(x.data(), v.data());
    // ONNX Runtime evaluates in the precision of the model, single (or half)
    // precision in most of them
    double const tol = (name.find("float16") != std::string::npos ? 2e-2
                        : name.find("double") != std::string::npos ||
                                name.find("matmul") != std::string::npos
                            ? 1e-12
                            : 2e-6);
    double err = 0.;
    for (size_t k = 0; k < nout; ++k)
      err = std::max(err, std::fabs(v[k] - y[k]) / (1. + std::fabs(y[k])));
    std::ostringstream os;
    os << name << ": values, relative error " << err;
    check(err <= tol, os.str());
  }

  std::cout << "\nModels refused\n";
  for (std::string const name :
       {"onnx_bad_leakyrelu.onnx", "onnx_bad_skip.onnx", "onnx_bad_twoinputs.onnx",
        "onnx_bad_size.onnx", "onnx_bad_external.onnx", "onnx_bad_truncated.onnx",
        "onnx_missing.onnx"})
  {
    mc::NETIO::t_Data data;
    std::string msg;
    bool const ok = mc::NETIO::read_onnx(name, data, msg);
    check(!ok && !msg.empty(), name + ": " + msg);
  }

  // A refused file leaves the network as it was
  mc::MLP<I> mlp2;
  mlp2.read_data("onnx_torch_tanh.onnx");
  size_t const uid = mlp2.uid();
  check(!mlp2.read_data("onnx_bad_skip.onnx") && mlp2.uid() == uid && mlp2.nin() == 3,
        "a refused file leaves the network as it was");

  // The network of an ONNX file has no Torch Module: EVALTORCH has no effect
  mlp2.options.EVALTORCH = true;
  std::vector<double> x = {0.1, -0.2, 0.3}, y1(2), y2(2);
  mlp2.eval(x.data(), y1.data());
  mlp2.options.EVALTORCH = false;
  mlp2.eval(x.data(), y2.data());
  check(y1 == y2, "no Torch Module with an ONNX file: EVALTORCH has no effect");

  std::cout << "\n=== Results: " << npass << "/" << npass + nfail << " passed ===\n";
  return nfail ? 1 : 0;
}
