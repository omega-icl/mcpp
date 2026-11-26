// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the EPL-2.0 with GPL-2.0-or-later as a Secondary License; see the LICENSE file.
// SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <sstream>

#include "netio.hpp"

namespace py = pybind11;

void
mc_netio(py::module_& m)
{
  typedef mc::NETIO::Dict Dict;

  py::class_<mc::NETIO> pyNETIO(m, "NETIO", R"doc(
Readers of network files, with no other library.

``NETIO.read_dict`` reads a network saved as a dictionary by ``torch.save``
once, into a :class:`NETIO.Dict` that tells its type, sizes and scalings,
and that :meth:`KAN.set_data` or :meth:`MLP.set_data` turns into a network.
:meth:`KAN.read_data` and :meth:`MLP.read_data` do both steps from the file
name, where nothing else is wanted of the file.

Examples
--------
>>> import pymcpp
>>> D = pymcpp.NETIO.read_dict("ackley2d_KAN-TITO-h5-g8x8.pt")
>>> D.type, D.widths
('KANmodel', [2, 5, 1])
>>> NN = pymcpp.KAN() if D.type == "KANmodel" else pymcpp.MLP()
>>> NN.set_data(D)
)doc");

  py::class_<Dict>(pyNETIO, "Dict", R"doc(
Network read from a dictionary saved by torch.save (see
:meth:`NETIO.read_dict`): its type, sizes and scalings. The layers are
turned into a network by :meth:`KAN.set_data` or :meth:`MLP.set_data`.
)doc")
      .def_readonly("type", &Dict::type, R"doc(
str: Type of the network, ``"KANmodel"`` (for a :class:`KAN`) or
``"ANNmodel"`` (for an :class:`MLP`).
)doc")
      .def_readonly("widths", &Dict::widths, R"doc(
list of int: Numbers of inputs, of nodes of each hidden layer, and of
outputs.
)doc")
      .def_readonly("xlower", &Dict::xlower, R"doc(
list of float: Lower values of the min-max scaling of the inputs -- with
``xupper``, the box of the raw inputs; empty where the inputs are not
scaled.
)doc")
      .def_readonly("xupper", &Dict::xupper, R"doc(
list of float: Upper values of the min-max scaling of the inputs; empty
where the inputs are not scaled.
)doc")
      .def_readonly("zlower", &Dict::zlower, R"doc(
list of float: Lower values of the min-max scaling of the outputs; empty
where the outputs are not scaled.
)doc")
      .def_readonly("zupper", &Dict::zupper, R"doc(
list of float: Upper values of the min-max scaling of the outputs; empty
where the outputs are not scaled.
)doc")
      .def("__repr__", [](Dict const& D)
           {
             std::ostringstream os;
             os << "NETIO.Dict(" << D.type << ", widths";
             for (size_t const n : D.widths) os << " " << n;
             os << ")";
             return os.str();
           });

  pyNETIO.def_static(
      "read_dict",
      [](std::string const& filename, int const activ)
      {
        Dict D;
        std::string msg;
        if (!mc::NETIO::read_dict(filename, D, msg, activ))
          throw std::runtime_error("NETIO.read_dict: " + filename + ": " + msg);
        return D;
      },
      py::arg("filename"), py::arg("activ") = -1, R"doc(
Read a network saved as a dictionary by torch.save.

The dictionary holds ``model_type`` (``"KANmodel"`` or ``"ANNmodel"``),
``config``, ``state_dict`` and ``normalization``; it is read with no other
library.

Parameters
----------
filename : str
    Path to the file.
activ : MLP.ACTIV_TYPE, optional
    Activation of the hidden layers of a dense network that does not name
    it. By default, it is the one that matches the bounds held in the file.

Returns
-------
D : NETIO.Dict
    The network read.

Raises
------
RuntimeError
    If the file cannot be read or holds anything else, with the reason.
)doc");
}
