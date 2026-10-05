// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

// Python bindings of the linear operators of ocbase.hpp: FFPartial, FFIntegral,
// FFEval. They are DAG operations (FFOp subclasses) and belong with
// FFGraph/FFVar in pymcpp; FFDom, the domain description consumed by FFModel,
// is bound by the cronos module.
//
// The overloads mirror the C++ ones.  Where C++ takes a monomial of directions
// (SMon<FFVar,lt_FFVar>), Python takes pymcpp.FFMon (ffmon.cpp: SMon<FFVar
// const*,lt_FFVar>), converted by value at the call.  Sides are ints, 0 = MINUS
// (tau^-, default) and 1 = PLUS (tau^+), so cronos.FFDom.MINUS / PLUS (IntEnum
// values) can be passed directly.
#include "ocbase.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace
{
typedef mc::SMon<mc::FFVar, mc::lt_FFVar>
    t_SMon;  // ocbase.hpp's direction monomial (keyed by value)
typedef mc::SMon<mc::FFVar const*, mc::lt_FFVar>
    FFMon;  // pymcpp.FFMon (keyed by pointer, ffmon.cpp)
typedef std::map<mc::FFVar, unsigned, mc::lt_FFVar> t_Ord;
typedef std::map<mc::FFVar, double, mc::lt_FFVar> t_Coord;
typedef std::map<mc::FFVar, int, mc::lt_FFVar> t_Side;

// pymcpp.FFMon -> the operators' monomial: the same directions and orders,
// re-keyed by value (the FFMon's pointers are dereferenced at call time, while
// the Python objects they point to are alive)
t_SMon
smon(FFMon const& mon)
{
  if (mon.expr.empty())
    throw std::invalid_argument("constant monomial: no direction given");
  t_Ord ord;
  for (auto const& [pv, k] : mon.expr) ord[*pv] = k;
  return t_SMon(ord);
}

t_Ord
orders(t_SMon const& mon)
{
  t_Ord o;
  for (auto const& [v, k] : mon.expr) o[v] = k;
  return o;
}
}  // namespace

void
mc_ocbase(py::module_& m)
{
  // A dict {FFVar: order} converts implicitly to pymcpp.FFMon, so every
  // operator taking a monomial of directions also takes the dict, e.g.
  // OpI(u, {t: 1, z: 1}). Requires mc_ffmon() to have run (pymcpp.cpp order).
  py::implicitly_convertible<std::map<mc::FFVar const*, unsigned, mc::lt_FFVar>,
                             FFMon>();

  // --- FFPartial ---
  py::class_<mc::FFPartial, mc::FFOp> pyFFPartial(m, "FFPartial", R"doc(
Partial derivative operator on the DAG.

``OpP(u, t)`` is the derivative of ``u`` along direction ``t``; ``OpP(u,
FFMon(t, 2))`` the second derivative; ``OpP(u, FFMon({t: 1, z: 1}))`` the
mixed derivative. Nested applications fold into a single node.

The operation is linear in its operand: forward differentiation
(``FFGraph.FAD``/``DFAD``) applies the operator to the tangent. Backward
differentiation is refused.

Example:

>>> from pymcpp import FFGraph, FFPartial
>>> G = FFGraph(); t = G.add_var("t"); u = G.add_var("u(t)")
>>> OpP = FFPartial()
>>> du = OpP(u, t)
)doc");
  pyFFPartial.def(py::init<>(), "Default constructor.")
      .def(
          "__call__",
          [](mc::FFPartial& self, mc::FFVar const& var, mc::FFVar const& dir)
          { return self(var, dir); }, py::arg("var"), py::arg("dir"), R"doc(
First derivative of ``var`` along ``dir``.

Parameters
----------
var : FFVar
    Operand.
dir : FFVar
    Direction (a variable with a domain declared on the model).

Returns
-------
FFVar
    The derivative node.
)doc")
      .def(
          "__call__",
          [](mc::FFPartial& self, mc::FFVar const& var, FFMon const& indep)
          { return self(var, smon(indep)); }, py::arg("var"), py::arg("indep"),
          R"doc(
Derivative of ``var`` given by a monomial of directions: ``FFMon({t: 1, z:
1})`` is the mixed derivative, ``FFMon(z, 2)`` the second derivative along
``z``.

Parameters
----------
var : FFVar
    Operand.
indep : FFMon
    Directions and their orders; a dict ``{direction: order}`` converts to
    ``FFMon``.
)doc")
      .def(
          "__call__",
          [](mc::FFPartial& self, std::vector<mc::FFVar> const& vars,
             mc::FFVar const& dir) { return self(vars, dir); },
          py::arg("vars"), py::arg("dir"), R"doc(
First derivatives of every variable in ``vars`` along ``dir``, one node each.
)doc")
      .def(
          "__call__",
          [](mc::FFPartial& self, std::vector<mc::FFVar> const& vars,
             FFMon const& indep) { return self(vars, smon(indep)); },
          py::arg("vars"), py::arg("indep"), R"doc(
Derivatives of every variable in ``vars`` given by a monomial of directions.
)doc")
      .def_property_readonly(
          "indep",
          [](mc::FFPartial const& self) { return orders(self.Indep()); },
          "Directions and orders of this operation, as ``{direction: "
          "order}``.");

  // --- FFIntegral ---
  py::class_<mc::FFIntegral, mc::FFOp> pyFFIntegral(m, "FFIntegral", R"doc(
Integral operator on the DAG.

``OpI(u, t)`` is the integral of ``u`` over the domain of ``t``; ``OpI(u,
FFMon({t: 1, z: 1}))`` the integral over both. Nested integrals over disjoint
directions fold into a single node. The domains are declared on the model
(``cronos.FFModel.add_domain``), not here.

The operation is linear in its operand: forward differentiation applies the
operator to the tangent.
)doc");
  pyFFIntegral.def(py::init<>(), "Default constructor.")
      .def(
          "__call__",
          [](mc::FFIntegral& self, mc::FFVar const& var, mc::FFVar const& dir)
          { return self(var, dir); }, py::arg("var"), py::arg("dir"), R"doc(
Integral of ``var`` over the domain of ``dir``.

Parameters
----------
var : FFVar
    Integrand.
dir : FFVar
    Direction integrated over.
)doc")
      .def(
          "__call__",
          [](mc::FFIntegral& self, mc::FFVar const& var, FFMon const& indep)
          { return self(var, smon(indep)); }, py::arg("var"), py::arg("indep"),
          R"doc(
Integral of ``var`` over the domains of the directions of a monomial, e.g.
``FFMon({t: 1, z: 1})`` for the integral over both ``t`` and ``z``; every
order must be 1.
)doc")
      .def(
          "__call__",
          [](mc::FFIntegral& self, std::vector<mc::FFVar> const& vars,
             mc::FFVar const& dir) { return self(vars, dir); },
          py::arg("vars"), py::arg("dir"),
          "Integrals of every variable in ``vars`` over the domain of ``dir``.")
      .def(
          "__call__",
          [](mc::FFIntegral& self, std::vector<mc::FFVar> const& vars,
             FFMon const& indep) { return self(vars, smon(indep)); },
          py::arg("vars"), py::arg("indep"),
          "Integrals of every variable in ``vars`` over the directions of a "
          "monomial.")
      .def_property_readonly(
          "indep",
          [](mc::FFIntegral const& self) { return orders(self.Indep()); },
          "Directions integrated over, as ``{direction: 1}``.");

  // --- FFEval ---
  py::class_<mc::FFEval, mc::FFOp> pyFFEval(m, "FFEval", R"doc(
Point-evaluation operator on the DAG.

``OpE(u, t, tau)`` is ``u`` at ``t = tau``; ``OpE(u, FFMon({t: 1, z: 1}), {t:
1., z: .5})`` is ``u`` at a point of several directions. At a coordinate where
the value may differ on either side (an element boundary of a discontinuous
input, a transition), ``side`` selects the limit: 0 = ``MINUS`` (tau^-,
default), 1 = ``PLUS`` (tau^+); the values ``cronos.FFDom.MINUS`` / ``PLUS``
may be passed directly.

The operation is linear in its operand: forward differentiation applies the
operator to the tangent.
)doc");
  pyFFEval.def(py::init<>(), "Default constructor.")
      .def(
          "__call__",
          [](mc::FFEval& self, mc::FFVar const& var, mc::FFVar const& dir,
             double coord, int side) { return self(var, dir, coord, side); },
          py::arg("var"), py::arg("dir"), py::arg("coord"), py::arg("side") = 0,
          R"doc(
Value of ``var`` at ``dir = coord``.

Parameters
----------
var : FFVar
    Operand.
dir : FFVar
    Direction evaluated.
coord : float
    Coordinate.
side : int, optional
    0 = MINUS (tau^-, default), 1 = PLUS (tau^+); ``cronos.FFDom.MINUS`` /
    ``PLUS`` may be passed.
)doc")
      .def(
          "__call__",
          [](mc::FFEval& self, mc::FFVar const& var, FFMon const& indep,
             t_Coord const& coord, t_Side const& side)
          { return self(var, smon(indep), coord, side); },
          py::arg("var"), py::arg("indep"), py::arg("coord"),
          py::arg("side") = t_Side(), R"doc(
Value of ``var`` at a point of several directions.

Parameters
----------
var : FFVar
    Operand.
indep : FFMon
    Directions evaluated, each of order 1, e.g. ``FFMon({t: 1, z: 1})`` or the
    dict ``{t: 1, z: 1}``.
coord : dict of FFVar to float
    Coordinate per direction.
side : dict of FFVar to int, optional
    Side per direction where it matters (0 = MINUS, default; 1 = PLUS).
)doc")
      .def(
          "__call__",
          [](mc::FFEval& self, std::vector<mc::FFVar> const& vars,
             mc::FFVar const& dir, double coord, int side)
          { return self(vars, dir, coord, side); },
          py::arg("vars"), py::arg("dir"), py::arg("coord"),
          py::arg("side") = 0,
          "Values of every variable in ``vars`` at ``dir = coord``.")
      .def(
          "__call__",
          [](mc::FFEval& self, std::vector<mc::FFVar> const& vars,
             FFMon const& indep, t_Coord const& coord, t_Side const& side)
          { return self(vars, smon(indep), coord, side); },
          py::arg("vars"), py::arg("indep"), py::arg("coord"),
          py::arg("side") = t_Side(),
          "Values of every variable in ``vars`` at a point of several "
          "directions.")
      .def_property_readonly(
          "coord", [](mc::FFEval const& self) { return self.Coord(); },
          "The evaluation point, ``{direction: coordinate}``.")
      .def_property_readonly(
          "side", [](mc::FFEval const& self) { return self.Side(); },
          "Sides stored as PLUS, ``{direction: 1}`` (MINUS is the default).");
}
