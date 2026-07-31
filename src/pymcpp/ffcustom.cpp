// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

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

#include "ffcustom.hpp"

namespace py = pybind11;

void
mc_ffcustom(py::module_& m)
{
  py::class_<mc::FFCustom<I>, mc::FFOp> pyFFCustom(m, "FFCustom", R"doc(
User-defined (black-box) operation backed by Python callbacks.

FFCustom lets you embed an arbitrary vector-valued operation
``y = f(x_0, ..., x_{n-1})`` as a single external node in an
:class:`FFGraph` DAG, with its behavior in each arithmetic supplied as
Python callables. After instantiating the operation, register one
evaluation callback per arithmetic you intend to use (``set_D_eval`` for
float, ``set_I_eval`` for Interval, ``set_MC_eval`` for McCormick, etc.),
then call the object with the participating DAG variables to insert the
operation into their DAG.

All evaluation callbacks share the same convention: they receive a single
argument, a list with the values of all operands in the respective
arithmetic (one entry per input variable, in the order the variables were
passed to ``__call__``), and must return a list with the values of all
results of the operation (one entry per dependent). Evaluating the DAG in
an arithmetic without a registered callback raises RuntimeError.

Symbolic differentiation of a DAG containing an FFCustom node (e.g. via
``FFGraph.fdiff``) requires a derivative operation registered with
:meth:`set_deriv`. Decomposition with :class:`SLift` needs no callback:
the results of the operation are automatically lifted as new auxiliary
variables.

Each distinct custom operation must be given a distinct integer identifier
`uid` when inserted in a DAG: two FFCustom nodes compare equal if and only
if they carry the same `uid`, so reusing an identifier for a different
function would silently merge the operations.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(2)
>>> def prod(x):
...     return [x[0] * x[1]]
>>> def dprod(x):
...     return [x[1], x[0]]  # d(x0*x1)/dx0, d(x0*x1)/dx1
>>> OpDY = pymcpp.FFCustom()
>>> OpDY.set_D_eval(dprod)
>>> OpY = pymcpp.FFCustom()
>>> OpY.set_D_eval(prod)
>>> OpY.set_I_eval(prod)  # interval product via overloaded operators
>>> OpY.set_deriv(OpDY, 0)
>>> Y = OpY(X, 1)
>>> DAG.eval([Y], X, [2., 3.])
[6.0]
)doc");

  pyFFCustom
 .def(
   py::init<>(),
   R"doc(
Construct a custom operation without any registered callbacks.

Register evaluation callbacks with the ``set_*_eval`` methods and
optionally a derivative operation with :meth:`set_deriv`, then call the
object with DAG variables to insert the operation into their DAG.
)doc"
 )
 .def(
   "__call__",
   []( mc::FFCustom<I>& self, std::vector<mc::FFVar> const& Var, int const uid )
   {
     return self( Var, uid );
   },
   py::arg("var"),
   py::arg("uid"),
   py::return_value_policy::reference_internal,
   R"doc(
Insert this operation with a single result into the DAG.

Adds one external operation node with the variables in `var` as operands
to the DAG they belong to.

Parameters
----------
var : list of FFVar
    Operand DAG variables; must all belong to the same DAG. Their order
    defines the order of the list passed to the evaluation callbacks.
uid : int
    Unique identifier of this operation within the DAG. Operations with
    equal `uid` are considered identical.

Returns
-------
dep : FFVar
    Dependent DAG variable representing the (scalar) operation result.
)doc"
 )
 .def(
   "__call__",
   []( mc::FFCustom<I>& self, size_t const nDep, std::vector<mc::FFVar> const& Var, int const uid )
   {
     auto ppDep = self( nDep, Var, uid );
     return std::vector<mc::FFVar*>( ppDep, ppDep+nDep );
   },
   py::arg("ndep"),
   py::arg("var"),
   py::arg("uid"),
   py::return_value_policy::reference_internal,
   R"doc(
Insert this operation with `ndep` results into the DAG.

Same as the single-result overload but for a vector-valued operation; the
evaluation callbacks must return lists of length `ndep`.

Parameters
----------
ndep : int
    Number of results (dependents) of the operation.
var : list of FFVar
    Operand DAG variables; must all belong to the same DAG.
uid : int
    Unique identifier of this operation within the DAG.

Returns
-------
dep : list of FFVar
    The `ndep` dependent DAG variables representing the operation results.
)doc"
 )
 .def(
   "__call__",
   []( mc::FFCustom<I>& self, size_t const iDep, size_t const nDep, std::vector<mc::FFVar> const& Var, int const uid )
   {
     return self( iDep, nDep, Var, uid );
   },
   py::arg("idep"),
   py::arg("ndep"),
   py::arg("var"),
   py::arg("uid"),
   py::return_value_policy::reference_internal,
   R"doc(
Insert this operation with `ndep` results and return result `idep` only.

Same as the multi-result overload but returns the single dependent with
index `idep` instead of the full list.

Parameters
----------
idep : int
    Index of the result to return (0 <= idep < ndep).
ndep : int
    Number of results (dependents) of the operation.
var : list of FFVar
    Operand DAG variables; must all belong to the same DAG.
uid : int
    Unique identifier of this operation within the DAG.

Returns
-------
dep : FFVar
    Dependent DAG variable representing result `idep`.
)doc"
 )
 .def(
   "set_D_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<double>( std::vector<double> const& )> const& DEval )
   {
     self.set_eval( DEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for float arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of float) -> list of float``,
    receiving the operand values and returning the result values (one
    entry per dependent of the operation).
)doc"
 )
 .def(
   "set_I_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<I>( std::vector<I> const& )> const& IEval )
   {
     self.set_eval( IEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for interval arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of Interval) -> list of
    Interval``, receiving operand enclosures and returning result
    enclosures. It must produce valid enclosures of the operation ranges
    for the bounds computed on the DAG to be rigorous.
)doc"
 )
 .def(
   "set_MC_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::McCormick<I>>( std::vector<mc::McCormick<I>> const& )> const& MCEval )
   {
     self.set_eval( MCEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for McCormick relaxation arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of McCormick) -> list of
    McCormick``, receiving operand relaxations and returning
    convex/concave relaxations of the results.
)doc"
 )
 .def(
   "set_SB_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::Specbnd<I>>( std::vector<mc::Specbnd<I>> const& )> const& SBEval )
   {
     self.set_eval( SBEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for spectral bound arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of Specbnd) -> list of Specbnd``.
)doc"
 )
 .def(
   "set_TM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::TVar<I>>( std::vector<mc::TVar<I>> const& )> const& TMEval )
   {
     self.set_eval( TMEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for Taylor model arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of TVar) -> list of TVar``.
)doc"
 )
 .def(
   "set_CM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::CVar<I>>( std::vector<mc::CVar<I>> const& )> const& CMEval )
   {
     self.set_eval( CMEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for Chebyshev model arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of CVar) -> list of CVar``.
)doc"
 )
 .def(
   "set_SCM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SCVar<I>>( std::vector<mc::SCVar<I>> const& )> const& SCMEval )
   {
     self.set_eval( SCMEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for sparse Chebyshev model arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of SCVar) -> list of SCVar``.
)doc"
 )
 .def(
   "set_PWCSM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SupVar<mc::PWCU>>( std::vector<mc::SupVar<mc::PWCU>> const& )> const& PWCSMEval )
   {
     self.set_eval( PWCSMEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for piecewise-constant superposition
model arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of PWCSVar) -> list of PWCSVar``.
)doc"
 )
 .def(
   "set_PWLSM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SupVar<mc::PWLU>>( std::vector<mc::SupVar<mc::PWLU>> const& )> const& PWLSMEval )
   {
     self.set_eval( PWLSMEval );
   },
   py::arg("f"),
   R"doc(
Register the evaluation callback for piecewise-linear superposition model
arithmetic.

Parameters
----------
f : callable
    Function with signature ``f(x: list of PWLSVar) -> list of PWLSVar``.
)doc"
 )
 .def(
   "set_deriv",
   []( mc::FFCustom<I>& self, mc::FFCustom<I> const& Deriv, int const uid )
   {
     self.set_deriv( Deriv, uid );
   },
   py::arg("deriv"),
   py::arg("uid"),
   R"doc(
Register another FFCustom operation as the derivative of this one.

The derivative operation `deriv` is copied and stored; when the DAG is
differentiated symbolically (e.g. ``FFGraph.fdiff`` or ``bdiff``), it is
inserted with the same operands as this operation and with
``nres * nvar`` results holding the Jacobian entries, where ``nres`` is
the number of results and ``nvar`` the number of operands of this
operation. Its evaluation callbacks must therefore return, for operand
values ``x``, the flat list ``[dF0/dx0, ..., dF0/dx{nvar-1}, dF1/dx0,
...]``, i.e. result index varying slowest (entry ``j*nvar + i`` is the
partial derivative of result ``j`` with respect to operand ``i``).

Without a registered derivative, symbolic differentiation of a DAG
containing this operation raises RuntimeError.

Parameters
----------
deriv : FFCustom
    Operation evaluating the Jacobian of this operation; copied
    internally.
uid : int
    Unique identifier given to the derivative operation when it is
    inserted in the DAG.
)doc"
 )
 .def(
   "uid",
   []( mc::FFCustom<I> const& self )
   {
     return self.type;
   },
   R"doc(
Return the internal identifier under which this operation is registered
in the DAG.

Returns
-------
uid : int
    Operation type identifier.
)doc"
 )
 .def(
   "name",
   &mc::FFCustom<I>::name,
   R"doc(
Return the name of this operation as displayed in the DAG, in the form
"Custom[<uid>]".
)doc"
 )
 .def(
   "__str__",
   []( mc::FFCustom<I> const& self )
   {
     std::ostringstream Oss;
     Oss << self;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFCustom<I> const& self )
   {
     std::ostringstream Oss;
     Oss << self;
     return Oss.str();
   }
 )
;
}
