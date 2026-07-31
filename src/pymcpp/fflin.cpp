#include <pybind11/operators.h>
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
    boost::numeric::interval_lib::rounded_transc_opp<double> >
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

#include "fflin.hpp"

namespace py = pybind11;

void
mc_fflin(py::module_& m)
{
  py::class_<mc::FFLin<I>, mc::FFOp> pyFFLin(m, "FFLin", R"doc(
Linear combination as a single external DAG operation.

FFLin defines the affine expression ``sum_k Coef[k]*Var[k] + Bias`` as one
external operation node in an :class:`FFGraph` DAG, instead of a chain of
elementary sum and product nodes. This keeps the DAG compact for
large-scale linear algebra (it is, for instance, the building block used
internally by :class:`MLP` for the neuron pre-activations) and enables
exact, single-pass propagation of the linear expression in every supported
arithmetic: float, Interval, McCormick, Specbnd, Taylor/Chebyshev models
(TVar, CVar, SCVar), superposition models (PWCSVar, PWLSVar), and
polyhedral relaxation. Reverse (constraint) propagation in interval and
polyhedral arithmetic as well as symbolic differentiation are supported,
so DAGs containing FFLin nodes can be evaluated, bounded, relaxed and
differentiated like any other DAG.

Instantiate the operation with the default constructor, then call it with
the participating variables and coefficients to insert the operation into
their DAG and obtain the resulting dependent.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_vars(3)
>>> Sum = pymcpp.FFLin()
>>> S1 = Sum(X, [1., 2., 3.], 0.5)  # 1*X0 + 2*X1 + 3*X2 + 0.5
>>> S2 = Sum(X, 1.)                 # X0 + X1 + X2
>>> DAG.eval([S1], X, [1., 1., 1.])
[6.5]
)doc");

  pyFFLin
      .def(py::init<>(), R"doc(
Construct an unattached linear-combination operation.

Call the resulting object with variables and coefficients to insert the
operation into a DAG.
)doc")
      .def(
          "__call__",
          [](mc::FFLin<I>& self, std::vector<mc::FFVar> const& Var,
             std::vector<double>& Coef, double const& Bias)
          { return self(Var, Coef, Bias); },
          py::arg("Var"), py::arg("Coef"), py::arg("Bias") = 0.,
          py::return_value_policy::reference_internal,
          R"doc(
Insert the linear combination ``sum_k Coef[k]*Var[k] + Bias`` in the DAG.

A single external operation node is added to the DAG that the variables in
`Var` belong to; the coefficients are deep-copied into the DAG.

Parameters
----------
Var : list of FFVar
    Participating DAG variables; must all belong to the same DAG.
Coef : list of float
    Weight coefficients, one per entry of `Var` (same length as `Var`).
Bias : float, optional
    Additive constant. Default is 0.

Returns
-------
res : FFVar
    Dependent DAG variable representing the linear combination.
)doc")
      .def(
          "__call__",
          [](mc::FFLin<I>& self, std::vector<mc::FFVar> const& Var,
             double const& Coef, double const& Bias)
          { return self(Var, Coef, Bias); },
          py::arg("Var"), py::arg("Coef"), py::arg("Bias") = 0.,
          py::return_value_policy::reference_internal,
          R"doc(
Insert the scaled sum ``Coef * sum_k Var[k] + Bias`` in the DAG.

Same as the first overload but applies a single scalar coefficient to
every variable.

Parameters
----------
Var : list of FFVar
    Participating DAG variables; must all belong to the same DAG.
Coef : float
    Common weight coefficient multiplying each variable.
Bias : float, optional
    Additive constant. Default is 0.

Returns
-------
res : FFVar
    Dependent DAG variable representing the linear combination.
)doc")
      .def(
          "Coef",
          [](mc::FFLin<I>& self) {
            return std::vector<double>(self.Coef(), self.Coef() + self.nCoef());
          },
          py::return_value_policy::move, R"doc(
Return the weight coefficients of this operation.

Returns
-------
coef : list of float
    New list with the weight coefficients; a single element if the
    operation was defined with a common scalar coefficient.
)doc")
      .def(
          "Bias", [](mc::FFLin<I>& self) { return self.Bias(); },
          R"doc(
Return the additive bias constant of this operation.

Returns
-------
bias : float
    Bias value.
)doc")
      // .def_readwrite(
      //   "type",
      //   &mc::FFLin<I>::type,
      //   "retrieve operation type"
      // )
      // .def_readwrite(
      //   "info",
      //   &mc::FFLin<I>::info,
      //   "retrieve operation id"
      // )
      // .def_readwrite(
      //   "varin",
      //   &mc::FFLin<I>::varin
      // )
      // .def_readwrite(
      //   "varout",
      //   &mc::FFLin<I>::varout
      // )
      .def("name", &mc::FFLin<I>::name, R"doc(
Return the name of this operation as displayed in the DAG (e.g.
"Lin[0x...]").
)doc")
      .def("__str__",
           [](mc::FFLin<I> const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFLin<I> const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           });
}
