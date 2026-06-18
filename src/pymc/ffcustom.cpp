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
  py::class_<mc::FFCustom<I>, mc::FFOp> pyFFCustom(m, "FFCustom");

  pyFFCustom
 .def(
   py::init<>(),
   "default constructor"
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
   "define custom operation with unique identifier in DAG"
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
   "define custom operation with unique identifier in DAG"
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
   "define custom operation with unique identifier in DAG"
 )
 .def(
   "set_D_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<double>( std::vector<double> const& )> const& DEval )
   {
     self.set_eval( DEval );
   },
   "set custom evaluation function in double arithmetic"
 )
 .def(
   "set_I_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<I>( std::vector<I> const& )> const& IEval )
   {
     self.set_eval( IEval );
   },
   "set custom evaluation function in interval arithmetic"
 )
 .def(
   "set_MC_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::McCormick<I>>( std::vector<mc::McCormick<I>> const& )> const& MCEval )
   {
     self.set_eval( MCEval );
   },
   "set custom evaluation function in McCormick arithmetic"
 )
 .def(
   "set_SB_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::Specbnd<I>>( std::vector<mc::Specbnd<I>> const& )> const& SBEval )
   {
     self.set_eval( SBEval );
   },
   "set custom evaluation function in Spectral bound arithmetic"
 )
 .def(
   "set_TM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::TVar<I>>( std::vector<mc::TVar<I>> const& )> const& TMEval )
   {
     self.set_eval( TMEval );
   },
   "set custom evaluation function in Taylor model arithmetic"
 )
 .def(
   "set_CM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::CVar<I>>( std::vector<mc::CVar<I>> const& )> const& CMEval )
   {
     self.set_eval( CMEval );
   },
   "set custom evaluation function in Chebyshev model arithmetic"
 )
 .def(
   "set_SCM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SCVar<I>>( std::vector<mc::SCVar<I>> const& )> const& SCMEval )
   {
     self.set_eval( SCMEval );
   },
   "set custom evaluation function in sparse Chebyshev model arithmetic"
 )
 .def(
   "set_PWCSM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SupVar<mc::PWCU>>( std::vector<mc::SupVar<mc::PWCU>> const& )> const& PWCSMEval )
   {
     self.set_eval( PWCSMEval );
   },
   "set custom evaluation function in piecewise-constant superposition model arithmetic"
 )
 .def(
   "set_PWLSM_eval",
   []( mc::FFCustom<I>& self, std::function<std::vector<mc::SupVar<mc::PWLU>>( std::vector<mc::SupVar<mc::PWLU>> const& )> const& PWLSMEval )
   {
     self.set_eval( PWLSMEval );
   },
   "set custom evaluation function in piecewise-linear superposition model arithmetic"
 )
 .def(
   "set_deriv",
   []( mc::FFCustom<I>& self, mc::FFCustom<I> const& Deriv, int const uid )
   {
     self.set_deriv( Deriv, uid );
   },
   "set custom differentiation function"
 )
 .def(
   "uid",
   []( mc::FFCustom<I> const& self )
   {
     return self.type;
   },
   "retreive operation unique identifier"
 )
 .def(
   "name",
   &mc::FFCustom<I>::name,
   "retreive operation name"
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
