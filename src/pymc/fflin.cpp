#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#ifdef MC__USE_PROFIL
 #include "mcprofil.hpp"
 typedef INTERVAL I;
#else
 #ifdef MC__USE_FILIB
  #include "mcfilib.hpp"
  typedef filib::interval<double,filib::native_switched,filib::i_mode_extended> I;
 #else
  #ifdef MC__USE_BOOST
   #include "mcboost.hpp"
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "fflin.hpp" 

namespace py = pybind11;


void mc_fflin( py::module_ &m )
{

py::class_<mc::FFLin<I>, mc::FFOp> pyFFLin( m, "FFLin" );

pyFFLin
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   "__call__",
   []( mc::FFLin<I>& self, std::vector<mc::FFVar> const& Var, std::vector<double>& Coef, double const& Bias )
   {
     return self( Var, Coef, Bias );
   },
   py::arg("Var"),
   py::arg("Coef"),
   py::arg("Bias") = 0.,
   py::return_value_policy::reference_internal,
   "define linear combination operation in DAG"
 )
 .def(
   "__call__",
   []( mc::FFLin<I>& self, std::vector<mc::FFVar> const& Var, double const& Coef, double const& Bias )
   {
     return self( Var, Coef, Bias );
   },
   py::arg("Var"),
   py::arg("Coef"),
   py::arg("Bias") = 0.,
   py::return_value_policy::reference_internal,
   "define linear combination operation in DAG"
 )
 .def(
   "Coef",
   []( mc::FFLin<I>& self )
     { return std::vector<double>( self.Coef(), self.Coef()+self.nCoef() ); },
   py::return_value_policy::move,
   "retrieve linear coefficients"
 )
 .def(
   "Bias",
   []( mc::FFLin<I>& self )
     { return self.Bias(); },
   "retrieve bias coefficient"
 )
// .def_readwrite(
//   "type",
//   &mc::FFLin<I>::type,
//   "retreive operation type"
// )
// .def_readwrite(
//   "info",
//   &mc::FFLin<I>::info,
//   "retreive operation id"
// )
// .def_readwrite(
//   "varin",
//   &mc::FFLin<I>::varin
// )
// .def_readwrite(
//   "varout",
//   &mc::FFLin<I>::varout
// )
 .def( "name",
   &mc::FFLin<I>::name,
   "retreive operation name"
 )
 .def(
   "__str__",
   []( mc::FFLin<I> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFLin<I> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
;

}

