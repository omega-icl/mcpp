#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "ffvect.hpp" 

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

namespace py = pybind11;


void mc_ffvect( py::module_ &m )
{

py::class_<mc::Vect> pyVect( m, "Vect" );

pyVect
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   py::init<mc::FFGraph*, std::vector<mc::FFVar> const&, std::vector<std::vector<mc::FFVar>> const&>(),
   "constructor passing current DAG, dependent and independent variables"
 )
 .def(
   py::init<mc::FFGraph*, std::vector<mc::FFVar> const&, std::vector<mc::FFVar> const&,
            std::vector<std::vector<mc::FFVar>> const&>(),
   "constructor passing current DAG, dependent, independent and constant variables"
 )
 .def(
   py::init<mc::Vect const&>(),
   "copy constructor"
 )
 .def( 
   "set",
   []( mc::Vect& self, mc::FFGraph* pDAG, std::vector<mc::FFVar> const& vVar,
       std::vector<std::vector<mc::FFVar>> const& vFun )
     { return self.set( pDAG, vVar, vFun ); },
   "set current DAG, dependent and independent variables"
 )
 .def( 
   "set",
   []( mc::Vect& self, mc::FFGraph* pDAG, std::vector<mc::FFVar> const& vVar,
       std::vector<mc::FFVar> const& vCst, std::vector<std::vector<mc::FFVar>> const& vFun )
     { return self.set( pDAG, vVar, vCst, vFun ); },
   "set current DAG, dependent, independent and constant variables"
 )
 .def_readwrite( 
   "options",
   &mc::Vect::options
 )
 .def_property_readonly(
   "DAG",
   []( mc::Vect& self )
   {
     return self.pDAG();
   },
   "Original DAG"
 )
 .def_property_readonly(
   "Var",
   []( mc::Vect& self )
   {
     return self.vVar();
   },
   "Participating variables"
 )
 .def_property_readonly(
   "Cst",
   []( mc::Vect& self )
   {
     return self.vCst();
   },
   "Participating variables"
 )
 .def_property_readonly(
   "Dep",
   []( mc::Vect& self )
   {
     return self.vFun();
   },
   "Vectorized dependents"
 )
;

py::class_<mc::Vect::Options> pyVectOptions( pyVect, "Vect.Options" );

pyVectOptions
 .def( py::init<>() )
 .def( py::init<mc::Vect::Options const&>() )
 .def_readwrite(
   "AUTODIFF",
   &mc::Vect::Options::AUTODIFF,
   "select forward or reverse automatic differentiation [Default: F]"
 )
;

py::enum_<mc::Vect::Options::AD>(pyVectOptions, "Vect.AD")
 .value("F", mc::Vect::Options::AD::F, "Forward differentiation")
 .value("B", mc::Vect::Options::AD::B, "Backward differentiation")
 .export_values()
;

py::class_<mc::FFVect<I>, mc::FFOp> pyFFVect( m, "FFVect" );

pyFFVect
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   "__call__",
   []( mc::FFVect<I>& self, mc::Vect* pFun )
   {
     auto pDep = self( pFun );
     return std::vector<mc::FFVar*>( pDep, pDep+pFun->nFun() );
   },
   py::return_value_policy::reference_internal,
   "define vector operation in DAG"
 )
 .def_readwrite(
   "type",
   &mc::FFVect<I>::type,
   "retreive operation type"
 )
 .def_readwrite(
   "info",
   &mc::FFVect<I>::info,
   "retreive operation id"
 )
 .def_readwrite(
   "varin",
   &mc::FFVect<I>::varin
 )
 .def_readwrite(
   "varout",
   &mc::FFVect<I>::varout
 )
 .def( "name",
   &mc::FFVect<I>::name,
   "retreive operation name"
 )
 .def(
   "__str__",
   []( mc::FFVect<I> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFVect<I> const& O )
   {
     std::ostringstream Oss;
     Oss << O;
     return Oss.str();
   }
 )
;

}

