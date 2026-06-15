// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "slift.hpp"

namespace py = pybind11;

void mc_slift( py::module_ &m )
{

// Typedefs
typedef mc::SLiftEnv SL;
typedef mc::SLiftVar SV;
typedef mc::FFGraph  FG;
typedef mc::FFVar    FV;
typedef mc::SMon<FV const*, mc::lt_FFVar>  SM;
typedef mc::SPoly<FV const*, mc::lt_FFVar> SP;

py::class_<SL> pySLiftEnv( m, "SLift" );
py::class_<SV> pySLiftVar( m, "SLiftVar" );

// Preserve the C++ class name as a Python alias, while keeping the existing
// shorter historical name.
m.attr( "SLiftEnv" ) = pySLiftEnv;

// Translate mc::SLiftEnv::Exceptions into a regular Python exception instead
// of pybind11's generic "unknown exception" fallback.
py::register_exception_translator( []( std::exception_ptr p ){
  try{
    if( p ) std::rethrow_exception( p );
  }
  catch( SL::Exceptions& e ){
    PyErr_SetString( PyExc_RuntimeError, e.what().c_str() );
  }
} );

pySLiftVar
 .def(
   py::init< double const >(),
   py::arg("d") = 0.,
   "Constructor as constant"
 )
 .def(
   py::init< SL*, FV const& >(),
   py::arg("env"),
   py::arg("x"),
   py::keep_alive<1,2>(), // SLiftVar stores a raw SLiftEnv*
   py::keep_alive<1,3>(), // sparse polynomial keys store FFVar const*
   "Constructor as DAG variable in environment"
 )
 .def(
   py::init<SV const&>(),
   "Copy constructor"
 )
 .def(
   "env",
   &SV::env,
   py::return_value_policy::reference_internal,
   "Associated SLift environment"
 )
 .def(
   "numer",
   &SV::numer,
   py::return_value_policy::reference_internal,
   "Get numerator sparse polynomial"
 )
 .def(
   "denom",
   &SV::denom,
   py::return_value_policy::reference_internal,
   "Get denominator sparse polynomial"
 )
 .def(
   "set",
   &SV::set,
   py::arg("env"),
   py::arg("x"),
   py::keep_alive<1,2>(), // SLiftVar stores a raw SLiftEnv*
   py::keep_alive<1,3>(), // sparse polynomial keys store FFVar const*
   py::return_value_policy::reference_internal,
   "Initialize variable in SLift environment"
 )
 .def(
   "__str__",
   []( SV const& v ){ std::ostringstream oss; oss << v; return oss.str(); }
 )
 .def(
   "__repr__",
   []( SV const& v ){ std::ostringstream oss; oss << v; return oss.str(); }
 )
 .def( + py::self )
 .def( py::self += py::self )
 .def( py::self += double() )
 .def( py::self -= py::self )
 .def( py::self -= double() )
 .def( py::self *= py::self )
 .def( py::self *= double() )
 .def( py::self /= py::self )
 .def( py::self /= double() )
 .def( py::self + py::self )
 .def( py::self + double() )
 .def( double() + py::self )
 .def( py::self - py::self )
 .def( py::self - double() )
 .def( double() - py::self )
 .def( py::self * py::self )
 .def( py::self * double() )
 .def( double() * py::self )
 .def( py::self / py::self )
 .def( py::self / double() )
 .def( double() / py::self )
;

m.def( "inv",   []( SV const& x ){ return mc::inv(x); } );
m.def( "sqr",   []( SV const& x ){ return mc::sqr(x); } );
m.def( "sqrt",  []( SV const& x ){ return mc::sqrt(x); } );
m.def( "exp",   []( SV const& x ){ return mc::exp(x); } );
m.def( "log",   []( SV const& x ){ return mc::log(x); } );
m.def( "xlog",  []( SV const& x ){ return mc::xlog(x); } );
m.def( "cos",   []( SV const& x ){ return mc::cos(x); } );
m.def( "sin",   []( SV const& x ){ return mc::sin(x); } );
m.def( "tan",   []( SV const& x ){ return mc::tan(x); } );
m.def( "acos",  []( SV const& x ){ return mc::acos(x); } );
m.def( "asin",  []( SV const& x ){ return mc::asin(x); } );
m.def( "atan",  []( SV const& x ){ return mc::atan(x); } );
m.def( "cosh",  []( SV const& x ){ return mc::cosh(x); } );
m.def( "sinh",  []( SV const& x ){ return mc::sinh(x); } );
m.def( "tanh",  []( SV const& x ){ return mc::tanh(x); } );
m.def( "fabs",  []( SV const& x ){ return mc::fabs(x); } );
m.def( "erf",   []( SV const& x ){ return mc::erf(x); } );
m.def( "fstep", []( SV const& x ){ return mc::fstep(x); } );
m.def( "pow",   []( SV const& x, int const n ){ return mc::pow(x,n); } );
m.def( "pow",   []( SV const& x, double const& r ){ return mc::pow(x,r); } );
m.def( "cheb",  []( SV const& x, unsigned const n ){ return mc::cheb(x,n); } );
m.def( "prod",
       []( std::vector<SV> const& x ){ return mc::prod( x.size(), x.data() ); },
       py::arg("x") );
m.def( "prod",
       []( unsigned n, std::vector<SV> const& x ){
         if( n > x.size() ) throw py::value_error( "prod: n exceeds len(x)" );
         return mc::prod( n, x.data() );
       },
       py::arg("n"), py::arg("x") );
m.def( "monom",
       []( std::vector<SV> const& x, std::vector<unsigned> const& k, bool const chebbasis ){
         if( k.size() != x.size() ) throw py::value_error( "monom: len(k) must match len(x)" );
         return mc::monom( x.size(), x.data(), k.data(), chebbasis );
       },
       py::arg("x"), py::arg("k"), py::arg("chebbasis") = false );
m.def( "monom",
       []( unsigned n, std::vector<SV> const& x, std::vector<unsigned> const& k, bool const chebbasis ){
         if( n > x.size() || n > k.size() ) throw py::value_error( "monom: n exceeds len(x) or len(k)" );
         return mc::monom( n, x.data(), k.data(), chebbasis );
       },
       py::arg("n"), py::arg("x"), py::arg("k"), py::arg("chebbasis") = false );
m.def( "max",   []( SV const& x, SV const& y ){ return mc::max(x,y); } );
m.def( "min",   []( SV const& x, SV const& y ){ return mc::min(x,y); } );
m.def( "lmtd",  []( SV const& x, SV const& y ){ return mc::lmtd(x,y); } );
m.def( "rlmtd", []( SV const& x, SV const& y ){ return mc::rlmtd(x,y); } );

// --- Exceptions Nested Class ---
py::class_<SL::Exceptions> pySLExceptions( pySLiftEnv, "Exceptions" );

py::enum_<SL::Exceptions::TYPE>( pySLExceptions, "TYPE" )
 .value( "DAGERR",   SL::Exceptions::DAGERR,   "Operation involving a factorable expression linked to a different DAG" )
 .value( "ENVERR",   SL::Exceptions::ENVERR,   "Operation between factorable expressions linked to different environments" )
 .value( "EXTERNAL", SL::Exceptions::EXTERNAL, "Invalid external operation" )
 .value( "INTERNAL", SL::Exceptions::INTERNAL, "Internal error" )
 .export_values();

pySLExceptions
 .def( "ierr", &SL::Exceptions::ierr, "Error flag" )
 .def( "what", &SL::Exceptions::what, "Error description" )
;

// --- Options Nested Struct ---
py::class_<SL::Options> pySLOptions( pySLiftEnv, "Options" );

pySLOptions
 .def( py::init<>(), "Default constructor" )
 .def( py::init< SL::Options const& >(), "Copy constructor" )
 .def( "reset", &SL::Options::reset, "Reset options to defaults" )
 .def_readwrite( "KEEPFACT", &SL::Options::KEEPFACT, "Whether to keep existing factorisations [Default: True]" )
 .def_readwrite( "NOAUXIL",  &SL::Options::NOAUXIL, "Whether to reformulate expressions without introducing auxiliary variables or lifting constraints [Default: False]" )
 .def_readwrite( "LIFTDIV",  &SL::Options::LIFTDIV,  "Whether to lift division terms using auxiliary variables [Default: False]" )
 .def_readwrite( "LIFTIPOW", &SL::Options::LIFTIPOW, "Whether to lift integral power terms using auxiliary variables [Default: False]" )
 .def_readwrite( "LOG2EXP",  &SL::Options::LOG2EXP,  "Whether to convert log univariate to exp univariate [Default: True]" )
 .def_readwrite( "SQRT2SQR", &SL::Options::SQRT2SQR, "Whether to convert sqrt univariate to sqr univariate [Default: True]" )
 .def_readwrite( "ACOS2COS", &SL::Options::ACOS2COS, "Whether to convert acos univariate to cos univariate [Default: True]" )
 .def_readwrite( "ASIN2SIN", &SL::Options::ASIN2SIN, "Whether to convert asin univariate to sin univariate [Default: True]" )
 .def_readwrite( "ATAN2TAN", &SL::Options::ATAN2TAN, "Whether to convert atan univariate to tan univariate [Default: False]" )
 .def_readwrite( "TAN2ATAN", &SL::Options::TAN2ATAN, "Whether to convert tan univariate to atan univariate [Default: True]" )
 .def_readwrite( "DISPFULL", &SL::Options::DISPFULL, "Whether to display a full subgraph of the lifted expressions [Default: False]" )
;

// --- SLiftEnv Main Class ---
pySLiftEnv
 .def(
   py::init<FG*>(),
   py::arg("dag") = nullptr,
   py::keep_alive<1,2>(), // SLiftEnv stores a raw FFGraph*
   "Constructor"
 )
 .def_readwrite(
   "options",
   &SL::options,
   "Options for lifting subexpressions"
 )
 .def_property_readonly(
   "dag",
   &SL::dag,
   py::return_value_policy::reference_internal,
   "Attached DAG"
 )
 .def(
   "set",
   &SL::set,
   py::arg("dag"),
   py::keep_alive<1,2>(), // SLiftEnv stores a raw FFGraph*
   "Attached DAG"
 )
 .def(
   "reset",
   &SL::reset, "Reset intermediate expressions"
 )
 .def(
   "process",
   py::overload_cast< std::vector<FV> const&, bool const >( &SL::process ),
   py::arg("vDep"),
   py::arg("add2dag") = true,
   "Process dependents in 'vDep'"
 )
 .def(
   "process",
   py::overload_cast< std::set<unsigned> const&, std::vector<FV> const&, bool const >( &SL::process ),
   py::arg("ndxDep"),
   py::arg("vDep"),
   py::arg("add2dag") = true,
   "Process selected dependents in 'vDep' index by 'ndxDep'"
 )
 // Accessors
 .def_property_readonly(
   "dep",
   &SL::Dep,
   py::return_value_policy::reference_internal,
   "List of subexpressions for initial dependents"
 )
 .def_property_readonly( 
   "poly",
   &SL::Poly,
   py::return_value_policy::reference_internal,
   "List of lifted polynomial constraints"
 )
 .def_property_readonly(
   "trans",
   &SL::Trans,
   py::return_value_policy::reference_internal,
   "List of lifted transcendental constraints"
 )
 .def_property_readonly(
   "var",
   &SL::Var,
   py::return_value_policy::reference_internal,
   "List of new (auxiliary) variables participating in lifted expressions"
 )
 .def_property_readonly(
   "aux",
   //&SL::Aux,
   []( SL& self )
   {
     auto const& aux_set = self.Aux();
     return std::vector<std::pair<mc::FFVar const*,mc::FFVar const*>>( aux_set.begin(), aux_set.end() );
   },
   py::return_value_policy::reference_internal,
   "List of tuples between initial DAG intermediates and new (auxiliary) variables"
 )
 .def_property_readonly(
   "OpLift", 
   &SL::OpLift,
   py::return_value_policy::reference_internal,
   "List of intermediate expressions for lifted operations (tuples: {Op, [Operands]})" 
 )
 .def_property_readonly(
   "AuxLift",
   &SL::AuxLift,
   py::return_value_policy::reference_internal,
   "List of intermediate expressions for lifted variables (tuples: {Var, SLiftVar})" 
 )
 .def(
   "insert_dag",
   py::overload_cast< SM const&, bool const, bool const >( &SL::insert_dag ),
   py::arg("mon"),
   py::arg("useprod") = false,
   py::arg("dagaux") = false,
   "Transcribe sparse monomial into DAG" 
 )
 .def(
   "insert_dag",
   py::overload_cast< SP const&, bool const, bool const >( &SL::insert_dag ),
   py::arg("poly"),
   py::arg("useprod") = false,
   py::arg("dagaux") = false,
   "Transcribe sparse polynomial into DAG" 
 )
 .def(
   "lift",
   []( SL& E, unsigned const nres, std::vector<SV> const& vVar )
   {
     std::vector<SV> vRes( nres );
     E.lift( nres, vRes.data(), vVar.size(), vVar.data() );
     return vRes;
   },
   py::arg("nres"),
   py::arg("vVar"),
   "Lift an external operation and return the lifted result variables"
 )
 .def(
   "__str__",
   []( SL const& E ){ std::ostringstream Ess; Ess << E; return Ess.str(); }
 )
 .def(
   "__repr__",
   []( SL const& E ){ std::ostringstream Ess; Ess << E; return Ess.str(); }
 )
;

}
