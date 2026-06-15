// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include <pybind11/pybind11.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

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

#include "ffdagext.hpp" 

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>);

void mc_ffdagext( py::module_ &m )
{

py::class_<mc::DAGEXT<I>> pyDAGEXT( m, "DAGEXT" );

pyDAGEXT
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   py::init<mc::FFGraph*, std::vector<mc::FFVar> const&, std::vector<mc::FFVar> const&>(),
   "data constructor"
 )
 .def(
   py::init<mc::DAGEXT<I> const&>(),
   "copy constructor"
 )
 .def_readwrite(
   "options",
   &mc::DAGEXT<I>::options
 )
 .def( 
   "set",
   []( mc::DAGEXT<I>& self,
       mc::FFGraph* dag,
       std::vector<mc::FFVar> const& varin,
       std::vector<mc::FFVar> const& varout )
     { self.set( dag, varin, varout ); },
   "set DAGEXT data"
 )
 .def_property_readonly(
   "dag",
   []( mc::DAGEXT<I> const& self ){ return self.dag(); },
   py::return_value_policy::reference_internal,
   "DAGEXT internal DAG"
 )
 .def_property_readonly(
   "nin",
   []( mc::DAGEXT<I> const& self ){ return self.nin(); },
   py::return_value_policy::reference_internal,
   "DAGEXT number of input variables"
 )
 .def_property_readonly(
   "varin",
   []( mc::DAGEXT<I> const& self ){ return self.varin(); },
   py::return_value_policy::reference_internal,
   "DAGEXT input variables"
 )
 .def_property_readonly(
   "nout",
   []( mc::DAGEXT<I> const& self ){ return self.nout(); },
   py::return_value_policy::reference_internal,
   "DAGEXT number of output variables"
 )
 .def_property_readonly(
   "varout",
   []( mc::DAGEXT<I> const& self ){ return self.varout(); },
   py::return_value_policy::reference_internal,
   "DAGEXT output variables"
 )
;

py::class_<mc::DAGEXT<I>::Options> pyDAGEXTOptions( pyDAGEXT, "Options" );//, py::module_local() );

pyDAGEXTOptions
 .def( py::init<>() )
 .def( py::init<mc::DAGEXT<I>::Options const&>() )
 .def( "reset", []( mc::DAGEXT<I>::Options& self ){ self.reset(); }, "reset DAGEXT options" )
 .def_readwrite( "AUTODIFF",  &mc::DAGEXT<I>::Options::AUTODIFF,  "Automatic differentiation mode [Default: F]" )
 .def_readwrite( "CPMAX",     &mc::DAGEXT<I>::Options::CPMAX,     "Maximum rounds of constraint propagation [Default: 1]" )
 .def_readwrite( "CPTHRES",   &mc::DAGEXT<I>::Options::CPTHRES,   "Threshold for repeating constraint propagation (minimum relative reduction in any variable) [Default: 1e4*DBL_EPSILON]" )
 .def_readwrite( "CPINF",     &mc::DAGEXT<I>::Options::CPINF,     "Infinite value for unbounded variables in constraint propagation [Default: 1e30]" )
;

py::enum_<mc::DAGEXT<I>::Options::AD_TYPE>(pyDAGEXTOptions, "AD_TYPE")
 .value("F", mc::DAGEXT<I>::Options::AD_TYPE::F, "Forward differentiation")
 .value("B", mc::DAGEXT<I>::Options::AD_TYPE::B, "Backward differentiation")
 .export_values()
;

py::class_<mc::FFDAGEXT<I>, mc::FFOp> pyFFDAGEXT( m, "FFDAGEXT" );

pyFFDAGEXT
 .def(
   py::init<bool const>(),
   py::arg("sparse")=true,
   "default constructor"
 )
 .def(
   py::init<mc::FFDAGEXT<I> const&>(),
   "copy constructor"
 )
 .def(
   "__call__",
   []( mc::FFDAGEXT<I>& self, std::vector<mc::FFVar> const& vVar, mc::DAGEXT<I>* pDAGEXT )
   {
     return self( vVar, pDAGEXT, mc::FFDAGEXT<I>::COPY );
   },
   py::return_value_policy::reference_internal,
   py::arg("var"),
   py::arg("dag"),
   "define DAGEXT operation in DAG"
 )
 .def(
   "__call__",
   []( mc::FFDAGEXT<I>& self, unsigned const idep, std::vector<mc::FFVar> const& vVar, mc::DAGEXT<I>* pDAGEXT )
   {
     return self( idep, vVar, pDAGEXT, mc::FFDAGEXT<I>::COPY );
   },
   py::return_value_policy::reference_internal,
   py::arg("dep"),
   py::arg("var"),
   py::arg("dag"),
   "define DAGEXT operation in DAG"
 )
 .def_property_readonly(
   "name",
   []( mc::FFDAGEXT<I> const& self ){ return self.name(); },
   py::return_value_policy::reference_internal,
   "DAGEXT operation name"
 )
 .def(
   "__str__",
   []( mc::FFDAGEXT<I> const& self )
   {
     std::ostringstream Oss;
     Oss << self;
     return Oss.str();
   }
 )
 .def(
   "__repr__",
   []( mc::FFDAGEXT<I> const& self )
   {
     std::ostringstream Oss;
     Oss << self;
     return Oss.str();
   }
 )
 .def_readwrite( 
   "options",
   &mc::FFDAGEXT<I>::options
 )
;

py::class_<mc::FFDAGEXT<I>::Options> pyFFDAGEXTOptions( pyFFDAGEXT, "Options" );

py::enum_<mc::FFDAGEXT<I>::Options::RELAX_TYPE>(pyFFDAGEXTOptions, "RELAX_TYPE")
 .value("INT",    mc::FFDAGEXT<I>::Options::RELAX_TYPE::INT,    "Interval bounds")
 .value("AUX",    mc::FFDAGEXT<I>::Options::RELAX_TYPE::AUX,    "Auxiliary variable polyhedral relaxation")
 .value("MC",     mc::FFDAGEXT<I>::Options::RELAX_TYPE::MC,     "McCormick relaxation with interval bounds")
 .value("SB",     mc::FFDAGEXT<I>::Options::RELAX_TYPE::SB,     "Spectral bounds")
 .value("SCM",    mc::FFDAGEXT<I>::Options::RELAX_TYPE::SCM,    "Sparse Chebyshev model relaxation")
 .value("PWCS",   mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWCS,   "Piecewise-constant superposition bounds")
 .value("PWLS",   mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWLS,   "Piecewise-linear superposition bounds")
 .export_values()
;

py::class_<std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>>(m, "RELAX_FFDAGEXT", py::module_local())
 .def(py::init<>())
 .def(
   "add",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s, mc::FFDAGEXT<I>::Options::RELAX_TYPE const &k )
     { s.insert(k); },
   "Add an element to the set")
 .def(
   "remove",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s, mc::FFDAGEXT<I>::Options::RELAX_TYPE const &k )
   { if( s.find(k) == s.end() ) throw py::key_error(); s.erase(k); },
   "Remove an element from the set")
 .def(
   "discard",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s, mc::FFDAGEXT<I>::Options::RELAX_TYPE const &k ){ s.erase(k); },
   "Remove an element if present")
 .def(
   "__contains__",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s, mc::FFDAGEXT<I>::Options::RELAX_TYPE const &k ){ return s.find(k) != s.end(); }
 )
 .def(
   "__len__",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s ){ return s.size(); }
 )
 .def(
   "__iter__", 
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s ){ return py::make_iterator(s.begin(), s.end()); },
   py::keep_alive<0, 1>()
 )
 .def(
   "__repr__",
   []( std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> &s )
   {
     std::ostringstream oss;
     oss << "{";
     bool first = true;
     for( auto const &v : s ){
       if( !first ) oss << ", ";
       switch( v ){
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::INT:    oss << "INT";    break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::AUX:    oss << "AUX";    break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::MC:     oss << "MC";     break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::SB:     oss << "SB";     break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::SCM:    oss << "SCM";    break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWCS:   oss << "PWCS";   break;
         case mc::FFDAGEXT<I>::Options::RELAX_TYPE::PWLS:   oss << "PWLS";   break;
       }
       first = false;
     }
     oss << "}";
     return oss.str();
   }
 );

pyFFDAGEXTOptions
 .def( py::init<>() )
 .def( py::init<mc::FFDAGEXT<I>::Options const&>() )
 .def( "reset", []( mc::FFDAGEXT<I>::Options& self ){ self.reset(); }, "reset options to defaults" )
 .def_property(
   "RELAX",
   []( mc::FFDAGEXT<I>::Options& self ) -> std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE>&
   { return self.RELAX; },
   []( mc::FFDAGEXT<I>::Options& self, py::object const& val )
   {
     // 1. Optimized path: Assignment from another RELAX_FFDAGEXT (C++ copy)
     if( py::isinstance< std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> >( val ) ){
       self.RELAX = val.cast<std::set<mc::FFDAGEXT<I>::Options::RELAX_TYPE> const&>();
       return;
     }
     // 2. Flexible path: Assignment from Python set/list/tuple
     else if( py::isinstance<py::iterable>(val) ){
       self.RELAX.clear();
       for( auto item : val )
         // Cast items to the Enum type and insert
         self.RELAX.insert( item.cast<mc::FFDAGEXT<I>::Options::RELAX_TYPE>() );
       return;
     }
     throw py::type_error( "Cannot assign to RELAX: Expected RELAX_FFDAGEXT or iterable of RELAX_TYPE enums" );
   },
   "Type of relaxations [Default: {MC}]"
 )
 .def_readwrite( "POLDEF",    &mc::FFDAGEXT<I>::Options::POLDEF,    "Whether options default to the parent polyhedral relaxation environment [Default: True]" )
 .def_readwrite( "POLIMG",    &mc::FFDAGEXT<I>::Options::POLIMG,    "Options for polyhedral relaxation environment" )
 .def_readwrite( "SBLIN",     &mc::FFDAGEXT<I>::Options::SBLIN,     "Number of linearization points of spectral bound relaxation [Default: 0 (NVAR+1 points)]" )
 .def_readwrite( "SBSEED",    &mc::FFDAGEXT<I>::Options::SBSEED,    "Random number generator seed in LHS sampler for linearization points of spectral bound relaxation [Default: 42]" )
 .def_readwrite( "SCMODORD",  &mc::FFDAGEXT<I>::Options::SCMODORD,  "Maximal degree of sparse Chebyshev model [Default: 3]" )
 .def_readwrite( "SCBERNORD", &mc::FFDAGEXT<I>::Options::SCBERNORD, "Degree of Bernstein basis conversion [Default: 0 (same as SCORD)]" )
 .def_readwrite( "SCMODEL",   &mc::FFDAGEXT<I>::Options::SCMODEL,   "Options for sparse Chebyshev model" )
 .def_readwrite( "PWCDIV",    &mc::FFDAGEXT<I>::Options::PWCDIV,    "Equipartition size in piecewise-constant superposition model [Default: 16]" )
 .def_readwrite( "PWCREL",    &mc::FFDAGEXT<I>::Options::PWCREL,    "Representation of piecewise-constant univariates - 0: continuous relaxation; 1: binary encoding [Default: 0]" )
 .def_readwrite( "PWCSLOPE",  &mc::FFDAGEXT<I>::Options::PWCSLOPE,  "Whether to append cuts from slopes in piecewise-constant superposition model relaxation [Default: True]" )
 .def_readwrite( "PWCSHADOW", &mc::FFDAGEXT<I>::Options::PWCSHADOW, "Whether to append cuts from shadow estimators in piecewise-constant superposition model relaxation [Default: True]" )
 .def_readwrite( "PWCSUP",    &mc::FFDAGEXT<I>::Options::PWCSUP,    "Options for superposition model with piecewise-constant univariate estimators" )
 .def_readwrite( "PWLINI",    &mc::FFDAGEXT<I>::Options::PWLINI,    "Initial partition size in piecewise-linear superposition model [Default: 16]" )
 .def_readwrite( "PWLMAX",    &mc::FFDAGEXT<I>::Options::PWLMAX,    "Maximal partition size in piecewise-linear superposition model [Default: 16]" )
 .def_readwrite( "PWLREL",    &mc::FFDAGEXT<I>::Options::PWLREL,    "Representation of piecewise-linear univariates - 0: continuous relaxation; 1: binary encoding; 2: SOS2 encoding [Default: 0]" )
 .def_readwrite( "PWLSHADOW", &mc::FFDAGEXT<I>::Options::PWLSHADOW, "Whether to append cuts from shadow estimators in piecewise-linear superposition model relaxation [Default: True]" )
 .def_readwrite( "PWLSUP",    &mc::FFDAGEXT<I>::Options::PWLSUP,    "Options for superposition model with adaptive piecewise-linear univariate estimators" )
;

}

