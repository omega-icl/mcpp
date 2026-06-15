// Copyright (C) Benoit Chachuat, Imperial College London.
// All Rights Reserved.
// This code is published under the Eclipse Public License.

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <sstream>

#include "ffinv.hpp"

namespace py = pybind11;

// Tells pybind11 to not copy this std::set, but wrap it directly
PYBIND11_MAKE_OPAQUE(std::set<mc::FFInv::Options::NLINV>);

void mc_ffinv( py::module_ &m )
{

typedef mc::FFInv FI;

py::class_<FI> pyFFInv( m, "FFInv" );
py::class_<FI::Exceptions> pyFIExceptions( pyFFInv, "Exceptions" );
py::class_<FI::Options> pyFIOptions( pyFFInv, "Options" );

py::enum_<FI::Exceptions::TYPE>( pyFIExceptions, "TYPE" )
 .value( "UNDEF",  FI::Exceptions::UNDEF,  "Undefined operation" )
 .value( "INTERN", FI::Exceptions::INTERN, "Internal error" )
 .export_values();

pyFIExceptions
 .def( "ierr", &FI::Exceptions::ierr, "Error flag" )
 .def( "what", &FI::Exceptions::what, "Error description" );

py::enum_<FI::Options::NLINV>( pyFIOptions, "NLINV", "Invertible nonlinear operations" )
 .value( "INV",  FI::Options::INV,  "Inverse" )
 .value( "SQRT", FI::Options::SQRT, "Square-root" )
 .value( "EXP",  FI::Options::EXP,  "Exponential" )
 .value( "LOG",  FI::Options::LOG,  "Logarithm" )
 .value( "IPOW", FI::Options::IPOW, "Integer power" )
 .value( "RPOW", FI::Options::RPOW, "Real power" )
 .export_values();

py::class_<std::set<FI::Options::NLINV>>(m, "NLINVSet")
 .def(py::init<>())
 .def(
   "add",
   [](std::set<FI::Options::NLINV> &s, FI::Options::NLINV const &k)
     { s.insert(k); },
   "Add an element to the set")
 .def(
   "remove",
   [](std::set<FI::Options::NLINV> &s, FI::Options::NLINV const &k)
   { if( s.find(k) == s.end() ) throw py::key_error(); s.erase(k); },
   "Remove an element from the set")
 .def(
   "discard",
   [](std::set<FI::Options::NLINV> &s, FI::Options::NLINV const &k){ s.erase(k); },
   "Remove an element if present")
 .def(
   "__contains__",
   [](std::set<FI::Options::NLINV> &s, FI::Options::NLINV const &k){ return s.find(k) != s.end(); }
 )
 .def(
   "__len__",
   [](std::set<FI::Options::NLINV> &s) { return s.size(); }
 )
 .def(
   "__iter__", 
   [](std::set<FI::Options::NLINV> &s){ return py::make_iterator(s.begin(), s.end()); },
   py::keep_alive<0, 1>()
 )
 .def(
   "__repr__",
   [](std::set<FI::Options::NLINV> &s)
   {
     std::ostringstream oss;
     oss << "{";
     bool first = true;
     for( auto const &v : s ){
       if( !first ) oss << ", ";
       switch( v ){
         case FI::Options::INV:  oss << "INV";  break;
         case FI::Options::SQRT: oss << "SQRT"; break;
         case FI::Options::EXP:  oss << "EXP";  break;
         case FI::Options::LOG:  oss << "LOG";  break;
         case FI::Options::IPOW: oss << "IPOW"; break;
         case FI::Options::RPOW: oss << "RPOW"; break;
       }
       // Casting to int for simple enum representation, or you can map back to string
       //oss << v;//static_cast<int>(v); 
       first = false;
     }
     oss << "}";
     return oss.str();
   }
 )
;

pyFIOptions
 .def( py::init<>() )
 .def( py::init<FI::Options const&>() )
 .def(
   "reset",
   []( FI::Options& self ){ self = FI::Options(); },
   "Reset options to default"
 )
 //.def_readwrite( "INVOP", &FI::Options::INVOP, "Set of allowed invertible operations" );
 .def_property(
   "INVOP",
   []( FI::Options& self ) -> std::set<FI::Options::NLINV>&
   { return self.INVOP; },
   []( FI::Options& self, py::object const& val )
   {
     // 1. Optimized path: Assignment from another NLINVSet (C++ copy)
     if( py::isinstance<std::set<FI::Options::NLINV>>(val) ){
       self.INVOP = val.cast<std::set<FI::Options::NLINV> const&>();
       return;
     }
     // 2. Flexible path: Assignment from Python set/list/tuple
     else if( py::isinstance<py::iterable>(val) ){
       self.INVOP.clear();
       for( auto item : val )
         // Cast items to the Enum type and insert
         self.INVOP.insert(item.cast<FI::Options::NLINV>());
       return;
     }
     throw py::type_error( "Cannot assign to INVOP: Expected NLINVSet or iterable of NLINV enums" );
   },
   py::return_value_policy::reference_internal,
   "Set of allowed invertible operations"
 )
;

py::enum_<FI::TYPE>( pyFFInv, "TYPE" )
 .value( "L", FI::TYPE::L, "Linear" )
 .value( "S", FI::TYPE::S, "Separably linear" )
 .value( "N", FI::TYPE::N, "Separably nonlinear" )
 .value( "U", FI::TYPE::U, "Undetermined" )
 .export_values()
;
        
// --- Main Class ---
pyFFInv
 // Constructors
 .def(
   py::init<>(),
   "Default constructor" 
 )
 .def(
   py::init<double const>(),
   py::arg("val"),
   "Constructor for real constant"
 )
 .def(
   py::init<FI const&>(),
   "Copy constructor"
 )

 // Modifiers
 .def(
   "indep",
   []( FI& self, int const ndx ){ return self.indep( ndx ); },
   py::arg("ndx"),
   "Initialize as independent variable with index 'ndx'"
 )

 // Accessors
 .def_readwrite_static( 
   "options",
   &FI::options
 )
 .def(
   "inv",
   []( FI const& self, int const ndx ){ return self.inv( ndx ); },
   py::arg("ndx"),
   "Determine invertibility of variable with index 'ndx'"
 )
 .def(
   "inv",
   []( FI const& self ){ return self.inv(); },
   "Dictionary of invertibilities"
 )
 .def(
   "update",
   &FI::update,
   "Update invertibility"
 )

// Operators
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

 // String representation
 .def( "__str__",  []( FI const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
 .def( "__repr__", []( FI const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
;

m.def( "sqr",   []( FI const& x ){ return mc::sqr(x); } );
m.def( "sqrt",  []( FI const& x ){ return mc::sqrt(x); } );
m.def( "exp",   []( FI const& x ){ return mc::exp(x); } );
m.def( "log",   []( FI const& x ){ return mc::log(x); } );
m.def( "xlog",  []( FI const& x ){ return mc::xlog(x); } );
m.def( "cos",   []( FI const& x ){ return mc::cos(x); } );
m.def( "sin",   []( FI const& x ){ return mc::sin(x); } );
m.def( "tan",   []( FI const& x ){ return mc::tan(x); } );
m.def( "acos",  []( FI const& x ){ return mc::acos(x); } );
m.def( "asin",  []( FI const& x ){ return mc::asin(x); } );
m.def( "atan",  []( FI const& x ){ return mc::atan(x); } );
m.def( "cosh",  []( FI const& x ){ return mc::cosh(x); } );
m.def( "sinh",  []( FI const& x ){ return mc::sinh(x); } );
m.def( "tanh",  []( FI const& x ){ return mc::tanh(x); } );
m.def( "fabs",  []( FI const& x ){ return mc::fabs(x); } );
m.def( "erf",   []( FI const& x ){ return mc::erf(x); } );
m.def( "fstep", []( FI const& x ){ return mc::fstep(x); } );
m.def( "bstep", []( FI const& x ){ return mc::bstep(x); } );
m.def( "pow",   []( FI const& x, int n ){ return mc::pow(x,n); } );
m.def( "pow",   []( FI const& x, double r ){ return mc::pow(x,r); } );
m.def( "pow",   []( FI const& x, FI const& y ){ return mc::pow(x,y); } );
m.def( "min",   []( FI const& x, FI const& y ){ return mc::min(x,y); } );
m.def( "max",   []( FI const& x, FI const& y ){ return mc::max(x,y); } );
m.def( "cheb",  []( FI const& x, unsigned n ){ return mc::cheb(x,n); } );
m.def( "prod",  []( unsigned n, std::vector<FI> const& x ){ return mc::prod(n, x.data()); } );
m.def( "monom", []( unsigned n, std::vector<FI> const& x, std::vector<unsigned> const& k )
                { if( x.size() != k.size() ) throw std::invalid_argument("Size mismatch between variables and exponents");
                  return mc::monom(n, x.data(), k.data()); });
}
