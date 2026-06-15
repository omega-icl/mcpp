#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <sstream>

#include "ffdep.hpp"

namespace py = pybind11;

void mc_ffdep( py::module_ &m )
{

typedef mc::FFDep FD;

py::class_<FD> pyFFDep( m, "FFDep" );
py::class_<FD::Exceptions> pyFDExceptions( pyFFDep, "Exceptions" );

py::enum_<FD::Exceptions::TYPE>( pyFDExceptions, "TYPE" )
 .value( "UNDEF",  FD::Exceptions::UNDEF,  "Undefined operation" )
 .value( "INTERN", FD::Exceptions::INTERN, "Internal error" )
 .export_values();

pyFDExceptions
 .def( "ierr", &FD::Exceptions::ierr, "Error flag" )
 .def( "what", &FD::Exceptions::what, "Error description" );

py::enum_<FD::TYPE>( pyFFDep, "TYPE" )
 .value( "L", FD::TYPE::L, "Linear" )
 .value( "Q", FD::TYPE::Q, "Quadratic" )
 .value( "P", FD::TYPE::P, "Polynomial" )
 .value( "R", FD::TYPE::R, "Rational" )
 .value( "N", FD::TYPE::N, "General nonlinear" )
 .export_values();

// --- Main Class ---
pyFFDep
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
   py::init<FD const&>(),
   "Copy constructor"
 )

 // Modifiers
 .def(
   "indep",
   []( FD& self, int const ndx ){ return self.indep( ndx ); },
   py::arg("ndx"),
   "Initialize as independent variable with index 'ndx'"
 )

 // Accessors
 .def(
   "dep",
   []( FD const& self, int const ndx ){ return self.dep( ndx ); },
   py::arg("ndx"),
   "Determine dependency on variable with index 'ndx'"
 )
 .def(
   "dep",
   []( FD const& self ){ return self.dep(); },
   "Dictionary of dependencies"
 )
 .def(
   "worst",
   &FD::worst,
   "Worst-case dependency"
 )
 .def(
   "update",
   &FD::update,
   "Update dependency"
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
 .def( "__str__",  []( FD const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
 .def( "__repr__", []( FD const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
;

// --- Free Functions ---
m.def( "sqr",   []( FD const& x ){ return mc::sqr(x); } );
m.def( "sqrt",  []( FD const& x ){ return mc::sqrt(x); } );
m.def( "exp",   []( FD const& x ){ return mc::exp(x); } );
m.def( "log",   []( FD const& x ){ return mc::log(x); } );
m.def( "xlog",  []( FD const& x ){ return mc::xlog(x); } );
m.def( "cos",   []( FD const& x ){ return mc::cos(x); } );
m.def( "sin",   []( FD const& x ){ return mc::sin(x); } );
m.def( "tan",   []( FD const& x ){ return mc::tan(x); } );
m.def( "acos",  []( FD const& x ){ return mc::acos(x); } );
m.def( "asin",  []( FD const& x ){ return mc::asin(x); } );
m.def( "atan",  []( FD const& x ){ return mc::atan(x); } );
m.def( "cosh",  []( FD const& x ){ return mc::cosh(x); } );
m.def( "sinh",  []( FD const& x ){ return mc::sinh(x); } );
m.def( "tanh",  []( FD const& x ){ return mc::tanh(x); } );
m.def( "fabs",  []( FD const& x ){ return mc::fabs(x); } );
m.def( "erf",   []( FD const& x ){ return mc::erf(x); } );
m.def( "fstep", []( FD const& x ){ return mc::fstep(x); } );
m.def( "bstep", []( FD const& x ){ return mc::bstep(x); } );
m.def( "pow",   []( FD const& x, int n ){ return mc::pow(x,n); } );
m.def( "pow",   []( FD const& x, double r ){ return mc::pow(x,r); } );
m.def( "pow",   []( FD const& x, FD const& y ){ return mc::pow(x,y); } );
m.def( "min",   []( FD const& x, FD const& y ){ return mc::min(x,y); } );
m.def( "max",   []( FD const& x, FD const& y ){ return mc::max(x,y); } );
m.def( "cheb",  []( FD const& x, unsigned n ){ return mc::cheb(x,n); } );
m.def( "prod",  []( unsigned n, std::vector<FD> const& x ){ return mc::prod(n, x.data()); } );
m.def( "monom", []( unsigned n, std::vector<FD> const& x, std::vector<unsigned> const& k )
                { if( x.size() != k.size() ) throw std::invalid_argument("Size mismatch between variables and exponents");
                  return mc::monom(n, x.data(), k.data()); });
}
