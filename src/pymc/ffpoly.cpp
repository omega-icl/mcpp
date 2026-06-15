#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "spoly.hpp"
#include "ffunc.hpp"

namespace py = pybind11;

void mc_ffpoly( py::module_ &m )
{

typedef mc::FFVar const*     KEY;
typedef mc::lt_FFVar         COMP;
typedef mc::SPoly<KEY, COMP> SP;
typedef mc::SMon<KEY, COMP>  SM;

py::class_<SP> pyFFPoly( m, "FFPoly" );

// --- Exceptions Nested Class ---
py::class_<SP::Exceptions> pySPExceptions( pyFFPoly, "Exceptions" );
py::enum_<SP::Exceptions::TYPE>( pySPExceptions, "TYPE" )
 .value( "DIVZERO",  SP::Exceptions::DIVZERO,  "Scalar division by zero" )
 .value( "DIVPOLY",  SP::Exceptions::DIVPOLY,  "Division between two polynomials" )
 .value( "INTERNAL", SP::Exceptions::INTERNAL, "Internal error" )
 .export_values();

pySPExceptions
 .def( "ierr", &SP::Exceptions::ierr, "Return error flag" )
 .def( "what", &SP::Exceptions::what, "Return error description" )
;

// --- Options Nested Struct ---
py::class_<SP::Options> pySPOptions( pyFFPoly, "Options" );

py::enum_<SP::Options::BASIS_TYPE>( pySPOptions, "BASIS_TYPE" )
 .value( "MONOM", SP::Options::MONOM, "Monomial basis" )
 .value( "CHEB",  SP::Options::CHEB,  "Chebyshev basis" )
 .export_values();

pySPOptions
 .def( py::init<>() )
 .def_readwrite( "BASIS",    &SP::Options::BASIS,    "Basis representation of sparse polynomial: 0-Monomials [Default], 1-Chebyshev" )
 .def_readwrite( "REMZERO",  &SP::Options::REMZERO,  "Whether to remove zeros entries from sparse polynomials [Default: True]" )
 .def_readwrite( "DISPLEN",  &SP::Options::DISPLEN,  "Number of digits in output stream [Default: 5]" )
 .def_readwrite( "DISPLINE", &SP::Options::DISPLINE, "Whether to display polynomial on single line [Default: True]" )
;

// --- FFPoly Main Class ---
pyFFPoly
 // Static Options
 .def_readwrite_static(
   "options",
   &SP::options,
   "Static options of sparse polynomials"
 )
 // Constructors
 .def(
   py::init<>(),
   "Default constructor"
 )
 .def(
   py::init<double const>(),
   py::arg("coef"),
   "Constructor for constant value"
 )
 .def(
   py::init<KEY const&, double const&>(),
   py::arg("var"),
   py::arg("coef") = 1.,
   "Constructor for variable 'var' and associated coefficient 'coef'"
 )
 .def(
   //py::init<std::pair<SM, double> const&>(),
   py::init( []( SM const& mon, double const& coef )
             { return new SP( { mon, coef } ); } ),
   py::arg("mon"),
   py::arg("coef") = 1.,
   "Constructor for monomial 'mon' and associated coefficient 'coef'"
 )
 .def(
   py::init<SP::t_poly const&>(),
   //py::init( []( std::vector< std::pair< SM, double > > const& coefmon )
   //          { return new SP( SP::t_poly( coefmon.cbegin(), coefmon.cend() ) ); } ),
   py::arg("coefmon"),
   "Constructor for dictionary {mon1: coef1, mon2: coef2, ...}"
 )
 .def(
   py::init<SP const&>(),
   "Copy constructor"
 )
 // Accessors
 .def(
   "maxord",
   &SP::maxord,
   "Maximal degree of monomial terms"
 )
 .def(
   "minord",
   py::overload_cast<>( &SP::minord, py::const_ ),
   "Minimal degree of monomial terms"
 )
 .def(
   "minord",
   py::overload_cast< KEY const& >( &SP::minord, py::const_ ),
   py::arg("x"),
   "Minimal degree of variable monomial terms"
 )
 .def_property_readonly(
   "nmon",
   &SP::nmon,
   "Total number of monomial terms"
 )
 .def( "coef",
   &SP::coef,
   py::arg("mon") = SM(),
 "Get coefficient of monomial" )
 .def(
   "nvar",
   &SP::nvar,
   "Total number of participating variables"
 )
 // Properties
 .def_property_readonly(
   "coefmon",
   py::overload_cast<>( &SP::mapmon, py::const_ ),
   "Dictionary of monomials and associated coefficients"
 )
 .def_property_readonly(
   "setvar",
   py::overload_cast<>( &SP::setvar, py::const_ ),
   "List of participating variables"
 )
 // Methods
 .def(
   "swap",
   &SP::swap,
   py::arg("spoly"),
   "Swap contents"
 )
 .def(
   "clean",
   &SP::clean,
   py::arg("tol") = 0.,
   "Remove entries below tolerance"
 )
 .def(
   "convert",
   &SP::convert,
   py::arg("basis"),
   "Convert to desired basis from options.BASIS"
 )
 .def(
   "ismultiple",
   &SP::ismultiple,
   py::arg("spoly"),
   "Check if polynomial is multiple of 'spoly'"
 )
 .def(
   "factor",
   &SP::factor,
   py::arg("x"),
   "Factor polynomial with respect to variable 'x'"
 )
 .def(
   "diff",
   &SP::diff,
   py::arg("x"),
   "Differentiate polynomial with respect to variable 'x'"
 )
 .def(
   "var",
   &SP::var,
   py::arg("x"),
   "Set equal to variable x"
 )
 .def_static(
   "display", 
   []( std::vector<std::pair<SM, double>> const& coefmon, int basis, int len, bool line )
   {
     SP::t_poly internal_map;
     for( auto const& [k,v] : coefmon ) internal_map.insert({k,v});
       return SP::display( internal_map, basis, len, line );
   },
   py::arg("coefmon"),
   py::arg("BASIS")   = SP::options.BASIS, 
   py::arg("DISPLEN") = SP::options.DISPLEN,
   py::arg("ONELINE") = SP::options.DISPLINE,
   "Display polynomial expression" 
 )
 .def(
   "eval",
   []( SP const& self, std::map<KEY,double,COMP> const& val )
   { return self.eval( val ); },
   py::arg("x"),
   "Evaluate polynomial (double arithmetic)" 
 )
 // Operators
 .def( py::self += py::self )
 .def( py::self += double() )
 .def( py::self += std::pair<SM,double>() )
 //.def( "add_mon", []( SP& self, SM const& m, double c ){ self += std::make_pair(m, c); return self; } )

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

 .def( "__eq__",  []( SP const& x, SP const& y ){ return x == y; },       py::is_operator())
 .def( "__ne__",  []( SP const& x, SP const& y ){ return !(x == y); },    py::is_operator())
 .def( "__pow__", []( SP const& s, unsigned n ){ return mc::pow(s, n); }, py::is_operator() )

 // String Representation
 .def( "__str__",  []( SP const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
 .def( "__repr__", []( SP const& self ){ std::ostringstream oss; oss << self; return oss.str(); } )
;

m.def( "sqr",  []( SP const& s ){ return mc::sqr(s); } );
m.def( "pow",  []( SP const& s, unsigned n ){ return mc::pow(s, n); } );
m.def( "cheb", []( SP const& s, unsigned n ){ return mc::cheb(s, n); } );
m.def( "prod", []( unsigned n, std::vector<SP> const& v ){ return mc::prod(n, v.data()); } );

}
