#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "smon.hpp"
#include "spoly.hpp"
#include "ffunc.hpp" 

namespace py = pybind11;

void mc_smon( py::module_ &m )
{

typedef unsigned              KEY;
typedef std::less<unsigned>   COMP;
typedef mc::SMon<KEY, COMP>   SMon;
typedef mc::SPoly<KEY, COMP>  SPoly;

py::class_<SMon> pySMon( m, "SMon" );

// --- Exceptions Nested Class ---
py::class_<SMon::Exceptions> pySMonExceptions( pySMon, "Exceptions" );
py::enum_<SMon::Exceptions::TYPE>( pySMonExceptions, "TYPE" )
 .value( "SUB",  SMon::Exceptions::SUB,  "Subtraction of a monomial that is not a proper subset" )
 .value( "DIV",  SMon::Exceptions::DIV,  "Division by a factor greater than the greatest common exponent" )
 .value( "CONV", SMon::Exceptions::CONV, "Conversion with different variable indexing failed" )
 .export_values();

pySMonExceptions
 .def( "ierr", &SMon::Exceptions::ierr, "Return error flag" )
 .def( "what", &SMon::Exceptions::what, "Return error description" )
;

// --- SMon Main Class ---
pySMon
 // Constructors
 .def( 
   py::init<>(),
   "Default constructor (constant 1)"
 )
 .def(
   py::init< KEY const&, unsigned const >(),
   py::arg("var"),
   py::arg("ord") = 1,
   "Constructor from variable and exponent" 
 )
 .def(
   py::init< std::map<KEY, unsigned, COMP> const& >(),
   py::arg("expr"),
   "Constructor from dictionary of variables and exponents" 
 )
 .def(
   py::init< SMon const& >(),
   "Copy constructor"
 )
// Accessors
 .def_readonly(
   "tord",
   &SMon::tord,
   "Monomial total order"
 )
 .def_property_readonly(
   "expr",
   //&SMon::expr,
   []( SMon& self )
   {
     auto const& expr = self.expr;
     return std::vector< std::pair< KEY, unsigned > >( expr.cbegin(), expr.cend() );
   },
   "List of variables and exponents tuples" 
 )
 // Methods
 .def( "exp", 
   &SMon::exp,
   py::arg("var"),
   "Exponent of variable 'var' in monomial" 
 )
 .def(
   "gcexp",
   &SMon::gcexp,
   "Greatest common exponent of terms"
 )
 .def(
   "lexp",
   &SMon::lexp,
   "Least exponent among all terms"
 )
 .def(
   "gexp",
   &SMon::gexp,
   "Greatest exponent among all terms"
 )
 .def(
   "hull",
   &SMon::hull,
   py::arg("mon"),
   "Union with other monomial"
 )
 .def(
   "inter",
   &SMon::inter,
   py::arg("mon"),
   "Test for intersection"
 )
 .def(
   "subset",
   &SMon::subset,
   py::arg("mon"),
   "Test for proper subset"
 )
 .def(
   "subseteq",
   &SMon::subseteq,
   py::arg("mon"),
   "Test for subset"
 )
 .def(
   "display",
   &SMon::display,
   py::arg("basis") = 0,
   "String representation with basis option: 0-Monomial, 1-Chebyshev"
 )
 // Operators
 .def(
   "__hash__",
   []( SMon const& self )
   {
     // Hash is computed based on 'tord' and the 'expr' map content.
     // We convert the map to a tuple of items (which is sorted by key in C++ map)
     // to ensure consistent hashing in Python.
     py::list expr_items;
     for( auto const& [var, ord] : self.expr ){
       expr_items.append( py::make_tuple(var, ord) );
     }
     // Return hash((tord, tuple(items)))
     return py::hash( py::make_tuple( self.tord, py::tuple(expr_items) ) );
   },
   "Compute hash of the monomial"
 )
 .def(
   "__getitem__", 
   []( SMon const& self, KEY const& var ){ return self[ var ]; },
   "Extract sub-monomial for variable" 
 )
 .def( py::self += py::self )
 .def( py::self -= py::self )
 .def( py::self +  py::self )
 .def( py::self -  py::self )
 .def( py::self == py::self )
 .def( py::self != py::self )
 .def( py::self *= unsigned() )
 .def( py::self /= unsigned() )
 .def( py::self *  unsigned() )
 .def( py::self /  unsigned() )

 // Comparison operators
 .def( "__eq__",  []( SMon const& Mon1, SMon const& Mon2 ){ return   Mon1 == Mon2;  }, py::is_operator() )
 .def( "__ne__",  []( SMon const& Mon1, SMon const& Mon2 ){ return !(Mon1 == Mon2); }, py::is_operator() )

 // String representation
 .def( "__str__",  []( SMon const& self ){ return self.display(SPoly::options.BASIS); } )
 .def( "__repr__", []( SMon const& self ){ return self.display(SPoly::options.BASIS); } )
;

}
