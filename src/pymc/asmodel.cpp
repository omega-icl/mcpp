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

#include "asmodel.hpp"
typedef mc::ASModel<I> ASM;
typedef mc::ASVar<I> ASV;

namespace py = pybind11;

void mc_asmodel(py::module &m) 
{

py::class_<ASM> pyASModel(m,"ASModel");
pyASModel
 .def(
   py::init<unsigned int const&, unsigned int const&>(),
   py::arg("nvar"),
   py::arg("ndiv")=1, 
   "constructor"
 )
 .def_readwrite_static(
   "options",
   &ASM::options
 )
 .def_property_readonly(
   "nvar",
   []( ASM& self ){ return self.nvar(); },
   "number of independent variables"
 )
 .def_property_readonly(
   "subvar",
   []( ASM& self ){ return self.psize(); },
   "subpartition of independent variables"
 )
 .def_property_readonly(
   "bndvar",
   []( ASM& self ){ return self.bndvar(); },
   "domain of independent variables"
 )
;

py::class_<ASV> pyASVar(m,"ASVar");
pyASVar
 .def(
   py::init<double const&>(),
   py::arg("val")=0., 
   "constructor for constant"
 )
 .def(
   py::init<ASM* const, unsigned int, I const&>(),
   "constructor for variable"
 )
 .def(
   py::init<ASV const&>(),
   "copy constructor"
 )
 //.def(
 //  py::init<ASV && var>(),
 //  "move constructor"
 //)
 .def(
   "set",
   []( ASV& self, ASM* const mod ){ return self.set( mod ); },
   "set model environment"
 )
 .def(
   "set",
   []( ASV& self, ASM* const mod, unsigned int ndx, I const& bnd ){ return self.set( mod, ndx, bnd ); },
   "set independent variable"
 )
 .def_property_readonly(
   "info",
   []( ASV& self ){ return self.get_ASVar(); },
   "information lfag: 1 - constant, 2 - linear, 3 - single superposition, 4 - shadowed superposition"
 )
 .def_property_readonly(
   "l",
   []( ASV& self ){ return self.lb(); },
   "lower bound of superposition"
 )
 .def_property_readonly(
   "u",
   []( ASV& self ){ return self.ub(); },
   "upper bound of superposition"
 )
 .def_property_readonly(
   "B",
   []( ASV& self ){ return self.bound(); },
   "bound of superposition"
 )
 .def(
   "R",
   []( ASV& self, std::vector<double> const& x )
     { if( self.mod() && x.size() != self.mod()->nvar() )
         throw( std::runtime_error( "Incorrect size" ) );
       return self.eval( x.data() );
     },
   "range of superposition at given point"
 )
 //.def( "copy", []( ASV const& self ){ return ASV( self ); } )
 //.def( "assign", py::overload_cast<double const&>( &ASV::operator= ), py::arg("val") )
 //.def( "assign", py::overload_cast<ASV const&>( &ASV::operator= ), py::arg("var") )
 .def(
   "__str__",
   []( ASV const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( ASV const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def( + py::self )
 .def( py::self += double() )
 .def( py::self += py::self )
 .def( double() + py::self )
 .def( py::self + double() )
 .def( py::self + py::self )
 .def( - py::self )
 .def( py::self -= double() )
 .def( py::self -= py::self )
 .def( double() - py::self )
 .def( py::self - double() )
 .def( py::self - py::self )
 .def( py::self *= double() )
 .def( py::self *= py::self )
 .def( double() * py::self )
 .def( py::self * double() )
 .def( py::self * py::self )
 .def( py::self /= double() )
 .def( py::self /= py::self )
 .def( double() / py::self )
 .def( py::self / double() )
 .def( py::self / py::self )
// .def( "__abs__", []( ASV const& m ){ return mc::Op<ASV>::abs(m); } )
// .def( "__pow__", []( ASV const& m, int const n ){ return mc::Op<ASV>::pow(m,n); } )
// .def( "__pow__", []( ASV const& m, double const& r ){ return mc::Op<ASV>::pow(m,r); } )
// .def( "__pow__", []( ASV const& m, ASV const& mm ){ return mc::Op<ASV>::pow(m,mm); } )
// .def( "__pow__", []( double const& r, ASV const& m ){ return mc::Op<ASV>::pow(r,m); } )
// .def( py::self == py::self ) 
// .def( py::self != py::self )
// .def( py::self <= py::self )
// .def( py::self >= py::self )
// .def( py::self < py::self )
// .def( py::self > py::self )
; 

py::class_<ASM::Options> pyASModelOptions( pyASModel, "Options");
pyASModelOptions
 .def( py::init<>() )
 .def_readwrite( "NSUB",      &ASM::Options::NSUB )
 .def_readwrite( "SHADOW_USE",  &ASM::Options::SHADOW_USE )
;

//m.def( "inv",    []( ASV const& x ){ return mc::inv(x); } );
//m.def( "sqr",    []( ASV const& x ){ return mc::sqr(x); } );
//m.def( "sqrt",   []( ASV const& x ){ return mc::sqrt(x); } );
//m.def( "exp",    []( ASV const& x ){ return mc::exp(x); } );
//m.def( "log",    []( ASV const& x ){ return mc::log(x); } );
//m.def( "cos",    []( ASV const& x ){ return mc::cos(x); } );
//m.def( "sin",    []( ASV const& x ){ return mc::sin(x); } );
//m.def( "tan",    []( ASV const& x ){ return mc::tan(x); } );
//m.def( "acos",   []( ASV const& x ){ return mc::acos(x); } );
//m.def( "asin",   []( ASV const& x ){ return mc::asin(x); } );
//m.def( "atan",   []( ASV const& x ){ return mc::atan(x); } );
//m.def( "cosh",   []( ASV const& x ){ return mc::cosh(x); } );
//m.def( "sinh",   []( ASV const& x ){ return mc::sinh(x); } );
//m.def( "tanh",   []( ASV const& x ){ return mc::tanh(x); } );
//m.def( "fabs",   []( ASV const& x ){ return mc::fabs(x); } );
m.def( "relu",   []( ASV const& x ){ return mc::relu(x); } );
//m.def( "xlog",   []( ASV const& x ){ return mc::xlog(x); } );
//m.def( "fstep",  []( ASV const& x ){ return mc::fstep(x); } );
//m.def( "bstep",  []( ASV const& x ){ return mc::bstep(x); } );
//m.def( "erf",    []( ASV const& x ){ return mc::erf(x); } );
//m.def( "erfc",   []( ASV const& x ){ return mc::erfc(x); } );
//m.def( "lmtd",   []( ASV const& x, ASV const& y ){ return mc::lmtd(x,y); } );
//m.def( "rlmtd",  []( ASV const& x, ASV const& y ){ return mc::rlmtd(x,y); } );
//m.def( "pow",    []( ASV const& x, int const n ){ return mc::pow(x,n); } );
//m.def( "pow",    []( ASV const& x, double const& r ){ return mc::pow(x,r); } );
//m.def( "pow",    []( ASV const& x, ASV const& y ){ return mc::pow(x,y); } );
//m.def( "pow",    []( double const& r, ASV const& y ){ return mc::pow(r,y); } );
//m.def( "cheb",   []( ASV const& x, unsigned const n ){ return mc::cheb(x,n); } );
//m.def( "hull",   []( ASV const& x, ASV const& y ){ return mc::hull(x,y); } );
m.def( "max",    []( ASV const& x, double const& y ){ return mc::max(x,y); } );
//m.def( "min",    []( ASV const& x, double const& y ){ return mc::min(x,y); } );
m.def( "max",    []( ASV const& x, ASV const& y ){ return mc::max(x,y); } );
//m.def( "min",    []( ASV const& x, ASV const& y ){ return mc::min(x,y); } );
//m.def( "inter",  []( ASV& z, ASV const& x, ASV const& y ){ return mc::inter(z,x,y); } );
//m.def( "ltcond", []( I const& z, ASV const& x, ASV const& y ){ return mc::ltcond(z,x,y); } );
//m.def( "ltcond", []( ASV const& z, ASV const& x, ASV const& y ){ return mc::ltcond(z,x,y); } );
//m.def( "gtcond", []( I const& z, ASV const& x, ASV const& y ){ return mc::ltcond(z,x,y); } );
//m.def( "gtcond", []( ASV const& z, ASV const& x, ASV const& y ){ return mc::ltcond(z,x,y); } );
}
