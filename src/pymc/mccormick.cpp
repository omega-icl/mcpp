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

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

namespace py = pybind11;

void mc_mccormick(py::module &m) 
{

py::class_<MC> pyMcCormick(m,"McCormick");
pyMcCormick
 .def(
   py::init<>()
 )
 .def(
   py::init<double const&>(),
   "constructor for constant value"
 )
 .def(
   py::init<I const&>(),
   "constructor for constant bounds"
 )
 .def(
   py::init<I const&, double const&>(),
   "constructor for range and equal convex and concave bounds"
 )
 .def(
   py::init<I const&, double const&, double const&>(),
   "constructor for range, convex and concave bounds"
 )
 .def(
   py::init<MC const&>(),
   "copy constructor"
 )
 .def_readwrite_static(
   "options",
   &MC::options
 )
 .def_property_readonly(
   "l",
   []( MC& self ){ return self.l(); },
   "variable lower range"
 )
 .def_property_readonly(
   "u",
   []( MC& self ){ return self.u(); },
   "variable upper range"
 )
 .def_property(
   "I",
   []( MC& self ){ return self.I(); },
   []( MC& self, I const& I ){ self.I() = I; },
   "variable range"
 )
 .def_property(
   "cv",
   []( MC& self ){ return self.cv(); },
   []( MC& self, double const& cv ){ self.cv( cv ); },
   "convex bound"
 )
 .def_property(
   "cc", 
   []( MC& self ){ return self.cc(); },
   []( MC& self, double const& cc ){ self.cc( cc ); },
   "concave bound"
 )
 .def(
   "c",
   []( MC& self, double const& c ){ return self.c( c ); },
   "set convex and concave bounds"
 )
 .def_property_readonly(
   "nsub",
   []( MC& self ){ return self.nsub(); },
   "subgradient size"
 )
 .def_property_readonly(
   "cvsub",
   []( MC& self )
     { return std::vector<double>( self.cvsub(), self.cvsub()+self.nsub() ); },
   py::return_value_policy::take_ownership,
   "convex bound subgradient"
 )
 .def_property_readonly(
   "ccsub",
   []( MC& self )
     { return std::vector<double>( self.ccsub(), self.ccsub()+self.nsub() ); },
   py::return_value_policy::take_ownership,
   "concave bound subgradient"
 )
 .def(
   "sub",
   static_cast< MC& (MC::*)( unsigned int const )>( &MC::sub ),
   py::arg( "size" ),
   "set subgradient size"
 )
 .def(
   "sub",
   static_cast< MC& (MC::*)( unsigned int const, unsigned int const )>( &MC::sub),
   py::arg("size"),
   py::arg("index"),
   "set subgradient size and values for variable index"
 )
 .def(
   "sub",
   []( MC& self, std::vector<double> const& cvsub, std::vector<double> const& ccsub )
     {
       if( cvsub.size() != ccsub.size() )
         throw( std::runtime_error( "Inconsistent subgradient sizes" ) );
       return self.sub( cvsub.size(), cvsub.data(), ccsub.data() );
     },
   py::return_value_policy::take_ownership,
   py::arg("cvsub"),
   py::arg("ccsub"), 
   "set subgradient values"
 )
 .def(
   "cut",
   static_cast<MC& (MC::*)()>( &MC::cut)
 )
 .def(
   "laff",
   []( MC& self, std::vector<double> const& val, std::vector<double> const& ref )
     {
       if( val.size() != self.nsub() || ref.size() != self.nsub() )
         throw( std::runtime_error( "Inconsistent sizes" ) );
       return self.laff( val.data(), ref.data() );
     },
   py::arg("val"),
   py::arg("ref"), 
   "affine underestimator value"
 )
 .def(
   "laff",
   []( MC& self, std::vector<I> const& rng, std::vector<double> const& ref )
     {
       if( rng.size() != self.nsub() || ref.size() != self.nsub() )
         throw( std::runtime_error( "Inconsistent sizes" ) );
       return self.laff( rng.data(), ref.data() );
     },
   py::arg("rng"),
   py::arg("ref"), 
   "affine underestimator bound"
 )
 .def(
   "uaff",
   []( MC& self, std::vector<double> const& val, std::vector<double> const& ref )
     {
       if( val.size() != self.nsub() || ref.size() != self.nsub() )
         throw( std::runtime_error( "Inconsistent sizes" ) );
       return self.uaff( val.data(), ref.data() );
     },
   py::arg("val"),
   py::arg("ref"), 
   "affine overestimator value"
 )
 .def(
   "uaff",
   []( MC& self, std::vector<I> const& rng, std::vector<double> const& ref )
     {
       if( rng.size() != self.nsub() || ref.size() != self.nsub() )
         throw( std::runtime_error( "Inconsistent sizes" ) );
       return self.uaff( rng.data(), ref.data() );
     },
   py::arg("rng"),
   py::arg("ref"), 
   "affine overestimator bound"
 )
 //.def( "copy", []( MC const& self ){ return MC( self ); } )
 //.def( "assign", py::overload_cast<double const&>( &MC::operator= ), py::arg("Value") )
 //.def( "assign", py::overload_cast<I const&>( &MC::operator= ), py::arg("Value") )
 //.def( "assign", py::overload_cast<MC const&>( &MC::operator= ), py::arg("Value") )
 .def(
   "__str__",
   []( MC const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( MC const& self ){ std::ostringstream os; os << self; return os.str(); }
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
 .def( "__abs__", []( MC const& m ){ return mc::Op<MC>::abs(m); } )
 .def( "__pow__", []( MC const& m, int const n ){ return mc::Op<MC>::pow(m,n); } )
 .def( "__pow__", []( MC const& m, double const& r ){ return mc::Op<MC>::pow(m,r); } )
 .def( "__pow__", []( MC const& m, MC const& mm ){ return mc::Op<MC>::pow(m,mm); } )
 .def( "__pow__", []( double const& r, MC const& m ){ return mc::Op<MC>::pow(r,m); } )
 .def( py::self == py::self ) 
 .def( py::self != py::self )
 .def( py::self <= py::self )
 .def( py::self >= py::self )
 .def( py::self < py::self )
 .def( py::self > py::self )
; 

py::class_<MC::Options> pyMcCormickOptions( pyMcCormick, "Options");
pyMcCormickOptions
 .def( py::init<>() )
 .def( py::init<MC::Options const&>() )
 .def_readwrite( "ENVEL_USE",      &MC::Options::ENVEL_USE )
 .def_readwrite( "ENVEL_MAXIT",    &MC::Options::ENVEL_MAXIT )
 .def_readwrite( "ENVEL_TOL",      &MC::Options::ENVEL_TOL )
 .def_readwrite( "MVCOMP_USE",     &MC::Options::MVCOMP_USE )
 .def_readwrite( "MVCOMP_TOL",     &MC::Options::MVCOMP_TOL )
 .def_readwrite( "DISPLAY_DIGITS", &MC::Options::DISPLAY_DIGITS )
;

m.def( "cut",    []( MC const& x ){ return mc::cut(x); } );
m.def( "inv",    []( MC const& x ){ return mc::inv(x); } );
m.def( "sqr",    []( MC const& x ){ return mc::sqr(x); } );
m.def( "sqrt",   []( MC const& x ){ return mc::sqrt(x); } );
m.def( "exp",    []( MC const& x ){ return mc::exp(x); } );
m.def( "log",    []( MC const& x ){ return mc::log(x); } );
m.def( "cos",    []( MC const& x ){ return mc::cos(x); } );
m.def( "sin",    []( MC const& x ){ return mc::sin(x); } );
m.def( "tan",    []( MC const& x ){ return mc::tan(x); } );
m.def( "acos",   []( MC const& x ){ return mc::acos(x); } );
m.def( "asin",   []( MC const& x ){ return mc::asin(x); } );
m.def( "atan",   []( MC const& x ){ return mc::atan(x); } );
m.def( "cosh",   []( MC const& x ){ return mc::cosh(x); } );
m.def( "sinh",   []( MC const& x ){ return mc::sinh(x); } );
m.def( "tanh",   []( MC const& x ){ return mc::tanh(x); } );
m.def( "fabs",   []( MC const& x ){ return mc::fabs(x); } );
m.def( "relu",   []( MC const& x ){ return mc::max(x,0.); } );
m.def( "xlog",   []( MC const& x ){ return mc::xlog(x); } );
m.def( "fstep",  []( MC const& x ){ return mc::fstep(x); } );
m.def( "bstep",  []( MC const& x ){ return mc::bstep(x); } );
m.def( "erf",    []( MC const& x ){ return mc::erf(x); } );
m.def( "erfc",   []( MC const& x ){ return mc::erfc(x); } );
m.def( "lmtd",   []( MC const& x, MC const& y ){ return mc::lmtd(x,y); } );
m.def( "rlmtd",  []( MC const& x, MC const& y ){ return mc::rlmtd(x,y); } );
m.def( "pow",    []( MC const& x, int const n ){ return mc::pow(x,n); } );
m.def( "pow",    []( MC const& x, double const& r ){ return mc::pow(x,r); } );
m.def( "pow",    []( MC const& x, MC const& y ){ return mc::pow(x,y); } );
m.def( "pow",    []( double const& r, MC const& y ){ return mc::pow(r,y); } );
m.def( "cheb",   []( MC const& x, unsigned const n ){ return mc::cheb(x,n); } );
m.def( "hull",   []( MC const& x, MC const& y ){ return mc::hull(x,y); } );
m.def( "max",    []( MC const& x, MC const& y ){ return mc::min(x,y); } );
m.def( "min",    []( MC const& x, MC const& y ){ return mc::max(x,y); } );
m.def( "inter",  []( MC& z, MC const& x, MC const& y ){ return mc::inter(z,x,y); } );
m.def( "ltcond", []( I const& z, MC const& x, MC const& y ){ return mc::ltcond(z,x,y); } );
m.def( "ltcond", []( MC const& z, MC const& x, MC const& y ){ return mc::ltcond(z,x,y); } );
m.def( "gtcond", []( I const& z, MC const& x, MC const& y ){ return mc::ltcond(z,x,y); } );
m.def( "gtcond", []( MC const& z, MC const& x, MC const& y ){ return mc::ltcond(z,x,y); } );
}
