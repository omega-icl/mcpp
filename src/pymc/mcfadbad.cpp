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

#include "mcfadbad.hpp"
#include "fadbad.h"
#include "fadiff.h"
#include "badiff.h"
#include "tadiff.h"
typedef fadbad::FTypeName<MC> FMC;
typedef fadbad::BTypeName<MC> BMC;
typedef fadbad::TTypeName<MC> TMC;


namespace py = pybind11;

void mc_fadbad_forward(py::module &m) 
{   
    py::class_<FMC> pyFadbadF(m,"FadbadF"); 
pyFadbadF
 .def(
    py::init<MC const&>(),
    "construct forward automatic differentiation"
 )
 .def(
    py::init<const double&>(),
    "construct forward automatic differentiation"
 )  
 .def(
    py::init<const FMC&>(),
    "copy constructor"
 )  
 .def_property_readonly(
    "value",
    []( FMC const& self ){ return self.val(); },
    "value"
 )
 .def_property_readonly(
    "size",
    []( FMC const& self ){ return self.size(); },
    "size"
 )
 .def(
    "__getitem__",
    []( FMC const& self, unsigned int i ){ return self.operator[](i); },
    "get ith derivative"
 )
 .def(
    "variable",
    []( FMC& self ){ return self.x(); },
    "return variable"
 )
 .def(
    "deriv",
    []( FMC& self, const unsigned int i ){ return self.deriv(i); },
    "get ith derivative"
 )  
 .def(
    "diff",
    []( FMC& self, unsigned int i, unsigned int N ){ return self.diff(i,N); },
    "calculate derivative"
 )
 .def(
    "depend",
    []( FMC const& self ){ return self.depend(); },
    "check if dependent"
 )
 .def(
    "setDepend",
    []( FMC& self, FMC const& val ){ self.setDepend(val); },
    "set dependent"
)
.def(
    "setDepend",
    []( FMC& self, FMC const& val1, FMC const& val2 ){ self.setDepend(val1,val2); },
    "set dependent"
)

 .def( + py::self )
 .def( py::self + double() )  
 .def( double() + py::self )
 .def( py::self + py::self )
 .def( py::self += double() )
 .def( py::self += py::self )   
 .def( - py::self ) 
 .def( py::self - double() )    
 .def( double() - py::self )
 .def( py::self - py::self )
 .def( py::self -= double() )
 .def( py::self -= py::self )
 .def( py::self * double() )    
 .def( double() * py::self )
 .def( py::self * py::self )
 .def( py::self *= double() )
 .def( py::self *= py::self )
 .def( py::self / double() )    
 .def( double() / py::self )
 .def( py::self / py::self )   
 .def( py::self /= double() )
 .def( py::self /= py::self )
 .def( py::self == py::self )
 .def( py::self != py::self )
 .def( py::self <  py::self )
 .def( py::self <= py::self )
 .def( py::self >  py::self )
 .def( py::self >= py::self )
 .def( py::self == double() )
 .def( double() == py::self )    
 .def( py::self != double() )
 .def( double() != py::self )
 .def( py::self <  double() )
 .def( double() <  py::self )
 .def( py::self <= double() )
 .def( double() <= py::self )    
 .def( py::self >  double() )
 .def( double() >  py::self )
 .def( py::self >= double() )
 .def( double() >= py::self )
 .def( "__pow__", []( double const& a, FMC const& m ){ return fadbad::pow(a,m); } )
 .def( "__pow__", []( int const& a, FMC const& m ){ return fadbad::pow(a,m); } )
 .def( "__pow__", []( FMC const& m, double const& a ){ return fadbad::pow(m,a); } )  
 .def( "__pow__", []( FMC const& m, int const& a ){ return fadbad::pow(m,a); } )   
 .def( "__pow__", []( FMC const& m1, FMC const& m2 ){ return fadbad::pow(m1,m2); } )   
 .def( "__pos__", []( FMC const& m ){ return operator+(m); } )
 .def( "__neg__", []( FMC const& m ){ return operator-(m); } )
;

m.def( "sqr",   []( FMC const& m ){ return fadbad::sqr(m); } );
m.def( "exp",   []( FMC const& m ){ return fadbad::exp(m); } );
m.def( "log",   []( FMC const& m ){ return fadbad::log(m); } );
m.def( "sqrt",  []( FMC const& m ){ return fadbad::sqrt(m); } );
m.def( "sin",   []( FMC const& m ){ return fadbad::sin(m); } );
m.def( "cos",   []( FMC const& m ){ return fadbad::cos(m); } );
m.def( "tan",   []( FMC const& m ){ return fadbad::tan(m); } );
// m.def( "cot",   []( FMC const& m ){ return fadbad::cot(m); } );    
m.def( "asin",  []( FMC const& m ){ return fadbad::asin(m); } );
m.def( "acos",  []( FMC const& m ){ return fadbad::acos(m); } );
m.def( "atan",  []( FMC const& m ){ return fadbad::atan(m); } );
m.def( "sinh",  []( FMC const& m ){ return fadbad::sinh(m); } );
m.def( "cosh",  []( FMC const& m ){ return fadbad::cosh(m); } );
m.def( "tanh",  []( FMC const& m ){ return fadbad::tanh(m); } );
// m.def( "coth",  []( FMC const& m ){ return fadbad::coth(m); } );
// cot and coth are not defined for the FMC
}
