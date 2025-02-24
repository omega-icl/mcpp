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


void mc_fadbad_backward(py::module &m) 
{   
    py::class_<BMC> pyFadbadB(m,"FadbadB");
pyFadbadB
 .def(
    py::init<MC const&>(),
    "construct backward automatic differentiation"
 )
 .def(
    py::init<const double&>(),
    "construct backward automatic differentiation"
 )
 .def(
    py::init<const BMC&>(),
    "copy constructor"
 )
 .def_property_readonly(
    "value",
    []( BMC const& self ){ return self.val(); },
    "value"
 ) 
 .def(
   "variable",
   []( BMC& self ){ return self.x(); },
   "return variable"
 )
 .def(
    "deriv",
    []( BMC& self, const unsigned int i ){ return self.deriv(i); },
    "get ith derivative"
 ) 
 .def(
    "diff",
    []( BMC self, const unsigned int i, const unsigned int N ){ return self.diff(i,N); },
    "calculate derivative"
 )

   .def( py::self += double() )
   .def( py::self += py::self )
   .def( py::self -= double() )
   .def( py::self -= py::self )
   .def( py::self *= double() )
   .def( py::self *= py::self )
   .def( py::self /= double() )
   .def( py::self /= py::self )

   .def( "__eq__", []( BMC const& m1, BMC const& m2 ){ return m1 == m2; } )
   .def( "__eq__", []( BMC const& m1, double const& m2 ){ return m1 == m2; } )
   .def( "__eq__", []( double const& m1, BMC const& m2 ){ return m1 == m2; } )
   .def( "__ne__", []( BMC const& m1, BMC const& m2 ){ return m1 != m2; } )
   .def( "__ne__", []( BMC const& m1, double const& m2 ){ return m1 != m2; } )
   .def( "__ne__", []( double const& m1, BMC const& m2 ){ return m1 != m2; } )
   .def( "__lt__", []( BMC const& m1, BMC const& m2 ){ return m1 < m2; } )
   .def( "__lt__", []( BMC const& m1, double const& m2 ){ return m1 < m2; } )
   .def( "__lt__", []( double const& m1, BMC const& m2 ){ return m1 < m2; } )
   .def( "__le__", []( BMC const& m1, BMC const& m2 ){ return m1 <= m2; } )
   .def( "__le__", []( BMC const& m1, double const& m2 ){ return m1 <= m2; } )
   .def( "__le__", []( double const& m1, BMC const& m2 ){ return m1 <= m2; } )
   .def( "__gt__", []( BMC const& m1, BMC const& m2 ){ return m1 > m2; } )
   .def( "__gt__", []( BMC const& m1, double const& m2 ){ return m1 > m2; } )
   .def( "__gt__", []( double const& m1, BMC const& m2 ){ return m1 > m2; } )
   .def( "__ge__", []( BMC const& m1, BMC const& m2 ){ return m1 >= m2; } )
   .def( "__ge__", []( BMC const& m1, double const& m2 ){ return m1 >= m2; } )
   .def( "__ge__", []( double const& m1, BMC const& m2 ){ return m1 >= m2; } )
   .def( "__pos__", []( BMC const& m ){ return operator+(m); } )
   .def( "__neg__", []( BMC const& m ){ return operator-(m); } )
   .def( "__pow__", []( BMC const& m1, BMC const& m2 ){ return fadbad::pow(m1,m2); } )
   .def( "__pow__", []( BMC const& m1, double const& m2 ){ return fadbad::pow(m1,m2); } )
   .def( "__pow__", []( double const& m1, BMC const& m2 ){ return fadbad::pow(m1,m2); } )
   .def( "__pow__", []( BMC const& m, int const& n ){ return fadbad::pow(m,n); } )
   .def( "__pow__", []( int const& n, BMC const& m ){ return fadbad::pow(n,m); } )

   .def( py::self + double() )   
   .def( double() + py::self )
   .def( py::self + py::self )
   .def( py::self - double() )
   .def( double() - py::self )
   .def( py::self - py::self )
   .def( py::self * double() )
   .def( double() * py::self )
   .def( py::self * py::self )
   .def( py::self / double() )
   .def( double() / py::self )
   .def( py::self / py::self )
   ;

m.def( "sqr",   []( BMC const& m ){ return fadbad::sqr(m); } );
m.def( "sqrt",  []( BMC const& m ){ return fadbad::sqrt(m); } );
m.def( "exp",   []( BMC const& m ){ return fadbad::exp(m); } );
m.def( "log",   []( BMC const& m ){ return fadbad::log(m); } );
m.def( "sin",   []( BMC const& m ){ return fadbad::sin(m); } );
m.def( "cos",   []( BMC const& m ){ return fadbad::cos(m); } );
m.def( "tan",   []( BMC const& m ){ return fadbad::tan(m); } );
// m.def( "cot",   []( BMC const& m ){ return fadbad::cot(m); } );
m.def( "asin",  []( BMC const& m ){ return fadbad::asin(m); } );
m.def( "acos",  []( BMC const& m ){ return fadbad::acos(m); } );
m.def( "atan",  []( BMC const& m ){ return fadbad::atan(m); } );
m.def( "sinh",  []( BMC const& m ){ return fadbad::sinh(m); } );
m.def( "cosh",  []( BMC const& m ){ return fadbad::cosh(m); } );
m.def( "tanh",  []( BMC const& m ){ return fadbad::tanh(m); } );
// m.def( "coth",  []( BMC const& m ){ return fadbad::coth(m); } );

;
}

void mc_fadbad_taylor(py::module &m) 
{   
    py::class_<TMC> pyFadbadT(m,"FadbadT");
pyFadbadT
 .def(
    py::init<MC const&>(),
    "construct taylor automatic differentiation"
 )
 .def(
    py::init<const double&>(),
    "construct taylor automatic differentiation"
 ) 
   .def(
      py::init<const TMC&>(),
      "copy constructor"
   )
 .def_property_readonly(
    "value",
    []( TMC const& self ){ return self.val(); },
    "value"
 ) 
.def_property_readonly(
   "size",
   []( TMC const& self ){ return self.length(); },
   "size" 
 )
.def(
   "__getitem__",
   []( TMC const& self, unsigned int i ){ return self[i]; },
   "get ith derivative"
 )
.def(
   "__setitem__",
   []( TMC& self, unsigned int i, double const& val ){ self[i] = val; },
   "set ith derivative"
)
 .def(
   "reset",
   []( TMC& self ){ return self.reset(); },
   "reset value"
 )
 .def(
   "eval",
   []( TMC& self, const unsigned int i ){ return self.eval(i); },
   "Taylor-expand f to degree i"
 )
//  .def(
//    "diff",
//    []( const TMC& self, const int i ){ return fadbad::diff(self,i); },
//    "differentiate f to order i"
//  )

 .def( py::self += double() )
 .def( py::self += py::self )
 .def( py::self -= double() )
 .def( py::self -= py::self )
 .def( py::self *= double() )
 .def( py::self *= py::self )
 .def( py::self /= double() )
 .def( py::self /= py::self )
 
 .def( "__eq__", []( TMC const& m1, TMC const& m2 ){ return m1 == m2; } )
 .def( "__eq__", []( TMC const& m1, double const& m2 ){ return m1 == m2; } )
 .def( "__eq__", []( double const& m1, TMC const& m2 ){ return m1 == m2; } )
 .def( "__ne__", []( TMC const& m1, TMC const& m2 ){ return m1 != m2; } )
 .def( "__ne__", []( TMC const& m1, double const& m2 ){ return m1 != m2; } )
 .def( "__ne__", []( double const& m1, TMC const& m2 ){ return m1 != m2; } )
 .def( "__lt__", []( TMC const& m1, TMC const& m2 ){ return m1 < m2; } )
 .def( "__lt__", []( TMC const& m1, double const& m2 ){ return m1 < m2; } )
 .def( "__lt__", []( double const& m1, TMC const& m2 ){ return m1 < m2; } )
 .def( "__le__", []( TMC const& m1, TMC const& m2 ){ return m1 <= m2; } )  
 .def( "__le__", []( TMC const& m1, double const& m2 ){ return m1 <= m2; } )
 .def( "__le__", []( double const& m1, TMC const& m2 ){ return m1 <= m2; } )
 .def( "__gt__", []( TMC const& m1, TMC const& m2 ){ return m1 > m2; } )
 .def( "__gt__", []( TMC const& m1, double const& m2 ){ return m1 > m2; } )
 .def( "__gt__", []( double const& m1, TMC const& m2 ){ return m1 > m2; } )
 .def( "__ge__", []( TMC const& m1, TMC const& m2 ){ return m1 >= m2; } )
 .def( "__ge__", []( TMC const& m1, double const& m2 ){ return m1 >= m2; } )
 .def( "__ge__", []( double const& m1, TMC const& m2 ){ return m1 >= m2; } )
 .def( "__pos__", []( TMC const& m ){ return operator+(m); } )
 .def( "__neg__", []( TMC const& m ){ return operator-(m); } )
 .def( "__pow__", []( TMC const& m1, TMC const& m2 ){ return fadbad::pow(m1,m2); } )
.def( "__pow__", []( TMC const& m1, double const& m2 ){ return fadbad::pow(m1,m2); } )
.def( "__pow__", []( double const& m1, TMC const& m2 ){ return fadbad::pow(m1,m2); } )
.def( "__pow__", []( TMC const& m, int const& n ){ return fadbad::pow(m,n); } )
.def( "__pow__", []( int const& m, TMC const& n ){ return fadbad::pow(m,n); } )

 .def( py::self + double() )
 .def( double() + py::self )
 .def( py::self + py::self )
 .def( py::self - double() )
 .def( double() - py::self )
 .def( py::self - py::self )
 .def( py::self * double() )
 .def( double() * py::self )
 .def( py::self * py::self )
 .def( py::self / double() )
 .def( double() / py::self )
 .def( py::self / py::self )
 ;

m.def( "sqr",   []( TMC const& m ){ return fadbad::sqr(m); } );
m.def( "sqrt",  []( TMC const& m ){ return fadbad::sqrt(m); } );
m.def( "exp",   []( TMC const& m ){ return fadbad::exp(m); } );
m.def( "log",   []( TMC const& m ){ return fadbad::log(m); } );
m.def( "sin",   []( TMC const& m ){ return fadbad::sin(m); } );
m.def( "cos",   []( TMC const& m ){ return fadbad::cos(m); } );
m.def( "tan",   []( TMC const& m ){ return fadbad::tan(m); } );
// m.def( "cot",   []( TMC const& m ){ return fadbad::cot(m); } );
m.def( "asin",  []( TMC const& m ){ return fadbad::asin(m); } );
m.def( "acos",  []( TMC const& m ){ return fadbad::acos(m); } );
m.def( "atan",  []( TMC const& m ){ return fadbad::atan(m); } );
m.def( "sinh",  []( TMC const& m ){ return fadbad::sinh(m); } );
m.def( "cosh",  []( TMC const& m ){ return fadbad::cosh(m); } );
m.def( "tanh",  []( TMC const& m ){ return fadbad::tanh(m); } );
// m.def( "coth",  []( TMC const& m ){ return fadbad::coth(m); } );

;
}