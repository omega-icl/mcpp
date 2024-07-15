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

namespace py = pybind11;

void mc_interval( py::module_ &m )
{

py::class_<I> pyInterval( m, "Interval" );
pyInterval
 .def( py::init<>() )
 .def( py::init<double const&>() )
 .def( py::init<double const&, double const&>() )
 .def( py::init<I const&>() )
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
 .def_readwrite_static( "options", &I::options )
 .def_property( "l", []( I& self ){ return mc::Op<I>::l(self); },
                     []( I& self, double const& l ){ self.l() = l; } )
 .def_property( "u", []( I& self ){ return mc::Op<I>::u(self); },
                     []( I& self, double const& u ){ self.u() = u; } )
#else
 .def_property_readonly( "l", []( I& self ){ return mc::Op<I>::l(self); } )
 .def_property_readonly( "u", []( I& self ){ return mc::Op<I>::u(self); } )
#endif
 .def( "__str__", []( I const& x ){ std::ostringstream Iss; Iss << x; return Iss.str(); } )
 .def( "__repr__", []( I const& x ){ std::ostringstream Iss; Iss << x; return Iss.str(); } )
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
 .def( "__abs__", []( I const& self ){ return mc::Op<I>::abs(self); } )
 .def( "__pow__", []( I const& self, int const n ){ return mc::Op<I>::pow(self,n); } )
 .def( "__pow__", []( I const& self, double const& r ){ return mc::Op<I>::pow(self,r); } )
 .def( "__pow__", []( I const& self, I const& J ){ return mc::Op<I>::pow(self,J); } )
;

#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
py::class_<I::Options> pyIntervalOptions( pyInterval, "Interval.Options" );
pyIntervalOptions
 .def( py::init<>() )
 .def( py::init<I::Options const&>() )
 .def_readwrite( "DISPLAY_DIGITS", &I::Options::DISPLAY_DIGITS )
;
#endif

m.def( "inv",   []( I const& x ){ return mc::Op<I>::inv(x); } );
m.def( "sqr",   []( I const& x ){ return mc::Op<I>::sqr(x); } );
m.def( "sqrt",  []( I const& x ){ return mc::Op<I>::sqrt(x); } );
m.def( "exp",   []( I const& x ){ return mc::Op<I>::exp(x); } );
m.def( "log",   []( I const& x ){ return mc::Op<I>::log(x); } );
m.def( "cos",   []( I const& x ){ return mc::Op<I>::cos(x); } );
m.def( "sin",   []( I const& x ){ return mc::Op<I>::sin(x); } );
m.def( "tan",   []( I const& x ){ return mc::Op<I>::tan(x); } );
m.def( "acos",  []( I const& x ){ return mc::Op<I>::acos(x); } );
m.def( "asin",  []( I const& x ){ return mc::Op<I>::asin(x); } );
m.def( "atan",  []( I const& x ){ return mc::Op<I>::atan(x); } );
m.def( "cosh",  []( I const& x ){ return mc::Op<I>::cosh(x); } );
m.def( "sinh",  []( I const& x ){ return mc::Op<I>::sinh(x); } );
m.def( "tanh",  []( I const& x ){ return mc::Op<I>::tanh(x); } );
m.def( "fabs",  []( I const& x ){ return mc::Op<I>::fabs(x); } );
m.def( "relu",  []( I const& x ){ return mc::Op<I>::max(x,0.); } );
m.def( "xlog",  []( I const& x ){ return mc::Op<I>::xlog(x); } );
m.def( "fstep", []( I const& x ){ return mc::Op<I>::xlog(x); } );
m.def( "bstep", []( I const& x ){ return mc::Op<I>::xlog(x); } );
m.def( "erf",   []( I const& x ){ return mc::Op<I>::erf(x); } );
m.def( "erfc",  []( I const& x ){ return mc::Op<I>::erfc(x); } );
#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
m.def( "lmtd",  []( mc::Interval const& I1, mc::Interval const& I2 ){ return mc::lmtd(I1,I2); } );
m.def( "rlmtd", []( mc::Interval const& I1, mc::Interval const& I2 ){ return mc::rlmtd(I1,I2); } );
m.def( "pow",   []( double const& r, I const& y ){ return mc::Op<I>::pow(r,y); } );
#endif
m.def( "cheb",  []( I const& x, unsigned const n ){ return mc::Op<I>::cheb(x,n); } );
m.def( "hull",  []( I const& x, I const& y ){ return mc::Op<I>::hull(x,y); } );
m.def( "max",   []( I const& x, I const& y ){ return mc::Op<I>::min(x,y); } );
m.def( "min",   []( I const& x, I const& y ){ return mc::Op<I>::max(x,y); } );
m.def( "inter",  []( I& z, I const& x, I const& y ){ return mc::Op<I>::inter(z,x,y); } );
}

