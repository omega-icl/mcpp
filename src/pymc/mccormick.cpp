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

namespace py = pybind11;

void mc_mccormick(py::module &m) 
{

py::class_<mc::McCormick<I>> pyMcCormick(m,"McCormick");
pyMcCormick
    
    // .def( py::init<>() )
    // .def( py::init<double const>(),py::arg("c")= 0.0 )
    // .def( py::init<I const&>(), py::arg("I")=I(0.0,1.0))
    .def( py::init<I const&, double const>(), py::arg("Interval")=I(0.0,1.0), py::arg("Value")= 0.0 )
    .def( py::init<I const&, double const, double const>(), py::arg("Interval")=I(0.0,1.0), py::arg("CV")= 0.0, py::arg("CC")= 0.0 )    
    .def( py::init<mc::McCormick<I> const&>(), py::arg("I")=I(0.0,1.0))

    .def_property_readonly("number_sub",static_cast<unsigned int (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::nsub))
    
    .def_property("Interval",static_cast<const I& (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::I), static_cast<void (mc::McCormick<I>::*)(const I&)>(&mc::McCormick<I>::I))  
    .def_property_readonly("LB",static_cast<double (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::l))
    .def_property_readonly("UB",static_cast<double (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::u))
    .def_property("CV", static_cast<double (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::cv), static_cast<void (mc::McCormick<I>::*)(const double&)>(&mc::McCormick<I>::cv))
    .def_property("CC",static_cast<double (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::cc), static_cast<void (mc::McCormick<I>::*)(const double&)>(&mc::McCormick<I>::cc))
    
    // .def_property_readonly("CV_sub",static_cast<const double* (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::cvsub))
    // .def_property_readonly("CC_sub",static_cast<const double* (mc::McCormick<I>::*)() const>(&mc::McCormick<I>::ccsub))

    .def_readwrite_static("options",&mc::McCormick<I>::options)
    
    .def("get_CV_sub", static_cast<double (mc::McCormick<I>::*)(const unsigned int) const>(&mc::McCormick<I>::cvsub),py::arg("i"))
    .def("get_CV_sub", [](const mc::McCormick<I>& self)
        {   
            pybind11::list list_res;
            const double* item = self.cvsub();
            for (unsigned int i = 0; i < self.nsub(); i++){
                list_res.append(item[i]);
            }
            return list_res;
        })

    .def("get_CC_sub", static_cast<double (mc::McCormick<I>::*)(const unsigned int) const>(&mc::McCormick<I>::ccsub),py::arg("i"))    
    .def("get_CC_sub", [](const mc::McCormick<I>& self)
        {   
            pybind11::list list_res;
            const double* item = self.ccsub();
            for (unsigned int i = 0; i < self.nsub(); i++){
                list_res.append(item[i]);
        }
            return list_res;
        })
    
    .def("set_value",&mc::McCormick<I>::c, py::arg("Value").noconvert())

    .def("set_subgradient",
            static_cast<mc::McCormick<I>& (mc::McCormick<I>::*)(const unsigned int)>(&mc::McCormick<I>::sub),
             py::arg("dim"))   
    .def("set_subgradient",
            static_cast<mc::McCormick<I>& (mc::McCormick<I>::*)(const unsigned int, const unsigned int)>(&mc::McCormick<I>::sub),
            py::arg("dim"), py::arg("i"))
    .def("set_subgradient",
            static_cast<mc::McCormick<I>& (mc::McCormick<I>::*)(const unsigned int, const double*, const double*)>(&mc::McCormick<I>::sub),py::arg("dim"), py::arg("cv_sub"), py::arg("cc_sub"))

    .def("cut", static_cast<mc::McCormick<I>& (mc::McCormick<I>::*)()>(&mc::McCormick<I>::cut))

    .def("laff", static_cast<double (mc::McCormick<I>::*)(const double*, const double*) const>(&mc::McCormick<I>::laff))
    .def("laff", static_cast<double (mc::McCormick<I>::*)(const I*, const double*) const>(&mc::McCormick<I>::laff))

    .def("uaff", static_cast<double (mc::McCormick<I>::*)(const double*, const double*) const>(&mc::McCormick<I>::uaff))
    .def("uaff", static_cast<double (mc::McCormick<I>::*)(const I*, const double*) const>(&mc::McCormick<I>::uaff))

    // .def("reset_subgradient", &mc::McCormick<I>::_sub_reset)
    
    .def("copy",[](const mc::McCormick<I>& self) {return mc::McCormick<I>(self);})

    .def("assign", py::overload_cast<const double>(&mc::McCormick<I>::operator=),py::arg("Value"))
    .def("assign", py::overload_cast<const I&>(&mc::McCormick<I>::operator=),py::arg("Value"))    
    .def("assign", py::overload_cast<const mc::McCormick<I>&>(&mc::McCormick<I>::operator=),py::arg("Value"))

    .def("__str__", [](const mc::McCormick<I>& self) 
        {
            std::ostringstream os;
            os << self;
            return os.str(); 
        } 
        )
    
    .def("__repr__", [](const mc::McCormick<I>& self) 
        {
            std::ostringstream os;
            os << self;
            return os.str(); 
        } 
        )
    
// overload operators
    .def(py::self + double())
    .def(double() + py::self)
    .def(py::self + py::self)
    .def(py::self += double())
    .def(py::self += mc::McCormick<I>())  
    .def(py::self - double())
    .def(double() - py::self)
    .def(py::self - py::self)
    .def(py::self -= double())
    .def(py::self -= mc::McCormick<I>())
    .def(double() * py::self)
    .def(py::self * double())
    .def(py::self * py::self)
    .def(py::self *= double())
    .def(py::self *= mc::McCormick<I>())
    .def(py::self / double())
    .def(double() / py::self)
    .def(py::self / py::self)
    .def(py::self /= double())
    .def(py::self /= mc::McCormick<I>())
    .def(py::self == mc::McCormick<I>()) 
    .def(py::self != mc::McCormick<I>())
    .def(py::self <= mc::McCormick<I>())
    .def(py::self >= mc::McCormick<I>())
    .def(py::self < mc::McCormick<I>())
    .def(py::self > mc::McCormick<I>())

    // .def("sum", []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::McCormick<I>::_sum1(MC1,MC2); } )

    .def("__pow__", []( mc::McCormick<I> const& MC, int const n ){ return mc::pow(MC,n); } )
    .def("__pow__", []( mc::McCormick<I> const& MC, double const n ){ return mc::pow(MC,n); } )
    .def("__pow__", []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::pow(MC1,MC2); } )
    .def("__pow__", []( double const a, mc::McCormick<I> const& MC ){ return mc::pow(a,MC); } )
; 

// mathematical functions
m.def( "inv",   []( mc::McCormick<I> const& MC ){ return mc::inv(MC); } );
m.def( "sqr",   []( mc::McCormick<I> const& MC ){ return mc::sqr(MC); } );
m.def( "sqrt",  []( mc::McCormick<I> const& MC ){ return mc::sqrt(MC); } );
m.def( "exp",   []( mc::McCormick<I> const& MC ){ return mc::exp(MC); } );
m.def( "log",   []( mc::McCormick<I> const& MC ){ return mc::log(MC); } );

m.def( "erfc",  []( mc::McCormick<I> const& MC ){ return mc::erfc(MC); } );
m.def( "erf",   []( mc::McCormick<I> const& MC ){ return mc::erf(MC); } );

m.def( "xlog",  []( mc::McCormick<I> const& MC ){ return mc::xlog(MC); } );
m.def( "arrh",  []( mc::McCormick<I> const& MC, double const k){ return mc::arrh(MC,k); } );

m.def( "pow",   []( mc::McCormick<I> const& MC, int const n ){ return mc::pow(MC,n); } );
m.def( "pow",   []( mc::McCormick<I> const& MC, double const a ){ return mc::pow(MC,a); } );
m.def( "pow",   []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::pow(MC1,MC2); } );
m.def( "pow",   []( double const a, mc::McCormick<I> const& MC ){ return mc::pow(a,MC); } );

m.def( "prod",   []( unsigned int const n, mc::McCormick<I> const* MC ){ return mc::prod(n,MC); } );  
m.def( "monom",  []( unsigned int const n, mc::McCormick<I> const* MC, unsigned const *k ){ return mc::monom(n,MC,k); } );
m.def( "cheb",   []( mc::McCormick<I> const& MC, unsigned const n ){ return mc::cheb(MC,n); } );
m.def( "fabs",  []( mc::McCormick<I> const& MC ){ return mc::fabs(MC); } );

m.def( "min",   []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::min(MC1,MC2); } );
m.def( "min",   []( mc::McCormick<I> const& MC, double const a ){ return mc::min(MC,a); } );
m.def( "min",   []( double const a, mc::McCormick<I> const& MC ){ return mc::min(a,MC); } );
m.def( "min",   []( unsigned int const n, mc::McCormick<I> const* MC ){ return mc::min(n,MC); } );
m.def( "max",   []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::max(MC1,MC2); } );
m.def( "max",   []( mc::McCormick<I> const& MC, double const a ){ return mc::max(MC,a); } );
m.def( "max",   []( double const a, mc::McCormick<I> const& MC ){ return mc::max(a,MC); } );
m.def( "max",   []( unsigned int const n, mc::McCormick<I> const* MC ){ return mc::max(n,MC); } );

m.def( "fstep", []( mc::McCormick<I> const& MC ){ return mc::fstep(MC); } );
m.def( "bstep", []( mc::McCormick<I> const& MC ){ return mc::bstep(MC); } );
m.def( "ltcond",[]( I const& I0, mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::ltcond(I0,MC1,MC2); } );
m.def( "ltcond",[]( mc::McCormick<I> const& MC0, mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::ltcond(MC0,MC1,MC2); } );
m.def( "gtcond",[]( I const& I0, mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::gtcond(I0,MC1,MC2); } );
m.def( "gtcond",[]( mc::McCormick<I> const& MC0, mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::gtcond(MC0,MC1,MC2); } );
m.def( "lmtd",  []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::lmtd(MC1,MC2); } );
m.def( "rlmtd", []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::rlmtd(MC1,MC2); } );
m.def( "erfc",  []( mc::McCormick<I> const& MC ){ return mc::erfc(MC); } );
m.def( "erf",   []( mc::McCormick<I> const& MC ){ return mc::erf(MC); } );

m.def( "cos",   []( mc::McCormick<I> const& MC ){ return mc::cos(MC); } );
m.def( "sin",   []( mc::McCormick<I> const& MC ){ return mc::sin(MC); } );
m.def( "tan",   []( mc::McCormick<I> const& MC ){ return mc::tan(MC); } );  
m.def( "asin",  []( mc::McCormick<I> const& MC ){ return mc::asin(MC); } );
m.def( "acos",  []( mc::McCormick<I> const& MC ){ return mc::acos(MC); } );
m.def( "atan",  []( mc::McCormick<I> const& MC ){ return mc::atan(MC); } );
m.def( "cosh",  []( mc::McCormick<I> const& MC ){ return mc::cosh(MC); } ); 
m.def( "sinh",  []( mc::McCormick<I> const& MC ){ return mc::sinh(MC); } );
m.def( "tanh",  []( mc::McCormick<I> const& MC ){ return mc::tanh(MC); } );

m.def( "hull",  []( mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2 ){ return mc::hull(MC1,MC2); } );
m.def( "inter", []( mc::McCormick<I>& MC0, mc::McCormick<I> const& MC1, mc::McCormick<I> const& MC2){return mc::inter(MC0,MC1,MC2); })

;


py::class_<mc::McCormick<I>::Options> pyMcCormickOptions( pyMcCormick, "McCormick.Options");
pyMcCormickOptions
 .def(py::init<>())  
 .def(py::init<mc::McCormick<I>::Options const&>())
 .def_readwrite("ENVEL_USE", &mc::McCormick<I>::Options::ENVEL_USE)
 .def_readwrite("ENVEL_MAXIT", &mc::McCormick<I>::Options::ENVEL_MAXIT)
 .def_readwrite("ENVEL_TOL", &mc::McCormick<I>::Options::ENVEL_TOL)
 .def_readwrite("MVCOMP_USE", &mc::McCormick<I>::Options::MVCOMP_USE)
 .def_readwrite("MVCOMP_TOL", &mc::McCormick<I>::Options::MVCOMP_TOL)
 .def_readwrite("DISPLAY_DIGITS", &mc::McCormick<I>::Options::DISPLAY_DIGITS)
;


py::class_<fadbad::Op<mc::McCormick<I>>> pyFadpadMC( pyMcCormick, "FadpadMC");
pyFadpadMC
 .def( py::init<>())
//  .def( py::init<mc::Op<fadbad::F<mc::McCormick<I>>> const&>())

;

}