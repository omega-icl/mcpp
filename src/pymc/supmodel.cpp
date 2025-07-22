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

#include "supmodel.hpp"

#include "pwlu.hpp"
typedef mc::PWLU PWLU;
typedef mc::SupModel<mc::PWLU> PWLSM;
typedef mc::SupVar<mc::PWLU> PWLSV;

#include "pwcu.hpp"
typedef mc::PWCU PWCU;
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;

namespace py = pybind11;

void mc_supmodel(py::module &m) 
{

py::class_<PWLSM>          pyPWLSModel( m, "PWLSModel" );
py::class_<PWLSM::Options> pyPWLSModelOptions( pyPWLSModel, "Options" );
py::class_<PWLSV>          pyPWLSVar(m,"PWLSVar");
py::class_<PWLU>           pyPWLU( m, "PWLU" );
py::class_<PWLU::Options>  pyPWLUOptions( pyPWLU, "Options" );

pyPWLSModel
 .def(
   py::init<size_t const&>(),
   py::arg("nvar"), 
   "constructor"
 )
 .def_readwrite(
   "options",
   &PWLSM::options
 )
 .def_property_readonly(
   "nvar",
   []( PWLSM& self ){ return self.nvar(); },
   "number of independent variables"
 )
 .def_property_readonly(
   "lbdvar",
   []( PWLSM& self ){ return self.lbdvar(); },
   "lower bounds of independent variables"
 )
 .def_property_readonly(
   "ubdvar",
   []( PWLSM& self ){ return self.ubdvar(); },
   "upper bounds of independent variables"
 )
 .def(
   "min",
   []( PWLSM& self, PWLSV& var, double const& minval )
     { return self.min( var, minval ); },
   "cut-off superposition relaxation from below at minval"
 )
 .def(
   "max",
   []( PWLSM& self, PWLSV& var, double const& maxval )
     { return self.max( var, maxval ); },
   "cut-off superposition relaxation from above at maxval"
 )
 .def(
   "uref",
   []( PWLSM& self, PWLSV& var, double const& lbd )
     { return self.min( var, lbd ); },
   "refine superposition underestimator to satisfy natural lower bound lbd"
 )
 .def(
   "oref",
   []( PWLSM& self, PWLSV& var, double const& ubd )
     { return self.oref( var, ubd ); },
   "refine superposition overestimator to satisfy natural upper bound ubd"
 )
;

pyPWLSModelOptions
 .def( py::init<>(), "default constructor" )
 .def( py::init<PWLSM::Options const&>(), "copy constructor" )
 .def_readwrite( "PROD_METH",      &PWLSM::Options::PROD_METH,      "reformulation method used for product terms" )
 .def_readwrite( "PROD_CUT",       &PWLSM::Options::PROD_CUT,       "wwether to cut superposition relaxations for product terms" )
 .def_readwrite( "SUM_TOL",        &PWLSM::Options::SUM_TOL,        "tolerance on range for univariate estimator propagation"  )
 .def_readwrite( "REF_WEIGHT",     &PWLSM::Options::REF_WEIGHT,     "Weight in the overloaded functions mc::min and mc::max" )
 .def_readwrite( "MAX_SUBDIV",     &PWLSM::Options::MAX_SUBDIV,     "Maximal number of subdivisions, for univariate estimators on adaptive grids only" )
 .def_readwrite( "USE_SHADOW",     &PWLSM::Options::USE_SHADOW,     "whether to enable shadow estimators" )
 .def_readwrite( "DISPLAY_SHADOW", &PWLSM::Options::DISPLAY_SHADOW, "whether to display shadow estimators" )
 .def_readwrite( "DISPLAY_DIGITS", &PWLSM::Options::DISPLAY_DIGITS, "number of digits displayed" )
 .def( "reset", []( PWLSM::Options& self ){ return self.reset(); }, "reset options" )
;

py::enum_<PWLSM::Options::PROD_REF>(pyPWLSModelOptions, "reformulation method for product term")
 .value("NONE",    PWLSM::Options::PROD_REF::NONE,    "DC decomposition w/o rescaling")
 .value("PARTIAL", PWLSM::Options::PROD_REF::PARTIAL, "DC decomposition w/ range rescaling")
 .value("FULL",    PWLSM::Options::PROD_REF::FULL,    "DC decomposition w/ range and midpoint rescaling")
 .value("LOG",     PWLSM::Options::PROD_REF::LOG,     "log-transform w/ range and midpoint rescaling")
 .export_values()
;

pyPWLSVar
 .def(
   py::init<double const&>(),
   py::arg("cst")=0., 
   "constructor for constant"
 )
 .def(
   py::init<PWLSM&>(),
   py::arg("mod"), 
   "default constructor"
 )
 .def(
   py::init<PWLSM&, unsigned int, I const&, size_t const>(),
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("bnd"), 
   py::arg("ndiv")=1, 
   "constructor for variable"
 )
 .def(
   py::init<PWLSM&, unsigned int, PWLU const&>(),
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("est"), 
   "constructor for variable"
 )
 .def(
   py::init<PWLSV const&>(),
   "copy constructor"
 )
 .def(
   "set",
   []( PWLSV& self, PWLSM& mod ){ return self.set( mod ); },
   "set model environment"
 )
 .def(
   "set",
   []( PWLSV& self, PWLSM& mod, unsigned int ndx, I const& bnd, size_t const ndiv )
     { return self.set( mod, ndx, bnd, ndiv ); },
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("bnd"), 
   py::arg("ndiv")=1, 
   "set independent variable"
 )
 .def(
   "set",
   []( PWLSV& self, PWLSM& mod, unsigned int ndx, PWLU const& est )
     { return self.set( mod, ndx, est ); },
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("est"), 
   "set independent variable"
 )
 .def(
   "set",
   []( PWLSV& self, double const& cst )
     { return self.set( cst ); },
   py::arg("cst"), 
   "set constant"
 )
 .def_property_readonly(
   "ndep",
   []( PWLSV& self ){ return self.ndep(); },
   "number of dependencies"
 )
 .def_property_readonly(
   "sdep",
   []( PWLSV& self ){ return self.sdep(); },
   "set of dependencies"
 )
 .def_property_readonly(
   "cst",
   []( PWLSV& self ){ return self.cst(); },
   "constant field"
 )
 .def(
   "uest",
   []( PWLSV& self, unsigned int opt ){ return self.uest( opt ); },
   py::arg("opt")=0, 
   "superposition underestimator"
 )
 .def(
   "oest",
   []( PWLSV& self, unsigned int opt ){ return self.oest( opt ); },
   py::arg("opt")=0, 
   "superposition overestimator"
 )
 .def(
   "l",
   []( PWLSV& self, unsigned int opt ){ return self.l( opt ); },
   py::arg("opt")=2, 
   "lower bound of superposition relaxation"
 )
 .def(
   "u",
   []( PWLSV& self, unsigned int opt ){ return self.u( opt ); },
   py::arg("opt")=2, 
   "upper bound of superposition relaxation"
 )
 .def(
   "uval",
   []( PWLSV& self, std::map<unsigned int,double> const& x, unsigned int opt ){ return self.uval( x, opt ); },
   py::arg("x"), 
   py::arg("opt")=2, 
   "evaluation of superposition underestimator"
 )
 .def(
   "oval",
   []( PWLSV& self, std::map<unsigned int,double> const& x, unsigned int opt ){ return self.oval( x, opt ); },
   py::arg("x"), 
   py::arg("opt")=2, 
   "evaluation of superposition overestimator"
 )
 //.def( "copy", []( PWLSV const& self ){ return PWLSV( self ); } )
 //.def( "assign", py::overload_cast<double const&>( &PWLSV::operator= ), py::arg("val") )
 //.def( "assign", py::overload_cast<PWLSV const&>( &PWLSV::operator= ), py::arg("var") )
 .def(
   "__str__",
   []( PWLSV const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( PWLSV const& self ){ std::ostringstream os; os << self; return os.str(); }
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
// .def( "__abs__", []( PWLSV const& m ){ return mc::Op<PWLSV>::abs(m); } )
 .def( "__pow__", []( PWLSV const& m, int const n ){ return mc::Op<PWLSV>::pow(m,n); } )
 .def( "__pow__", []( PWLSV const& m, double const& r ){ return mc::Op<PWLSV>::pow(m,r); } )
 .def( "__pow__", []( PWLSV const& m, PWLSV const& mm ){ return mc::Op<PWLSV>::pow(m,mm); } )
 .def( "__pow__", []( double const& r, PWLSV const& m ){ return mc::Op<PWLSV>::pow(r,m); } )
// .def( py::self == py::self ) 
// .def( py::self != py::self )
// .def( py::self <= py::self )
// .def( py::self >= py::self )
// .def( py::self < py::self )
// .def( py::self > py::self )
; 

m.def( "inv",    []( PWLSV const& x ){ return mc::inv(x); } );
m.def( "sqr",    []( PWLSV const& x ){ return mc::sqr(x); } );
m.def( "sqrt",   []( PWLSV const& x ){ return mc::sqrt(x); } );
m.def( "exp",    []( PWLSV const& x ){ return mc::exp(x); } );
m.def( "log",    []( PWLSV const& x ){ return mc::log(x); } );
//m.def( "cos",    []( PWLSV const& x ){ return mc::cos(x); } );
//m.def( "sin",    []( PWLSV const& x ){ return mc::sin(x); } );
//m.def( "tan",    []( PWLSV const& x ){ return mc::tan(x); } );
//m.def( "acos",   []( PWLSV const& x ){ return mc::acos(x); } );
//m.def( "asin",   []( PWLSV const& x ){ return mc::asin(x); } );
//m.def( "atan",   []( PWLSV const& x ){ return mc::atan(x); } );
//m.def( "cosh",   []( PWLSV const& x ){ return mc::cosh(x); } );
//m.def( "sinh",   []( PWLSV const& x ){ return mc::sinh(x); } );
m.def( "tanh",   []( PWLSV const& x ){ return mc::tanh(x); } );
m.def( "fabs",   []( PWLSV const& x ){ return mc::fabs(x); } );
m.def( "relu",   []( PWLSV const& x ){ return mc::relu(x); } );
m.def( "xlog",   []( PWLSV const& x ){ return mc::xlog(x); } );
//m.def( "erf",    []( PWLSV const& x ){ return mc::erf(x); } );
//m.def( "erfc",   []( PWLSV const& x ){ return mc::erfc(x); } );
m.def( "pow",    []( PWLSV const& x, int const n ){ return mc::pow(x,n); } );
m.def( "pow",    []( PWLSV const& x, double const& r ){ return mc::pow(x,r); } );
m.def( "pow",    []( PWLSV const& x, PWLSV const& y ){ return mc::pow(x,y); } );
m.def( "pow",    []( double const& r, PWLSV const& y ){ return mc::pow(r,y); } );
m.def( "cheb",   []( PWLSV const& x, unsigned const n ){ return mc::cheb(x,n); } );
m.def( "max",    []( PWLSV const& x, double const& y ){ return mc::max(x,y); } );
m.def( "min",    []( PWLSV const& x, double const& y ){ return mc::min(x,y); } );
m.def( "max",    []( PWLSV const& x, PWLSV const& y ){ return mc::max(x,y); } );
m.def( "min",    []( PWLSV const& x, PWLSV const& y ){ return mc::min(x,y); } );

pyPWLU
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   py::init<double const&, double const&, double const&>(),
   py::arg("xL"),
   py::arg("xU"),
   py::arg("y"),
   "constructor of constant"
 )
 .def(
   py::init<double const&, double const&, size_t const>(),
   py::arg("xL"),
   py::arg("xU"),
   py::arg("N")=1,
   "constructor of independent variable"
 )
 .def(
   py::init<double const&, std::vector<double> const&>(),
   py::arg("xL"),
   py::arg("dx"),
   "constructor of variable"
 )
 .def(
   py::init<PWLU const&>(),
   "copy constructor"
 )
 .def_readwrite_static(
   "options",
   &PWLU::options
 )
 .def(
   "set",
   []( PWLU& self, double const& y )
     { return self.set( y ); },
   py::arg("y"),
   "set constant"
 )
 .def(
   "set",
   []( PWLU& self, double const& xL, double const& xU, double const& y )
     { return self.set( xL, xU, y ); },
   py::arg("xL"),
   py::arg("xU"),
   py::arg("y"),
   "set constant"
 )
 .def(
   "set",
   []( PWLU& self, double const& xL, double const& xU, size_t const N )
     { return self.set( xL, xU, N ); },
   py::arg("xL"),
   py::arg("xU"),
   py::arg("N")=1,
   "set independent variable"
 )
 .def(
   "set",
   []( PWLU& self, double const& xL, std::vector<double> const& dx )
     { return self.set( xL, dx ); },
   py::arg("xL"),
   py::arg("dx"),
   "set independent variable"
 )
 .def(
   "insert",
   []( PWLU& self, std::vector<double> const& dx )
     { return self.insert( dx ); },
   py::arg("dx"),
   "insert breakpoints"
 )
 .def(
   "insert",
   []( PWLU& self, double const& x )
     { return self.insert( x ); },
   py::arg("x"),
   "insert breakpoint"
 )
 .def_property_readonly(
   "l",
   []( PWLU& self ){ return self.l(); },
   "lower bound of univariate estimator"
 )
 .def_property_readonly(
   "u",
   []( PWLU& self ){ return self.u(); },
   "upper bound of univariate estimator"
 )
 .def_property_readonly(
   "w",
   []( PWLU& self ){ return self.w(); },
   "width of univariate estimator"
 )
 .def(
   "lval",
   []( PWLU& self, double const& x ){ return self.l( x ); },
   "evaluate univariate estimator (lower range)"
 )
 .def(
   "uval",
   []( PWLU& self, double const& x ){ return self.l( x ); },
   "evaluate univariate estimator (upper range)"
 )
 .def_property(
   "xL",
   []( PWLU& self ){ return self.xL(); },
   []( PWLU& self, double const& x ){ self.xL() = x; },
   "variable lower range"
 )
 .def_property(
   "yL",
   []( PWLU& self ){ return self.yL(); },
   []( PWLU& self, double const& y ){ self.yL() = y; },
   "estimator value at lower range"
 )
 .def_property(
   "dx",
   []( PWLU& self ){ return self.dx(); },
   []( PWLU& self, std::vector<double> const& dx ){ self.dx() = dx; },
   "variable partition steps"
 )
 .def_property(
   "dy",
   []( PWLU& self ){ return self.dy(); },
   []( PWLU& self, std::vector<double> const& dy ){ self.dy() = dy; },
   "estimator slopes on partition"
 )
 .def(
   "__str__",
   []( PWLU const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( PWLU const& self ){ std::ostringstream os; os << self; return os.str(); }
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
 .def( double() * py::self )
 .def( py::self * double() )
 .def( py::self /= double() )
 .def( py::self / double() )
 .def(
   "min",
   []( PWLU& self, double const& cutoff )
     { return self.min( cutoff ); },
   "cut-off estimator from below"
 )
 .def(
   "max",
   []( PWLU& self, double const& cutoff )
     { return self.max( cutoff ); },
   "cut-off estimator from above"
 )
 .def(
   "clean",
   []( PWLU& self, bool const under )
     { return self.clean( under ); },
   "remove partitions with step smaller than set tolerances (options.BKPTATOl, options.BKPTRTOL)"
 )
 .def(
   "merge",
   []( PWLU& self, bool const under )
     { return self.clean( under ); },
   "merge partitions with identical slope within set tolerances (options.BKPTATOl, options.BKPTRTOL)"
 )
 .def(
   "reduce",
   []( PWLU& self, bool const under, size_t const nseg )
     { return self.clean( under ); },
   "reduce partition size to within nseg segments through successive relaxation"
 )
;

pyPWLUOptions
 .def( py::init<>(), "default constructor" )
 .def( py::init<PWLU::Options const&>(), "copy constructor" )
 .def_readwrite( "BKPTATOL",   &PWLU::Options::BKPTATOL,   "absolute tolerance in managing breakpoints" )
 .def_readwrite( "BKPTRTOL",   &PWLU::Options::BKPTRTOL,   "relative tolerance in managing breakpoints" )
 .def_readwrite( "REDUCEMETH", &PWLU::Options::REDUCEMETH, "Method for breakpoint reduction" )
 .def_readwrite( "DISPNUM",    &PWLU::Options::DISPNUM,    "number of digits displayed" )
 .def( "reset", []( PWLU::Options& self ){ return self.reset(); }, "reset options" )
;



py::class_<PWCSM>          pyPWCSModel( m, "PWCSModel" );
py::class_<PWCSM::Options> pyPWCSModelOptions( pyPWCSModel, "Options" );
py::class_<PWCSV>          pyPWCSVar(m,"PWCSVar");
py::class_<PWCU>           pyPWCU( m, "PWCU" );
py::class_<PWCU::Options>  pyPWCUOptions( pyPWCU, "Options" );

pyPWCSModel
 .def(
   py::init<size_t const&>(),
   py::arg("nvar"), 
   "constructor"
 )
 .def_readwrite(
   "options",
   &PWCSM::options
 )
 .def_property_readonly(
   "nvar",
   []( PWCSM& self ){ return self.nvar(); },
   "number of independent variables"
 )
 .def_property_readonly(
   "lbdvar",
   []( PWCSM& self ){ return self.lbdvar(); },
   "lower bounds of independent variables"
 )
 .def_property_readonly(
   "ubdvar",
   []( PWCSM& self ){ return self.ubdvar(); },
   "upper bounds of independent variables"
 )
 .def(
   "min",
   []( PWCSM& self, PWCSV& var, double const& minval )
     { return self.min( var, minval ); },
   "cut-off superposition relaxation from below at minval"
 )
 .def(
   "max",
   []( PWCSM& self, PWCSV& var, double const& maxval )
     { return self.max( var, maxval ); },
   "cut-off superposition relaxation from above at maxval"
 )
 .def(
   "uref",
   []( PWCSM& self, PWCSV& var, double const& lbd )
     { return self.min( var, lbd ); },
   "refine superposition underestimator to satisfy natural lower bound lbd"
 )
 .def(
   "oref",
   []( PWCSM& self, PWCSV& var, double const& ubd )
     { return self.oref( var, ubd ); },
   "refine superposition overestimator to satisfy natural upper bound ubd"
 )
;

pyPWCSModelOptions
 .def( py::init<>(), "default constructor" )
 .def( py::init<PWCSM::Options const&>(), "copy constructor" )
 .def_readwrite( "PROD_METH",      &PWCSM::Options::PROD_METH,      "reformulation method used for product terms" )
 .def_readwrite( "PROD_CUT",       &PWCSM::Options::PROD_CUT,       "wwether to cut superposition relaxations for product terms" )
 .def_readwrite( "SUM_TOL",        &PWCSM::Options::SUM_TOL,        "tolerance on range for univariate estimator propagation"  )
 .def_readwrite( "REF_WEIGHT",     &PWCSM::Options::REF_WEIGHT,     "Weight in the overloaded functions mc::min and mc::max" )
 .def_readwrite( "MAX_SUBDIV",     &PWCSM::Options::MAX_SUBDIV,     "Maximal number of subdivisions, for univariate estimators on adaptive grids only" )
 .def_readwrite( "USE_SHADOW",     &PWCSM::Options::USE_SHADOW,     "whether to enable shadow estimators" )
 .def_readwrite( "DISPLAY_SHADOW", &PWCSM::Options::DISPLAY_SHADOW, "whether to display shadow estimators" )
 .def_readwrite( "DISPLAY_DIGITS", &PWCSM::Options::DISPLAY_DIGITS, "number of digits displayed" )
 .def( "reset", []( PWCSM::Options& self ){ return self.reset(); }, "reset options" )
;

py::enum_<PWCSM::Options::PROD_REF>(pyPWCSModelOptions, "reformulation method for product term")
 .value("NONE",    PWCSM::Options::PROD_REF::NONE,    "DC decomposition w/o rescaling")
 .value("PARTIAL", PWCSM::Options::PROD_REF::PARTIAL, "DC decomposition w/ range rescaling")
 .value("FULL",    PWCSM::Options::PROD_REF::FULL,    "DC decomposition w/ range and midpoint rescaling")
 .value("LOG",     PWCSM::Options::PROD_REF::LOG,     "log-transform w/ range and midpoint rescaling")
 .export_values()
;

pyPWCSVar
 .def(
   py::init<double const&>(),
   py::arg("cst")=0., 
   "constructor for constant"
 )
 .def(
   py::init<PWCSM&>(),
   py::arg("mod"), 
   "default constructor"
 )
 .def(
   py::init<PWCSM&, unsigned int, I const&, size_t const>(),
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("bnd"), 
   py::arg("ndiv")=1, 
   "constructor for variable"
 )
 .def(
   py::init<PWCSM&, unsigned int, PWCU const&>(),
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("est"), 
   "constructor for variable"
 )
 .def(
   py::init<PWCSV const&>(),
   "copy constructor"
 )
 .def(
   "set",
   []( PWCSV& self, PWCSM& mod ){ return self.set( mod ); },
   "set model environment"
 )
 .def(
   "set",
   []( PWCSV& self, PWCSM& mod, unsigned int ndx, I const& bnd, size_t const ndiv )
     { return self.set( mod, ndx, bnd, ndiv ); },
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("bnd"), 
   py::arg("ndiv")=1, 
   "set independent variable"
 )
 .def(
   "set",
   []( PWCSV& self, PWCSM& mod, unsigned int ndx, PWCU const& est )
     { return self.set( mod, ndx, est ); },
   py::arg("mod"), 
   py::arg("ndx"), 
   py::arg("est"), 
   "set independent variable"
 )
 .def(
   "set",
   []( PWCSV& self, double const& cst )
     { return self.set( cst ); },
   py::arg("cst"), 
   "set constant"
 )
 .def_property_readonly(
   "ndep",
   []( PWCSV& self ){ return self.ndep(); },
   "number of dependencies"
 )
 .def_property_readonly(
   "sdep",
   []( PWCSV& self ){ return self.sdep(); },
   "set of dependencies"
 )
 .def_property_readonly(
   "cst",
   []( PWCSV& self ){ return self.cst(); },
   "constant field"
 )
 .def(
   "uest",
   []( PWCSV& self, unsigned int opt ){ return self.uest( opt ); },
   py::arg("opt")=0, 
   "superposition underestimator"
 )
 .def(
   "oest",
   []( PWCSV& self, unsigned int opt ){ return self.oest( opt ); },
   py::arg("opt")=0, 
   "superposition overestimator"
 )
 .def(
   "l",
   []( PWCSV& self, unsigned int opt ){ return self.l( opt ); },
   py::arg("opt")=2, 
   "lower bound of superposition relaxation"
 )
 .def(
   "u",
   []( PWCSV& self, unsigned int opt ){ return self.u( opt ); },
   py::arg("opt")=2, 
   "upper bound of superposition relaxation"
 )
 .def(
   "uval",
   []( PWCSV& self, std::map<unsigned int,double> const& x, unsigned int opt ){ return self.uval( x, opt ); },
   py::arg("x"), 
   py::arg("opt")=2, 
   "evaluation of superposition underestimator"
 )
 .def(
   "oval",
   []( PWCSV& self, std::map<unsigned int,double> const& x, unsigned int opt ){ return self.oval( x, opt ); },
   py::arg("x"), 
   py::arg("opt")=2, 
   "evaluation of superposition overestimator"
 )
 //.def( "copy", []( PWCSV const& self ){ return PWCSV( self ); } )
 //.def( "assign", py::overload_cast<double const&>( &PWCSV::operator= ), py::arg("val") )
 //.def( "assign", py::overload_cast<PWCSV const&>( &PWCSV::operator= ), py::arg("var") )
 .def(
   "__str__",
   []( PWCSV const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( PWCSV const& self ){ std::ostringstream os; os << self; return os.str(); }
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
// .def( "__abs__", []( PWCSV const& m ){ return mc::Op<PWCSV>::abs(m); } )
 .def( "__pow__", []( PWCSV const& m, int const n ){ return mc::Op<PWCSV>::pow(m,n); } )
 .def( "__pow__", []( PWCSV const& m, double const& r ){ return mc::Op<PWCSV>::pow(m,r); } )
 .def( "__pow__", []( PWCSV const& m, PWCSV const& mm ){ return mc::Op<PWCSV>::pow(m,mm); } )
 .def( "__pow__", []( double const& r, PWCSV const& m ){ return mc::Op<PWCSV>::pow(r,m); } )
// .def( py::self == py::self ) 
// .def( py::self != py::self )
// .def( py::self <= py::self )
// .def( py::self >= py::self )
// .def( py::self < py::self )
// .def( py::self > py::self )
; 

m.def( "inv",    []( PWCSV const& x ){ return mc::inv(x); } );
m.def( "sqr",    []( PWCSV const& x ){ return mc::sqr(x); } );
m.def( "sqrt",   []( PWCSV const& x ){ return mc::sqrt(x); } );
m.def( "exp",    []( PWCSV const& x ){ return mc::exp(x); } );
m.def( "log",    []( PWCSV const& x ){ return mc::log(x); } );
//m.def( "cos",    []( PWCSV const& x ){ return mc::cos(x); } );
//m.def( "sin",    []( PWCSV const& x ){ return mc::sin(x); } );
//m.def( "tan",    []( PWCSV const& x ){ return mc::tan(x); } );
//m.def( "acos",   []( PWCSV const& x ){ return mc::acos(x); } );
//m.def( "asin",   []( PWCSV const& x ){ return mc::asin(x); } );
//m.def( "atan",   []( PWCSV const& x ){ return mc::atan(x); } );
//m.def( "cosh",   []( PWCSV const& x ){ return mc::cosh(x); } );
//m.def( "sinh",   []( PWCSV const& x ){ return mc::sinh(x); } );
m.def( "tanh",   []( PWCSV const& x ){ return mc::tanh(x); } );
m.def( "fabs",   []( PWCSV const& x ){ return mc::fabs(x); } );
m.def( "relu",   []( PWCSV const& x ){ return mc::relu(x); } );
m.def( "xlog",   []( PWCSV const& x ){ return mc::xlog(x); } );
//m.def( "erf",    []( PWCSV const& x ){ return mc::erf(x); } );
//m.def( "erfc",   []( PWCSV const& x ){ return mc::erfc(x); } );
m.def( "pow",    []( PWCSV const& x, int const n ){ return mc::pow(x,n); } );
m.def( "pow",    []( PWCSV const& x, double const& r ){ return mc::pow(x,r); } );
m.def( "pow",    []( PWCSV const& x, PWCSV const& y ){ return mc::pow(x,y); } );
m.def( "pow",    []( double const& r, PWCSV const& y ){ return mc::pow(r,y); } );
m.def( "cheb",   []( PWCSV const& x, unsigned const n ){ return mc::cheb(x,n); } );
m.def( "max",    []( PWCSV const& x, double const& y ){ return mc::max(x,y); } );
m.def( "min",    []( PWCSV const& x, double const& y ){ return mc::min(x,y); } );
m.def( "max",    []( PWCSV const& x, PWCSV const& y ){ return mc::max(x,y); } );
m.def( "min",    []( PWCSV const& x, PWCSV const& y ){ return mc::min(x,y); } );

pyPWCU
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   py::init<double const&, double const&, double const&, size_t const&>(),
   py::arg("xL"),
   py::arg("xU"),
   py::arg("y"),
   py::arg("N"),
   "constructor of constant"
 )
 .def(
   py::init<double const&, double const&, size_t const>(),
   py::arg("xL"),
   py::arg("xU"),
   py::arg("N"),
   "constructor of variable"
 )
 .def(
   py::init<PWCU const&>(),
   "copy constructor"
 )
 .def_readwrite_static(
   "options",
   &PWCU::options
 )
 .def(
   "set",
   []( PWCU& self, double const& y )
     { return self.set( y ); },
   py::arg("y"),
   "set constant"
 )
 .def(
   "set",
   []( PWCU& self, double const& xL, double const& xU, double const& y, size_t const& N )
     { return self.set( xL, xU, y, N ); },
   py::arg("xL"),
   py::arg("xU"),
   py::arg("y"),
   py::arg("N"),
   "set constant"
 )
 .def(
   "set",
   []( PWCU& self, double const& xL, double const& xU, size_t const& N )
     { return self.set( xL, xU, N ); },
   py::arg("xL"),
   py::arg("xU"),
   py::arg("N"),
   "set variable"
 )
 .def_property_readonly(
   "l",
   []( PWCU& self ){ return self.l(); },
   "lower bound of univariate estimator"
 )
 .def_property_readonly(
   "u",
   []( PWCU& self ){ return self.u(); },
   "upper bound of univariate estimator"
 )
 .def_property_readonly(
   "w",
   []( PWCU& self ){ return self.w(); },
   "width of univariate estimator"
 )
 .def(
   "lval",
   []( PWCU& self, double const& x ){ return self.l( x ); },
   "evaluate univariate estimator (lower range)"
 )
 .def(
   "uval",
   []( PWCU& self, double const& x ){ return self.l( x ); },
   "evaluate univariate estimator (upper range)"
 )
 .def_property_readonly(
   "n",
   []( PWCU& self ){ return self.size(); },
   "partition size"
 )
 .def_property(
   "xL",
   []( PWCU& self ){ return self.xL(); },
   []( PWCU& self, double const& x ){ self.xL() = x; },
   "variable lower range"
 )
 .def_property(
   "xU",
   []( PWCU& self ){ return self.xU(); },
   []( PWCU& self, double const& x ){ self.xU() = x; },
   "variable lower range"
 )
 .def_property(
   "yL",
   []( PWCU& self ){ return self.yL(); },
   []( PWCU& self, std::vector<double> const& y ){ self.yL() = y; },
   "underestimator value on partition"
 )
 .def_property(
   "yU",
   []( PWCU& self ){ return self.yU(); },
   []( PWCU& self, std::vector<double> const& y ){ self.yU() = y; },
   "overestimator value on partition"
 )
 .def(
   "__str__",
   []( PWCU const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( PWCU const& self ){ std::ostringstream os; os << self; return os.str(); }
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
 .def( double() * py::self )
 .def( py::self * double() )
 .def( py::self /= double() )
 .def( py::self / double() )
 .def(
   "min",
   []( PWCU& self, double const& cutoff )
     { return self.min( cutoff ); },
   "cut-off estimator from below"
 )
 .def(
   "max",
   []( PWCU& self, double const& cutoff )
     { return self.max( cutoff ); },
   "cut-off estimator from above"
 )
;

pyPWCUOptions
 .def( py::init<>(), "default constructor" )
 .def( py::init<PWCU::Options const&>(), "copy constructor" )
 .def_readwrite( "BKPTATOL",   &PWCU::Options::BKPTATOL,   "absolute tolerance in managing breakpoints" )
 .def_readwrite( "BKPTRTOL",   &PWCU::Options::BKPTRTOL,   "relative tolerance in managing breakpoints" )
 .def_readwrite( "DISPNUM",    &PWCU::Options::DISPNUM,    "number of digits displayed" )
 .def( "reset", []( PWCU::Options& self ){ return self.reset(); }, "reset options" )
;
}
