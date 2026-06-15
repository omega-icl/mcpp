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
   typedef boost::numeric::interval_lib::save_state<
             boost::numeric::interval_lib::rounded_transc_opp<double>
           > T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "specbnd.hpp"

typedef mc::Specbnd<I> SB;

namespace py = pybind11;

void mc_specbnd(py::module &m)
{

py::class_<SB> pySpecbnd(m, "Specbnd");

pySpecbnd
 // Constructors
 .def(
   py::init<>(),
   "default constructor"
 )
 .def(
   py::init<double const &>(),
   "constructor for constant value"
 )
 .def(
   py::init<I const &>(),
   "constructor for constant interval"
 )
 .def(
   py::init<I const &, size_t, size_t>(),
   py::arg("B"),
   py::arg("i"),
   py::arg("n"),
   "constructor for variable with range B, index i and n independent variables"
 )
 .def(
   py::init<SB const&>(),
   "copy constructor"
 )
 // Modifiers
 .def(
   "set",
   []( SB& self, I const& B, size_t i, size_t n )
     -> SB&
     { return self.set(B, i, n); },
   py::arg("B"),
   py::arg("i"),
   py::arg("n"),
   py::return_value_policy::reference_internal,
   "set variable with range B and index i of n independent variables"
 )
 .def(
   "dep",
   []( SB &self, size_t i, size_t n )
     -> SB&
     { return self.dep(i, n); },
   py::arg("i"),
   py::arg("n"),
   py::return_value_policy::reference_internal,
   "set variable index i of n independent variables"
    )
  // Accessors
  .def_property_readonly(
    "n",
    []( SB const &self )
      { return self.n(); },
    "number of independent variables"
  )
  .def_property_readonly(
    "I",
    [](SB const &self) { return self.I(); },
    "function range"
  )
  .def_property_readonly(
    "FI",
    []( SB const &self )
      { return self.FI(); },
    "gradient range"
  )
  .def_property_readonly(
    "SI",
    []( SB const &self )
      { return self.SI(); },
    "spectral range"
  )
 .def_static(
   "spectrum",
   []( std::vector<double> const& hess )
     { size_t N = std::sqrt( hess.size() ); assert( N*N == hess.size() );
       return SB::spectrum( N, hess.data() ); },
   py::arg("hess"),
   "compute spectrum of dense symmetric real matrix"
 )
 .def_static(
   "spectrum",
   []( size_t ndim, std::vector<double> const& hess, std::vector<unsigned> const& row,
       std::vector<unsigned> const& col )
     { size_t NNZ = row.size(); assert( !NNZ || NNZ == hess.size() );
       return SB::spectrum( ndim, hess.data(), NNZ,
              row.empty()? nullptr: row.data(),
              col.empty()? nullptr: col.data() ); },
   py::arg("ndim"),
   py::arg("hess"),
   py::arg("row")=std::vector<unsigned>(),
   py::arg("col")=std::vector<unsigned>(),
   "compute spectrum of sparse symmetric real matrix"
 )
 .def_static(
   "spectral_bound_gershgorin",
   []( std::vector<I> const& hess, std::vector<double> const& scal )
     { size_t N = std::sqrt( hess.size() ); assert( N*N == hess.size() );
       return SB::spectral_bound_gershgorin( N, hess.data(), 0, nullptr, nullptr,
              scal.empty()? nullptr: scal.data() ); },
   py::arg("hess"),
   py::arg("scal")=std::vector<double>(),
   "compute spectral bound of dense symmetric interval matrix using Gershgorin's method"
 )
 .def_static(
   "spectral_bound_gershgorin",
   []( size_t ndim, std::vector<I> const& hess, std::vector<int> const& row,
       std::vector<int> const& col, std::vector<double> const& scal )
     { size_t NNZ = row.size(); assert( !NNZ || NNZ == hess.size() );
       return SB::spectral_bound_gershgorin( ndim, hess.data(), NNZ, 
              row.empty()? nullptr: reinterpret_cast<unsigned const*>(row.data()),
              col.empty()? nullptr: reinterpret_cast<unsigned const*>(col.data()),
              scal.empty()? nullptr: scal.data() ); },
   py::arg("ndim"),
   py::arg("hess"),
   py::arg("row")=std::vector<int>(),
   py::arg("col")=std::vector<int>(),
   py::arg("scal")=std::vector<double>(),
   "compute spectral bound of sparse symmetric interval matrix using Gershgorin's method"
 )
 .def_static(
   "spectral_bound_rohn",
   []( std::vector<I> const& hess )
     { size_t N = std::sqrt( hess.size() ); assert( N*N == hess.size() );
       return SB::spectral_bound_rohn( N, hess.data() ); },
   py::arg("hess"),
   "compute spectral bound of dense symmetric interval matrix using Rohn's method"
 )
 .def_static(
   "spectral_bound_rohn",
   []( size_t ndim, std::vector<I> const& hess, std::vector<int> const& row,
       std::vector<int> const& col )
     { size_t NNZ = row.size(); assert( !NNZ || NNZ == hess.size() );
       return SB::spectral_bound_rohn( ndim, hess.data(), NNZ,
              row.empty()? nullptr: reinterpret_cast<unsigned const*>(row.data()),
              col.empty()? nullptr: reinterpret_cast<unsigned const*>(col.data()) ); },
   py::arg("hess"),
   py::arg("ndim"),
   py::arg("row")=std::vector<int>(),
   py::arg("col")=std::vector<int>(),
   "compute spectral bound of symmetric interval matrix using Rohn's method"
 )
 .def_static(
   "spectral_bound_hertz",
   []( std::vector<I> const& hess )
     { size_t N = std::sqrt( hess.size() ); assert( N*N == hess.size() );
       return SB::spectral_bound_hertz( N, hess.data() ); },
   py::arg("hess"),
   "compute spectral bound of dense symmetric interval matrix using Hertz's method"
 )
 .def_static(
   "spectral_bound_hertz",
   []( size_t ndim, std::vector<I> const& hess, std::vector<int> const& row,
       std::vector<int> const& col )
     { size_t NNZ = row.size(); assert( !NNZ || NNZ == hess.size() );
       return SB::spectral_bound_hertz( ndim, hess.data(), NNZ,
              row.empty()? nullptr: reinterpret_cast<unsigned const*>(row.data()),
              col.empty()? nullptr: reinterpret_cast<unsigned const*>(col.data()) ); },
   py::arg("hess"),
   py::arg("ndim"),
   py::arg("row")=std::vector<int>(),
   py::arg("col")=std::vector<int>(),
   "compute spectral bound of symmetric interval matrix using Hertz's method"
 )
 .def(
   "__str__",
   []( SB const& self ){ std::ostringstream os; os << self; return os.str(); }
 )
 .def(
   "__repr__",
   []( SB const& self ){ std::ostringstream os; os << self; return os.str(); }
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
 .def( double() / py::self )
 .def( py::self / double() )
 .def( py::self / py::self )
 .def( "__abs__", []( SB const& m ){ return mc::Op<SB>::abs(m); } )
 .def( "__pow__", []( SB const& m, int const n ){ return mc::pow(m,n); } )
 .def( "__pow__", []( SB const& m, double const& r ){ return mc::pow(m,r); } )
 .def( "__pow__", []( SB const& m, SB const& mm ){ return mc::pow(m,mm); } )
 .def( "__pow__", []( double const& r, SB const& m ){ return mc::pow(r,m); } )
;

m.def( "inv",    []( SB const& x ){ return mc::inv(x); } );
m.def( "sqr",    []( SB const& x ){ return mc::sqr(x); } );
m.def( "sqrt",   []( SB const& x ){ return mc::sqrt(x); } );
m.def( "exp",    []( SB const& x ){ return mc::exp(x); } );
m.def( "log",    []( SB const& x ){ return mc::log(x); } );
m.def( "xlog",   []( SB const& x ){ return mc::xlog(x); } );
m.def( "cos",    []( SB const& x ){ return mc::cos(x); } );
m.def( "sin",    []( SB const& x ){ return mc::sin(x); } );
m.def( "tan",    []( SB const& x ){ return mc::tan(x); } );
m.def( "acos",   []( SB const& x ){ return mc::acos(x); } );
m.def( "asin",   []( SB const& x ){ return mc::asin(x); } );
m.def( "atan",   []( SB const& x ){ return mc::atan(x); } );
m.def( "cosh",   []( SB const& x ){ return mc::cosh(x); } );
m.def( "sinh",   []( SB const& x ){ return mc::sinh(x); } );
m.def( "tanh",   []( SB const& x ){ return mc::tanh(x); } );
m.def( "erf",    []( SB const& x ){ return mc::erf(x); } );
m.def( "erfc",   []( SB const& x ){ return mc::erfc(x); } );
m.def( "pow",    []( SB const& x, int const n ){ return mc::pow(x,n); } );
m.def( "pow",    []( SB const& x, double const& r ){ return mc::pow(x,r); } );
m.def( "pow",    []( SB const& x, SB const& y ){ return mc::pow(x,y); } );
m.def( "pow",    []( double const& r, SB const& y ){ return mc::pow(r,y); } );
m.def( "cheb",   []( SB const& x, unsigned const n ){ return mc::cheb(x,n); } );

// Options class
//py::class_<SB::Options> pySpecbndOptions(pySpecbnd, "Options");
//pySpecbndOptions
// .def(py::init<>())
// .def(py::init<SB::Options const &>())
// .def_readwrite(
//   "HESSBND",
//   &SB::Options::HESSBND,
//   "strategy for computing spectral bounds in interval Hessian matrix"
// )
//;

// Enum for HESSBND strategy
//py::enum_<SB::Options::HESSBND_STRATEGY>(pySpecbndOptions, "HESSBND_STRATEGY")
// .value("GERSHGORIN", SB::Options::GERSHGORIN)
// .value("HERTZROHN",  SB::Options::HERTZROHN)
// .export_values()
//;
}

