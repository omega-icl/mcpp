#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "mcfunc.hpp"

namespace py = pybind11;

void mc_mcfunc( py::module_ &m )
{
m.def( "inv",   []( double const& x ){ return 1/x; } );
m.def( "sqr",   []( double const& x ){ return mc::sqr(x); } );
m.def( "sqrt",  []( double const& x ){ return std::sqrt(x); } );
m.def( "exp",   []( double const& x ){ return std::exp(x); } );
m.def( "log",   []( double const& x ){ return std::log(x); } );
m.def( "cos",   []( double const& x ){ return std::cos(x); } );
m.def( "sin",   []( double const& x ){ return std::sin(x); } );
m.def( "tan",   []( double const& x ){ return std::tan(x); } );
m.def( "acos",  []( double const& x ){ return std::acos(x); } );
m.def( "asin",  []( double const& x ){ return std::asin(x); } );
m.def( "atan",  []( double const& x ){ return std::atan(x); } );
m.def( "cosh",  []( double const& x ){ return std::cosh(x); } );
m.def( "sinh",  []( double const& x ){ return std::sinh(x); } );
m.def( "tanh",  []( double const& x ){ return std::tanh(x); } );
m.def( "fabs",  []( double const& x ){ return std::fabs(x); } );
m.def( "relu",  []( double const& x ){ return mc::relu(x); } );
m.def( "xlog",  []( double const& x ){ return mc::xlog(x); } );
m.def( "fstep", []( double const& x ){ return mc::fstep(x); } );
m.def( "bstep", []( double const& x ){ return mc::bstep(x); } );
m.def( "erf",   []( double const& x ){ return std::erf(x); } );
m.def( "erfc",  []( double const& x ){ return std::erfc(x); } );
m.def( "lmtd",  []( double const& x, double const& y ){ return mc::lmtd(x,y); } );
m.def( "rlmtd", []( double const& x, double const& y ){ return mc::rlmtd(x,y); } );
m.def( "pow",   []( double const& x, int const n ){ return std::pow(x,n); } );
m.def( "pow",   []( double const& x, double const& r ){ return std::pow(x,r); } );
m.def( "cheb",  []( double const& x, unsigned const n ){ return mc::cheb(x,n); } );
m.def( "max",   []( double const& x, double const& y ){ return mc::max(x,y); } );
m.def( "min",   []( double const& x, double const& y ){ return mc::min(x,y); } );

m.def( "machprec", [](){ return mc::machprec(); } );
}

