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

#include <chrono>
#include <fstream>
#include "ffunc.hpp" 

namespace py = pybind11;

void mc_ffunc( py::module_ &m )
{

py::class_<mc::FFNum> pyFFNum( m, "FFNum" );
pyFFNum
 .def( py::init<int const>(), "constructor for an integer constant and default constructor", py::arg("i")=0 )
 .def( py::init<double const&>(), "constructor for a real constant" )
 .def( py::init<mc::FFNum const&>(), "copy constructor for DAG constant" )
 .def_property_readonly( "val", &mc::FFNum::val, "retreive DAG constant" )
 .def( "__str__", []( mc::FFNum const& V ){ std::ostringstream Vss; Vss << V.val(); return Vss.str(); } )
 .def( "__repr__", []( mc::FFNum const& V ){ std::ostringstream Vss; Vss << V.val(); return Vss.str(); } )
       ;
 
py::class_<mc::FFVar> pyFFVar( m, "FFVar" );
pyFFVar
 .def( py::init<int const>(), "constructor for an integer constant and default constructor", py::arg("i")=0 )
 .def( py::init<double const&>(), "constructor for a real constant" )
 .def( py::init<mc::FFBase*, std::string const&>(), "constructor for DAG variable", py::arg("dag"), py::arg("name")="" )
 .def( py::init<mc::FFVar const&>(), "copy constructor for DAG variable" )
 .def( "set", py::overload_cast<mc::FFBase*, std::string const&>(&mc::FFVar::set), "attach variable to DAG", py::arg("dag"), py::arg("name")="" )
 .def( "set", py::overload_cast<int const>(&mc::FFVar::set, py::const_), "set variable to constant integer", py::arg("i")=0 )
 .def( "set", py::overload_cast<double const&>(&mc::FFVar::set, py::const_), "set variable to constant real" )
 .def( "set", py::overload_cast<std::string const&>(&mc::FFVar::set, py::const_), "set variable name" )
 .def( "unset", &mc::FFVar::unset, "unset constness" )
 .def( "num", &mc::FFVar::num, "retreive constant value" )
 .def( "cst", &mc::FFVar::cst, "retreive constness" )
 .def( "dag", &mc::FFVar::dag, "retreive DAG" )
 .def_property_readonly( "opdef", py::overload_cast<>(&mc::FFVar::opdef, py::const_), "retreive DAG defining operation" )
 .def_property_readonly( "id", py::overload_cast<>(&mc::FFVar::id, py::const_), "retreive identifier" )
 .def( "__str__", []( mc::FFVar const& V ){ std::ostringstream Vss; Vss << V; return Vss.str(); } )
 .def( "__repr__", []( mc::FFVar const& V ){ std::ostringstream Vss; Vss << V; return Vss.str(); } )
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
 .def( "__pow__", []( mc::FFVar const& V, int const n ){ return mc::pow(V,n); } )
 .def( "__pow__", []( mc::FFVar const& V, double const& r ){ return mc::pow(V,r); } )
 .def( "__pow__", []( mc::FFVar const& V, mc::FFVar const& W ){ return mc::pow(V,W); } )
;

m.def( "inv",   []( mc::FFVar const& V ){ return mc::inv(V); } );
m.def( "sqr",   []( mc::FFVar const& V ){ return mc::sqr(V); } );
m.def( "sqrt",  []( mc::FFVar const& V ){ return mc::sqrt(V); } );
m.def( "exp",   []( mc::FFVar const& V ){ return mc::exp(V); } );
m.def( "log",   []( mc::FFVar const& V ){ return mc::log(V); } );
m.def( "cos",   []( mc::FFVar const& V ){ return mc::cos(V); } );
m.def( "sin",   []( mc::FFVar const& V ){ return mc::sin(V); } );
m.def( "tan",   []( mc::FFVar const& V ){ return mc::tan(V); } );
m.def( "acos",  []( mc::FFVar const& V ){ return mc::acos(V); } );
m.def( "asin",  []( mc::FFVar const& V ){ return mc::asin(V); } );
m.def( "atan",  []( mc::FFVar const& V ){ return mc::atan(V); } );
m.def( "cosh",  []( mc::FFVar const& V ){ return mc::cosh(V); } );
m.def( "sinh",  []( mc::FFVar const& V ){ return mc::sinh(V); } );
m.def( "tanh",  []( mc::FFVar const& V ){ return mc::tanh(V); } );
m.def( "fabs",  []( mc::FFVar const& V ){ return mc::fabs(V); } );
m.def( "relu",  []( mc::FFVar const& V ){ return mc::max(V,0.); } );
m.def( "xlog",  []( mc::FFVar const& V ){ return mc::xlog(V); } );
m.def( "fstep", []( mc::FFVar const& V ){ return mc::xlog(V); } );
m.def( "bstep", []( mc::FFVar const& V ){ return mc::xlog(V); } );
m.def( "erf",   []( mc::FFVar const& V ){ return mc::erf(V); } );
m.def( "erfc",  []( mc::FFVar const& V ){ return mc::erfc(V); } );
m.def( "pow",   []( I const& x, int const n ){ return mc::Op<I>::pow(x,n); } );
m.def( "pow",   []( I const& x, double const& r ){ return mc::Op<I>::pow(x,r); } );
m.def( "pow",   []( I const& x, I const& y ){ return mc::Op<I>::pow(x,y); } );
m.def( "pow",   []( double const& r, I const& y ){ return mc::Op<I>::pow(r,y); } );
m.def( "pow",   []( mc::FFVar const& V, int const n ){ return mc::pow(V,n); } );
m.def( "pow",   []( mc::FFVar const& V, double const& r ){ return mc::pow(V,r); } );
m.def( "pow",   []( mc::FFVar const& V, mc::FFVar const& W ){ return mc::pow(V,W); } );
m.def( "pow",   []( double const& r, mc::FFVar const& W ){ return mc::pow(r,W); } );
m.def( "cheb",  []( mc::FFVar const& V, unsigned const n ){ return mc::cheb(V,n); } );
m.def( "max",   []( mc::FFVar const& V, mc::FFVar const& W ){ return mc::min(V,W); } );
m.def( "min",   []( mc::FFVar const& V, mc::FFVar const& W ){ return mc::max(V,W); } );
m.def( "inter", []( mc::FFVar const& V, mc::FFVar const& W ){ return mc::inter(V,W); } );

py::enum_<mc::FFVar::TYPE>(pyFFVar, "TYPE")
 .value("VAR",   mc::FFVar::TYPE::VAR)
 .value("AUX",   mc::FFVar::TYPE::AUX)
 .value("CINT",  mc::FFVar::TYPE::CINT)
 .value("CREAL", mc::FFVar::TYPE::CREAL)
 .export_values()
;

py::class_<mc::FFOp> pyFFOp( m, "FFOp" );
pyFFOp
 .def_readwrite( "type", &mc::FFOp::type, "retreive operation type" )
 .def_readwrite( "varin", &mc::FFOp::varin )
 .def_readwrite( "varout", &mc::FFOp::varout )
 .def( "name", &mc::FFOp::name, "retreive operation name" )
 .def( "__str__", []( mc::FFOp const& O ){ std::ostringstream Oss; Oss << O; return Oss.str(); } )
 .def( "__repr__", []( mc::FFOp const& O ){ std::ostringstream Oss; Oss << O; return Oss.str(); } )
;

py::enum_<mc::FFOp::TYPE>(pyFFOp, "TYPE")
 .value("CNST",   mc::FFOp::TYPE::CNST)
 .value("VAR",    mc::FFOp::TYPE::VAR)
 .value("PLUS",   mc::FFOp::TYPE::PLUS)
 .value("SHIFT",  mc::FFOp::TYPE::SHIFT)
 .value("NEG",    mc::FFOp::TYPE::NEG)
 .value("MINUS",  mc::FFOp::TYPE::MINUS)
 .value("TIMES",  mc::FFOp::TYPE::TIMES)
 .value("SCALE",  mc::FFOp::TYPE::SCALE)
 .value("DIV",    mc::FFOp::TYPE::DIV)
 .value("INV",    mc::FFOp::TYPE::INV)
 .value("PROD",   mc::FFOp::TYPE::PROD)
 .value("IPOW",   mc::FFOp::TYPE::IPOW)
 .value("DPOW",   mc::FFOp::TYPE::DPOW)
 .value("CHEB",   mc::FFOp::TYPE::CHEB)
 .value("SQR",    mc::FFOp::TYPE::SQR)
 .value("SQRT",   mc::FFOp::TYPE::SQRT)
 .value("EXP",    mc::FFOp::TYPE::EXP)
 .value("LOG",    mc::FFOp::TYPE::LOG)
 .value("XLOG",   mc::FFOp::TYPE::XLOG)
 .value("SIN",    mc::FFOp::TYPE::SIN)
 .value("COS",    mc::FFOp::TYPE::COS)
 .value("TAN",    mc::FFOp::TYPE::TAN)
 .value("ASIN",   mc::FFOp::TYPE::ASIN)
 .value("ACOS",   mc::FFOp::TYPE::ACOS)
 .value("ATAN",   mc::FFOp::TYPE::ATAN)
 .value("SINH",   mc::FFOp::TYPE::SINH)
 .value("COSH",   mc::FFOp::TYPE::COSH)
 .value("TANH",   mc::FFOp::TYPE::TANH)
 .value("ERF",    mc::FFOp::TYPE::ERF)
 .value("FABS",   mc::FFOp::TYPE::FABS)
 .value("FSTEP",  mc::FFOp::TYPE::FSTEP)
 .value("MINF",   mc::FFOp::TYPE::MINF)
 .value("MAXF",   mc::FFOp::TYPE::MAXF)
 .value("INTER",  mc::FFOp::TYPE::INTER)
 .value("EXTERN", mc::FFOp::TYPE::EXTERN)
 .export_values()
;

py::class_<mc::FFSubgraph> pyFFSubgraph( m, "FFSubgraph" );
pyFFSubgraph
 .def( py::init<>() )
 .def( py::init< mc::FFSubgraph const&>() )
 .def( "clear", &mc::FFSubgraph::clear, "clear subgraph" )
 .def_readonly( "len_tap", &mc::FFSubgraph::len_tap )
 .def_readonly( "len_wrk", &mc::FFSubgraph::len_wrk )
;

py::class_<mc::FFBase> pyFFBase( m, "FFBase" );
pyFFBase
 .def( py::init<>() )
 .def( "clear",
       &mc::FFBase::clear,
       "clear graph"
 )
 .def( "subgraph",
       []( mc::FFBase& G, std::vector<mc::FFVar const*> const& V )
       {
         return G.subgraph( V );
       },
       "create subgraph"
 )
 .def( "output",
       []( mc::FFBase& G, std::vector<mc::FFVar const*> const& V, std::string const& S )
       {
         return mc::FFBase::output( G.subgraph( V ), S );
       },
       py::arg("vDep"),
       py::arg("str") = "",
       "output subgraph"
 )
 .def( "output",
       []( mc::FFBase const& G, mc::FFSubgraph const& SG, std::string const& S )
       {
         mc::FFBase::output( SG, S );
       },
       py::arg("sgDep"),
       py::arg("str") = "",
       "output subgraph"
 )
 .def( "dot_script",
       []( mc::FFBase const& G, std::vector<mc::FFVar const*> const& V, std::string const& fname )
       {
         if( fname == "" ) return G.dot_script( V );
         std::ofstream ofs( fname );
         return G.dot_script( V, ofs );
       },
       "output dot script"
 )
 .def( "__str__",
       []( mc::FFBase const& G )
       {
         std::ostringstream Gss;
         Gss << G;
         return Gss.str();
       }
 )
 .def( "__repr__",
       []( mc::FFBase const& G )
       {
         std::ostringstream Gss;
         Gss << G;
         return Gss.str();
       }
 )
;

py::class_<mc::FFGraph, mc::FFBase> pyFFGraph( m, "FFGraph" );
pyFFGraph
 .def( py::init<>() )
 .def_readwrite( 
   "options",
   &mc::FFGraph::options
 )
 .def( "fdiff",
       []( mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep )
       {
         return G.SDFAD( vDep, vIndep, std::vector<mc::FFVar const*>() );
       },
       py::return_value_policy::reference_internal,
       "apply forward differentiation"
     )
 .def( "fdiff",
       []( mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep,
           std::vector<mc::FFVar const*> const& vDir )
       {
         return G.SDFAD( vDep, vIndep, vDir ); 
       },
       py::return_value_policy::reference_internal,
       "apply directional forward differentiation"
     )
 .def( "bdiff",
       []( mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vIndep )
       {
         return G.SDBAD( vDep, std::vector<mc::FFVar const*>(), vIndep );
       },
       py::return_value_policy::reference_internal,
       "apply backward differentiation"
     )
 .def( "bdiff",
       []( mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep, std::vector<mc::FFVar const*> const& vDir,
           std::vector<mc::FFVar const*> const& vIndep )
       {
         return G.SDBAD( vDep, vDir, vIndep );
       },
       py::return_value_policy::reference_internal,
       "apply backward differentiation"
     )
 .def( "tdiff",
       []( mc::FFGraph& G, unsigned int const ordermax, std::vector<mc::FFVar const*> const& vDep,
           std::vector<mc::FFVar const*> const& vVar, mc::FFVar const* const pIndep )
       {
         return G.TAD( ordermax, vDep, vVar, pIndep );
       },
       py::return_value_policy::reference_internal,
       "apply Taylor expansion"
     )
 .def( "compose",
       []( mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDepOut, std::vector< std::pair<mc::FFVar const*,
           mc::FFVar const*> > const& vDepIn )
       {
         return G.compose( vDepOut, vDepIn );
       },
       "apply compostion"
     )
 .def( "__str__",
       []( mc::FFGraph const& G ){ std::ostringstream Gss; Gss << G; return Gss.str(); }
     )
 .def( "__repr__",
       []( mc::FFGraph const& G ){ std::ostringstream Gss; Gss << G; return Gss.str(); }
     )
 .def( "eval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar, std::vector<double> const& DVar )
       {
         size_t const nDep = vDep.size();
         std::vector<double> DDep( nDep );
         G.eval( SgDep, vDep, DDep, vVar, DVar );
         return DDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in double arithmetic"
     )
 .def( "eval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep, std::vector<mc::FFVar> const& vVar,
           std::vector<double> const& DVar )
       {
         size_t const nDep = vDep.size();
         std::vector<double> DDep( nDep );
         G.eval( vDep, DDep, vVar, DVar );
         return DDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in double arithmetic"
     )
 .def( "eval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar, std::vector<I> const& IVar )
       {
         size_t const nDep = vDep.size();
         std::vector<I> IDep( nDep );
         G.eval( SgDep, vDep, IDep, vVar, IVar );
         return IDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in interval arithmetic"
     )
 .def( "eval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep, std::vector<mc::FFVar> const& vVar,
           std::vector<I> const& IVar )
       {
         size_t const nDep = vDep.size();
         std::vector<I> IDep( nDep );
         G.eval( vDep, IDep, vVar, IVar );
         return IDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in interval arithmetic"
     )
 .def( "eval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar, std::vector<MC> const& MCVar )
       {
         size_t const nDep = vDep.size();
         std::vector<MC> MCDep( nDep );
         G.eval( SgDep, vDep, MCDep, vVar, MCVar );
         return MCDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in McCormick arithmetic"
     )
 .def( "eval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep, std::vector<mc::FFVar> const& vVar,
           std::vector<MC> const& MCVar )
       {
         size_t const nDep = vDep.size();
         std::vector<MC> MCDep( nDep );
         G.eval( vDep, MCDep, vVar, MCVar );
         return MCDep;
       },
       py::return_value_policy::take_ownership,
       "evaluate subgraph in McCormick arithmetic"
     )
 .def( "reval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<I>& IDep, std::vector<mc::FFVar> const& vVar, std::vector<I>& IVar,
           I const& InfVal, unsigned const MaxPass, double const& ThresPass )
       {
         G.reval( SgDep, vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass );
         return std::make_pair(IVar,IDep);
       },
       py::arg("SgDep"),
       py::arg("vDep"),
       py::arg("IDep"),
       py::arg("vVar"),
       py::arg("IVar"),
       py::arg("InfVal"),
       py::arg("MaxPass")   = 5,
       py::arg("ThresPass") = 0e0,
       //py::return_value_policy::take_ownership,
       "propagate dependent ranges through subgraph in interval arithmetic"
     )
 .def( "reval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep, std::vector<I>& IDep,
           std::vector<mc::FFVar> const& vVar, std::vector<I>& IVar, I const& InfVal,
           unsigned const MaxPass, double const& ThresPass )
       {
         G.reval( vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass );
         return std::make_pair(IVar,IDep);
       },
       py::arg("vDep"),
       py::arg("IDep"),
       py::arg("vVar"),
       py::arg("IVar"),
       py::arg("InfVal"),
       py::arg("MaxPass")   = 5,
       py::arg("ThresPass") = 0e0,
       //py::return_value_policy::take_ownership,
       "propagate dependent ranges through subgraph in interval arithmetic"
     )
 .def( "veval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<double>>& DVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = DVar1.size();
         std::vector<std::vector<double>> DDep( nSam, std::vector<double>( nDep ) );
         G.veval( SgDep, vDep, DDep, vVar1, DVar1, vVar2, DVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return DDep;
       },
       py::arg("SgDep"),
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("DVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("DVar2") = std::vector<double>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in double arithmetic for multiple scenarios"
     )
 .def( "veval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<double>>& DVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = DVar1.size();
         std::vector<std::vector<double>> DDep( nSam, std::vector<double>( nDep ) );
         G.veval( vDep, DDep, vVar1, DVar1, vVar2, DVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return DDep;
       },
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("DVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("DVar2") = std::vector<double>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in double arithmetic for multiple scenarios"
     )
 .def( "veval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<I>>& IVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = IVar1.size();
         std::vector<std::vector<I>> IDep( nSam, std::vector<I>( nDep ) );
         G.veval( SgDep, vDep, IDep, vVar1, IVar1, vVar2, IVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return IDep;
       },
       py::arg("SgDep"),
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("IVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("IVar2") = std::vector<I>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in interval arithmetic for multiple scenarios"
     )
 .def( "veval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<I>>& IVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = IVar1.size();
         std::vector<std::vector<I>> IDep( nSam, std::vector<I>( nDep ) );
         G.veval( vDep, IDep, vVar1, IVar1, vVar2, IVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return IDep;
       },
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("IVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("IVar2") = std::vector<I>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in interval arithmetic for multiple scenarios"
     )
 .def( "veval",
       []( mc::FFGraph& G, mc::FFSubgraph& SgDep, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<MC>>& MCVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = MCVar1.size();
         std::vector<std::vector<MC>> MCDep( nSam, std::vector<MC>( nDep ) );
         G.veval( SgDep, vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return MCDep;
       },
       py::arg("SgDep"),
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("MCVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("MCVar2") = std::vector<MC>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in McCormick arithmetic for multiple scenarios"
     )
 .def( "veval",
       []( mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
           std::vector<mc::FFVar> const& vVar1, std::vector<std::vector<MC>>& MCVar1,
           std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2, bool const walltime )
       {
         auto starttime = std::chrono::system_clock::now();
         size_t const nDep = vDep.size(), nSam = MCVar1.size();
         std::vector<std::vector<MC>> MCDep( nSam, std::vector<MC>( nDep ) );
         G.veval( vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2 );
         if( walltime ){
           auto wtime = std::chrono::duration_cast<std::chrono::microseconds>( std::chrono::system_clock::now() - starttime );
           std::cerr << "vectorized DAG evaluation on " << G.options.MAXTHREAD << " threads: " << wtime.count()*1e-6 << " sec\n";
         }
         return MCDep;
       },
       py::arg("vDep"),
       py::arg("vVar1"),
       py::arg("MCVar1"),
       py::arg("vVar2") = std::vector<mc::FFVar>(),
       py::arg("MCVar2") = std::vector<MC>(),
       py::arg("walltime") = false,
       py::return_value_policy::take_ownership,
       "evaluate subgraph in McCormick arithmetic for multiple scenarios"
     )
;

py::class_<mc::FFGraph::Options> pyFFGraphOptions( pyFFGraph, "FFGraph.Options" );
pyFFGraphOptions
 .def( py::init<>() )
 .def( py::init<mc::FFGraph::Options const&>() )
 .def_readwrite( "DETECTSIGNOM",   &mc::FFGraph::Options::DETECTSIGNOM,   "Whether to detect signomial terms as exp(d.log(x)) and handle them as x^d signomial terms [Default: True]" )
 .def_readwrite( "USEMOVE",        &mc::FFGraph::Options::USEMOVE,        "Whether to enable the move semantic during DAG evaluation [Default: False]" )
 .def_readwrite( "MAXTHREAD",      &mc::FFGraph::Options::MAXTHREAD,      "Maximum number of threads in vectorized DAG evaluation [Default: 1]" )
;
}

