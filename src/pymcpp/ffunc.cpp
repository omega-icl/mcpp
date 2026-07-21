#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#ifdef MC__USE_PROFIL
#include "mcprofil.hpp"
typedef INTERVAL I;
#else
#ifdef MC__USE_FILIB
#include "mcfilib.hpp"
typedef filib::interval<double, filib::native_switched, filib::i_mode_extended>
    I;
#else
#ifdef MC__USE_BOOST
#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<
    boost::numeric::interval_lib::rounded_transc_opp<double>>
    T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round, T_boost_check>
    T_boost_policy;
typedef boost::numeric::interval<double, T_boost_policy> I;
#else
#include "interval.hpp"
typedef mc::Interval I;
#endif
#endif
#endif

#include "mccormick.hpp"
typedef mc::McCormick<I> MC;

#include "specbnd.hpp"
typedef mc::Specbnd<I> SB;

#include "tmodel.hpp"
typedef mc::TVar<I> T;

#include "cmodel.hpp"
typedef mc::CVar<I> C;

#include "scmodel.hpp"
typedef mc::SCVar<I> SC;

#include "pwlu.hpp"
#include "supmodel.hpp"
typedef mc::SupVar<mc::PWLU> PWLS;
#include "pwcu.hpp"
typedef mc::SupVar<mc::PWCU> PWCS;

#include "polimage.hpp"
typedef mc::PolVar<I> PV;

#include "ellimage.hpp"
typedef mc::EllVar<I> EV;

#include "ffdep.hpp"
typedef mc::FFDep FD;

#include "ffinv.hpp"
typedef mc::FFInv FI;

#include <chrono>
#include <fstream>

#include "ffexpr.hpp"
#include "ffunc.hpp"

namespace py = pybind11;

void
mc_ffunc(py::module_& m)
{
  py::class_<mc::FFNum> pyFFNum(m, "FFNum");
  pyFFNum
      .def(py::init<int const>(),
           "constructor for an integer constant and default constructor",
           py::arg("i") = 0)
      .def(py::init<double const&>(), "constructor for a real constant")
      .def(py::init<mc::FFNum const&>(), "copy constructor for DAG constant")
      .def_property_readonly("val", &mc::FFNum::val, "retreive DAG constant")
      .def("__str__",
           [](mc::FFNum const& V)
           {
             std::ostringstream Vss;
             Vss << V.val();
             return Vss.str();
           })
      .def("__repr__",
           [](mc::FFNum const& V)
           {
             std::ostringstream Vss;
             Vss << V.val();
             return Vss.str();
           });

  py::class_<mc::FFVar> pyFFVar(m, "FFVar");
  pyFFVar
      .def(py::init<int const>(),
           "constructor for an integer constant and default constructor",
           py::arg("i") = 0)
      .def(py::init<double const&>(), "constructor for a real constant")
      .def(py::init<mc::FFBase*, std::string const&>(),
           "constructor for DAG variable", py::arg("dag"), py::arg("name") = "")
      .def(py::init<mc::FFVar const&>(), "copy constructor for DAG variable")
      .def("set",
           py::overload_cast<mc::FFBase*, std::string const&>(&mc::FFVar::set),
           "attach variable to DAG", py::arg("dag"), py::arg("name") = "")
      .def("set", py::overload_cast<int const>(&mc::FFVar::set, py::const_),
           "set variable to constant integer", py::arg("i") = 0)
      .def("set", py::overload_cast<double const&>(&mc::FFVar::set, py::const_),
           "set variable to constant real")
      .def("set",
           py::overload_cast<std::string const&>(&mc::FFVar::set, py::const_),
           "set variable name")
      .def("unset", &mc::FFVar::unset, "unset constness")
      .def("num", &mc::FFVar::num, "retreive constant value")
      .def("cst", &mc::FFVar::cst, "retreive constness")
      .def("dag", &mc::FFVar::dag, "retreive DAG")
      .def("str",
           [](mc::FFVar const& V) { return mc::FFExpr::dep(V).ostr().str(); })
      .def_property_readonly("opdef",
                             py::overload_cast<>(&mc::FFVar::opdef, py::const_),
                             "retreive DAG defining operation")
      .def_property_readonly("id",
                             py::overload_cast<>(&mc::FFVar::id, py::const_),
                             "retreive identifier")
      .def("__str__",
           [](mc::FFVar const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
           })
      .def("__repr__",
           [](mc::FFVar const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
           })
      .def(+py::self)
      .def(py::self += double())
      .def(py::self += py::self)
      .def(double() + py::self)
      .def(py::self + double())
      .def(py::self + py::self)
      .def(-py::self)
      .def(py::self -= double())
      .def(py::self -= py::self)
      .def(double() - py::self)
      .def(py::self - double())
      .def(py::self - py::self)
      .def(py::self *= double())
      .def(py::self *= py::self)
      .def(double() * py::self)
      .def(py::self * double())
      .def(py::self * py::self)
      .def(py::self /= double())
      .def(py::self /= py::self)
      .def(double() / py::self)
      .def(py::self / double())
      .def(py::self / py::self)
      .def("__pow__",
           [](mc::FFVar const& V, int const n) { return mc::pow(V, n); })
      .def("__pow__",
           [](mc::FFVar const& V, double const& r) { return mc::pow(V, r); })
      .def("__pow__", [](mc::FFVar const& V, mc::FFVar const& W)
           { return mc::pow(V, W); });

  m.def("inv", [](mc::FFVar const& V) { return mc::inv(V); });
  m.def("sqr", [](mc::FFVar const& V) { return mc::sqr(V); });
  m.def("sqrt", [](mc::FFVar const& V) { return mc::sqrt(V); });
  m.def("exp", [](mc::FFVar const& V) { return mc::exp(V); });
  m.def("log", [](mc::FFVar const& V) { return mc::log(V); });
  m.def("cos", [](mc::FFVar const& V) { return mc::cos(V); });
  m.def("sin", [](mc::FFVar const& V) { return mc::sin(V); });
  m.def("tan", [](mc::FFVar const& V) { return mc::tan(V); });
  m.def("acos", [](mc::FFVar const& V) { return mc::acos(V); });
  m.def("asin", [](mc::FFVar const& V) { return mc::asin(V); });
  m.def("atan", [](mc::FFVar const& V) { return mc::atan(V); });
  m.def("cosh", [](mc::FFVar const& V) { return mc::cosh(V); });
  m.def("sinh", [](mc::FFVar const& V) { return mc::sinh(V); });
  m.def("tanh", [](mc::FFVar const& V) { return mc::tanh(V); });
  m.def("fabs", [](mc::FFVar const& V) { return mc::fabs(V); });
  m.def("relu", [](mc::FFVar const& V) { return mc::max(V, 0.); });
  m.def("xlog", [](mc::FFVar const& V) { return mc::xlog(V); });
  m.def("fstep", [](mc::FFVar const& V) { return mc::fstep(V); });
  m.def("bstep", [](mc::FFVar const& V) { return mc::bstep(V); });
  m.def("erf", [](mc::FFVar const& V) { return mc::erf(V); });
  m.def("erfc", [](mc::FFVar const& V) { return mc::erfc(V); });
  m.def("pow", [](I const& x, int const n) { return mc::Op<I>::pow(x, n); });
  m.def("pow",
        [](I const& x, double const& r) { return mc::Op<I>::pow(x, r); });
  m.def("pow", [](I const& x, I const& y) { return mc::Op<I>::pow(x, y); });
  m.def("pow",
        [](double const& r, I const& y) { return mc::Op<I>::pow(r, y); });
  m.def("pow", [](mc::FFVar const& V, int const n) { return mc::pow(V, n); });
  m.def("pow",
        [](mc::FFVar const& V, double const& r) { return mc::pow(V, r); });
  m.def("pow",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::pow(V, W); });
  m.def("pow",
        [](double const& r, mc::FFVar const& W) { return mc::pow(r, W); });
  m.def("cheb",
        [](mc::FFVar const& V, unsigned const n) { return mc::cheb(V, n); });
  m.def("max",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::max(V, W); });
  m.def("min",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::min(V, W); });
  m.def("inter",
        [](mc::FFVar const& V, mc::FFVar const& W) { return mc::inter(V, W); });

  py::enum_<mc::FFVar::TYPE>(pyFFVar, "TYPE")
      .value("VAR", mc::FFVar::TYPE::VAR)
      .value("AUX", mc::FFVar::TYPE::AUX)
      .value("CINT", mc::FFVar::TYPE::CINT)
      .value("CREAL", mc::FFVar::TYPE::CREAL)
      .export_values();

  py::class_<mc::FFOp> pyFFOp(m, "FFOp");
  pyFFOp.def_readwrite("type", &mc::FFOp::type, "retreive operation type")
      .def_readwrite("varin", &mc::FFOp::varin)
      .def_readwrite("varout", &mc::FFOp::varout)
      .def("name", &mc::FFOp::name, "retreive operation name")
      .def("__str__",
           [](mc::FFOp const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           })
      .def("__repr__",
           [](mc::FFOp const& O)
           {
             std::ostringstream Oss;
             Oss << O;
             return Oss.str();
           });

  py::enum_<mc::FFOp::TYPE>(pyFFOp, "TYPE")
      .value("CNST", mc::FFOp::TYPE::CNST)
      .value("VAR", mc::FFOp::TYPE::VAR)
      .value("PLUS", mc::FFOp::TYPE::PLUS)
      .value("SHIFT", mc::FFOp::TYPE::SHIFT)
      .value("NEG", mc::FFOp::TYPE::NEG)
      .value("MINUS", mc::FFOp::TYPE::MINUS)
      .value("TIMES", mc::FFOp::TYPE::TIMES)
      .value("SCALE", mc::FFOp::TYPE::SCALE)
      .value("DIV", mc::FFOp::TYPE::DIV)
      .value("INV", mc::FFOp::TYPE::INV)
      .value("PROD", mc::FFOp::TYPE::PROD)
      .value("IPOW", mc::FFOp::TYPE::IPOW)
      .value("DPOW", mc::FFOp::TYPE::DPOW)
      .value("CHEB", mc::FFOp::TYPE::CHEB)
      .value("SQR", mc::FFOp::TYPE::SQR)
      .value("SQRT", mc::FFOp::TYPE::SQRT)
      .value("EXP", mc::FFOp::TYPE::EXP)
      .value("LOG", mc::FFOp::TYPE::LOG)
      .value("XLOG", mc::FFOp::TYPE::XLOG)
      .value("SIN", mc::FFOp::TYPE::SIN)
      .value("COS", mc::FFOp::TYPE::COS)
      .value("TAN", mc::FFOp::TYPE::TAN)
      .value("ASIN", mc::FFOp::TYPE::ASIN)
      .value("ACOS", mc::FFOp::TYPE::ACOS)
      .value("ATAN", mc::FFOp::TYPE::ATAN)
      .value("SINH", mc::FFOp::TYPE::SINH)
      .value("COSH", mc::FFOp::TYPE::COSH)
      .value("TANH", mc::FFOp::TYPE::TANH)
      .value("ERF", mc::FFOp::TYPE::ERF)
      .value("FABS", mc::FFOp::TYPE::FABS)
      .value("FSTEP", mc::FFOp::TYPE::FSTEP)
      .value("MINF", mc::FFOp::TYPE::MINF)
      .value("MAXF", mc::FFOp::TYPE::MAXF)
      .value("INTER", mc::FFOp::TYPE::INTER)
      .value("EXTERN", mc::FFOp::TYPE::EXTERN)
      .export_values();

  py::class_<mc::FFSubgraph> pyFFSubgraph(m, "FFSubgraph");
  pyFFSubgraph.def(py::init<>())
      .def(py::init<mc::FFSubgraph const&>())
      .def("clear", &mc::FFSubgraph::clear, "clear subgraph")
      .def_readonly("len_tap", &mc::FFSubgraph::len_tap)
      .def_readonly("len_wrk", &mc::FFSubgraph::len_wrk);

  py::class_<mc::FFBase> pyFFBase(m, "FFBase");
  pyFFBase.def(py::init<>())
      .def(
          "add_var", [](mc::FFBase& G, std::string const& name)
          { return G.add_var(name); }, py::arg("name") = "",
          "add variable to graph")
      .def(
          "add_vars", [](mc::FFBase& G, size_t dim, std::string const& name)
          { return G.add_vars(dim, name); }, py::arg("dim"),
          py::arg("name") = "", "add variables to graph")
      .def("clear", &mc::FFBase::clear, "clear graph")
      .def(
          "subgraph", [](mc::FFBase& G, std::vector<mc::FFVar const*> const& V)
          { return G.subgraph(V); }, "create subgraph")
      .def(
          "output",
          [](mc::FFBase& G, std::vector<mc::FFVar const*> const& V,
             std::string const& S) { mc::FFBase::output(G.subgraph(V), S); },
          py::arg("vDep"), py::arg("str") = "", "output subgraph")
      .def(
          "output",
          [](mc::FFBase const& G, mc::FFSubgraph const& SG,
             std::string const& S) { mc::FFBase::output(SG, S); },
          py::arg("sgDep"), py::arg("str") = "", "output subgraph")
      .def(
          "dot_script",
          [](mc::FFBase const& G, std::vector<mc::FFVar const*> const& V,
             std::string const& fname)
          {
            if (fname == "") return G.dot_script(V);
            std::ofstream ofs(fname);
            return G.dot_script(V, ofs);
          },
          "output dot script")
      .def("__str__",
           [](mc::FFBase const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def("__repr__",
           [](mc::FFBase const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           });

  py::class_<mc::FFGraph, mc::FFBase> pyFFGraph(m, "FFGraph");
  py::class_<mc::FFGraph::Options> pyFFGraphOptions(pyFFGraph, "Options");

  pyFFGraph.def(py::init<>())
      .def_readwrite("options", &mc::FFGraph::options)
      .def(
          "fdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDFAD(vDep, vIndep, std::vector<mc::FFVar const*>()); },
          py::return_value_policy::reference_internal,
          "apply forward differentiation")
      .def(
          "fdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep,
             std::vector<mc::FFVar const*> const& vDir)
          { return G.SDFAD(vDep, vIndep, vDir); },
          py::return_value_policy::reference_internal,
          "apply directional forward differentiation")
      .def(
          "bdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDBAD(vDep, std::vector<mc::FFVar const*>(), vIndep); },
          py::return_value_policy::reference_internal,
          "apply backward differentiation")
      .def(
          "bdiff",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vDir,
             std::vector<mc::FFVar const*> const& vIndep)
          { return G.SDBAD(vDep, vDir, vIndep); },
          py::return_value_policy::reference_internal,
          "apply backward differentiation")
      .def(
          "tdiff",
          [](mc::FFGraph& G, unsigned int const ordermax,
             std::vector<mc::FFVar const*> const& vDep,
             std::vector<mc::FFVar const*> const& vVar,
             mc::FFVar const* const pIndep)
          { return G.TAD(ordermax, vDep, vVar, pIndep); },
          py::return_value_policy::reference_internal, "apply Taylor expansion")
      .def(
          "compose",
          [](mc::FFGraph& G, std::vector<mc::FFVar const*> const& vDepOut,
             std::vector<std::pair<mc::FFVar const*, mc::FFVar const*>> const&
                 vDepIn) { return G.compose(vDepOut, vDepIn); },
          py::arg("vDepOut"), py::arg("vDepIn"),
          py::return_value_policy::reference_internal, "apply composition")
      .def(
          "insert",
          [](mc::FFGraph& G, mc::FFGraph& dag,
             std::vector<mc::FFVar> const& vDepIn,
             std::vector<mc::FFVar> vDepOut)
          {
            G.insert(&dag, vDepIn, vDepOut);
            return vDepOut;
          },
          py::arg("dag"), py::arg("vDepIn"),
          py::arg("vDepOut") = std::vector<mc::FFVar>(),
          "insert dependents from dag into current graph")
      .def(
          "substitute",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDepOut,
             std::vector<mc::FFVar> const& vAuxTarg,
             std::vector<mc::FFVar> const& vAuxSubst)
          { return G.substitute(vDepOut, vAuxTarg, vAuxSubst); },
          py::arg("vDepOut"), py::arg("vAuxTarg"), py::arg("vAuxSubst"),
          "substitute variables or auxiliaries in dependents")
      .def("__str__",
           [](mc::FFGraph const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def("__repr__",
           [](mc::FFGraph const& G)
           {
             std::ostringstream Gss;
             Gss << G;
             return Gss.str();
           })
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<double> const& DVar)
          {
            size_t const nDep = vDep.size();
            std::vector<double> DDep(nDep);
            G.eval(SgDep, vDep, DDep, vVar, DVar);
            return DDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in double arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<double> const& DVar)
          {
            size_t const nDep = vDep.size();
            std::vector<double> DDep(nDep);
            G.eval(vDep, DDep, vVar, DVar);
            return DDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in double arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I> const& IVar)
          {
            size_t const nDep = vDep.size();
            std::vector<I> IDep(nDep);
            G.eval(SgDep, vDep, IDep, vVar, IVar);
            return IDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in interval arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I> const& IVar)
          {
            size_t const nDep = vDep.size();
            std::vector<I> IDep(nDep);
            G.eval(vDep, IDep, vVar, IVar);
            return IDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in interval arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<MC> const& MCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<MC> MCDep(nDep);
            G.eval(SgDep, vDep, MCDep, vVar, MCVar);
            return MCDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in McCormick relaxation arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<MC> const& MCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<MC> MCDep(nDep);
            G.eval(vDep, MCDep, vVar, MCVar);
            return MCDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in McCormick relaxation arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SB> const& SBVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SB> SBDep(nDep);
            G.eval(SgDep, vDep, SBDep, vVar, SBVar);
            return SBDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in spectral bound arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SB> const& SBVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SB> SBDep(nDep);
            G.eval(vDep, SBDep, vVar, SBVar);
            return SBDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in spectral bound arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<T> const& TVar)
          {
            size_t const nDep = vDep.size();
            std::vector<T> TDep(nDep);
            G.eval(SgDep, vDep, TDep, vVar, TVar);
            return TDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<T> const& TVar)
          {
            size_t const nDep = vDep.size();
            std::vector<T> TDep(nDep);
            G.eval(vDep, TDep, vVar, TVar);
            return TDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<C> const& CVar)
          {
            size_t const nDep = vDep.size();
            std::vector<C> CDep(nDep);
            G.eval(SgDep, vDep, CDep, vVar, CVar);
            return CDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in dense Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<C> const& CVar)
          {
            size_t const nDep = vDep.size();
            std::vector<C> CDep(nDep);
            G.eval(vDep, CDep, vVar, CVar);
            return CDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in dense Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SC> const& SCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SC> SCDep(nDep);
            G.eval(SgDep, vDep, SCDep, vVar, SCVar);
            return SCDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<SC> const& SCVar)
          {
            size_t const nDep = vDep.size();
            std::vector<SC> SCDep(nDep);
            G.eval(vDep, SCDep, vVar, SCVar);
            return SCDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWLS> const& PWLSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWLS> PWLSDep(nDep);
            G.eval(SgDep, vDep, PWLSDep, vVar, PWLSVar);
            return PWLSDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise linear univariates")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWLS> const& PWLSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWLS> PWLSDep(nDep);
            G.eval(vDep, PWLSDep, vVar, PWLSVar);
            return PWLSDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise linear univariates")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWCS> const& PWCSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWCS> PWCSDep(nDep);
            G.eval(SgDep, vDep, PWCSDep, vVar, PWCSVar);
            return PWCSDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise constant univariates")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar,
             std::vector<PWCS> const& PWCSVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PWCS> PWCSDep(nDep);
            G.eval(vDep, PWCSDep, vVar, PWCSVar);
            return PWCSDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise constant univariates")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<PV> const& PVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PV> PVDep(nDep);
            G.eval(SgDep, vDep, PVDep, vVar, PVVar);
            return PVDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<PV> const& PVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<PV> PVDep(nDep);
            G.eval(vDep, PVDep, vVar, PVVar);
            return PVDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<EV> const& EVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<EV> EVDep(nDep);
            G.eval(SgDep, vDep, EVDep, vVar, EVVar);
            return EVDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<EV> const& EVVar)
          {
            size_t const nDep = vDep.size();
            std::vector<EV> EVDep(nDep);
            G.eval(vDep, EVDep, vVar, EVVar);
            return EVDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FD> const& FDVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FD> FDDep(nDep);
            G.eval(SgDep, vDep, FDDep, vVar, FDVar);
            return FDDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FD> const& FDVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FD> FDDep(nDep);
            G.eval(vDep, FDDep, vVar, FDVar);
            return FDDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FI> const& FIVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FI> FIDep(nDep);
            G.eval(SgDep, vDep, FIDep, vVar, FIVar);
            return FIDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "eval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar, std::vector<FI> const& FIVar)
          {
            size_t const nDep = vDep.size();
            std::vector<FI> FIDep(nDep);
            G.eval(vDep, FIDep, vVar, FIVar);
            return FIDep;
          },
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic")
      .def(
          "reval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep, std::vector<I>& IDep,
             std::vector<mc::FFVar> const& vVar, std::vector<I>& IVar,
             I const& InfVal, unsigned const MaxPass, double const& ThresPass)
          {
            G.reval(SgDep, vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass);
            return std::make_pair(IVar, IDep);
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("IDep"), py::arg("vVar"),
          py::arg("IVar"), py::arg("InfVal"), py::arg("MaxPass") = 5,
          py::arg("ThresPass") = 0e0,
          // py::return_value_policy::take_ownership,
          "propagate dependent ranges through subgraph in interval arithmetic")
      .def(
          "reval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<I>& IDep, std::vector<mc::FFVar> const& vVar,
             std::vector<I>& IVar, I const& InfVal, unsigned const MaxPass,
             double const& ThresPass)
          {
            G.reval(vDep, IDep, vVar, IVar, InfVal, MaxPass, ThresPass);
            return std::make_pair(IVar, IDep);
          },
          py::arg("vDep"), py::arg("IDep"), py::arg("vVar"), py::arg("IVar"),
          py::arg("InfVal"), py::arg("MaxPass") = 5, py::arg("ThresPass") = 0e0,
          // py::return_value_policy::take_ownership,
          "propagate dependent ranges through subgraph in interval arithmetic")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<double>>& DVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = DVar1.size();
            std::vector<std::vector<double>> DDep(nSam,
                                                  std::vector<double>(nDep));
            G.veval(SgDep, vDep, DDep, vVar1, DVar1, vVar2, DVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return DDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("DVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("DVar2") = std::vector<double>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in double arithmetic for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<double>>& DVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<double>& DVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = DVar1.size();
            std::vector<std::vector<double>> DDep(nSam,
                                                  std::vector<double>(nDep));
            G.veval(vDep, DDep, vVar1, DVar1, vVar2, DVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return DDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("DVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("DVar2") = std::vector<double>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in double arithmetic for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<I>>& IVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = IVar1.size();
            std::vector<std::vector<I>> IDep(nSam, std::vector<I>(nDep));
            G.veval(SgDep, vDep, IDep, vVar1, IVar1, vVar2, IVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return IDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("IVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("IVar2") = std::vector<I>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in interval arithmetic for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<I>>& IVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<I>& IVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = IVar1.size();
            std::vector<std::vector<I>> IDep(nSam, std::vector<I>(nDep));
            G.veval(vDep, IDep, vVar1, IVar1, vVar2, IVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return IDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("IVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("IVar2") = std::vector<I>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in interval arithmetic for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<MC>>& MCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = MCVar1.size();
            std::vector<std::vector<MC>> MCDep(nSam, std::vector<MC>(nDep));
            G.veval(SgDep, vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return MCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("MCVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("MCVar2") = std::vector<MC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in McCormick relaxation arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<MC>>& MCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<MC>& MCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = MCVar1.size();
            std::vector<std::vector<MC>> MCDep(nSam, std::vector<MC>(nDep));
            G.veval(vDep, MCDep, vVar1, MCVar1, vVar2, MCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return MCDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("MCVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("MCVar2") = std::vector<MC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in McCormick relaxation arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SB>>& SBVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SB>& SBVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SBVar1.size();
            std::vector<std::vector<SB>> SBDep(nSam, std::vector<SB>(nDep));
            G.veval(SgDep, vDep, SBDep, vVar1, SBVar1, vVar2, SBVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SBDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("SBVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("SBVar2") = std::vector<SB>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in spectral bound arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SB>>& SBVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SB>& SBVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SBVar1.size();
            std::vector<std::vector<SB>> SBDep(nSam, std::vector<SB>(nDep));
            G.veval(vDep, SBDep, vVar1, SBVar1, vVar2, SBVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SBDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("SBVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("SBVar2") = std::vector<SB>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in spectral bound arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<T>>& TVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<T>& TVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = TVar1.size();
            std::vector<std::vector<T>> TDep(nSam, std::vector<T>(nDep));
            G.veval(SgDep, vDep, TDep, vVar1, TVar1, vVar2, TVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return TDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("TVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("TVar2") = std::vector<T>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<T>>& TVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<T>& TVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = TVar1.size();
            std::vector<std::vector<T>> TDep(nSam, std::vector<T>(nDep));
            G.veval(vDep, TDep, vVar1, TVar1, vVar2, TVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return TDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("TVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("TVar2") = std::vector<T>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<C>>& CVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<C>& CVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = CVar1.size();
            std::vector<std::vector<C>> CDep(nSam, std::vector<C>(nDep));
            G.veval(SgDep, vDep, CDep, vVar1, CVar1, vVar2, CVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return CDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"), py::arg("CVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("CVar2") = std::vector<C>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in dense Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<C>>& CVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<C>& CVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = CVar1.size();
            std::vector<std::vector<C>> CDep(nSam, std::vector<C>(nDep));
            G.veval(vDep, CDep, vVar1, CVar1, vVar2, CVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return CDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("CVar1"),
          py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("CVar2") = std::vector<C>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in dense Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SC>>& SCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SC>& SCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SCVar1.size();
            std::vector<std::vector<SC>> SCDep(nSam, std::vector<SC>(nDep));
            G.veval(SgDep, vDep, SCDep, vVar1, SCVar1, vVar2, SCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SCDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("SCVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("SCVar2") = std::vector<SC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<SC>>& SCVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<SC>& SCVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = SCVar1.size();
            std::vector<std::vector<SC>> SCDep(nSam, std::vector<SC>(nDep));
            G.veval(vDep, SCDep, vVar1, SCVar1, vVar2, SCVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return SCDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("SCVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("SCVar2") = std::vector<SC>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in sparse Chebyshev model arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWLS>>& PWLSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWLS>& PWLSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWLSVar1.size();
            std::vector<std::vector<PWLS>> PWLSDep(nSam,
                                                   std::vector<PWLS>(nDep));
            G.veval(SgDep, vDep, PWLSDep, vVar1, PWLSVar1, vVar2, PWLSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWLSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PWLSVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PWLSVar2") = std::vector<PWLS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise linear univariates for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWLS>>& PWLSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWLS>& PWLSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWLSVar1.size();
            std::vector<std::vector<PWLS>> PWLSDep(nSam,
                                                   std::vector<PWLS>(nDep));
            G.veval(vDep, PWLSDep, vVar1, PWLSVar1, vVar2, PWLSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWLSDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PWLSVar1"),
          py::arg("vVar2")    = std::vector<mc::FFVar>(),
          py::arg("PWLSVar2") = std::vector<PWLS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise linear univariates for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWCS>>& PWCSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWCS>& PWCSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWCSVar1.size();
            std::vector<std::vector<PWCS>> PWCSDep(nSam,
                                                   std::vector<PWCS>(nDep));
            G.veval(SgDep, vDep, PWCSDep, vVar1, PWCSVar1, vVar2, PWCSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWCSDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PWCSVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PWCSVar2") = std::vector<PWCS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise constant univariates for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PWCS>>& PWCSVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PWCS>& PWCSVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PWCSVar1.size();
            std::vector<std::vector<PWCS>> PWCSDep(nSam,
                                                   std::vector<PWCS>(nDep));
            G.veval(vDep, PWCSDep, vVar1, PWCSVar1, vVar2, PWCSVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PWCSDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PWCSVar1"),
          py::arg("vVar2")    = std::vector<mc::FFVar>(),
          py::arg("PWCSVar2") = std::vector<PWCS>(),
          py::arg("walltime") = false, py::return_value_policy::take_ownership,
          "evaluate subgraph in superposition relaxation arithmetic with "
          "piecewise constant univariates for multiple scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PV>>& PVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PV>& PVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PVVar1.size();
            std::vector<std::vector<PV>> PVDep(nSam, std::vector<PV>(nDep));
            G.veval(SgDep, vDep, PVDep, vVar1, PVVar1, vVar2, PVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("PVVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("PVVar2") = std::vector<PV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<PV>>& PVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<PV>& PVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = PVVar1.size();
            std::vector<std::vector<PV>> PVDep(nSam, std::vector<PV>(nDep));
            G.veval(vDep, PVDep, vVar1, PVVar1, vVar2, PVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return PVDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("PVVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("PVVar2") = std::vector<PV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<EV>>& EVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<EV>& EVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = EVVar1.size();
            std::vector<std::vector<EV>> EVDep(nSam, std::vector<EV>(nDep));
            G.veval(SgDep, vDep, EVDep, vVar1, EVVar1, vVar2, EVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return EVDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("EVVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("EVVar2") = std::vector<EV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<EV>>& EVVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<EV>& EVVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = EVVar1.size();
            std::vector<std::vector<EV>> EVDep(nSam, std::vector<EV>(nDep));
            G.veval(vDep, EVDep, vVar1, EVVar1, vVar2, EVVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return EVDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("EVVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("EVVar2") = std::vector<EV>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FD>>& FDVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FD>& FDVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FDVar1.size();
            std::vector<std::vector<FD>> FDDep(nSam, std::vector<FD>(nDep));
            G.veval(SgDep, vDep, FDDep, vVar1, FDVar1, vVar2, FDVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FDDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("FDVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("FDVar2") = std::vector<FD>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FD>>& FDVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FD>& FDVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FDVar1.size();
            std::vector<std::vector<FD>> FDDep(nSam, std::vector<FD>(nDep));
            G.veval(vDep, FDDep, vVar1, FDVar1, vVar2, FDVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FDDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("FDVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("FDVar2") = std::vector<FD>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, mc::FFSubgraph& SgDep,
             std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FI>>& FIVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FI>& FIVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FIVar1.size();
            std::vector<std::vector<FI>> FIDep(nSam, std::vector<FI>(nDep));
            G.veval(SgDep, vDep, FIDep, vVar1, FIVar1, vVar2, FIVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FIDep;
          },
          py::arg("SgDep"), py::arg("vDep"), py::arg("vVar1"),
          py::arg("FIVar1"), py::arg("vVar2") = std::vector<mc::FFVar>(),
          py::arg("FIVar2") = std::vector<FI>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios")
      .def(
          "veval",
          [](mc::FFGraph& G, std::vector<mc::FFVar> const& vDep,
             std::vector<mc::FFVar> const& vVar1,
             std::vector<std::vector<FI>>& FIVar1,
             std::vector<mc::FFVar> const& vVar2, std::vector<FI>& FIVar2,
             bool const walltime)
          {
            auto starttime    = std::chrono::system_clock::now();
            size_t const nDep = vDep.size(), nSam = FIVar1.size();
            std::vector<std::vector<FI>> FIDep(nSam, std::vector<FI>(nDep));
            G.veval(vDep, FIDep, vVar1, FIVar1, vVar2, FIVar2);
            if (walltime)
            {
              auto wtime =
                  std::chrono::duration_cast<std::chrono::microseconds>(
                      std::chrono::system_clock::now() - starttime);
              std::cerr << "vectorized DAG evaluation on "
                        << G.options.MAXTHREAD
                        << " threads: " << wtime.count() * 1e-6 << " sec\n";
            }
            return FIDep;
          },
          py::arg("vDep"), py::arg("vVar1"), py::arg("FIVar1"),
          py::arg("vVar2")  = std::vector<mc::FFVar>(),
          py::arg("FIVar2") = std::vector<FI>(), py::arg("walltime") = false,
          py::return_value_policy::take_ownership,
          "evaluate subgraph in polyhedral image arithmetic for multiple "
          "scenarios");

  pyFFGraphOptions.def(py::init<>())
      .def(py::init<mc::FFGraph::Options const&>())
      .def_readwrite("DETECTSIGNOM", &mc::FFGraph::Options::DETECTSIGNOM,
                     "Whether to detect signomial terms as exp(d.log(x)) and "
                     "handle them as x^d signomial terms [Default: True]")
      .def_readwrite("CHEBRECURS", &mc::FFGraph::Options::CHEBRECURS,
                     "Whether to intersect Chebyshev variables with their "
                     "recursive expressions [Default: False]")
      .def_readwrite("USEMOVE", &mc::FFGraph::Options::USEMOVE,
                     "Whether to enable the move semantic during DAG "
                     "evaluation [Default: False]")
      .def_readwrite("MAXTHREAD", &mc::FFGraph::Options::MAXTHREAD,
                     "Maximum number of threads in vectorized DAG evaluation "
                     "[Default: 1]");
}
