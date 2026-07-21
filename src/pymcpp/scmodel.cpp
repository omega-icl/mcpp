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
    boost::numeric::interval_lib::rounded_transc_opp<double> >
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

#include "scmodel.hpp"

namespace py = pybind11;

void
mc_scmodel(py::module& m)
{
  typedef mc::SCModel<I> SCM;
  typedef mc::SCVar<I> SCV;

  py::class_<SCM> pySCModel(m, "SCModel");

  py::class_<SCM::Options> pySCModelOptions(pySCModel, "Options");

  py::enum_<SCM::Options::BOUNDER>(pySCModelOptions, "BOUNDER")
      .value("NAIVE", SCM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", SCM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .export_values();

  pySCModel
      // Constructors
      .def(py::init<unsigned const>(), py::arg("maxord") = 3,
           "constructor for sparse Chebyshev model with maximal order 'maxord'")
      // Accessors
      .def_readwrite("options", &SCM::options)
      .def_property_readonly(
          "maxord", [](SCM const& self) { return self.maxord(); },
          "maximal degree of model")
      .def_property_readonly(
          "setvar", [](SCM const& self) { return self.setvar(); },
          "model variables")
      .def_property_readonly(
          "bndvar", [](SCM const& self) { return self.bndvar(); },
          "bounds on model variables")
      .def_property_readonly(
          "refvar", [](SCM const& self) { return self.refvar(); },
          "reference values of model variables")
      .def_property_readonly(
          "scalvar", [](SCM const& self) { return self.scalvar(); },
          "scaling factors of model variables")
      .def_property_readonly(
          "setaux", [](SCM const& self) { return self.setvar(); },
          "auxiliary variables")
      .def(
          "append_aux", [](SCM& self) { return self.append_aux(); },
          "append new auxiliary variable")
      .def(
          "reset_aux", [](SCM& self) { self.reset_aux(); },
          "reset auxiliary variables");

  py::class_<SCV> pySCVar(m, "SCVar");

  pySCVar
      // Constructors
      .def(py::init<double const&, SCM*>(), py::arg("cst"),
           py::arg("env") = nullptr, "constructor for constant scalar")
      .def(py::init<I const&, SCM*>(), py::arg("bnd"), py::arg("env") = nullptr,
           "constructor for constant bound")
      .def(py::init<SCM*, unsigned const&, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"),
           "constructor for variable with range 'rng' and index 'ndx' in model "
           "'env'")
      .def(py::init<SCV const&>(), "copy constructor")
      // Modifiers
      .def(
          "set", [](SCV& self, SCM* env) -> SCV& { return self.set(env); },
          py::return_value_policy::reference_internal, "set model environment")
      .def(
          "set", [](SCV& self, I const& rem) -> SCV& { return self.set(rem); },
          py::return_value_policy::reference_internal, "set remainder bound")
      .def(
          "set", [](SCV& self, SCV::t_poly const& coefmon, I const& rem) -> SCV&
          { return self.set(coefmon); }, py::arg("coefmon") = 0.,
          py::arg("rem") = 0., py::return_value_policy::reference_internal,
          "set sparse polynomial coefficients and monomials")
      .def(
          "set",
          [](SCV& self, SCM* env, unsigned const& ndx, I const& rng) -> SCV&
          { return self.set(env, ndx, rng); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal,
          "set variable with range 'rng' and index 'ndx' in model 'env'")
      // Accessors
      .def("env", &SCV::env, "model environment")
      .def_property_readonly(
          "nord", [](SCV const& self) { return self.nord(); }, "current degree")
      .def_property_readonly(
          "nvar", [](SCV const& self) { return self.nvar(); },
          "number of variable dependencies")
      .def_property_readonly(
          "nmon", [](SCV const& self) { return self.nmon(); },
          "number of monomials")
      .def_property_readonly(
          "coefmon", py::overload_cast<>(&SCV::coefmon, py::const_),
          "dictionary of monomials and associated coefficients")
      .def_property_readonly("ndxvar",
                             py::overload_cast<>(&SCV::ndxvar, py::const_),
                             //[]( SCV const& self )
                             //  { return self.ndxvar(); },
                             "list of participating variables")
      .def_property_readonly(
          "bound", [](SCV const& self) { return self.bound(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def_property_readonly(
          "B", [](SCV const& self) { return self.B(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def(
          "bndpol", [](SCV const& self) { return self.bndpol(); },
          "bound on multivariate polynomial using default bounder")
      .def(
          "bndpol", [](SCV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          "bound on multivariate polynomial using bounder 'type'")
      .def("bndord", &SCV::bndord, py::arg("minord"),
           "compute bound on all monomials with (total) order no less than "
           "'minord'")
      .def(
          "P", [](SCV const& self, std::map<unsigned, double> const& x)
          { return self.P(x); }, py::arg("x"),
          "evaluate polynomial part at 'x'")
      .def(
          "P", [](SCV const& self) { return self.P(); },
          "cancel remainder and return new model with same polynpomial part")
      .def_property_readonly("R", &SCV::R, "remainder bound")
      .def("center", &SCV::center, "center remainder term in model")
      .def("C", &SCV::C, "center remainder term in model")
      .def("constant", &SCV::constant, py::arg("reset") = false,
           "get coefficient of constant term and reset to zero if 'reset' is "
           "True")
      .def("linear", &SCV::linear, py::arg("id"), py::arg("reset") = false,
           "get coefficient of linear term for variable 'id' and reset to zero "
           "if 'reset' is True")
      .def(
          "lift",
          [](SCV& self, SCM& scm, double const& atol, double const& rtol)
          { return self.lift(&scm, atol, rtol); }, py::arg("scm"),
          py::arg("atol"), py::arg("rtol"),
          "lift model by appending new auxiliary variable from the remainder "
          "term, if the remainder has magnitude greater than ATOL/RTOL")
      .def(
          "project", [](SCV& self, bool const reset)
          { return self.project(reset); }, py::arg("reset") = false,
          "project model by removing monomials with any auxiliary variables "
          "participating")
      .def(
          "project", [](SCV& self, unsigned const& id)
          { return self.project(id); }, py::arg("id"),
          "project model by removing monomials with variable 'id' "
          "participating")
      .def(
          "scale", [](SCV& self, unsigned const& id, I const& dom)
          { return self.scale(id, dom); }, py::arg("id"), py::arg("dom"),
          "Scale coefficients in Chebyshev variable for the modified domain "
          "dom of variable id")
      .def(
          "rescale", [](SCV& self, std::map<unsigned, I> const& dom)
          { return self.scale(dom); }, py::arg("dom"),
          "rescale monomial coefficients in model over new domain 'dom'")
      .def("unscale", &SCV::unscale,
           "unscale monomial coefficients in model over original variable "
           "domain")
      .def("simplify", &SCV::simplify, py::arg("atol") = 0e0,
           py::arg("rtol") = 0e0, py::arg("tord") = -1,
           "simplify model by appending coefficients with magnitude less than "
           "ATOL/RTOL or order greater no less than TORD to the remainder term")
      .def(
          "to_monomial", [](SCV const& self, bool const scaled)
          { return self.to_monomial(scaled); }, py::arg("scaled") = false,
          "express coefficient map in monomial basis representation")
      .def(
          "to_monomial",
          [](SCV const& self, bool const scaled, double const& atol,
             double const& rtol, int const tord)
          { return self.to_monomial(scaled, atol, rtol, tord); },
          py::arg("scaled"), py::arg("atol"), py::arg("rtol"),
          py::arg("tord") = -1,
          "express coefficient map in monomial basis representation after "
          "removing terms with coefficients less than ATOL/RTOL or order "
          "greater than or equal to ORD, and return bound on removed terms")
      .def("__str__",
           [](SCV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](SCV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
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
      .def(double() / py::self)
      .def(py::self / double())
      .def(py::self / py::self)
      .def("__abs__", [](SCV const& m) { return mc::Op<SCV>::abs(m); })
      .def("__pow__", [](SCV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](SCV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](SCV const& m, SCV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, SCV const& m) { return mc::pow(r, m); });

  m.def("inv", [](SCV const& x) { return mc::inv(x); });
  m.def("sqr", [](SCV const& x) { return mc::sqr(x); });
  m.def("sqrt", [](SCV const& x) { return mc::sqrt(x); });
  m.def("exp", [](SCV const& x) { return mc::exp(x); });
  m.def("log", [](SCV const& x) { return mc::log(x); });
  m.def("xlog", [](SCV const& x) { return mc::xlog(x); });
  m.def("cos", [](SCV const& x) { return mc::cos(x); });
  m.def("sin", [](SCV const& x) { return mc::sin(x); });
  m.def("tan", [](SCV const& x) { return mc::tan(x); });
  m.def("acos", [](SCV const& x) { return mc::acos(x); });
  m.def("asin", [](SCV const& x) { return mc::asin(x); });
  m.def("atan", [](SCV const& x) { return mc::atan(x); });
  m.def("cosh", [](SCV const& x) { return mc::cosh(x); });
  m.def("sinh", [](SCV const& x) { return mc::sinh(x); });
  m.def("tanh", [](SCV const& x) { return mc::tanh(x); });
  m.def("erf", [](SCV const& x) { return mc::erf(x); });
  m.def("erfc", [](SCV const& x) { return mc::erfc(x); });
  m.def("pow", [](SCV const& x, int const n) { return mc::pow(x, n); });
  m.def("pow", [](SCV const& x, double const& r) { return mc::pow(x, r); });
  m.def("pow", [](SCV const& x, SCV const& y) { return mc::pow(x, y); });
  m.def("pow", [](double const& r, SCV const& y) { return mc::pow(r, y); });
  m.def("cheb", [](SCV const& x, unsigned const n) { return mc::cheb(x, n); });
  m.def("fabs", [](SCV const& x) { return mc::fabs(x); });
  m.def("hull", [](SCV const& x, SCV const& y) { return mc::hull(x, y); });
  m.def("inter", [](SCV& xy, SCV const& x, SCV const& y)
        { return mc::inter(xy, x, y); });

  // Nested class Options

  pySCModelOptions.def(py::init<>())
      .def(py::init<SCM::Options const&>())
      .def(
          "reset", [](SCM::Options& self) { self = SCM::Options(); },
          "Reset options to default")
      .def_readwrite("BASIS", &SCM::Options::BASIS,
                     "Basis representation of the monomials: 0-Monomial basis; "
                     "1-Chebyshev basis [Default: 1 (CHEB)]")
      .def_readwrite(
          "LIFT_USE", &SCM::Options::LIFT_USE,
          "Whether to lift the remainder term in nonlinear operations by "
          "introducing auxiliary variables in the model [Default: False]")
      .def_readwrite("LIFT_ATOL", &SCM::Options::LIFT_ATOL,
                     "Absolute tolerance for lifting the reminder term - only "
                     "if LIFT_USE == true [Default: 1e-10]")
      .def_readwrite("LIFT_RTOL", &SCM::Options::LIFT_RTOL,
                     "Relative tolerance for lifting the reminder term - only "
                     "if LIFT_USE == true [Default: 1e-3]")
      .def_readwrite(
          "REMEZ_USE", &SCM::Options::REMEZ_USE,
          "Whether to use the Remez algorithm for computing a minimax "
          "approximation for univariate terms [Default: True]")
      .def_readwrite(
          "REMEZ_MAXIT", &SCM::Options::REMEZ_MAXIT,
          "Maximal number of iterations in Remez algorithm for computing a "
          "minimax approximation for univariate terms [Default: 10]")
      .def_readwrite(
          "REMEZ_TOL", &SCM::Options::REMEZ_TOL,
          "Stopping tolerance in Remez algorithm for computing a minimax "
          "approximation for univariate terms [Default: 1e-5]")
      .def_readwrite("REMEZ_MIG", &SCM::Options::REMEZ_MIG,
                     "Threshold for interval width below which Remez algorithm "
                     "is not used [Default: 1e-10]")
      .def_readwrite(
          "INTERP_EXTRA", &SCM::Options::INTERP_EXTRA,
          "Extra terms in Chebyshev interpolation of univariates: 0-Chebyshev "
          "interpolation of order NORD; extra terms allow approximation of "
          "Chebyshev truncated series [Default: 0]")
      .def_readwrite(
          "INTERP_THRES", &SCM::Options::INTERP_THRES,
          "Threshold for coefficient values in Chebyshev expansion for "
          "bounding of transcendental univariates [Default: 1e2*machprec()]")
      .def_readwrite("BOUNDER_TYPE", &SCM::Options::BOUNDER_TYPE,
                     "Chebyshev model range bounder [Default: LSB]")
      .def_readwrite("MIG_USE", &SCM::Options::MIG_USE,
                     "Whether to simplify the monomial terms with small "
                     "magnitude in the model [Default: False]")
      .def_readwrite("MIG_ATOL", &SCM::Options::MIG_ATOL,
                     "Absolute tolerance for simplifying monomial terms - only "
                     "if MIG_USE == true [Default: 0]")
      .def_readwrite("MIG_RTOL", &SCM::Options::MIG_RTOL,
                     "Relative tolerance for simplifying monomial terms - only "
                     "if MIG_USE == true [Default: machprec()]")
      .def_readwrite("MIXED_IA", &SCM::Options::MIXED_IA,
                     "Whether to intersect internal bounds with underlying "
                     "bounds in the templated arithmetics [Default: False]")
      .def_readwrite(
          "REF_POLY", &SCM::Options::REF_POLY,
          "Scalar in [0,1] related to the choice of the polynomial part in the "
          "overloaded functions mc::inter and mc::hull [Default: 0]")
      .def_readwrite("DISPLAY_DIGITS", &SCM::Options::DISPLAY_DIGITS,
                     "Number of digits in output stream for Chebyshev model "
                     "coefficients [Default: 7]");


  py::enum_<SCM::Options::MONBASIS>(pySCModelOptions, "MONBASIS")
      .value("MONOM", SCM::Options::MONBASIS::MONOM, "Monomial basis")
      .value("CHEB", SCM::Options::MONBASIS::CHEB, "Chebyshev basis")
      .export_values();

  // Nested class Exceptions
  py::class_<SCM::Exceptions> pySCModelExceptions(pySCModel, "Exceptions");

  py::enum_<SCM::Exceptions::TYPE>(pySCModelExceptions, "TYPE")
      .value("DIV", SCM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", SCM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", SCM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", SCM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("DPOW", SCM::Exceptions::TYPE::DPOW,
             "Real power operation with negative numbers in range")
      .value("TAN", SCM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", SCM::Exceptions::TYPE::ACOS,
             "Cosine inverse operation with range outside [-1,1]")
      .value("ASIN", SCM::Exceptions::TYPE::ASIN,
             "Sine inverse operation with range outside [-1,1]")
      .value("COMPOSE", SCM::Exceptions::TYPE::COMPOSE,
             "Failed to compose Chebyshev variable")
      .value("INIT", SCM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", SCM::Exceptions::TYPE::INCON,
             "Inconsistent bounds with template parameter arithmetic")
      .value("MODEL", SCM::Exceptions::TYPE::MODEL,
             "Operation between variables linked to different models")
      .value("INTERNAL", SCM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", SCM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented")
      .export_values();

  pySCModelExceptions.def(py::init<SCM::Exceptions::TYPE>())
      .def("ierr", &SCM::Exceptions::ierr, "Error flag")
      .def("what", &SCM::Exceptions::what, "Error description");
}
