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

#include "sicmodel.hpp"

namespace py = pybind11;

void
mc_sicmodel(py::module& m)
{
  typedef mc::SICModel<I> SICM;
  typedef mc::SICVar<I> SICV;

  py::class_<SICM> pySICModel(m, "SICModel");

  py::class_<SICM::Options> pySICModelOptions(pySICModel, "Options");

  py::enum_<SICM::Options::BOUNDER>(pySICModelOptions, "BOUNDER")
      .value("NAIVE", SICM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", SICM::Options::BOUNDER::LSB,
             "Lin & Stadtherr range bounder")
      .value("BERNSTEIN", SICM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .export_values();

  pySICModel
      // Constructors
      .def(py::init<unsigned const>(), py::arg("maxord") = 3,
           "constructor for sparse interval Chebyshev model with maximal order "
           "'maxord'")
      // Accessors
      .def_readwrite("options", &SICM::options)
      .def_property_readonly(
          "maxord", [](SICM const& self) { return self.maxord(); },
          "maximal degree of model")
      .def_property_readonly(
          "setvar", [](SICM const& self) { return self.setvar(); },
          "model variables")
      .def_property_readonly(
          "bndvar", [](SICM const& self) { return self.bndvar(); },
          "bounds on model variables")
      .def_property_readonly(
          "refvar", [](SICM const& self) { return self.refvar(); },
          "reference values of model variables")
      .def_property_readonly(
          "scalvar", [](SICM const& self) { return self.scalvar(); },
          "scaling factors of model variables")
      .def(
          "append_aux", [](SICM& self) { return self.append_aux(); },
          "append new auxiliary variable")
      .def(
          "reset_aux", [](SICM& self) { self.reset_aux(); },
          "reset auxiliary variables");

  py::class_<SICV> pySICVar(m, "SICVar");

  pySICVar
      // Constructors
      .def(py::init<double const&, SICM*>(), py::arg("cst"),
           py::arg("env") = nullptr, "constructor for constant scalar")
      .def(py::init<I const&, SICM*>(), py::arg("bnd"),
           py::arg("env") = nullptr, "constructor for constant bound")
      .def(py::init<SICM*, unsigned const&, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"),
           "constructor for variable with range 'rng' and index 'ndx' in model "
           "'env'")
      .def(py::init<SICV const&>(), "copy constructor")
      // Modifiers
      .def(
          "set", [](SICV& self, SICM* env) -> SICV& { return self.set(env); },
          py::return_value_policy::reference_internal, "set model environment")
      .def(
          "set",
          [](SICV& self, I const& bnd) -> SICV& { return self.set(bnd); },
          py::return_value_policy::reference_internal,
          "set as interval-valued constant")
      .def(
          "set",
          [](SICV& self, SICM* env, unsigned const& ndx, I const& rng) -> SICV&
          { return self.set(env, ndx, rng); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal,
          "set variable with range 'rng' and index 'ndx' in model 'env'")
      // Accessors
      .def("env", &SICV::env, "model environment")
      .def_property_readonly(
          "nord", [](SICV const& self) { return self.nord(); },
          "current degree")
      .def_property_readonly(
          "nvar", [](SICV const& self) { return self.nvar(); },
          "number of variable dependencies")
      .def_property_readonly(
          "nmon", [](SICV const& self) { return self.nmon(); },
          "number of monomials")
      .def_property_readonly(
          "coefmon", py::overload_cast<>(&SICV::coefmon, py::const_),
          "dictionary of monomials and associated interval coefficients")
      .def_property_readonly("ndxvar",
                             py::overload_cast<>(&SICV::ndxvar, py::const_),
                             "list of participating variables")
      .def_property_readonly(
          "bound", [](SICV const& self) { return self.bound(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def_property_readonly(
          "B", [](SICV const& self) { return self.B(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def(
          "bndpol", [](SICV const& self) { return self.bndpol(); },
          "bound on the interval-valued polynomial using default bounder")
      .def("bndord", &SICV::bndord, py::arg("minord"),
           "compute bound on all monomials with (total) order no less than "
           "'minord'")
      .def(
          "P", [](SICV const& self, std::map<unsigned, double> const& x)
          { return self.P(x); }, py::arg("x"),
          "evaluate mid-point polynomial part at 'x'")
      .def(
          "IP", [](SICV const& self, std::map<unsigned, double> const& x)
          { return self.IP(x); }, py::arg("x"),
          "evaluate interval-valued polynomial at 'x'")
      .def(
          "P", [](SICV const& self) { return self.P(); },
          "return new model with interval range on the constant coefficient "
          "only")
      .def("constant", &SICV::constant, py::arg("reset") = false,
           "get interval coefficient of constant term and reset to zero if "
           "'reset' is True")
      .def("linear", &SICV::linear, py::arg("id"), py::arg("reset") = false,
           "get interval coefficient of linear term for variable 'id' and "
           "reset to zero if 'reset' is True")
      .def(
          "lift",
          [](SICV& self, SICM& sicm, double const& atol, double const& rtol)
          { return self.lift(&sicm, atol, rtol); }, py::arg("sicm"),
          py::arg("atol"), py::arg("rtol"),
          "lift model by appending new auxiliary variable from the "
          "uncertainty, if its magnitude is greater than ATOL/RTOL")
      .def(
          "project", [](SICV& self, bool const reset)
          { return self.project(reset); }, py::arg("reset") = false,
          "project model by removing monomials with any auxiliary variables "
          "participating")
      .def(
          "project", [](SICV& self, unsigned const& id)
          { return self.project(id); }, py::arg("id"),
          "project model by removing monomials with variable 'id' "
          "participating")
      .def(
          "scale", [](SICV& self, unsigned const& id, I const& dom)
          { return self.scale(id, dom); }, py::arg("id"), py::arg("dom"),
          "scale coefficients for the modified domain 'dom' of variable 'id'")
      .def(
          "rescale", [](SICV& self, std::map<unsigned, I> const& dom)
          { return self.scale(dom); }, py::arg("dom"),
          "rescale monomial coefficients in model over new domain 'dom'")
      .def("unscale", &SICV::unscale,
           "unscale monomial coefficients in model over original variable "
           "domain")
      .def("simplify", &SICV::simplify, py::arg("atol") = 0e0,
           py::arg("rtol") = 0e0, py::arg("tord") = -1,
           "simplify model by appending coefficients with magnitude less than "
           "ATOL/RTOL or order no less than TORD to the constant term")
      .def(
          "to_monomial", [](SICV const& self, bool const scaled)
          { return self.to_monomial(scaled); }, py::arg("scaled") = false,
          "express coefficient map in monomial basis representation")
      .def(
          "to_monomial",
          [](SICV const& self, bool const scaled, double const& atol,
             double const& rtol, int const tord)
          { return self.to_monomial(scaled, atol, rtol, tord); },
          py::arg("scaled"), py::arg("atol"), py::arg("rtol"),
          py::arg("tord") = -1,
          "express coefficient map in monomial basis after removing terms with "
          "coefficients less than ATOL/RTOL or order no less than ORD, and "
          "return bound on removed terms")
      .def("__str__",
           [](SICV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](SICV const& self)
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
      .def("__abs__", [](SICV const& m) { return mc::Op<SICV>::abs(m); })
      .def("__pow__", [](SICV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](SICV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__",
           [](SICV const& m, SICV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, SICV const& m) { return mc::pow(r, m); });

  m.def("inv", [](SICV const& x) { return mc::inv(x); });
  m.def("sqr", [](SICV const& x) { return mc::sqr(x); });
  m.def("sqrt", [](SICV const& x) { return mc::sqrt(x); });
  m.def("exp", [](SICV const& x) { return mc::exp(x); });
  m.def("log", [](SICV const& x) { return mc::log(x); });
  m.def("xlog", [](SICV const& x) { return mc::xlog(x); });
  m.def("cos", [](SICV const& x) { return mc::cos(x); });
  m.def("sin", [](SICV const& x) { return mc::sin(x); });
  m.def("tan", [](SICV const& x) { return mc::tan(x); });
  m.def("acos", [](SICV const& x) { return mc::acos(x); });
  m.def("asin", [](SICV const& x) { return mc::asin(x); });
  m.def("atan", [](SICV const& x) { return mc::atan(x); });
  m.def("cosh", [](SICV const& x) { return mc::cosh(x); });
  m.def("sinh", [](SICV const& x) { return mc::sinh(x); });
  m.def("tanh", [](SICV const& x) { return mc::tanh(x); });
  m.def("erf", [](SICV const& x) { return mc::erf(x); });
  m.def("erfc", [](SICV const& x) { return mc::erfc(x); });
  m.def("pow", [](SICV const& x, int const n) { return mc::pow(x, n); });
  m.def("pow", [](SICV const& x, double const& r) { return mc::pow(x, r); });
  m.def("pow", [](SICV const& x, SICV const& y) { return mc::pow(x, y); });
  m.def("pow", [](double const& r, SICV const& y) { return mc::pow(r, y); });
  m.def("cheb", [](SICV const& x, unsigned const n) { return mc::cheb(x, n); });
  m.def("fabs", [](SICV const& x) { return mc::fabs(x); });
  m.def("hull", [](SICV const& x, SICV const& y) { return mc::hull(x, y); });
  m.def("inter", [](SICV& xy, SICV const& x, SICV const& y)
        { return mc::inter(xy, x, y); });

  // Nested class Options

  pySICModelOptions.def(py::init<>())
      .def(py::init<SICM::Options const&>())
      .def(
          "reset", [](SICM::Options& self) { self.reset(); },
          "Reset options to default")
      .def_readwrite("BASIS", &SICM::Options::BASIS,
                     "Basis representation of the monomials: 0-Monomial basis; "
                     "1-Chebyshev basis [Default: 1 (CHEB)]")
      .def_readwrite("HOT_SPLIT", &SICM::Options::HOT_SPLIT,
                     "Strategy for allocating higher-order terms to monomials "
                     "in sparse product terms [Default: FULL]")
      .def_readwrite("PRODBND_SPLIT", &SICM::Options::PRODBND_SPLIT,
                     "Whether to distribute uncertainty in product with "
                     "interval bound [Default: False]")
      .def_readwrite("LIFT_USE", &SICM::Options::LIFT_USE,
                     "Whether to lift the uncertainty in nonlinear operations "
                     "by introducing auxiliary variables [Default: False]")
      .def_readwrite("LIFT_ATOL", &SICM::Options::LIFT_ATOL,
                     "Absolute tolerance for lifting - only if LIFT_USE == "
                     "true [Default: 1e-10]")
      .def_readwrite("LIFT_RTOL", &SICM::Options::LIFT_RTOL,
                     "Relative tolerance for lifting - only if LIFT_USE == "
                     "true [Default: 1e-3]")
      .def_readwrite("REMEZ_USE", &SICM::Options::REMEZ_USE,
                     "Whether to use the Remez algorithm for univariate "
                     "minimax approximation [Default: True]")
      .def_readwrite("REMEZ_MAXIT", &SICM::Options::REMEZ_MAXIT,
                     "Maximal number of Remez iterations [Default: 10]")
      .def_readwrite("REMEZ_TOL", &SICM::Options::REMEZ_TOL,
                     "Stopping tolerance in Remez algorithm [Default: 1e-5]")
      .def_readwrite("REMEZ_MIG", &SICM::Options::REMEZ_MIG,
                     "Threshold for interval width below which Remez is not "
                     "used [Default: 1e-10]")
      .def_readwrite(
          "INTERP_EXTRA", &SICM::Options::INTERP_EXTRA,
          "Extra terms in Chebyshev interpolation of univariates [Default: 0]")
      .def_readwrite("INTERP_THRES", &SICM::Options::INTERP_THRES,
                     "Threshold for coefficient values in Chebyshev "
                     "interpolation of univariates [Default: 1e2*machprec()]")
      .def_readwrite("BOUNDER_TYPE", &SICM::Options::BOUNDER_TYPE,
                     "Chebyshev model range bounder [Default: LSB]")
      .def_readwrite("BERNSTEIN_ORDER", &SICM::Options::BERNSTEIN_ORDER,
                     "Degree of the Bernstein basis when BOUNDER_TYPE is set "
                     "to BERNSTEIN [Default: 0]")
      .def_readwrite("MIG_USE", &SICM::Options::MIG_USE,
                     "Whether to simplify monomial terms with small magnitude "
                     "[Default: False]")
      .def_readwrite("MIG_ATOL", &SICM::Options::MIG_ATOL,
                     "Absolute tolerance for simplifying monomial terms - only "
                     "if MIG_USE == true [Default: 0]")
      .def_readwrite("MIG_RTOL", &SICM::Options::MIG_RTOL,
                     "Relative tolerance for simplifying monomial terms - only "
                     "if MIG_USE == true [Default: machprec()]")
      .def_readwrite("MIXED_IA", &SICM::Options::MIXED_IA,
                     "Whether to intersect internal bounds with underlying "
                     "bounds in the templated arithmetics [Default: False]")
      .def_readwrite("REF_POLY", &SICM::Options::REF_POLY,
                     "Scalar in [0,1] related to the choice of the polynomial "
                     "part in mc::inter and mc::hull [Default: 0]")
      .def_readwrite("DISPLAY_DIGITS", &SICM::Options::DISPLAY_DIGITS,
                     "Number of digits in output stream for Chebyshev model "
                     "coefficients [Default: 7]");


  py::enum_<SICM::Options::ALLOCATION>(pySICModelOptions, "ALLOCATION")
      .value("NONE", SICM::Options::ALLOCATION::NONE,
             "No split of HOT; allocate to constant coefficient")
      .value("SIMPLE", SICM::Options::ALLOCATION::SIMPLE,
             "Split HOT between current term and new variable")
      .value("FULL", SICM::Options::ALLOCATION::FULL,
             "Split HOT among all variables")
      .export_values();

  py::enum_<SICM::Options::MONBASIS>(pySICModelOptions, "MONBASIS")
      .value("MONOM", SICM::Options::MONBASIS::MONOM, "Monomial basis")
      .value("CHEB", SICM::Options::MONBASIS::CHEB, "Chebyshev basis")
      .export_values();

  // Nested class Exceptions
  py::class_<SICM::Exceptions> pySICModelExceptions(pySICModel, "Exceptions");

  py::enum_<SICM::Exceptions::TYPE>(pySICModelExceptions, "TYPE")
      .value("DIV", SICM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", SICM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", SICM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", SICM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("DPOW", SICM::Exceptions::TYPE::DPOW,
             "Real power operation with negative numbers in range")
      .value("TAN", SICM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", SICM::Exceptions::TYPE::ACOS,
             "Cosine inverse operation with range outside [-1,1]")
      .value("ASIN", SICM::Exceptions::TYPE::ASIN,
             "Sine inverse operation with range outside [-1,1]")
      .value("COMPOSE", SICM::Exceptions::TYPE::COMPOSE,
             "Failed to compose Chebyshev variable")
      .value("INIT", SICM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", SICM::Exceptions::TYPE::INCON,
             "Inconsistent bounds with template parameter arithmetic")
      .value("MODEL", SICM::Exceptions::TYPE::MODEL,
             "Operation between variables linked to different models")
      .value("INTERNAL", SICM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", SICM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented")
      .export_values();

  pySICModelExceptions.def(py::init<SICM::Exceptions::TYPE>())
      .def("ierr", &SICM::Exceptions::ierr, "Error flag")
      .def("what", &SICM::Exceptions::what, "Error description");
}
