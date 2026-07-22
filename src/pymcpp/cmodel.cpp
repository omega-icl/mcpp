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

#include "cmodel.hpp"

namespace py = pybind11;

void
mc_cmodel(py::module& m)
{
  typedef mc::CModel<I> CM;
  typedef mc::CVar<I> CV;

  py::class_<CM> pyCModel(m, "CModel");

  py::class_<CM::Options> pyCModelOptions(pyCModel, "Options");

  py::enum_<CM::Options::BOUNDER>(pyCModelOptions, "BOUNDER")
      .value("NAIVE", CM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", CM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .value("EIGEN", CM::Options::BOUNDER::EIGEN,
             "Eigenvalue decomposition-based bounder")
      .value("BERNSTEIN", CM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .value("HYBRID", CM::Options::BOUNDER::HYBRID,
             "Hybrid LSB + EIGEN range bounder")
      .export_values();

  pyCModel
      // Constructors
      .def(py::init<unsigned const, unsigned const, bool const>(),
           py::arg("nvar"), py::arg("maxord"), py::arg("sparse") = false,
           "constructor for dense Chebyshev model")
      // Accessors
      .def_readwrite("options", &CM::options)
      // Inherited from PolyModel
      .def_property_readonly("sparse", &CM::sparse,
                             "whether or not sparse operations are enabled")
      .def_property_readonly("nvar", &CM::nvar, "number of variables in model")
      .def_property_readonly("nord", &CM::nord, "maximal degree of model")
      .def_property_readonly("nmon", &CM::nmon,
                             "number of monomial terms in polynomial model")
      .def(
          "expmon",
          [](CM const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), "Get variable exponents in monomial term i")
      .def(
          "loc_expmon",
          [](CM const& self, std::vector<unsigned> const& iexp)
          {
            if (iexp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.loc_expmon(iexp.data());
          },
          py::arg("iexp"), "Get index of monomial term with variable exponents")
      .def_property_readonly(
          "posord",
          [](CM const& self)
          {
            const unsigned* p = self.posord();
            return std::vector<unsigned>(p, p + self.nord() + 2);
          },
          "Indices of first monomial term of order i")
      .def("get_binom", &CM::get_binom, py::arg("n"), py::arg("k"),
           "Get binomial coefficient n choose k")
      // CModel specific
      .def("bndvar", &CM::bndvar, py::arg("i"), "bounds on model variable i")
      .def("refvar", &CM::refvar, py::arg("i"),
           "reference value of model variable i")
      .def("scalvar", &CM::scalvar, py::arg("i"),
           "scaling factor of model variable i")
      //.def(
      //  "reset",
      //  &CM::reset,
      //  "reset the bounds on Chebyshev basis functions"
      //)
      .def(
          "lift", [](CM const& self, unsigned const n, std::vector<I> const& X)
          { return self.lift(n, X.data()); }, py::arg("n"), py::arg("X"),
          py::return_value_policy::take_ownership,
          "lift Chebyshev model by n extra dimensions");

  py::class_<CV> pyCVar(m, "CVar");

  pyCVar
      // Constructors
      .def(py::init<double const>(), py::arg("cst") = 0.,
           "constructor for constant scalar")
      .def(py::init<I const&>(), py::arg("bnd"),
           "constructor for constant bound")
      .def(py::init<CM*, unsigned const, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"),
           "constructor for variable with range 'rng' and index 'ndx' in model "
           "'mod'")
      .def(py::init<CV const&>(), "copy constructor")
      // Modifiers
      .def(
          "set", [](CV& self, CM* env, unsigned const ix, I const& X) -> CV&
          { return self.set(env, ix, X); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal,
          "set variable with range 'rng' and index 'ndx' in model 'env'")
      .def(
          "set", [](CV& self, CM* env, bool const reset) -> CV&
          { return self.set(env, reset); }, py::arg("env"),
          py::arg("reset") = false, py::return_value_policy::reference_internal,
          "set model environment")
      // Accessors
      .def("env", &CV::env, "model environment")
      // Inherited from PolyVar
      .def_property_readonly("nord", &CV::nord, "current degree")
      .def_property_readonly("nvar", &CV::nvar,
                             "number of variable dependencies")
      .def_property_readonly("nmon", &CV::nmon, "number of monomials")
      .def(
          "coefmon",
          [](CV const& self)
          {
            // Return all coefficients as a list/vector
            // Note: dense CVar stores them contiguously
            return std::vector<double>(self.coefmon().second,
                                       self.coefmon().second + self.nmon());
          },
          "polynomial coefficients (all)")
      .def(
          "coefmon",
          [](CV const& self, std::vector<unsigned> const& exp)
          {
            if (exp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.coefmon(exp.data());
          },
          py::arg("exp"),
          "Get coefficient in monomial term with variable exponents 'exp'")
      .def(
          "expmon",
          [](CV const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), "Get monomial exponents for term i")
      .def_property_readonly(
          "ndxmon", &CV::ndxmon,
          "Indices of nonzero terms in sparse representation")
      .def_property_readonly(
          "bound", [](CV const& self) { return self.bound(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def_property_readonly(
          "B", [](CV const& self) { return self.B(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def(
          "bndpol", [](CV const& self) { return self.bndpol(); },
          "Retrieve bound on multivariate polynomial using default bounder")
      .def(
          "bndpol", [](CV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          "Retrieve bound on multivariate polynomial using bounder 'type'")
      .def("bndord", &CV::bndord, py::arg("minord"),
           "compute bound on all monomials with (total) order no less than "
           "'minord'")
      .def(
          "P", [](CV const& self, std::vector<double> const& x)
          { return self.P(x.data()); }, py::arg("x"),
          "evaluate polynomial part at 'x'")
      .def(
          "P", [](CV const& self) { return self.P(); },
          "cancel remainder and return new model with same polynomial part")
      .def_property_readonly(
          "R", [](CV const& self) { return self.R(); }, "remainder bound")
      .def("center", &CV::center, "center remainder term in model")
      .def("C", &CV::C, "center remainder term in model")
      .def("constant", &CV::constant, py::arg("reset") = false,
           "get coefficient of constant term and reset to zero if 'reset' is "
           "True")
      .def(
          "linear", [](CV& self, unsigned const ivar, bool const reset)
          { return self.linear(ivar, reset); }, py::arg("id"),
          py::arg("reset") = false,
          "get coefficient of linear term for variable 'id' and reset to zero "
          "if 'reset' is True")
      //.def(
      //  "scale",
      //  []( CV const& self, std::vector<I> const& X )
      //    { return self.scale( X.data() ); },
      //  py::arg("dom"),
      //  "Scale coefficients in Chebyshev variable for the reduced variable
      //  dom"
      //)
      .def("__str__",
           [](CV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](CV const& self)
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
      .def("__abs__", [](CV const& m) { return mc::Op<CV>::abs(m); })
      .def("__pow__",
           [](CV const& m, int const n) { return mc::Op<CV>::pow(m, n); })
      .def("__pow__",
           [](CV const& m, double const& r) { return mc::Op<CV>::pow(m, r); })
      .def("__pow__",
           [](CV const& m, CV const& mm) { return mc::Op<CV>::pow(m, mm); })
      .def("__pow__",
           [](double const& r, CV const& m) { return mc::Op<CV>::pow(r, m); });

  m.def("inv", [](CV const& x) { return mc::inv(x); });
  m.def("sqr", [](CV const& x) { return mc::sqr(x); });
  m.def("sqrt", [](CV const& x) { return mc::sqrt(x); });
  m.def("exp", [](CV const& x) { return mc::exp(x); });
  m.def("log", [](CV const& x) { return mc::log(x); });
  m.def("xlog", [](CV const& x) { return mc::xlog(x); });
  m.def("cos", [](CV const& x) { return mc::cos(x); });
  m.def("sin", [](CV const& x) { return mc::sin(x); });
  m.def("tan", [](CV const& x) { return mc::tan(x); });
  m.def("acos", [](CV const& x) { return mc::acos(x); });
  m.def("asin", [](CV const& x) { return mc::asin(x); });
  m.def("atan", [](CV const& x) { return mc::atan(x); });
  m.def("cosh", [](CV const& x) { return mc::cosh(x); });
  m.def("sinh", [](CV const& x) { return mc::sinh(x); });
  m.def("tanh", [](CV const& x) { return mc::tanh(x); });
  m.def("pow", [](CV const& x, int const n) { return mc::pow(x, n); });
  m.def("pow", [](CV const& x, double const& r) { return mc::pow(x, r); });
  m.def("pow", [](CV const& x, CV const& y) { return mc::pow(x, y); });
  m.def("pow", [](double const& r, CV const& y) { return mc::pow(r, y); });
  m.def("cheb", [](CV const& x, unsigned const n) { return mc::cheb(x, n); });
  m.def("fabs", [](CV const& x) { return mc::fabs(x); });
  m.def("hull", [](CV const& x, CV const& y) { return mc::hull(x, y); });
  m.def("inter",
        [](CV& xy, CV const& x, CV const& y) { return mc::inter(xy, x, y); });

  // Nested class Options


  pyCModelOptions.def(py::init<>())
      .def(py::init<CM::Options const&>())
      .def("reset", &CM::Options::reset, "Reset options to defaults")
      .def_readwrite(
          "INTERP_EXTRA", &CM::Options::INTERP_EXTRA,
          "Extra terms in Chebyshev interpolation of univariates: 0-Chebyshev "
          "interpolation of order NORD; extra terms allow approximation of "
          "Chebyshev truncated series [Default: 0]")
      .def_readwrite(
          "INTERP_THRES", &CM::Options::INTERP_THRES,
          "Threshold for coefficient values in Chebyshev expansion for "
          "bounding of transcendental univariates [Default: 1e2*DBL_EPSILON]")
      .def_readwrite("BOUNDER_TYPE", &CM::Options::BOUNDER_TYPE,
                     "Chebyshev model range bounder [Default: LSB]")
      .def_readwrite("BERNSTEIN_ORDER", &CM::Options::BERNSTEIN_ORDER,
                     "Degree of Bernstein polynomial for Chebyshev model range "
                     "bounding, when BOUNDER_TYPE = BERNSTEIN is selected "
                     "[Default: 0 (same as model order)]")
      .def_readwrite("MIXED_IA", &CM::Options::MIXED_IA,
                     "Whether to intersect internal bounds with underlying "
                     "bounds in the templated arithmetics [Default: False]")
      .def_readwrite(
          "REF_POLY", &CM::Options::REF_POLY,
          "Scalar in [0,1] related to the choice of the polynomial part in the "
          "overloaded functions mc::inter and mc::hull [Default: 0]")
      .def_readwrite("DISPLAY_DIGITS", &CM::Options::DISPLAY_DIGITS,
                     "Number of digits in output stream for Chebyshev model "
                     "coefficients [Default: 7]");


  // Nested class Exceptions
  py::class_<CM::Exceptions> pyCModelExceptions(pyCModel, "Exceptions");

  py::enum_<CM::Exceptions::TYPE>(pyCModelExceptions, "TYPE")
      .value("DIV", CM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", CM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", CM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", CM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("TAN", CM::Exceptions::TYPE::TAN,
             "Tangent operation with (k+1/2)·PI in range")
      .value("ACOS", CM::Exceptions::TYPE::ACOS,
             "Sine/Cosine inverse operation with range outside [-1,1]")
      .value("EIGEN", CM::Exceptions::TYPE::EIGEN,
             "Failed to compute eigenvalue decomposition in range bounder "
             "CModel::Options::EIGEN")
      .value("INIT", CM::Exceptions::TYPE::INIT,
             "Failed to construct Chebyshev variable")
      .value("INCON", CM::Exceptions::TYPE::INCON,
             "Chebyshev model bound does not intersect with bound in template "
             "parameter arithmetic")
      .value("CMODEL", CM::Exceptions::TYPE::CMODEL,
             "Operation between Chebyshev variables linked to different "
             "Chebyshev models")
      .value("INTERNAL", CM::Exceptions::TYPE::INTERNAL, "Internal error")
      .value("UNDEF", CM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented in module")
      .export_values();

  pyCModelExceptions.def(py::init<CM::Exceptions::TYPE>())
      .def("ierr", &CM::Exceptions::ierr, "Error flag")
      .def("what", &CM::Exceptions::what, "Error description");
}
