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

#include "tmodel.hpp"

namespace py = pybind11;

void
mc_tmodel(py::module& m)
{
  typedef mc::TModel<I> TM;
  typedef mc::TVar<I> TV;

  py::class_<TM> pyTModel(m, "TModel");

  py::class_<TM::Options> pyTModelOptions(pyTModel, "Options");

  py::enum_<TM::Options::BOUNDER>(pyTModelOptions, "BOUNDER")
      .value("NAIVE", TM::Options::BOUNDER::NAIVE,
             "Naive polynomial range bounder")
      .value("LSB", TM::Options::BOUNDER::LSB, "Lin & Stadtherr range bounder")
      .value("EIGEN", TM::Options::BOUNDER::EIGEN,
             "Eigenvalue decomposition-based bounder")
      .value("BERNSTEIN", TM::Options::BOUNDER::BERNSTEIN,
             "Bernstein range bounder")
      .value("HYBRID", TM::Options::BOUNDER::HYBRID,
             "Hybrid LSB + EIGEN range bounder")
      .export_values();

  pyTModel
      // Constructors
      .def(py::init<unsigned const, unsigned const>(), py::arg("nvar"),
           py::arg("maxord"), "constructor for Taylor model")
      // Accessors
      .def_readwrite("options", &TM::options)
      // Inherited from PolyModel
      .def_property_readonly("sparse", &TM::sparse,
                             "whether or not sparse operations are enabled")
      .def_property_readonly("nvar", &TM::nvar, "number of variables in model")
      .def_property_readonly("nord", &TM::nord, "maximal degree of model")
      .def_property_readonly("nmon", &TM::nmon,
                             "number of monomial terms in polynomial model")
      .def(
          "expmon",
          [](TM const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), "Get variable exponents in monomial term i")
      .def(
          "loc_expmon",
          [](TM const& self, std::vector<unsigned> const& iexp)
          {
            if (iexp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.loc_expmon(iexp.data());
          },
          py::arg("iexp"), "Get index of monomial term with variable exponents")
      .def_property_readonly(
          "posord",
          [](TM const& self)
          {
            const unsigned* p = self.posord();
            return std::vector<unsigned>(p, p + self.nord() + 2);
          },
          "Indices of first monomial term of order i")
      .def("get_binom", &TM::get_binom, py::arg("n"), py::arg("k"),
           "Get binomial coefficient n choose k")
      // TModel specific and inherited from PolyModel
      .def(
          "bndvar",
          [](TM const& self, unsigned const i)
          {
            if (i >= self.nvar()) throw py::index_error();
            return self.bndvar()[i];
          },
          py::arg("i"), "bounds on model variable i")
      .def("refvar", &TM::refvar, py::arg("i"),
           "reference value of model variable i")
      .def("scalvar", &TM::scalvar, py::arg("i"),
           "scaling factor of model variable i")
      .def("reset", &TM::reset, "reset the bounds on powers of variable ranges")
      .def(
          "lift", [](TM const& self, unsigned const n, std::vector<I> const& X)
          { return self.lift(n, X.data()); }, py::arg("n"), py::arg("X"),
          py::return_value_policy::take_ownership,
          "lift Taylor model by n extra dimensions");

  py::class_<TV> pyTVar(m, "TVar");

  pyTVar
      // Constructors
      .def(py::init<double const>(), py::arg("cst") = 0.,
           "constructor for constant scalar")
      .def(py::init<I const&>(), py::arg("bnd"),
           "constructor for constant bound")
      .def(py::init<TM*, unsigned const, I const&>(), py::arg("env"),
           py::arg("ndx"), py::arg("rng"),
           "constructor for variable with range 'rng' and index 'ndx' in model "
           "'env' (reference point at midpoint)")
      .def(py::init<TM*, unsigned const, I const&, double const>(),
           py::arg("env"), py::arg("ndx"), py::arg("rng"), py::arg("ref"),
           "constructor for variable with range 'rng', index 'ndx' and "
           "reference point 'ref' in model 'env'")
      .def(py::init<TV const&>(), "copy constructor")
      // Modifiers
      .def(
          "set",
          [](TV& self, TM* env, unsigned const ix, I const& X,
             double const& ref) -> TV& { return self.set(env, ix, X, ref); },
          py::arg("env"), py::arg("ndx"), py::arg("rng"), py::arg("ref"),
          py::return_value_policy::reference_internal,
          "set variable with range 'rng', index 'ndx' and reference point "
          "'ref' in model 'env'")
      .def(
          "set", [](TV& self, TM* env, unsigned const ix, I const& X) -> TV&
          { return self.set(env, ix, X); }, py::arg("env"), py::arg("ndx"),
          py::arg("rng"), py::return_value_policy::reference_internal,
          "set variable with range 'rng' and index 'ndx' in model 'env'")
      .def(
          "set", [](TV& self, TM* env, bool const reset) -> TV&
          { return self.set(env, reset); }, py::arg("env"),
          py::arg("reset") = false, py::return_value_policy::reference_internal,
          "set model environment")
      // Accessors
      .def("env", &TV::env, "model environment")
      // Inherited from PolyVar
      .def_property_readonly("nord", &TV::nord, "current degree")
      .def_property_readonly("nvar", &TV::nvar,
                             "number of variable dependencies")
      .def_property_readonly("nmon", &TV::nmon, "number of monomials")
      .def(
          "coefmon",
          [](TV const& self)
          {
            // Return all coefficients as a list/vector
            // Note: dense TVar stores them contiguously
            return std::vector<double>(self.coefmon().second,
                                       self.coefmon().second + self.nmon());
          },
          "polynomial coefficients (all)")
      .def(
          "coefmon",
          [](TV const& self, std::vector<unsigned> const& exp)
          {
            if (exp.size() != self.nvar())
              throw std::invalid_argument("Size mismatch");
            return self.coefmon(exp.data());
          },
          py::arg("exp"),
          "Get coefficient in monomial term with variable exponents 'exp'")
      .def(
          "expmon",
          [](TV const& self, unsigned const i)
          {
            if (i >= self.nmon()) throw py::index_error();
            const unsigned* p = self.expmon(i);
            return std::vector<unsigned>(p, p + self.nvar());
          },
          py::arg("i"), "Get monomial exponents for term i")
      .def_property_readonly("ndxmon", &TV::ndxmon, "Indices of nonzero terms")
      .def_property_readonly(
          "bound", [](TV const& self) { return self.bound(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def_property_readonly(
          "B", [](TV const& self) { return self.B(); },
          "bound on model using selected bounder in option 'BOUNDER_TYPE'")
      .def(
          "bndpol", [](TV const& self) { return self.bndpol(); },
          "Retrieve bound on multivariate polynomial using default bounder")
      .def(
          "bndpol", [](TV const& self, int const type)
          { return self.bndpol(type); }, py::arg("type"),
          "Retrieve bound on multivariate polynomial using bounder 'type'")
      .def("bndord", &TV::bndord, py::arg("minord"),
           "compute bound on all monomials with (total) order no less than "
           "'minord'")
      .def(
          "P", [](TV const& self, std::vector<double> const& x)
          { return self.P(x.data()); }, py::arg("x"),
          "evaluate polynomial part at 'x'")
      .def(
          "P", [](TV const& self) { return self.P(); },
          "cancel remainder and return new model with same polynomial part")
      .def_property_readonly(
          "R", [](TV const& self) { return self.R(); }, "remainder bound")
      .def("center", &TV::center, "center remainder term in model")
      .def("C", &TV::C, "center remainder term in model")
      .def("constant", &TV::constant, py::arg("reset") = false,
           "get coefficient of constant term and reset to zero if 'reset' is "
           "True")
      .def(
          "linear", [](TV& self, unsigned const ivar, bool const reset)
          { return self.linear(ivar, reset); }, py::arg("id"),
          py::arg("reset") = false,
          "get coefficient of linear term for variable 'id' and reset to zero "
          "if 'reset' is True")
      .def(
          "linear",
          [](TV const& self)
          {
            double* p = self.linear();
            std::vector<double> res(p, p + self.nvar());
            delete[] p;
            return res;
          },
          "get coefficients of linear terms")
      .def("__str__",
           [](TV const& self)
           {
             std::ostringstream os;
             os << self;
             return os.str();
           })
      .def("__repr__",
           [](TV const& self)
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
      .def("__abs__", [](TV const& m) { return mc::Op<TV>::abs(m); })
      .def("__pow__", [](TV const& m, int const n) { return mc::pow(m, n); })
      .def("__pow__",
           [](TV const& m, double const& r) { return mc::pow(m, r); })
      .def("__pow__", [](TV const& m, TV const& mm) { return mc::pow(m, mm); })
      .def("__pow__",
           [](double const& r, TV const& m) { return mc::pow(r, m); });

  m.def("inv", [](TV const& x) { return mc::inv(x); });
  m.def("sqr", [](TV const& x) { return mc::sqr(x); });
  m.def("sqrt", [](TV const& x) { return mc::sqrt(x); });
  m.def("exp", [](TV const& x) { return mc::exp(x); });
  m.def("log", [](TV const& x) { return mc::log(x); });
  m.def("xlog", [](TV const& x) { return mc::xlog(x); });
  m.def("cos", [](TV const& x) { return mc::cos(x); });
  m.def("sin", [](TV const& x) { return mc::sin(x); });
  m.def("tan", [](TV const& x) { return mc::tan(x); });
  m.def("acos", [](TV const& x) { return mc::acos(x); });
  m.def("asin", [](TV const& x) { return mc::asin(x); });
  m.def("atan", [](TV const& x) { return mc::atan(x); });
  m.def("cosh", [](TV const& x) { return mc::cosh(x); });
  m.def("sinh", [](TV const& x) { return mc::sinh(x); });
  m.def("tanh", [](TV const& x) { return mc::tanh(x); });
  m.def("pow", [](TV const& x, int const n) { return mc::pow(x, n); });
  m.def("pow", [](TV const& x, double const& r) { return mc::pow(x, r); });
  m.def("pow", [](TV const& x, TV const& y) { return mc::pow(x, y); });
  m.def("pow", [](double const& r, TV const& y) { return mc::pow(r, y); });
  m.def("cheb", [](TV const& x, unsigned const n) { return mc::cheb(x, n); });
  m.def("hull", [](TV const& x, TV const& y) { return mc::hull(x, y); });
  m.def("inter",
        [](TV& xy, TV const& x, TV const& y) { return mc::inter(xy, x, y); });


  pyTModelOptions.def(py::init<>())
      .def(py::init<TM::Options const&>())
      .def("reset", &TM::Options::reset, "Reset options to defaults")
      .def_readwrite("BOUNDER_TYPE", &TM::Options::BOUNDER_TYPE,
                     "Taylor model range bounder [Default: LSB]")
      .def_readwrite("BERNSTEIN_ORDER", &TM::Options::BERNSTEIN_ORDER,
                     "Order of Bernstein polynomial for Taylor model range "
                     "bounding, when BOUNDER_TYPE = BERNSTEIN is selected "
                     "[Default: 0 (same as model order)]")
      .def_readwrite("SCALE_VARIABLES", &TM::Options::SCALE_VARIABLES,
                     "Whether to scale the variable ranges to [-1,1] "
                     "internally [Default: False]")
      .def_readwrite("CENTER_REMAINDER", &TM::Options::CENTER_REMAINDER,
                     "Whether to center the remainder term during Taylor model "
                     "propagation [Default: False]")
      .def_readwrite(
          "REF_MIDPOINT", &TM::Options::REF_MIDPOINT,
          "Whether to take the midpoint of the inner range as the reference in "
          "the outer composition with a univariate function [Default: False]")
      .def_readwrite(
          "REF_POLY", &TM::Options::REF_POLY,
          "Scalar in [0,1] related to the choice of the polynomial part in the "
          "overloaded functions mc::inter and mc::hull [Default: 0]")
      .def_readwrite("BERNSTEIN_USE", &TM::Options::BERNSTEIN_USE,
                     "Whether to compute a Berstein model of the outer "
                     "function in a univariate composition [Default: False]")
      .def_readwrite("BERNSTEIN_OPT", &TM::Options::BERNSTEIN_OPT,
                     "Whether to compute exact remainder bounds for Berstein "
                     "models of convex/concave univariates [Default: True]")
      .def_readwrite(
          "BERNSTEIN_MAXIT", &TM::Options::BERNSTEIN_MAXIT,
          "Maximum number of iterations for determination of the exact "
          "remainder bounds in a Berstein model [Default: 100]")
      .def_readwrite("BERNSTEIN_TOL", &TM::Options::BERNSTEIN_TOL,
                     "Termination tolerance for determination of the exact "
                     "remainder bounds in a Berstein model [Default: 1e-10]")
      .def_readwrite("DISPLAY_DIGITS", &TM::Options::DISPLAY_DIGITS,
                     "Number of digits in output stream for Taylor model "
                     "coefficients [Default: 5]");


  // Nested class Exceptions
  py::class_<TM::Exceptions> pyTModelExceptions(pyTModel, "Exceptions");

  py::enum_<TM::Exceptions::TYPE>(pyTModelExceptions, "TYPE")
      .value("DIV", TM::Exceptions::TYPE::DIV, "Division by zero scalar")
      .value("INV", TM::Exceptions::TYPE::INV,
             "Inverse operation with zero in range")
      .value("LOG", TM::Exceptions::TYPE::LOG,
             "Log operation with non-positive numbers in range")
      .value("SQRT", TM::Exceptions::TYPE::SQRT,
             "Square-root operation with negative numbers in range")
      .value("ASIN", TM::Exceptions::TYPE::ASIN,
             "Sine/Cosine inverse operation with range outside [-1,1]")
      .value("EIGEN", TM::Exceptions::TYPE::EIGEN,
             "Failed to compute eigenvalue decomposition in range bounder "
             "TModel::Options::EIGEN")
      .value("BERNSTEIN", TM::Exceptions::TYPE::BERNSTEIN,
             "Failed to compute the maximum gap between a univariate term and "
             "its Bernstein model")
      .value("INIT", TM::Exceptions::TYPE::INIT,
             "Failed to construct Taylor variable")
      .value("INCON", TM::Exceptions::TYPE::INCON,
             "Taylor model bound does not intersect with bound in template "
             "parameter arithmetic")
      .value("TMODEL", TM::Exceptions::TYPE::TMODEL,
             "Operation between Taylor variables linked to different Taylor "
             "models")
      .value("UNDEF", TM::Exceptions::TYPE::UNDEF,
             "Feature not yet implemented in mc::TModel")
      .export_values();

  pyTModelExceptions.def(py::init<TM::Exceptions::TYPE>())
      .def("ierr", &TM::Exceptions::ierr, "Error flag")
      .def("what", &TM::Exceptions::what, "Error description");
}
