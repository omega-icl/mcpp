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

#include "polimage.hpp"

typedef mc::PolImg<I> PI;
typedef mc::PolVar<I> PV;
typedef mc::PolCut<I> PC;

namespace py = pybind11;

void
mc_polimage(py::module_& m)
{
  py::class_<PV> pyPolVar(m, "PolVar");

  py::enum_<PV::TYPE>(pyPolVar, "TYPE")
      .value("VARCONT", PV::TYPE::VARCONT, "DAG continuous variable")
      .value("VARINT", PV::TYPE::VARINT, "DAG integer variable")
      .value("AUXCONT", PV::TYPE::AUXCONT, "Auxiliary continuous variable")
      .value("AUXINT", PV::TYPE::AUXINT, "Auxiliary integer variable")
      .value("AUXCST", PV::TYPE::AUXCST, "Auxiliary constant")
      .export_values();

  pyPolVar
      // --- Static Constants ---
      .def_property_static(
          "VARCONTNAME", [](py::object) { return PV::VARCONTNAME; },
          [](py::object, std::string const& str) { PV::VARCONTNAME = str; })
      .def_property_static(
          "VARINTNAME", [](py::object) { return PV::VARINTNAME; },
          [](py::object, std::string const& str) { PV::VARINTNAME = str; })
      .def_property_static(
          "AUXCONTNAME", [](py::object) { return PV::AUXCONTNAME; },
          [](py::object, std::string const& str) { PV::AUXCONTNAME = str; })
      .def_property_static(
          "AUXINTNAME", [](py::object) { return PV::AUXINTNAME; },
          [](py::object, std::string const& str) { PV::AUXINTNAME = str; })
      .def_property_static(
          "AUXCSTNAME", [](py::object) { return PV::AUXCSTNAME; },
          [](py::object, std::string const& str) { PV::AUXCSTNAME = str; })

      // --- Constructors ---
      .def(py::init<double const&>(), py::arg("d") = 0.,
           "Constructor for a constant value")
      .def(py::init<int const>(), py::arg("n"),
           "Constructor for a constant integer value")
      .def(py::init<PI*, mc::FFVar const&, I const&, bool const>(),
           py::arg("img"), py::arg("var"), py::arg("rng") = 0.,
           py::arg("cnt") = true,
           "Constructor for DAG variable 'var' with range 'rng' in polytope "
           "image 'img'")
      .def(py::init<PI*, I const&, bool const>(), py::arg("img"),
           py::arg("rng") = 0., py::arg("cnt") = true,
           "Constructor for auxiliary variable with range 'rng' in polytope "
           "image 'img'")
      .def(py::init<PV const&>(), "Copy constructor")
      // --- Modifiers ---
      .def("update", py::overload_cast<I const&>(&PV::update), py::arg("rng"),
           py::return_value_policy::reference_internal, "Update variable range")
      .def("update", py::overload_cast<bool const>(&PV::update), py::arg("cnt"),
           py::return_value_policy::reference_internal,
           "Update variable type (continuous/discrete)")
      .def("set",
           py::overload_cast<PI*, mc::FFVar const&, I const&, bool const>(
               &PV::set),
           py::arg("img"), py::arg("var"), py::arg("rng") = 0.,
           py::arg("cnt") = true, py::return_value_policy::reference_internal,
           "Set as DAG variable in polytope image")
      .def("set", py::overload_cast<PI*, I const&, bool const>(&PV::set),
           py::arg("img"), py::arg("rng") = 0., py::arg("cnt") = true,
           py::return_value_policy::reference_internal,
           "Set as auxiliary variable in polytope image")
      .def("add_breakpt", &PV::add_breakpt, py::arg("bkpt"),
           "Add break-point to variable")
      .def("reset_subdiv", &PV::reset_subdiv, "Reset variable subdivision")
      .def("reset_cuts", &PV::reset_cuts, "Reset cuts flag")
      .def("set_cuts", &PV::set_cuts, "Set cuts flag")
      .def("create_subdiv", &PV::create_subdiv, py::arg("lbd"), py::arg("ubd"),
           py::arg("reset") = false, "Create variable subdivision")
      .def("sos2_subdiv", &PV::SOS2_subdiv, py::arg("op") = nullptr,
           py::arg("reset") = false,
           "Set SOS2 variable subdivision using SOS2 encoding")
      .def("bin_subdiv", &PV::BIN_subdiv, py::arg("op") = nullptr,
           py::arg("reset") = false,
           "Set variable subdivision using binary encoding")
      .def("cnt_subdiv", &PV::CONT_subdiv, py::arg("pOp") = nullptr,
           py::arg("reset") = false,
           "Set variable subdivision using relaxed continuous encoding")
      // --- Accessors ---
      .def_property_readonly("name", &PV::name, "Variable name")
      .def_property_readonly(
          "cnt", [](PV const& self) { return !self.discr(); },
          "Check if variable is continuous")
      .def_property_readonly("cst", &PV::cst, "Check if variable is constant")
      .def_property_readonly("rng", &PV::range, "Variable range")
      .def_property_readonly("id", &PV::id, "Variable unique identifier")
      .def_property_readonly("img", &PV::image,
                             py::return_value_policy::reference_internal,
                             "Polyhedral image")
      .def_property_readonly("var", py::overload_cast<>(&PV::var),
                             py::return_value_policy::reference_internal,
                             "Associated DAG variable")
      .def_property_readonly("breakpts", &PV::breakpts, "Set of breakpoints")
      .def_property_readonly("subdiv", &PV::subdiv,
                             "Variable subdivision (pairs of points/variables)")
      .def("has_cuts", &PV::has_cuts, "Check if cuts have been generated")
      // --- Operators ---
      /*
       //.def( py::self = py::self )
       //.def( py::self = double() )
       //.def( py::self = int() )
       .def( py::self += py::self )
       .def( py::self -= py::self )
       .def( py::self *= py::self )
       .def( py::self /= py::self )
       //
       .def( + py::self )
       .def( py::self + py::self )
       .def( py::self + double() )
       .def( - py::self )
       .def( py::self - py::self )
       .def( py::self * py::self )
       .def( py::self * double() )
       .def( py::self / py::self )
       .def( double() / py::self )
       //
       .def( "__pow__", []( PV const& x, int const n ){ return mc::pow(x,n); },
       py::is_operator() ) .def( "__pow__", []( PV const& x, double const& r ){
       return mc::pow(x,r); },    py::is_operator() )
       //.def( "__pow__", []( PV const& x, PV const& y ){ return mc::pow(x,y);
       },        py::is_operator() ) .def( "__le__",  []( PV const& x, PV const&
       y ){ return mc::Op<PV>::le(x,y); }, py::is_operator()) .def( "__lt__",
       []( PV const& x, PV const& y ){ return mc::Op<PV>::lt(x,y); },
       py::is_operator()) .def( "__ge__",  []( PV const& x, PV const& y ){
       return mc::Op<PV>::ge(x,y); }, py::is_operator()) .def( "__gt__",  []( PV
       const& x, PV const& y ){ return mc::Op<PV>::gt(x,y); },
       py::is_operator())
      */
      .def(
          "__eq__", [](PV const& x, PV const& y) { return x.id() == y.id(); },
          py::is_operator())
      .def(
          "__ne__", [](PV const& x, PV const& y) { return x.id() != y.id(); },
          py::is_operator())
      .def("__str__",
           [](PV const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
           })
      .def("__repr__",
           [](PV const& V)
           {
             std::ostringstream Vss;
             Vss << V;
             return Vss.str();
           });
  /*
  m.def( "abs",   []( PV const& x ){ return mc::Op<PV>::abs(x); } );
  m.def( "mid",   []( PV const& x ){ return mc::Op<PV>::mid(x); } );
  m.def( "diam",  []( PV const& x ){ return mc::Op<PV>::diam(x); } );
  m.def( "inv",   []( PV const& x ){ return mc::inv(x); } );
  m.def( "sqr",   []( PV const& x ){ return mc::sqr(x); } );
  m.def( "sqrt",  []( PV const& x ){ return mc::sqrt(x); } );
  m.def( "exp",   []( PV const& x ){ return mc::exp(x); } );
  m.def( "log",   []( PV const& x ){ return mc::log(x); } );
  m.def( "cos",   []( PV const& x ){ return mc::cos(x); } );
  m.def( "sin",   []( PV const& x ){ return mc::sin(x); } );
  m.def( "tan",   []( PV const& x ){ return mc::tan(x); } );
  m.def( "acos",  []( PV const& x ){ return mc::acos(x); } );
  m.def( "asin",  []( PV const& x ){ return mc::asin(x); } );
  m.def( "atan",  []( PV const& x ){ return mc::atan(x); } );
  m.def( "cosh",  []( PV const& x ){ return mc::cosh(x); } );
  m.def( "sinh",  []( PV const& x ){ return mc::sinh(x); } );
  m.def( "tanh",  []( PV const& x ){ return mc::tanh(x); } );
  m.def( "fabs",  []( PV const& x ){ return mc::fabs(x); } );
  //m.def( "relu",  []( PV const& x ){ return mc::max(x,0.); } );
  m.def( "xlog",  []( PV const& x ){ return mc::xlog(x); } );
  m.def( "fstep", []( PV const& x ){ return mc::fstep(x); } );
  //m.def( "bstep", []( PV const& x ){ return mc::bstep(x); } );
  m.def( "erf",   []( PV const& x ){ return mc::erf(x); } );
  //m.def( "erfc",  []( PV const& x ){ return mc::erfc(x); } );
  m.def( "pow",   []( PV const& x, int const n ){ return mc::pow(x,n); } );
  m.def( "pow",   []( PV const& x, double const& r ){ return mc::pow(x,r); } );
  //m.def( "pow",   []( PV const& x, PV const& y ){ return mc::pow(x,y); } );
  //m.def( "pow",   []( double const& r, PV const& y ){ return
  mc::exp(y*std::log(r)); } ); m.def( "cheb",  []( PV const& x, unsigned const n
  ){ return mc::cheb(x,n); } ); m.def( "max",   []( PV const& x, PV const& y ){
  return mc::max(x,y); } ); m.def( "min",   []( PV const& x, PV const& y ){
  return mc::min(x,y); } ); m.def( "hull",  []( PV const& x, PV const& y ){
  return mc::Op<PV>::hull(x,y); } ); m.def( "inter",  []( PV& z, PV const& x, PV
  const& y ){ return mc::Op<PV>::inter(z,x,y); } );
  */
  py::class_<PC> pyPolCut(m, "PolCut");

  // Enumeration for PolCut::TYPE
  py::enum_<PC::TYPE>(pyPolCut, "TYPE")
      .value("EQ", PC::TYPE::EQ, "Equality constraint Ax=b")
      .value("LE", PC::TYPE::LE, "Inequality constraint Ax<=b")
      .value("GE", PC::TYPE::GE, "Inequality constraint Ax>=b")
      .value("SOS1", PC::TYPE::SOS1, "SOS1-type constraint")
      .value("SOS2", PC::TYPE::SOS2, "SOS2-type constraint")
      .value("NLIN", PC::TYPE::NLIN, "Nonlinear constraint y=f(x)")
      .export_values();

  pyPolCut
      // --- Constructors ---
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&>(), py::arg("op"),
           py::arg("type"), py::arg("b"),
           "Constructor for cut without participating variable")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), "Constructor for cut with 1 linear variable")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"),
           "Constructor for cut with 2 linear variables")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&, PV const&,
                    double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"), py::arg("x3"),
           py::arg("a3"), "Constructor for cut with 3 linear variables")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&, PV const&,
                    double const&, PV const&, double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"), py::arg("x3"),
           py::arg("a3"), py::arg("x4"), py::arg("a4"),
           "Constructor for cut with 4 linear variables")
      .def(py::init(
               [](mc::FFOp const* op, PC::TYPE type, double b,
                  std::vector<PV> const& x, std::vector<double> const& a)
               {
                 if (x.size() != a.size())
                   throw std::invalid_argument(
                       "Size mismatch between variables and coefficients");
                 return new PC(op, type, b, x.size(), x.data(), a.data());
               }),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x"),
           py::arg("a"),
           // py::return_value_policy::take_ownership,
           "Constructor for cut with linear variable list and coefficients")
      .def(py::init([](mc::FFOp const* op, PC::TYPE type, double b,
                       std::vector<PV> const& x, double const& a0)
                    { return new PC(op, type, b, x.size(), x.data(), a0); }),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x"),
           py::arg("a0"),
           "Constructor for cut with linear variable list and a shared "
           "coefficient")
      .def(py::init<mc::FFOp const*, PV const&>(), py::arg("op"), py::arg("x1"),
           "Constructor for cut with 1 nonlinear variables")
      .def(py::init<mc::FFOp const*, PV const&, PV const&>(), py::arg("op"),
           py::arg("x1"), py::arg("x2"),
           "Constructor for cut with 2 nonlinear variables")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, double const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("b"),
           "Constructor for cut with 2 nonlinear variables and 1 constant")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, PV const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("x3"),
           "Constructor for cut with 3 nonlinear variables")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, PV const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("x3"),
           "Constructor for cut with nonlinear variable list")
      .def(py::init(
               [](mc::FFOp const* op, PV const& X1, std::vector<PV> const& X)
               { return new PC(op, X1, X.size(), X.data()); }),
           py::arg("op"), py::arg("x1"), py::arg("x"),
           "Constructor for cut with nonlinear variable list")
      // --- Modifiers (Append) ---
      .def("append", py::overload_cast<PV const&>(&PC::append), py::arg("x1"),
           py::return_value_policy::reference_internal,
           "Append variable dependency 'x1' to cut")
      .def("append", py::overload_cast<PV const&, double const&>(&PC::append),
           py::arg("x1"), py::arg("a1"),
           py::return_value_policy::reference_internal,
           "Append linear term 'a1·x1' to cut")
      .def(
          "append",
          [](PC& self, std::vector<PV> const& x,
             std::vector<double> const& a) -> PC&
          {
            if (!a.empty() && x.size() != a.size())
              throw std::invalid_argument("Size mismatch");
            return self.append(x.size(), x.data(),
                               a.empty() ? nullptr : a.data());
          },
          py::arg("x"), py::arg("a") = std::vector<double>(),
          py::return_value_policy::reference_internal,
          "Append linear terms 'a[0]·x[0]+a[1]·x[1]+···' to cut")
      .def("append",
           py::overload_cast<PV const&, PV const&, double const&>(&PC::append),
           py::arg("x1"), py::arg("x2"), py::arg("a12"),
           py::return_value_policy::reference_internal,
           "Append quadratic term 'a12·x1·x2' to cut")
      .def(
          "append",
          [](PC& self, std::vector<PV> const& x1, std::vector<PV> const& x2,
             std::vector<double> const& a12) -> PC&
          {
            if (x1.size() != x2.size() || x1.size() != a12.size())
              throw std::invalid_argument("Size mismatch");
            return self.append(x1.size(), x1.data(), x2.data(), a12.data());
          },
          py::arg("x1"), py::arg("x2"), py::arg("a12"),
          py::return_value_policy::reference_internal,
          "Append quadratic terms 'a12[0]·x1[0]·x2[0]+a12[1]·x1[1]·x2[1]+···' "
          "to cut")
      // --- Accessors ---
      .def_property_readonly("type", &PC::type, "Type of cut")
      .def_property_readonly("nvar", &PC::nvar, "Number of linear terms")
      .def_property_readonly("nqvar", &PC::nqvar, "Number of quadratic terms")
      .def_property_readonly("op", &PC::op, "Associated operation in DAG")
      .def_property_readonly(
          "coef",
          [](PC const& self) {
            return std::vector<double>(self.coef(), self.coef() + self.nvar());
          },
          "Coefficients in linear terms")
      .def_property_readonly(
          "var", [](PC const& self)
          { return std::vector<PV>(self.var(), self.var() + self.nvar()); },
          "Variables in linear terms")
      .def_property_readonly(
          "qcoef",
          [](PC const& self) {
            return std::vector<double>(self.qcoef(),
                                       self.qcoef() + self.nqvar());
          },
          "Coefficients in quadratic terms")
      .def_property_readonly(
          "qvar1",
          [](PC const& self) {
            return std::vector<PV>(self.qvar1(), self.qvar1() + self.nqvar());
          },
          "Left variable operands in quadratic terms")
      .def_property_readonly(
          "qvar2",
          [](PC const& self) {
            return std::vector<PV>(self.qvar2(), self.qvar2() + self.nqvar());
          },
          "Right variable operands in quadratic terms")
      .def_property(
          "rhs", py::overload_cast<>(&PC::rhs, py::const_),
          [](PC& self, double v) { self.rhs() = v; },
          "Right-hand side constant")
      // --- String Representation ---
      .def("__str__",
           [](PC const& C)
           {
             std::ostringstream Css;
             Css << C;
             return Css.str();
           })
      .def("__repr__",
           [](PC const& C)
           {
             std::ostringstream Css;
             Css << C;
             return Css.str();
           });

  py::class_<PI> pyPolImg(m, "PolImg");

  pyPolImg.def(py::init<>())
      .def_readwrite("options", &PI::options,
                     "options of polyhedral image propagation")
      .def_property_readonly("vars", py::overload_cast<>(&PI::Vars),
                             py::return_value_policy::reference_internal,
                             "dictionary of DAG variables")
      .def_property_readonly("aux", py::overload_cast<>(&PI::Aux),
                             py::return_value_policy::reference_internal,
                             "list of auxiliary variables")
      .def_property_readonly(
          "cuts",
          //&PI::Cuts,
          [](PI const& self)
          {
            auto const& cuts_set = self.Cuts();
            return std::vector<PC*>(cuts_set.begin(), cuts_set.end());
          },
          // py::return_value_policy::reference_internal,
          "list of cuts")
      .def("reset", &PI::reset, "reset polyheral image")
      .def("reset_cuts", &PI::reset_cuts, "erase all cuts and aux variables")
      .def("erase_cuts", py::overload_cast<>(&PI::erase_cuts), "erase all cuts")
      .def(
          "generate_cuts",
          [](PI& self, std::vector<PV> const& vdep, bool const reset)
          { self.generate_cuts(vdep, reset); }, py::arg("vdep"),
          py::arg("reset") = false,
          "append relaxation cuts for all 'vdep' dependents to the polyhedral "
          "image")
      .def(
          "generate_cuts",
          [](PI& self, std::set<unsigned> const& ndxdep,
             std::vector<PV> const& vdep, bool const reset)
          { self.generate_cuts(ndxdep, vdep, reset); },
          py::arg("ndxdep"), py::arg("vdep"), py::arg("reset") = false,
          "append relaxation cuts for selected 'vdep' dependents indexed by "
          "'ndxdep' to the polyhedral image")
      .def("__str__",
           [](PI const& P)
           {
             std::ostringstream Pss;
             Pss << P;
             return Pss.str();
           })
      .def("__repr__",
           [](PI const& P)
           {
             std::ostringstream Pss;
             Pss << P;
             return Pss.str();
           });

  py::class_<PI::Options> pyPolImgOptions(pyPolImg, "Options");

  pyPolImgOptions.def(py::init<>(), "Default constructor")
      .def(py::init<PI::Options const&>(), "Copy constructor")
      .def("reset", &PI::Options::reset, "Reset options to defaults")
      .def_readwrite("AGGREG_LQ", &PI::Options::AGGREG_LQ,
                     "Whether or not to aggregate linear expressions in cuts "
                     "[Default: True]")
      .def_readwrite("ROOT_USE", &PI::Options::ROOT_USE,
                     "Whether or not to use root search to construct envelopes "
                     "of univariate terms [Default: True]")
      .def_readwrite(
          "ROOT_MAXIT", &PI::Options::ROOT_MAXIT,
          "Maximal number of iterations in root search [Default: 100]")
      .def_readwrite("ROOT_TOL", &PI::Options::ROOT_TOL,
                     "Termination tolerance in root search [Default: 1e-10]")
      .def_readwrite("SANDWICH_ATOL", &PI::Options::SANDWICH_ATOL,
                     "Absolute tolerance in outer-approximation of univariate "
                     "terms [Default: 1e-10]")
      .def_readwrite("SANDWICH_RTOL", &PI::Options::SANDWICH_RTOL,
                     "Relative tolerance in outer-approximation of univariate "
                     "terms [Default: 1e-3]")
      .def_readwrite("SANDWICH_MAXCUT", &PI::Options::SANDWICH_MAXCUT,
                     "Maximal number of cuts in outer approximation of "
                     "univariate terms [Default: 5]")
      .def_readwrite("SANDWICH_RULE", &PI::Options::SANDWICH_RULE,
                     "Rule for outer-approximation of nonlinear convex/concave "
                     "terms [Default: MAXERR]")
      .def_readwrite("FRACTIONAL_ATOL", &PI::Options::FRACTIONAL_ATOL,
                     "Absolute tolerance to prevent division by zero")
      .def_readwrite("FRACTIONAL_RTOL", &PI::Options::FRACTIONAL_RTOL,
                     "Relative tolerance to prevent division by zero")
      .def_readwrite("BREAKPOINT_TYPE", &PI::Options::BREAKPOINT_TYPE,
                     "Rule for piecewise linear cuts of nonlinear "
                     "convex/concave terms [Default: NONE]")
      .def_readwrite("BREAKPOINT_ATOL", &PI::Options::BREAKPOINT_ATOL,
                     "Absolute tolerance in adding breakpoints in piecewise "
                     "linear cuts [Default: 1e-5]")
      .def_readwrite("BREAKPOINT_RTOL", &PI::Options::BREAKPOINT_RTOL,
                     "Relative tolerance in adding breakpoints in piecewise "
                     "linear cuts [Default: 1e-5]")
      .def_readwrite("ALLOW_QUAD", &PI::Options::ALLOW_QUAD,
                     "Whether or not to retain quadratic terms in polyhedral "
                     "relaxation [Default: False]")
      .def_readwrite("ALLOW_NLIN", &PI::Options::ALLOW_NLIN,
                     "Set of nonlinear terms to retain in polyhedral "
                     "relaxation [Default: {}]")
      .def_readwrite("ALLOW_DISJ", &PI::Options::ALLOW_DISJ,
                     "Set of disjunctive terms to retain in polyhedral "
                     "relaxation [Default: {}]");

  py::enum_<PI::Options::SANDWICH>(pyPolImgOptions, "SANDWICH_TYPE")
      .value("BISECT", PI::Options::SANDWICH::BISECT, "Range bisection")
      .value("MAXERR", PI::Options::SANDWICH::MAXERR, "Maximum error rule")
      .export_values();

  py::enum_<PI::Options::REFINE>(pyPolImgOptions, "REFINE_TYPE")
      .value("NONE", PI::Options::REFINE::NONE,
             "No semi-linear cuts (use secant approximation)")
      .value("CONT", PI::Options::REFINE::CONT,
             "Semilinear cuts with linear relaxed (continuous) reformulation")
      .value("BIN", PI::Options::REFINE::BIN,
             "Semilinear cuts with linear binary reformulation")
      .value("SOS2", PI::Options::REFINE::SOS2,
             "Semilinear cuts with SOS2 reformulation")
      .export_values();
}
