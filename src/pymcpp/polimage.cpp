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
  py::class_<PV> pyPolVar(m, "PolVar", R"doc(
Variable participating in a polyhedral image (`PolImg`).

A polyhedral image encloses the graph of a factorable function by
lifting it into a higher-dimensional space: every intermediate
operation in the function's DAG receives its own variable, and each
nonlinear operation is enclosed by linear cuts relating these
variables. A `PolVar` is one coordinate of that lifted space. It can
be:

- a DAG variable, attached to an `FFVar` of an `FFGraph` together
  with an interval range (the domain over which relaxation cuts are
  valid);
- an auxiliary variable, appended automatically to the image for each
  intermediate DAG operation when the DAG is evaluated in `PolVar`
  arithmetic via `FFGraph.eval`;
- a constant.

DAG variables are constructed explicitly by the user. The variable
range determines the tightness of the relaxation cuts subsequently
generated with `PolImg.generate_cuts`; both continuous and integer
(discrete) variables are supported via the `cnt` flag.

The `==` and `!=` operators compare variable identifiers (type and
index within the image, see `id`), not numerical values.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_var("X")
>>> IMG = pymcpp.PolImg()
>>> XP = pymcpp.PolVar(IMG, X, pymcpp.Interval(-2, 1), True)
>>> XP.name, XP.cnt
('X', True)
)doc");

  py::enum_<PV::TYPE>(pyPolVar, "TYPE", R"doc(
Kind of a `PolVar` within a polyhedral image.

DAG variables (VARCONT, VARINT) are attached to an `FFVar` of the
DAG; auxiliary variables (AUXCONT, AUXINT) are introduced by the
lifting of intermediate operations; AUXCST denotes a constant.
)doc")
      .value("VARCONT", PV::TYPE::VARCONT, "DAG continuous variable")
      .value("VARINT", PV::TYPE::VARINT, "DAG integer variable")
      .value("AUXCONT", PV::TYPE::AUXCONT, "Auxiliary continuous variable")
      .value("AUXINT", PV::TYPE::AUXINT, "Auxiliary integer variable")
      .value("AUXCST", PV::TYPE::AUXCST, "Auxiliary constant")
      .export_values();

  py::class_<PI> pyPolImg(m, "PolImg", R"doc(
Polyhedral relaxation environment for factorable functions.

`PolImg` constructs a polyhedral enclosure of the image set
``{ g(x) : xL <= x <= xU }`` of a factorable function ``g`` defined
in an `FFGraph` DAG. The enclosure is built in three steps:

1. Decomposition/lifting: each unary or binary operation in the DAG
   of ``g`` introduces an auxiliary variable ``v_k`` together with a
   defining relation, e.g. ``v_k = v_i*v_j`` or ``v_k = phi(v_i)``.
   Common subexpressions are shared, which tightens the relaxation.
2. Relaxation: each nonlinear relation is replaced by its polyhedral
   envelope over the variable ranges, e.g. the four McCormick
   inequalities for a bilinear term, or secant/tangent enclosures for
   convex, concave and convexo-concave univariates.
3. Outer approximation: the remaining nonlinear convex/concave parts
   are outer-approximated by supporting hyperplanes at finitely many
   points, selected iteratively by a sandwich algorithm whose rule,
   tolerances and maximal number of cuts are set in `PolImg.Options`.

The result is a finite set of cuts (`PolCut` objects, mostly linear
equalities/inequalities) between the lifted variables, suitable for
embedding in an LP or MILP model whose optimum bounds the function.
Options also allow keeping quadratic (`ALLOW_QUAD`), nonlinear
(`ALLOW_NLIN`) or disjunctive (`ALLOW_DISJ`) terms as such, and
adding piecewise-linear cuts at variable breakpoints
(`BREAKPOINT_TYPE`), which yields tighter MILP relaxations.

Typical workflow: create the image, attach `PolVar` variables for
the DAG independents with their interval ranges, evaluate the DAG in
`PolVar` arithmetic to lift the intermediates, then call
`generate_cuts` and read off `vars`, `aux` and `cuts`.

Examples
--------
>>> import pymcpp
>>> DAG = pymcpp.FFGraph()
>>> X = DAG.add_var("X")
>>> Y = DAG.add_var("Y")
>>> F = X * (pymcpp.exp(X) - Y) ** 2
>>> IMG = pymcpp.PolImg()
>>> XP = pymcpp.PolVar(IMG, X, pymcpp.Interval(-2, 1), True)
>>> YP = pymcpp.PolVar(IMG, Y, pymcpp.Interval(-1, 2), True)
>>> [FP] = DAG.eval([F], [X, Y], [XP, YP])  # lift intermediates
>>> IMG.generate_cuts([FP])                 # build relaxation cuts
>>> len(IMG.cuts) > 0
True
)doc");

  pyPolVar
      // --- Static Constants ---
      .def_property_static(
          "VARCONTNAME", [](py::object) { return PV::VARCONTNAME; },
          [](py::object, std::string const& str) { PV::VARCONTNAME = str; },
          R"doc(
Class-wide name prefix for unnamed continuous DAG variables.

Used by `PolVar.name` when the underlying `FFVar` has no user-defined
name. Default is "V".
)doc")
      .def_property_static(
          "VARINTNAME", [](py::object) { return PV::VARINTNAME; },
          [](py::object, std::string const& str) { PV::VARINTNAME = str; },
          R"doc(
Class-wide name prefix for unnamed integer DAG variables.

Default is "Y".
)doc")
      .def_property_static(
          "AUXCONTNAME", [](py::object) { return PV::AUXCONTNAME; },
          [](py::object, std::string const& str) { PV::AUXCONTNAME = str; },
          R"doc(
Class-wide name prefix for continuous auxiliary variables.

Default is "W".
)doc")
      .def_property_static(
          "AUXINTNAME", [](py::object) { return PV::AUXINTNAME; },
          [](py::object, std::string const& str) { PV::AUXINTNAME = str; },
          R"doc(
Class-wide name prefix for integer auxiliary variables.

Default is "Z".
)doc")
      .def_property_static(
          "AUXCSTNAME", [](py::object) { return PV::AUXCSTNAME; },
          [](py::object, std::string const& str) { PV::AUXCSTNAME = str; },
          R"doc(
Class-wide name prefix for auxiliary constants.

Default is "C".
)doc")

      // --- Constructors ---
      .def(py::init<double const&>(), py::arg("d") = 0., R"doc(
Construct a constant `PolVar` holding the value ``d``.

The variable has type ``AUXCST`` and belongs to no polyhedral image.

Parameters
----------
d : float, optional
    Constant value. Default is 0.
)doc")
      .def(py::init<int const>(), py::arg("n"), R"doc(
Construct a constant `PolVar` holding the integer value ``n``.

Parameters
----------
n : int
    Constant value.
)doc")
      .def(py::init<PI*, mc::FFVar const&, I const&, bool const>(),
           py::arg("img"), py::arg("var"), py::arg("rng") = 0.,
           py::arg("cnt") = true, R"doc(
Construct a DAG variable in a polyhedral image.

Registers the DAG variable ``var`` in image ``img`` with interval
range ``rng``. If ``var`` is already registered in ``img``, its range
and continuity flag are updated instead.

Parameters
----------
img : PolImg
    Polyhedral image the variable participates in.
var : FFVar
    Underlying DAG variable.
rng : Interval, optional
    Range of the variable, over which relaxation cuts are valid.
    Default is the degenerate interval [0,0].
cnt : bool, optional
    True for a continuous variable, False for an integer (discrete)
    variable. Default is True.
)doc")
      .def(py::init<PI*, I const&, bool const>(), py::arg("img"),
           py::arg("rng") = 0., py::arg("cnt") = true, R"doc(
Construct an auxiliary variable in a polyhedral image.

Appends a new auxiliary variable (not attached to any DAG variable)
to image ``img`` with range ``rng``.

Parameters
----------
img : PolImg
    Polyhedral image the variable is appended to.
rng : Interval, optional
    Range of the variable. Default is the degenerate interval [0,0].
cnt : bool, optional
    True for continuous, False for integer. Default is True.
)doc")
      .def(py::init<PV const&>(), R"doc(
Copy constructor.
)doc")
      // --- Modifiers ---
      .def("update", py::overload_cast<I const&>(&PV::update), py::arg("rng"),
           py::return_value_policy::reference_internal, R"doc(
Update the range of the variable in its polyhedral image.

Also resets the variable subdivision and its cuts flag. No effect if
the variable does not belong to an image.

Parameters
----------
rng : Interval
    New variable range.

Returns
-------
PolVar
    This variable (for chaining).
)doc")
      .def("update", py::overload_cast<bool const>(&PV::update), py::arg("cnt"),
           py::return_value_policy::reference_internal, R"doc(
Update the continuity (continuous/integer) of the variable.

Same as above but changes the variable type instead of the range.

Parameters
----------
cnt : bool
    True for continuous, False for integer.

Returns
-------
PolVar
    This variable (for chaining).
)doc")
      .def("set",
           py::overload_cast<PI*, mc::FFVar const&, I const&, bool const>(
               &PV::set),
           py::arg("img"), py::arg("var"), py::arg("rng") = 0.,
           py::arg("cnt") = true, py::return_value_policy::reference_internal,
           R"doc(
Re-initialize this object as DAG variable ``var`` in image ``img``.

Equivalent to the DAG-variable constructor; existing breakpoints,
subdivision and cuts flag are cleared.

Parameters
----------
img : PolImg
    Polyhedral image the variable participates in.
var : FFVar
    Underlying DAG variable.
rng : Interval, optional
    Variable range. Default is the degenerate interval [0,0].
cnt : bool, optional
    True for continuous, False for integer. Default is True.

Returns
-------
PolVar
    This variable (for chaining).
)doc")
      .def("set", py::overload_cast<PI*, I const&, bool const>(&PV::set),
           py::arg("img"), py::arg("rng") = 0., py::arg("cnt") = true,
           py::return_value_policy::reference_internal, R"doc(
Re-initialize this object as a new auxiliary variable in ``img``.

Parameters
----------
img : PolImg
    Polyhedral image the variable is appended to.
rng : Interval, optional
    Variable range. Default is the degenerate interval [0,0].
cnt : bool, optional
    True for continuous, False for integer. Default is True.

Returns
-------
PolVar
    This variable (for chaining).
)doc")
      .def("add_breakpt", &PV::add_breakpt, py::arg("bkpt"), R"doc(
Add a breakpoint to the variable.

Breakpoints subdivide the variable range for piecewise-linear
(semilinear) cuts, controlled by `PolImg.Options.BREAKPOINT_TYPE`.
Points outside the variable range, or closer to an existing
breakpoint or range bound than the tolerances
`PolImg.Options.BREAKPOINT_ATOL`/`BREAKPOINT_RTOL`, are ignored.

Parameters
----------
bkpt : float
    Breakpoint location.
)doc")
      .def("reset_subdiv", &PV::reset_subdiv, R"doc(
Clear the current subdivision (points and auxiliary variables).
)doc")
      .def("reset_cuts", &PV::reset_cuts, R"doc(
Mark the operation associated with this variable as not yet cut.

A subsequent `PolImg.generate_cuts` call will generate cuts for it
again.
)doc")
      .def("set_cuts", &PV::set_cuts, R"doc(
Mark the operation associated with this variable as already cut.

`PolImg.generate_cuts` will skip cut generation for it.
)doc")
      .def("create_subdiv", &PV::create_subdiv, py::arg("lbd"), py::arg("ubd"),
           py::arg("reset") = false, R"doc(
Create the subdivision of interval [lbd, ubd] at the breakpoints.

Collects ``lbd``, the breakpoints lying strictly inside the interval,
and ``ubd`` into an ordered list of subdivision points.

Parameters
----------
lbd : float
    Lower bound of the subdivided interval.
ubd : float
    Upper bound of the subdivided interval.
reset : bool, optional
    Whether to discard any previously created subdivision first.
    Default is False.

Returns
-------
list of float
    The subdivision points (empty if the variable has no image).
)doc")
      .def("sos2_subdiv", &PV::SOS2_subdiv, py::arg("op") = nullptr,
           py::arg("reset") = false, R"doc(
Encode the variable subdivision with an SOS2 constraint.

Appends auxiliary weighting variables and cuts to the image so that
the variable is expressed as a convex combination of its subdivision
points with an SOS2 restriction (at most two adjacent nonzero
weights). Requires a subdivision created with `create_subdiv`.

Parameters
----------
op : FFOp, optional
    DAG operation the generated cuts are attributed to. Default is
    None.
reset : bool, optional
    Whether to rebuild the encoding from scratch. Default is False.

Returns
-------
list of PolVar
    The auxiliary weighting variables introduced.
)doc")
      .def("bin_subdiv", &PV::BIN_subdiv, py::arg("op") = nullptr,
           py::arg("reset") = false, R"doc(
Encode the variable subdivision with binary variables.

Same purpose as `sos2_subdiv` but uses an incremental (delta-method)
formulation with auxiliary binary variables selecting the active
subinterval.

Parameters
----------
op : FFOp, optional
    DAG operation the generated cuts are attributed to. Default is
    None.
reset : bool, optional
    Whether to rebuild the encoding from scratch. Default is False.

Returns
-------
list of PolVar
    The auxiliary (continuous increment) variables introduced; the
    binary selection variables are appended to the image as well.
)doc")
      .def("cnt_subdiv", &PV::CONT_subdiv, py::arg("pOp") = nullptr,
           py::arg("reset") = false, R"doc(
Encode the variable subdivision with a relaxed continuous encoding.

Same as `bin_subdiv` but with the binary selection variables relaxed
to continuous variables in [0,1], giving an LP-representable (weaker)
subdivision.

Parameters
----------
pOp : FFOp, optional
    DAG operation the generated cuts are attributed to. Default is
    None.
reset : bool, optional
    Whether to rebuild the encoding from scratch. Default is False.

Returns
-------
list of PolVar
    The auxiliary (continuous increment) variables introduced.
)doc")
      // --- Accessors ---
      .def_property_readonly("name", &PV::name, R"doc(
str: Name of the variable.

The name of the underlying DAG variable if user-defined; otherwise a
default name made of the type prefix (see `VARCONTNAME` etc.) and the
variable index.
)doc")
      .def_property_readonly(
          "cnt", [](PV const& self) { return !self.discr(); }, R"doc(
bool: Whether the variable is continuous (False for integer).
)doc")
      .def_property_readonly("cst", &PV::cst, R"doc(
bool: Whether the variable is a constant (type AUXCST).
)doc")
      .def_property_readonly("rng", &PV::range, R"doc(
Interval: Range of the variable.

For constants, the degenerate interval at the constant value.
)doc")
      .def_property_readonly("id", &PV::id, R"doc(
tuple of (PolVar.TYPE, int): Unique identifier of the variable.

The first element is the variable kind, the second its index within
the variables of that kind in the image.
)doc")
      .def_property_readonly("img", &PV::image,
                             py::return_value_policy::reference_internal,
                             R"doc(
PolImg: Polyhedral image the variable belongs to, or None for
constants.
)doc")
      .def_property_readonly("var", py::overload_cast<>(&PV::var),
                             py::return_value_policy::reference_internal,
                             R"doc(
FFVar: Underlying DAG variable.

Only valid for DAG variables (types VARCONT/VARINT); accessing it on
an auxiliary or constant is undefined behavior.
)doc")
      .def_property_readonly("breakpts", &PV::breakpts, R"doc(
set of float: Breakpoints added to the variable with `add_breakpt`.
)doc")
      .def_property_readonly("subdiv", &PV::subdiv, R"doc(
tuple of (list of float, list of PolVar): Current subdivision.

Returns
-------
points : list of float
    Subdivision points created by `create_subdiv`.
vars : list of PolVar
    Auxiliary variables introduced by `sos2_subdiv`, `bin_subdiv` or
    `cnt_subdiv` (empty before any encoding is applied).
)doc")
      .def("has_cuts", &PV::has_cuts, R"doc(
Return whether cuts were already generated for this variable.

Returns
-------
bool
    True if the operation associated with this variable has been
    processed by `PolImg.generate_cuts`.
)doc")
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
          py::is_operator(),
          "Compare variable identifiers (type and index), not values.")
      .def(
          "__ne__", [](PV const& x, PV const& y) { return x.id() != y.id(); },
          py::is_operator(),
          "Compare variable identifiers (type and index), not values.")
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
  py::class_<PC> pyPolCut(m, "PolCut", R"doc(
Single cut in a polyhedral image.

A cut of type ``EQ``, ``LE`` or ``GE`` is a constraint between
polyhedral-image variables of the form::

    sum_k coef[k]*var[k]
      + sum_k qcoef[k]*qvar1[k]*qvar2[k]   {=, <=, >=}   rhs

Quadratic terms only occur when `PolImg.Options.ALLOW_QUAD` is
enabled. Two further cut types encode non-algebraic restrictions:
``SOS1``/``SOS2`` cuts declare that the variables in `var` (with
weights `coef`) form a special-ordered set of type 1 or 2, and
``NLIN`` cuts record a nonlinear relation ``y = f(x)`` between the
listed variables, kept verbatim when the corresponding operation is
included in `PolImg.Options.ALLOW_NLIN` or `ALLOW_DISJ`. The DAG
operation a cut originates from is available as `op`.

Cuts are normally created by `PolImg.generate_cuts` and retrieved
through `PolImg.cuts`; the constructors below allow building cuts
manually.

Examples
--------
>>> for cut in IMG.cuts:  # doctest: +SKIP
...     print(cut.type, list(zip(cut.coef, cut.var)), cut.rhs)
)doc");

  // Enumeration for PolCut::TYPE
  py::enum_<PC::TYPE>(pyPolCut, "TYPE", R"doc(
Type of a `PolCut` constraint.
)doc")
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
           py::arg("type"), py::arg("b"), R"doc(
Construct a cut without any participating variable.

Terms can be added afterwards with `append`.

Parameters
----------
op : FFOp
    DAG operation the cut is attributed to (may be None).
type : PolCut.TYPE
    Cut type (relation to the right-hand side).
b : float
    Right-hand side constant.
)doc")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), R"doc(
Construct the cut a1*x1 {=, <=, >=} b with one linear term.
)doc")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"), R"doc(
Construct the cut a1*x1 + a2*x2 {=, <=, >=} b with two linear terms.
)doc")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&, PV const&,
                    double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"), py::arg("x3"),
           py::arg("a3"), R"doc(
Construct the cut a1*x1 + a2*x2 + a3*x3 {=, <=, >=} b.
)doc")
      .def(py::init<mc::FFOp const*, PC::TYPE, double const&, PV const&,
                    double const&, PV const&, double const&, PV const&,
                    double const&, PV const&, double const&>(),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x1"),
           py::arg("a1"), py::arg("x2"), py::arg("a2"), py::arg("x3"),
           py::arg("a3"), py::arg("x4"), py::arg("a4"), R"doc(
Construct the cut a1*x1 + a2*x2 + a3*x3 + a4*x4 {=, <=, >=} b.
)doc")
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
           R"doc(
Construct the cut sum_k a[k]*x[k] {=, <=, >=} b from parallel lists.

Parameters
----------
op : FFOp
    DAG operation the cut is attributed to (may be None).
type : PolCut.TYPE
    Cut type.
b : float
    Right-hand side constant.
x : list of PolVar
    Participating variables.
a : list of float
    Coefficients, one per variable in ``x``.

Raises
------
ValueError
    If ``x`` and ``a`` have different lengths.
)doc")
      .def(py::init([](mc::FFOp const* op, PC::TYPE type, double b,
                       std::vector<PV> const& x, double const& a0)
                    { return new PC(op, type, b, x.size(), x.data(), a0); }),
           py::arg("op"), py::arg("type"), py::arg("b"), py::arg("x"),
           py::arg("a0"), R"doc(
Construct the cut a0*(x[0] + x[1] + ...) {=, <=, >=} b, with the same
coefficient ``a0`` on every variable.
)doc")
      .def(py::init<mc::FFOp const*, PV const&>(), py::arg("op"), py::arg("x1"),
           R"doc(
Construct a nonlinear (NLIN) cut with one participating variable.

Records the dependency of ``x1`` on operation ``op``.
)doc")
      .def(py::init<mc::FFOp const*, PV const&, PV const&>(), py::arg("op"),
           py::arg("x1"), py::arg("x2"), R"doc(
Construct a nonlinear (NLIN) cut with two participating variables.
)doc")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, double const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("b"), R"doc(
Construct a nonlinear (NLIN) cut with two variables and a constant.
)doc")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, PV const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("x3"), R"doc(
Construct a nonlinear (NLIN) cut with three participating variables.
)doc")
      .def(py::init<mc::FFOp const*, PV const&, PV const&, PV const&>(),
           py::arg("op"), py::arg("x1"), py::arg("x2"), py::arg("x3"), R"doc(
Same as the previous overload (three participating variables).
)doc")
      .def(py::init(
               [](mc::FFOp const* op, PV const& X1, std::vector<PV> const& X)
               { return new PC(op, X1, X.size(), X.data()); }),
           py::arg("op"), py::arg("x1"), py::arg("x"), R"doc(
Construct a nonlinear (NLIN) cut relating ``x1`` and the variable
list ``x``.
)doc")
      // --- Modifiers (Append) ---
      .def("append", py::overload_cast<PV const&>(&PC::append), py::arg("x1"),
           py::return_value_policy::reference_internal, R"doc(
Append variable ``x1`` as a dependency with zero coefficient.

Returns
-------
PolCut
    This cut (for chaining).
)doc")
      .def("append", py::overload_cast<PV const&, double const&>(&PC::append),
           py::arg("x1"), py::arg("a1"),
           py::return_value_policy::reference_internal, R"doc(
Append the linear term a1*x1 to the cut.

Returns
-------
PolCut
    This cut (for chaining).
)doc")
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
          py::return_value_policy::reference_internal, R"doc(
Append the linear terms a[0]*x[0] + a[1]*x[1] + ... to the cut.

Parameters
----------
x : list of PolVar
    Variables to append.
a : list of float, optional
    Coefficients, one per variable. If omitted or empty, the
    variables are appended with zero coefficients (pure
    dependencies).

Returns
-------
PolCut
    This cut (for chaining).

Raises
------
ValueError
    If ``a`` is nonempty and its length differs from that of ``x``.
)doc")
      .def("append",
           py::overload_cast<PV const&, PV const&, double const&>(&PC::append),
           py::arg("x1"), py::arg("x2"), py::arg("a12"),
           py::return_value_policy::reference_internal, R"doc(
Append the quadratic term a12*x1*x2 to the cut.

Returns
-------
PolCut
    This cut (for chaining).
)doc")
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
          py::return_value_policy::reference_internal, R"doc(
Append the quadratic terms a12[k]*x1[k]*x2[k], k = 0, 1, ..., to the
cut.

Parameters
----------
x1 : list of PolVar
    Left factors.
x2 : list of PolVar
    Right factors (same length as ``x1``).
a12 : list of float
    Coefficients (same length as ``x1``).

Returns
-------
PolCut
    This cut (for chaining).

Raises
------
ValueError
    If the three lists have different lengths.
)doc")
      // --- Accessors ---
      .def_property_readonly("type", &PC::type, R"doc(
PolCut.TYPE: Type of the cut (relation to the right-hand side).
)doc")
      .def_property_readonly("nvar", &PC::nvar, R"doc(
int: Number of linear terms in the cut.
)doc")
      .def_property_readonly("nqvar", &PC::nqvar, R"doc(
int: Number of quadratic terms in the cut.
)doc")
      .def_property_readonly("op", &PC::op, R"doc(
FFOp: DAG operation the cut originates from, or None.
)doc")
      .def_property_readonly(
          "coef",
          [](PC const& self) {
            return std::vector<double>(self.coef(), self.coef() + self.nvar());
          },
          R"doc(
list of float: Coefficients of the linear terms (length `nvar`).

Position k corresponds to variable ``var[k]``.
)doc")
      .def_property_readonly(
          "var", [](PC const& self)
          { return std::vector<PV>(self.var(), self.var() + self.nvar()); },
          R"doc(
list of PolVar: Variables of the linear terms (length `nvar`).
)doc")
      .def_property_readonly(
          "qcoef",
          [](PC const& self) {
            return std::vector<double>(self.qcoef(),
                                       self.qcoef() + self.nqvar());
          },
          R"doc(
list of float: Coefficients of the quadratic terms (length `nqvar`).

Position k corresponds to the product ``qvar1[k]*qvar2[k]``.
)doc")
      .def_property_readonly(
          "qvar1",
          [](PC const& self) {
            return std::vector<PV>(self.qvar1(), self.qvar1() + self.nqvar());
          },
          R"doc(
list of PolVar: Left factors of the quadratic terms (length `nqvar`).
)doc")
      .def_property_readonly(
          "qvar2",
          [](PC const& self) {
            return std::vector<PV>(self.qvar2(), self.qvar2() + self.nqvar());
          },
          R"doc(
list of PolVar: Right factors of the quadratic terms (length
`nqvar`).
)doc")
      .def_property(
          "rhs", py::overload_cast<>(&PC::rhs, py::const_),
          [](PC& self, double v) { self.rhs() = v; },
          R"doc(
float: Right-hand side constant of the cut (readable and writable).
)doc")
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


  py::class_<PI::Options> pyPolImgOptions(pyPolImg, "Options", R"doc(
Options controlling polyhedral cut generation in a `PolImg`.

An instance is available as the `PolImg.options` attribute of every
image; modify its fields before calling `PolImg.generate_cuts`.

Examples
--------
>>> IMG = pymcpp.PolImg()
>>> IMG.options.SANDWICH_MAXCUT = 7
>>> IMG.options.BREAKPOINT_TYPE = pymcpp.PolImg.Options.REFINE_TYPE.SOS2
)doc");

  py::enum_<PI::Options::SANDWICH>(pyPolImgOptions, "SANDWICH_TYPE", R"doc(
Rule for selecting new linearization points in the sandwich
algorithm that outer-approximates convex/concave univariate terms.
)doc")
      .value("BISECT", PI::Options::SANDWICH::BISECT, "Range bisection")
      .value("MAXERR", PI::Options::SANDWICH::MAXERR, "Maximum error rule")
      .export_values();

  py::enum_<PI::Options::REFINE>(pyPolImgOptions, "REFINE_TYPE", R"doc(
Reformulation used for piecewise-linear (semilinear) cuts at
variable breakpoints (see `PolVar.add_breakpt`).
)doc")
      .value("NONE", PI::Options::REFINE::NONE,
             "No semi-linear cuts (use secant approximation)")
      .value("CONT", PI::Options::REFINE::CONT,
             "Semilinear cuts with linear relaxed (continuous) reformulation")
      .value("BIN", PI::Options::REFINE::BIN,
             "Semilinear cuts with linear binary reformulation")
      .value("SOS2", PI::Options::REFINE::SOS2,
             "Semilinear cuts with SOS2 reformulation")
      .export_values();

  pyPolImg.def(py::init<>(), R"doc(
Construct an empty polyhedral image.
)doc")
      .def_readwrite("options", &PI::options, R"doc(
PolImg.Options: Options controlling cut generation for this image.
)doc")
      .def_property_readonly("vars", py::overload_cast<>(&PI::Vars),
                             py::return_value_policy::reference_internal,
                             R"doc(
dict: DAG variables registered in the image.

Maps the identifier of each underlying `FFVar` (a tuple
``(FFVar.TYPE, int)``) to the corresponding `PolVar`.
)doc")
      .def_property_readonly("aux", py::overload_cast<>(&PI::Aux),
                             py::return_value_policy::reference_internal,
                             R"doc(
list of PolVar: Auxiliary variables appended to the image, in
creation order.
)doc")
      .def_property_readonly(
          "cuts",
          //&PI::Cuts,
          [](PI const& self)
          {
            auto const& cuts_set = self.Cuts();
            return std::vector<PC*>(cuts_set.begin(), cuts_set.end());
          },
          // py::return_value_policy::reference_internal,
          R"doc(
list of PolCut: Cuts currently in the image.

A new list is returned on each access; it is empty until
`generate_cuts` is called. See `PolCut` for the structure of each
cut and the notebook example for exporting cuts to an LP/MILP
solver.
)doc")
      .def("reset", &PI::reset, R"doc(
Clear the image entirely: all DAG variables, auxiliary variables and
cuts are erased.
)doc")
      .def("reset_cuts", &PI::reset_cuts, R"doc(
Erase all cuts and auxiliary variables.

The registered DAG variables are kept, but their subdivisions and
cuts flags are reset, so cuts can be regenerated from scratch.
)doc")
      .def("erase_cuts", py::overload_cast<>(&PI::erase_cuts), R"doc(
Erase all cuts, keeping every variable (DAG and auxiliary).
)doc")
      .def(
          "generate_cuts",
          [](PI& self, std::vector<PV> const& vdep, bool const reset)
          { self.generate_cuts(vdep, reset); }, py::arg("vdep"),
          py::arg("reset") = false, R"doc(
Append relaxation cuts for the dependents ``vdep`` to the image.

Traverses the operations defining the dependents backward through
the DAG and generates the polyhedral relaxation cuts of every
operation not yet processed (see `PolVar.has_cuts`). The dependents
must have been obtained beforehand by evaluating the DAG in `PolVar`
arithmetic (`FFGraph.eval`).

Parameters
----------
vdep : list of PolVar
    Dependent variables (function outputs) to relax.
reset : bool, optional
    If True, call `reset_cuts` first, discarding all existing cuts
    and auxiliaries. Default is False.
)doc")
      .def(
          "generate_cuts",
          [](PI& self, std::set<unsigned> const& ndxdep,
             std::vector<PV> const& vdep, bool const reset)
          { self.generate_cuts(ndxdep, vdep, reset); },
          py::arg("ndxdep"), py::arg("vdep"), py::arg("reset") = false, R"doc(
Append relaxation cuts for selected dependents only.

Same as above, but restricted to the dependents ``vdep[i]`` with
index ``i`` in ``ndxdep``.

Parameters
----------
ndxdep : set of int
    Indices into ``vdep`` selecting which dependents to relax.
vdep : list of PolVar
    Dependent variables (function outputs).
reset : bool, optional
    If True, call `reset_cuts` first. Default is False.
)doc")
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


  pyPolImgOptions.def(py::init<>(), R"doc(
Construct an option set with default values.
)doc")
      .def(py::init<PI::Options const&>(), R"doc(
Copy constructor.
)doc")
      .def("reset", &PI::Options::reset, R"doc(
Reset all options to their default values.
)doc")
      .def_readwrite("AGGREG_LQ", &PI::Options::AGGREG_LQ, R"doc(
bool: Whether to aggregate linear subexpressions into single cuts
instead of introducing one auxiliary variable per operation.
Default is False.
)doc")
      .def_readwrite("ROOT_USE", &PI::Options::ROOT_USE, R"doc(
bool: Whether to use root search for the junction points of the
envelopes of convexo-concave univariate terms (e.g. odd powers,
sin/cos, tanh). If False, weaker interval-based enclosures are used.
Default is True.
)doc")
      .def_readwrite(
          "ROOT_MAXIT", &PI::Options::ROOT_MAXIT, R"doc(
int: Maximal number of iterations in the envelope root search
(0 means no limit). Default is 100.
)doc")
      .def_readwrite("ROOT_TOL", &PI::Options::ROOT_TOL, R"doc(
float: Termination tolerance of the envelope root search.
Default is 1e-10.
)doc")
      .def_readwrite("SANDWICH_ATOL", &PI::Options::SANDWICH_ATOL, R"doc(
float: Absolute tolerance on the maximal outer-approximation error
in the sandwich algorithm for convex/concave univariate terms.
Default is 1e-10.
)doc")
      .def_readwrite("SANDWICH_RTOL", &PI::Options::SANDWICH_RTOL, R"doc(
float: Relative tolerance on the maximal outer-approximation error
in the sandwich algorithm. Default is 1e-3.
)doc")
      .def_readwrite("SANDWICH_MAXCUT", &PI::Options::SANDWICH_MAXCUT, R"doc(
int: Maximal number of cuts generated per convex/concave nonlinear
constraint by the sandwich algorithm. Default is 5.
)doc")
      .def_readwrite("SANDWICH_RULE", &PI::Options::SANDWICH_RULE, R"doc(
PolImg.Options.SANDWICH_TYPE: Rule for choosing new linearization
points in the sandwich algorithm: BISECT (interval bisection) or
MAXERR (subdivide at the maximum-error point).
Default is SANDWICH_TYPE.MAXERR.
)doc")
      .def_readwrite("FRACTIONAL_ATOL", &PI::Options::FRACTIONAL_ATOL, R"doc(
float: Absolute tolerance to prevent division by zero when relaxing
fractional terms. Default is machine epsilon (about 2.2e-16).
)doc")
      .def_readwrite("FRACTIONAL_RTOL", &PI::Options::FRACTIONAL_RTOL, R"doc(
float: Relative tolerance to prevent division by zero when relaxing
fractional terms. Default is machine epsilon (about 2.2e-16).
)doc")
      .def_readwrite("BREAKPOINT_TYPE", &PI::Options::BREAKPOINT_TYPE, R"doc(
PolImg.Options.REFINE_TYPE: Reformulation used for piecewise-linear
cuts at variable breakpoints: NONE (no semilinear cuts), CONT
(relaxed continuous encoding), BIN (binary encoding) or SOS2 (SOS2
encoding). Only takes effect for variables with breakpoints (see
`PolVar.add_breakpt`). Default is REFINE_TYPE.BIN.
)doc")
      .def_readwrite("BREAKPOINT_ATOL", &PI::Options::BREAKPOINT_ATOL, R"doc(
float: Absolute tolerance below which a new breakpoint is considered
duplicate of an existing point and discarded. Default is 1e-5.
)doc")
      .def_readwrite("BREAKPOINT_RTOL", &PI::Options::BREAKPOINT_RTOL, R"doc(
float: Relative tolerance below which a new breakpoint is considered
duplicate of an existing point and discarded. Default is 1e-5.
)doc")
      .def_readwrite("ALLOW_QUAD", &PI::Options::ALLOW_QUAD, R"doc(
bool: Whether to keep quadratic terms (squares and bilinear
products) as quadratic cuts instead of relaxing them polyhedrally.
The resulting model is a (potentially nonconvex) QP/QCP rather than
an LP. Default is False.
)doc")
      .def_readwrite("ALLOW_NLIN", &PI::Options::ALLOW_NLIN, R"doc(
set of FFOp.TYPE: Nonlinear operations to keep verbatim as NLIN cuts
(e.g. {FFOp.TYPE.EXP}) instead of relaxing them polyhedrally, for
use with solvers supporting these constraints natively.
Default is the empty set.
)doc")
      .def_readwrite("ALLOW_DISJ", &PI::Options::ALLOW_DISJ, R"doc(
set of FFOp.TYPE: Disjunctive operations (min, max, abs, fstep) to
keep verbatim as NLIN cuts instead of relaxing them polyhedrally.
Default is the empty set.
)doc");


}
