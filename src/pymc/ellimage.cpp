#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

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

#include "ellimage.hpp"

namespace py = pybind11;

namespace
{
using ELL = mc::Ellipsoid;
using EI  = mc::EllImg<I>;
using EV  = mc::EllVar<I>;
using HP  = std::pair<arma::vec,double>;

arma::vec
_to_arma_vec( py::handle obj )
{
  auto a = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure( obj );
  if( !a ) throw py::type_error( "expected a one-dimensional array-like object" );
  py::buffer_info info = a.request();

  arma::uword n = 0;
  if( info.ndim == 1 ){
    n = static_cast<arma::uword>( info.shape[0] );
  }
  else if( info.ndim == 2 && ( info.shape[0] == 1 || info.shape[1] == 1 ) ){
    n = static_cast<arma::uword>( info.shape[0] * info.shape[1] );
  }
  else{
    throw py::value_error( "expected a vector with shape (n,), (n,1), or (1,n)" );
  }

  arma::vec v( n );
  const double* data = static_cast<const double*>( info.ptr );
  for( arma::uword i=0; i<n; ++i ) v(i) = data[i];
  return v;
}

arma::mat
_to_arma_mat( py::handle obj )
{
  auto a = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure( obj );
  if( !a ) throw py::type_error( "expected a two-dimensional array-like object" );
  py::buffer_info info = a.request();

  if( info.ndim != 2 ) throw py::value_error( "expected a matrix with shape (m,n)" );
  const arma::uword nr = static_cast<arma::uword>( info.shape[0] );
  const arma::uword nc = static_cast<arma::uword>( info.shape[1] );

  arma::mat M( nr, nc );
  const double* data = static_cast<const double*>( info.ptr );
  for( arma::uword i=0; i<nr; ++i )
    for( arma::uword j=0; j<nc; ++j )
      M(i,j) = data[i*nc+j];
  return M;
}

py::array_t<double>
_from_arma_vec( arma::vec const& v )
{
  py::array_t<double> out( { static_cast<py::ssize_t>( v.n_elem ) } );
  auto buf = out.mutable_unchecked<1>();
  for( arma::uword i=0; i<v.n_elem; ++i ) buf( static_cast<py::ssize_t>( i ) ) = v(i);
  return out;
}

py::array_t<double>
_from_arma_mat( arma::mat const& M )
{
  py::array_t<double> out( { static_cast<py::ssize_t>( M.n_rows ), static_cast<py::ssize_t>( M.n_cols ) } );
  auto buf = out.mutable_unchecked<2>();
  for( arma::uword i=0; i<M.n_rows; ++i )
    for( arma::uword j=0; j<M.n_cols; ++j )
      buf( static_cast<py::ssize_t>( i ), static_cast<py::ssize_t>( j ) ) = M(i,j);
  return out;
}

py::array_t<double>
_from_arma_spmat( arma::sp_mat const& S )
{
  return _from_arma_mat( arma::mat( S ) );
}

HP
_to_hp( py::handle obj )
{
  py::tuple t = py::reinterpret_borrow<py::tuple>( obj );
  if( !py::isinstance<py::tuple>( obj ) || t.size() != 2 )
    throw py::type_error( "expected a halfspace/hyperplane as (normal, rhs)" );
  return HP( _to_arma_vec( t[0] ), py::cast<double>( t[1] ) );
}

std::vector<HP>
_to_hps( py::iterable obj )
{
  std::vector<HP> out;
  for( py::handle item : obj ) out.push_back( _to_hp( item ) );
  return out;
}

std::vector<EV>
_to_ev_vector( py::iterable obj )
{
  std::vector<EV> out;
  for( py::handle item : obj ) out.push_back( py::cast<EV>( item ) );
  return out;
}

std::vector<arma::vec>
_to_vec_vector( py::iterable obj )
{
  std::vector<arma::vec> out;
  for( py::handle item : obj ) out.push_back( _to_arma_vec( item ) );
  return out;
}

std::string
_to_string( ELL const& E )
{
  std::ostringstream os;
  os << E;
  return os.str();
}

std::string
_to_string( EI const& E )
{
  std::ostringstream os;
  os << E;
  return os.str();
}

std::string
_to_string( EV const& V )
{
  std::ostringstream os;
  os << V;
  return os.str();
}
} // namespace

void
mc_ellimage( py::module &m )
{

// Surface the C++ ellipsoid/ellimage exceptions -- which are NOT derived from
// std::exception -- as informative Python exceptions, instead of the opaque
// "Caught an unknown exception!" that pybind's default handler would raise.
static py::exception<ELL::Exceptions> pyEllipsoidError( m, "EllipsoidError" );
static py::exception<EI::Exceptions>  pyEllImgError   ( m, "EllImgError" );
py::register_exception_translator( []( std::exception_ptr p ){
  try{ if( p ) std::rethrow_exception( p ); }
  catch( EI::Exceptions&  e ){ py::set_error( pyEllImgError,    e.what().c_str() ); }
  catch( ELL::Exceptions& e ){ py::set_error( pyEllipsoidError, e.what().c_str() ); }
} );


py::class_<ELL> pyEllipsoid( m, "Ellipsoid" );

pyEllipsoid
 .def( py::init<>() )
 .def(
   py::init( []( py::handle Q, py::object c ){
     return ELL( _to_arma_mat( Q ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   } ),
   py::arg("Q"),
   py::arg("c") = py::none(),
   "constructor from dense shape matrix Q and optional center c"
 )
 .def(
   py::init( []( unsigned const n, std::vector<double> const& Qlt, py::object c ){
     if( Qlt.size() != static_cast<size_t>( n*(n+1)/2 ) )
       throw py::value_error( "Qlt must contain n*(n+1)/2 lower-triangular entries" );
     std::vector<double> cv;
     const double* cp = nullptr;
     if( !c.is_none() ){
       arma::vec cav = _to_arma_vec( c );
       if( cav.n_elem != n ) throw py::value_error( "center c must have length n" );
       cv.resize( n );
       for( unsigned i=0; i<n; ++i ) cv[i] = cav(i);
       cp = cv.data();
     }
     return ELL( n, Qlt.data(), cp );
   } ),
   py::arg("n"),
   py::arg("Qlt"),
   py::arg("c") = py::none(),
   "constructor from lower-triangular shape entries stored columnwise"
 )
 .def_static(
   "from_radius",
   []( py::handle r, py::object c ){
     return ELL( _to_arma_vec( r ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("r"),
   py::arg("c") = py::none(),
   "construct an axis-aligned ellipsoid enclosing the box c +/- r"
 )
 .def_static(
   "from_box",
   []( std::vector<I> const& box ){
     return ELL( static_cast<unsigned>( box.size() ), box.data() );
   },
   py::arg("box"),
   "construct an ellipsoid enclosing an interval vector"
 )
 .def( py::init<ELL const&>(), "copy constructor" )
 .def_readwrite_static(
   "options",
   &ELL::options,
   "global options for ellipsoidal calculus"
 )
 .def(
   "unitball",
   &ELL::unitball,
   py::arg("n"),
   py::return_value_policy::reference_internal,
   "reset to the n-dimensional unit ball"
 )
 .def(
   "set",
   []( ELL& self, py::handle Q, py::object c ) -> ELL& {
     return self.set( _to_arma_mat( Q ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("Q"),
   py::arg("c") = py::none(),
   py::return_value_policy::reference_internal,
   "set dense shape matrix Q and optional center c"
 )
 .def(
   "set_lower",
   []( ELL& self, unsigned const n, std::vector<double> const& Qlt, py::object c ) -> ELL& {
     if( Qlt.size() != static_cast<size_t>( n*(n+1)/2 ) )
       throw py::value_error( "Qlt must contain n*(n+1)/2 lower-triangular entries" );
     std::vector<double> cv;
     const double* cp = nullptr;
     if( !c.is_none() ){
       arma::vec cav = _to_arma_vec( c );
       if( cav.n_elem != n ) throw py::value_error( "center c must have length n" );
       cv.resize( n );
       for( unsigned i=0; i<n; ++i ) cv[i] = cav(i);
       cp = cv.data();
     }
     return self.set( n, Qlt.data(), cp );
   },
   py::arg("n"),
   py::arg("Qlt"),
   py::arg("c") = py::none(),
   py::return_value_policy::reference_internal,
   "set lower-triangular shape entries stored columnwise"
 )
 .def(
   "set_radius",
   []( ELL& self, py::handle r, py::object c ) -> ELL& {
     return self.set( _to_arma_vec( r ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("r"),
   py::arg("c") = py::none(),
   py::return_value_policy::reference_internal,
   "set an axis-aligned ellipsoid enclosing the box c +/- r"
 )
 .def(
   "set_box",
   []( ELL& self, std::vector<I> const& box ) -> ELL& {
     return self.set( static_cast<unsigned>( box.size() ), box.data() );
   },
   py::arg("box"),
   py::return_value_policy::reference_internal,
   "set an ellipsoid enclosing an interval vector"
 )
 .def( "reset",     &ELL::reset,     py::return_value_policy::reference_internal, "reset ellipsoid" )
 .def( "reset_aux", &ELL::reset_aux, py::return_value_policy::reference_internal, "clear cached decompositions" )
 .def(
   "extend",
   []( ELL& self, py::handle Qi, double const& ci ) -> ELL& {
     return self.extend( arma::rowvec( _to_arma_vec( Qi ).t() ), ci );
   },
   py::arg("Qi"),
   py::arg("ci") = 0.,
   py::return_value_policy::reference_internal,
   "append one dimension using lower-triangular row Qi and center entry ci"
 )
 .def( "O", &ELL::O, "return the centred ellipsoid with the same shape matrix" )
 .def_property_readonly( "n", &ELL::n, "dimension of the ellipsoid" )
 .def_property_readonly( "c", []( ELL const& self ){ return _from_arma_vec( self.c() ); }, "center vector" )
 .def_property_readonly( "Q", []( ELL const& self ){ return _from_arma_mat( self.Q() ); }, "shape matrix" )
 .def( "ci", py::overload_cast<unsigned>( &ELL::c, py::const_ ), py::arg("i"), "center coefficient" )
 .def( "Qij", py::overload_cast<unsigned,unsigned>( &ELL::Q, py::const_ ), py::arg("i"), py::arg("j"), "shape matrix coefficient" )
 .def(
   "set_ci",
   []( ELL& self, unsigned const i, double const& v ) -> ELL& {
     if( i >= self.n() ) throw py::index_error( "center index out of range" );
     self.c()(i) = v;
     return self;
   },
   py::arg("i"),
   py::arg("value"),
   py::return_value_policy::reference_internal,
   "set center coefficient"
 )
 .def(
   "set_Qij",
   []( ELL& self, unsigned const i, unsigned const j, double const& v ) -> ELL& {
     if( i >= self.n() || j >= self.n() ) throw py::index_error( "shape index out of range" );
     self.Q(i,j) = v;
     self.Q(j,i) = v;
     self.reset_aux();
     return self;
   },
   py::arg("i"),
   py::arg("j"),
   py::arg("value"),
   py::return_value_policy::reference_internal,
   "set symmetric shape matrix coefficient and clear cached decompositions"
 )
 .def( "trQ",   &ELL::trQ,   "trace of the shape matrix" )
 .def( "psdQ",  &ELL::psdQ,  "test whether the shape matrix is positive semidefinite" )
 .def( "rankQ", &ELL::rankQ, "rank of the shape matrix with current tolerances" )
 .def( "eigQ",  []( ELL& self ){ auto const& e = self.eigQ(); return py::make_tuple( _from_arma_vec(e.first), _from_arma_mat(e.second) ); }, "eigenvalue decomposition of Q" )
 .def( "sqrtQ", []( ELL& self, bool const complete ){ return _from_arma_mat( self.sqrtQ( complete ) ); }, py::arg("complete")=false, "matrix square root of Q" )
 .def( "svdQ",  []( ELL& self ){ auto const& s = self.svdQ(); return py::make_tuple( _from_arma_vec(s.first), _from_arma_mat(s.second.first), _from_arma_mat(s.second.second) ); }, "singular value decomposition of Q" )
 .def( "regQ",  []( ELL& self ){ return _from_arma_mat( self.regQ() ); }, "regularise and return Q" )
 .def( "invQ",  []( ELL& self ){ return _from_arma_mat( self.invQ() ); }, "inverse or pseudo-inverse of Q" )
 .def( "invQij", []( ELL& self, unsigned const i, unsigned const j ){ return self.invQ(i,j); }, py::arg("i"), py::arg("j"), "coefficient of inverse shape matrix" )
 .def( "align", []( ELL const& self, py::handle v, py::handle x ){ return _from_arma_mat( self.align( _to_arma_vec(v), _to_arma_vec(x) ) ); }, py::arg("v"), py::arg("x"), "orthogonal map rotating x toward v" )
 .def( "l", &ELL::l, py::arg("i"), "coordinate lower bound" )
 .def( "u", &ELL::u, py::arg("i"), "coordinate upper bound" )
 .def( "r", &ELL::r, py::arg("i"), "coordinate radius" )
 .def(
   "translate",
   []( ELL& self, py::handle d ) -> ELL& { self += _to_arma_vec( d ); return self; },
   py::arg("d"),
   py::return_value_policy::reference_internal,
   "translate center by d"
 )
 .def(
   "__iadd__",
   []( ELL& self, py::handle d ) -> ELL& { self += _to_arma_vec( d ); return self; },
   py::is_operator(),
   py::return_value_policy::reference_internal
 )
 .def(
   "__isub__",
   []( ELL& self, py::handle d ) -> ELL& { self -= _to_arma_vec( d ); return self; },
   py::is_operator(),
   py::return_value_policy::reference_internal
 )
 .def( "__str__",  static_cast<std::string(*)(ELL const&)>(&_to_string) )
 .def( "__repr__", static_cast<std::string(*)(ELL const&)>(&_to_string) )
;

// Nested Ellipsoid::Options
py::class_<ELL::Options> pyEllipsoidOptions( pyEllipsoid, "Options" );
pyEllipsoidOptions
 .def( py::init<>() )
 .def( py::init<ELL::Options const&>() )
 .def( "reset", []( ELL::Options& self ){ self = ELL::Options(); }, "reset options to defaults" )
 .def_readwrite( "PSDCHK",     &ELL::Options::PSDCHK,     "whether to check positive semidefiniteness of shape matrices [Default: false]" )
 .def_readwrite( "PSDTOL",     &ELL::Options::PSDTOL,     "absolute tolerance for PSD checks [Default: 1e2*machprec()]" )
 .def_readwrite( "RKTOLA",     &ELL::Options::RKTOLA,     "absolute rank/regularisation tolerance [Default: machprec()]" )
 .def_readwrite( "RKTOLR",     &ELL::Options::RKTOLR,     "relative rank/regularisation tolerance [Default: machprec()]" )
 .def_readwrite( "ROOTTOL",    &ELL::Options::ROOTTOL,    "absolute stopping tolerance for root finding [Default: 1e-10]" )
 .def_readwrite( "ROOTSECANT", &ELL::Options::ROOTSECANT, "whether to use secant root finding before golden-section fallback [Default: false]" )
 .def_readwrite( "ROOTMAXIT",  &ELL::Options::ROOTMAXIT,  "maximum root-finding iterations; 0 means no maximum [Default: 0]" )
;

// Nested Ellipsoid::Exceptions
py::class_<ELL::Exceptions> pyEllipsoidExceptions( pyEllipsoid, "Exceptions" );
pyEllipsoidExceptions
 .def( py::init<ELL::Exceptions::TYPE>() )
 .def( "ierr", &ELL::Exceptions::ierr, "error flag" )
 .def( "what", &ELL::Exceptions::what, "error description" )
;

py::enum_<ELL::Exceptions::TYPE>( pyEllipsoidExceptions, "TYPE" )
 .value( "NONPSD", ELL::Exceptions::TYPE::NONPSD, "non-positive-semidefinite shape matrix" )
 .value( "LAPACK", ELL::Exceptions::TYPE::LAPACK, "linear algebra routine failed" )
 .value( "ROOT",   ELL::Exceptions::TYPE::ROOT,   "root-finding routine failed" )
 .export_values()
;

py::class_<EI, ELL> pyEllImg( m, "EllImg" );

pyEllImg
 .def( py::init<>() )
 .def(
   py::init( []( py::handle Q, py::object c ){
     return EI( _to_arma_mat( Q ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   } ),
   py::arg("Q"),
   py::arg("c") = py::none(),
   "constructor from dense shape matrix Q and optional center c"
 )
 .def(
   py::init( []( unsigned const n, std::vector<double> const& Qlt, py::object c ){
     if( Qlt.size() != static_cast<size_t>( n*(n+1)/2 ) )
       throw py::value_error( "Qlt must contain n*(n+1)/2 lower-triangular entries" );
     std::vector<double> cv;
     const double* cp = nullptr;
     if( !c.is_none() ){
       arma::vec cav = _to_arma_vec( c );
       if( cav.n_elem != n ) throw py::value_error( "center c must have length n" );
       cv.resize( n );
       for( unsigned i=0; i<n; ++i ) cv[i] = cav(i);
       cp = cv.data();
     }
     return EI( n, Qlt.data(), cp );
   } ),
   py::arg("n"),
   py::arg("Qlt"),
   py::arg("c") = py::none(),
   "constructor from lower-triangular shape entries stored columnwise"
 )
 .def_static(
   "from_radius",
   []( py::handle r, py::object c ){
     return EI( _to_arma_vec( r ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("r"),
   py::arg("c") = py::none(),
   "construct an axis-aligned ellipsoidal image enclosing the box c +/- r"
 )
 .def( py::init<EI const&>(), "copy constructor" )
 .def_readwrite_static(
   "options",
   &EI::options,
   "global options for ellipsoidal image propagation"
 )
 .def_property_readonly( "Q_lift", []( EI& self ){ return _from_arma_spmat( self.Q_lift() ); }, "dense copy of the lifted shape matrix" )
 .def_property_readonly( "c_lift", []( EI& self ){ return _from_arma_vec( self.c_lift() ); }, "lifted center vector" )
 .def_property_readonly( "qdim", []( EI& self ){ return static_cast<long>( self.c_lift().n_elem ); }, "current lifted dimension" )
 .def(
   "set",
   []( EI& self, py::handle Q, py::object c ) -> EI& {
     return self.set( _to_arma_mat( Q ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("Q"),
   py::arg("c") = py::none(),
   py::return_value_policy::reference_internal,
   "set dense shape matrix Q and optional center c"
 )
 .def(
   "set_radius",
   []( EI& self, py::handle r, py::object c ) -> EI& {
     return self.set( _to_arma_vec( r ), c.is_none()? arma::vec(): _to_arma_vec( c ) );
   },
   py::arg("r"),
   py::arg("c") = py::none(),
   py::return_value_policy::reference_internal,
   "set an axis-aligned ellipsoidal image enclosing the box c +/- r"
 )
 .def( "reset", []( EI& self ) -> EI& { return self.reset(); }, py::return_value_policy::reference_internal, "reset lifted image to the underlying ellipsoid" )
 .def(
   "get",
   []( EI& self, py::iterable vars ){
     std::vector<EV> v = _to_ev_vector( vars );
     return self.get( static_cast<unsigned>( v.size() ), v.data() );
   },
   py::arg("vars"),
   "project the lifted ellipsoid onto the supplied EllVar variables"
 )
 .def( "__str__",  static_cast<std::string(*)(EI const&)>(&_to_string) )
 .def( "__repr__", static_cast<std::string(*)(EI const&)>(&_to_string) )
;

// Nested EllImg::Options
py::class_<EI::Options> pyEllImgOptions( pyEllImg, "Options" );
pyEllImgOptions
 .def( py::init<>() )
 .def( py::init<EI::Options const&>() )
 .def( "reset", []( EI::Options& self ){ self = EI::Options(); }, "reset options to defaults" )
 .def_readwrite( "PREALLOC",   &EI::Options::PREALLOC,   "number of rows to preallocate in the lifted shape matrix and center [Default: 0]" )
 .def_readwrite( "MINK_TOL",   &EI::Options::MINK_TOL,   "tolerance in Minkowski sums [Default: 1e-10]" )
 .def_readwrite( "REMEZ_USE",  &EI::Options::REMEZ_USE,  "whether to use degree-1 Remez linearisation for univariate terms [Default: true]" )
 .def_readwrite( "REMEZ_MAXIT",&EI::Options::REMEZ_MAXIT,"maximum Remez iterations [Default: 5]" )
 .def_readwrite( "REMEZ_TOL",  &EI::Options::REMEZ_TOL,  "Remez stopping tolerance [Default: 1e-5]" )
 .def_readwrite( "REMEZ_MIG",  &EI::Options::REMEZ_MIG,  "minimum interval diameter for invoking Remez [Default: 1e-10]" )
 .def_readwrite( "DCPROD_USE", &EI::Options::DCPROD_USE, "whether to use DC decomposition to lift bilinear terms [Default: false]" )
;

// Nested EllImg::Exceptions
py::class_<EI::Exceptions> pyEllImgExceptions( pyEllImg, "Exceptions" );
pyEllImgExceptions
 .def( py::init<EI::Exceptions::TYPE>() )
 .def( "ierr", &EI::Exceptions::ierr, "error flag" )
 .def( "what", &EI::Exceptions::what, "error description" )
;

py::enum_<EI::Exceptions::TYPE>( pyEllImgExceptions, "TYPE" )
 .value( "DIV",   EI::Exceptions::TYPE::DIV,   "division by zero scalar" )
 .value( "INV",   EI::Exceptions::TYPE::INV,   "inverse operation with zero in domain" )
 .value( "LOG",   EI::Exceptions::TYPE::LOG,   "log operation with non-positive numbers in domain" )
 .value( "SQRT",  EI::Exceptions::TYPE::SQRT,  "square-root operation with negative numbers in domain" )
 .value( "TAN",   EI::Exceptions::TYPE::TAN,   "tangent operation with zero in cosine domain" )
 .value( "ACOS",  EI::Exceptions::TYPE::ACOS,  "inverse cosine operation with domain outside [-1,1]" )
 .value( "ASIN",  EI::Exceptions::TYPE::ASIN,  "inverse sine operation with domain outside [-1,1]" )
 .value( "INIT",  EI::Exceptions::TYPE::INIT,  "failed to construct ellipsoidal variable" )
 .value( "EIMG",  EI::Exceptions::TYPE::EIMG,  "variables belong to different ellipsoidal images" )
 .value( "UNDEF", EI::Exceptions::TYPE::UNDEF, "feature not yet implemented" )
 .export_values()
;

py::class_<EV> pyEllVar( m, "EllVar" );

pyEllVar
 .def( py::init<>() )
 .def( py::init<EV const&>(), "copy constructor" )
 .def( py::init<double const&>(), py::arg("cst"), "constructor for constant scalar" )
 .def( py::init<I const&>(), py::arg("range"), "constructor for interval/range" )
 .def( py::init<double const&, double const&>(), py::arg("l"), py::arg("u"), "constructor from lower and upper bounds" )
 .def( py::init<EI&, unsigned const>(), py::arg("image"), py::arg("index"), py::keep_alive<1,2>(), "constructor for image variable" )
 .def( py::init<EI&, unsigned const, I const&>(), py::arg("image"), py::arg("index"), py::arg("range"), py::keep_alive<1,2>(), "constructor for image variable with tailored range" )
 .def(
   "set",
   []( EV& self, EI& img, unsigned const i ) -> EV& { return self.set( img, i ); },
   py::arg("image"),
   py::arg("index"),
   py::return_value_policy::reference_internal,
   py::keep_alive<1,2>(),
   "set variable in ellipsoidal image environment"
 )
 .def(
   "set",
   []( EV& self, EI& img, unsigned const i, I const& range ) -> EV& { return self.set( img, i, range ); },
   py::arg("image"),
   py::arg("index"),
   py::arg("range"),
   py::return_value_policy::reference_internal,
   py::keep_alive<1,2>(),
   "set variable in ellipsoidal image environment with tailored range"
 )
 .def_property_readonly( "range", &EV::range, "interval/range enclosure of the variable" )
 .def_property_readonly( "image", &EV::image, py::return_value_policy::reference_internal, "associated ellipsoidal image, or None for constants/ranges" )
 .def_property_readonly( "index", &EV::index, "row index in the lifted image, or -1 for constants/ranges" )
 .def( "__str__",  static_cast<std::string(*)(EV const&)>(&_to_string) )
 .def( "__repr__", static_cast<std::string(*)(EV const&)>(&_to_string) )
 // compound assignment: in place, returns self -> no new image reference
 .def( py::self += double() )
 .def( py::self += I() )
 .def( py::self += py::self )
 .def( py::self -= double() )
 .def( py::self -= I() )
 .def( py::self -= py::self )
 .def( py::self *= double() )
 .def( py::self *= I() )
 .def( py::self *= py::self )
 .def( py::self /= double() )
 .def( py::self /= I() )
 .def( py::self /= py::self )
 // unary and binary operators returning a NEW EllVar: tie the result's
 // lifetime to its EllVar operand(s) so the underlying EllImg stays alive
 .def( "__pos__", []( EV const& a ){ return +a; }, py::keep_alive<0,1>() )
 .def( "__neg__", []( EV const& a ){ return -a; }, py::keep_alive<0,1>() )
 .def( "__add__",      []( EV const& a, EV const& b ){ return a + b; }, py::keep_alive<0,1>(), py::keep_alive<0,2>() )
 .def( "__add__",      []( EV const& a, double b ){ return a + b; }, py::keep_alive<0,1>() )
 .def( "__add__",      []( EV const& a, I const& b ){ return a + b; }, py::keep_alive<0,1>() )
 .def( "__radd__",     []( EV const& a, double b ){ return b + a; }, py::keep_alive<0,1>() )
 .def( "__radd__",     []( EV const& a, I const& b ){ return b + a; }, py::keep_alive<0,1>() )
 .def( "__sub__",      []( EV const& a, EV const& b ){ return a - b; }, py::keep_alive<0,1>(), py::keep_alive<0,2>() )
 .def( "__sub__",      []( EV const& a, double b ){ return a - b; }, py::keep_alive<0,1>() )
 .def( "__sub__",      []( EV const& a, I const& b ){ return a - b; }, py::keep_alive<0,1>() )
 .def( "__rsub__",     []( EV const& a, double b ){ return b - a; }, py::keep_alive<0,1>() )
 .def( "__rsub__",     []( EV const& a, I const& b ){ return b - a; }, py::keep_alive<0,1>() )
 .def( "__mul__",      []( EV const& a, EV const& b ){ return a * b; }, py::keep_alive<0,1>(), py::keep_alive<0,2>() )
 .def( "__mul__",      []( EV const& a, double b ){ return a * b; }, py::keep_alive<0,1>() )
 .def( "__mul__",      []( EV const& a, I const& b ){ return a * b; }, py::keep_alive<0,1>() )
 .def( "__rmul__",     []( EV const& a, double b ){ return b * a; }, py::keep_alive<0,1>() )
 .def( "__rmul__",     []( EV const& a, I const& b ){ return b * a; }, py::keep_alive<0,1>() )
 .def( "__truediv__",  []( EV const& a, EV const& b ){ return a / b; }, py::keep_alive<0,1>(), py::keep_alive<0,2>() )
 .def( "__truediv__",  []( EV const& a, double b ){ return a / b; }, py::keep_alive<0,1>() )
 .def( "__truediv__",  []( EV const& a, I const& b ){ return a / b; }, py::keep_alive<0,1>() )
 .def( "__rtruediv__", []( EV const& a, double b ){ return b / a; }, py::keep_alive<0,1>() )
 .def( "__rtruediv__", []( EV const& a, I const& b ){ return b / a; }, py::keep_alive<0,1>() )
 .def( "__pow__", []( EV const& x, int const n ){ return mc::pow( x, n ); }, py::keep_alive<0,1>() )
;

m.def( "inv",  []( EV const& x ){ return mc::inv( x ); }, py::keep_alive<0,1>() );
m.def( "sqr",  []( EV const& x ){ return mc::sqr( x ); }, py::keep_alive<0,1>() );
m.def( "sqrt", []( EV const& x ){ return mc::sqrt( x ); }, py::keep_alive<0,1>() );
m.def( "exp",  []( EV const& x ){ return mc::exp( x ); }, py::keep_alive<0,1>() );
m.def( "log",  []( EV const& x ){ return mc::log( x ); }, py::keep_alive<0,1>() );
m.def( "xlog", []( EV const& x ){ return mc::xlog( x ); }, py::keep_alive<0,1>() );
m.def( "cos",  []( EV const& x ){ return mc::cos( x ); }, py::keep_alive<0,1>() );
m.def( "sin",  []( EV const& x ){ return mc::sin( x ); }, py::keep_alive<0,1>() );
m.def( "tan",  []( EV const& x ){ return mc::tan( x ); }, py::keep_alive<0,1>() );
m.def( "acos", []( EV const& x ){ return mc::acos( x ); }, py::keep_alive<0,1>() );
m.def( "asin", []( EV const& x ){ return mc::asin( x ); }, py::keep_alive<0,1>() );
m.def( "atan", []( EV const& x ){ return mc::atan( x ); }, py::keep_alive<0,1>() );
m.def( "cosh", []( EV const& x ){ return mc::cosh( x ); }, py::keep_alive<0,1>() );
m.def( "sinh", []( EV const& x ){ return mc::sinh( x ); }, py::keep_alive<0,1>() );
m.def( "tanh", []( EV const& x ){ return mc::tanh( x ); }, py::keep_alive<0,1>() );
m.def( "erf",  []( EV const& x ){ return mc::erf( x ); }, py::keep_alive<0,1>() );
m.def( "erfc", []( EV const& x ){ return mc::erfc( x ); }, py::keep_alive<0,1>() );
m.def( "pow",  []( EV const& x, int const n ){ return mc::pow( x, n ); }, py::keep_alive<0,1>() );
m.def( "cheb", []( EV const& x, unsigned const n ){ return mc::cheb( x, n ); }, py::keep_alive<0,1>() );

m.def( "ell_unitball", []( unsigned const n ){ return mc::ell_unitball( n ); }, py::arg("n") );
m.def( "mtimes", []( ELL const& E, py::handle A, py::object b ){ return mc::mtimes( E, _to_arma_mat( A ), b.is_none()? arma::vec(): _to_arma_vec( b ) ); }, py::arg("E"), py::arg("A"), py::arg("b")=py::none() );
m.def( "minksum_ea", []( ELL const& E1, ELL const& E2, double const eps ){ return mc::minksum_ea( E1, E2, eps ); }, py::arg("E1"), py::arg("E2"), py::arg("eps")=mc::machprec() );
m.def(
  "minksum_ea",
  []( std::vector<ELL> const& E, py::object D ) -> py::object {
    auto a = py::array_t<double, py::array::c_style | py::array::forcecast>::ensure( D );
    if( a && a.request().ndim == 1 )
      return py::cast( mc::minksum_ea( E, _to_arma_vec( D ) ) );
    return py::cast( mc::minksum_ea( E, _to_vec_vector( D ) ) );
  },
  py::arg("ellipsoids"),
  py::arg("direction_or_directions"),
  "external approximation of a Minkowski sum of ellipsoids for one direction, or for an iterable of directions"
);
m.def( "minksum_box", []( ELL const& E, py::handle r, py::object c, double const tol, double const eps ){ return mc::minksum_ea( E, std::make_pair( _to_arma_vec( r ), c.is_none()? arma::vec(): _to_arma_vec( c ) ), tol, eps ); }, py::arg("E"), py::arg("r"), py::arg("c")=py::none(), py::arg("tol")=1e-10, py::arg("eps")=mc::machprec() );
m.def( "minksum_interval", []( ELL const& E, std::vector<I> const& box, double const tol, double const eps ){ return mc::minksum_ea( E, box.data(), tol, eps ); }, py::arg("E"), py::arg("box"), py::arg("tol")=1e-10, py::arg("eps")=mc::machprec() );
m.def( "inv", []( ELL const& E ){ return mc::inv( E ); } );
m.def( "inv", []( std::vector<ELL> const& E ){ return mc::inv( E ); } );
m.def( "dist", []( ELL const& E, py::tuple hp ){ return mc::dist( E, _to_hp( hp ) ); }, py::arg("E"), py::arg("hp") );
m.def( "dist", []( ELL const& E1, ELL const& E2 ){ return mc::dist( E1, E2 ); }, py::arg("E1"), py::arg("E2") );
m.def( "hpintersection", []( ELL const& E, py::tuple hp ){ return mc::hpintersection( E, _to_hp( hp ) ); }, py::arg("E"), py::arg("hp") );
m.def( "hpintersection", []( ELL const& E, py::iterable hp ){ return mc::hpintersection( E, _to_hps( hp ) ); }, py::arg("E"), py::arg("hp") );
m.def( "intersection_ea", []( ELL const& E1, ELL const& E2, double const tol ){ return mc::intersection_ea( E1, E2, tol ); }, py::arg("E1"), py::arg("E2"), py::arg("tol")=mc::machprec() );
m.def( "intersection_ea", []( ELL const& E, py::tuple hp, double const tol ){ return mc::intersection_ea( E, _to_hp( hp ), tol ); }, py::arg("E"), py::arg("hp"), py::arg("tol")=mc::machprec() );
m.def( "ellintersection_ia", []( py::iterable hp, double const tol, unsigned const maxit ){ return mc::ellintersection_ia( _to_hps( hp ), tol, maxit ); }, py::arg("hp"), py::arg("tol")=1e-4, py::arg("maxit")=100 );

}
