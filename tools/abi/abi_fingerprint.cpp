// abi_fingerprint.cpp -- the binary layout of the MC++ classes that modules built on pymcpp share with it.
//
// A module such as CRONOS's `cronos` is compiled against MC++'s headers and then receives live FFVar, FFGraph and
// operation objects created by the separately compiled `pymcpp`.  That is sound only if both agree on these classes'
// layout and on pybind11's internals.  MC++'s rule (INSTALL.md, "Releasing"): a PATCH release (x.y.z -> x.y.z+1)
// leaves this fingerprint unchanged; a change needs a new MINOR release.  tools/abi/check_abi.py enforces it in CI.
//
// Prints the layout-affecting configuration, one "name sizeof alignof polymorphic" line per class, then the pybind11
// version and internals version.  A layout is a function of the MC++ version AND of the configuration: build this
// program exactly as the wheels are built (CMake target `abi_fingerprint`, MC++'s own definitions).  The macros that
// change shared layouts: MC__USE_THREADLOCAL (FFBase/FFGraph: _curOp a thread_local static instead of a member),
// MC__USE_THREAD (Vect), MC__USE_ARMADILLO (CModel, TModel, specbnd) and the interval library (the template
// operations' argument).
// What it cannot see: a virtual function reordered without a size change, or a member's type changed to one of the
// same size -- the rule still forbids those; the fingerprint catches the common breaks.
#include <cstdio>
#include <type_traits>
#include <pybind11/pybind11.h>
// The interval type pymcpp's wheels instantiate the template operations with (MC_INTERVAL_LIBRARY=BOOST, the default;
// src/pymcpp/*.cpp).  Kept in step with them by hand: change both together.
#ifndef MC__USE_BOOST
# error "abi_fingerprint: built for the Boost interval library, the wheels' (MC_INTERVAL_LIBRARY=BOOST)"
#endif
#include "mcboost.hpp"
typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
typedef boost::numeric::interval_lib::policies<T_boost_round, T_boost_check> T_boost_policy;
typedef boost::numeric::interval<double, T_boost_policy> I;
#include "ffunc.hpp"
#include "ffexpr.hpp"
#include "ffcustom.hpp"
#include "ffdagext.hpp"
#include "fflin.hpp"
#include "ffmlp.hpp"
#include "ffvect.hpp"

#define FP( T ) std::printf( "%-22s %4zu %2zu %d\n", #T, sizeof( T ), alignof( T ), (int)std::is_polymorphic<T>::value )

int main(){
  std::printf( "config THREADLOCAL=%d THREAD=%d ARMADILLO=%d interval=%s\n",
#ifdef MC__USE_THREADLOCAL
    1,
#else
    0,
#endif
#ifdef MC__USE_THREAD
    1,
#else
    0,
#endif
#ifdef MC__USE_ARMADILLO
    1,
#else
    0,
#endif
    "BOOST" );
  FP( mc::FFBase );   FP( mc::FFGraph );  FP( mc::FFGraph::Options );  FP( mc::FFVar );  FP( mc::FFNum );
  FP( mc::FFOp );     FP( mc::FFSubgraph );
  FP( mc::FFPartial ); FP( mc::FFEval );  FP( mc::FFIntegral );
  FP( mc::FFCustom<I> ); FP( mc::FFDAGEXT<I> );  FP( mc::FFLin<I> );  FP( mc::FFMLP<I> );  FP( mc::FFVect<I> );
  std::printf( "pybind11 %d.%d.%d internals %d\n", PYBIND11_VERSION_MAJOR, PYBIND11_VERSION_MINOR,
               PYBIND11_VERSION_PATCH, PYBIND11_INTERNALS_VERSION );
  return 0;
}
