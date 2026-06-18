#include <iostream>

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

#include "ellimage.hpp"
typedef mc::EllImg<I> EI;
typedef mc::EllVar<I> EV;

inline void
require(bool ok, const char* msg)
{
  if (!ok)
  {
    std::cerr << "FAIL: " << msg << "\n";
    std::exit(1);
  }
}

static bool
psd(const arma::mat& Q, double tol = 1e-8)
{
  arma::vec eig;
  arma::eig_sym(eig, Q);
  return eig.min() >= -tol;
}

int
main()
{
  // Reproduce the documented example with the current default EllImg options.
  EI::Options defaults;
  EI::options = defaults;

  arma::vec cx(2);
  arma::mat Qx(2, 2, arma::fill::zeros);
  cx(0)    = 3.;
  Qx(0, 0) = 5.;
  cx(1)    = 4.;
  Qx(1, 0) = 4.;
  Qx(0, 1) = 4.;
  Qx(1, 1) = 5.;

  EI Ex(Qx, cx);
  EV X1(Ex, 0);
  EV X2(Ex, 1);

  EV F[2] = {log(X1) + mc::sqr(X2), sin(X1) - cos(X2)};

  arma::vec qlift = Ex.c_lift();
  arma::mat Qlift(Ex.Q_lift());
  require(qlift.n_elem == 8, "documented example lifted dimension");
  require(Qlift.n_rows == 8 && Qlift.n_cols == 8,
          "documented example lifted shape size");
  require(psd(Qlift), "documented example lifted shape is PSD");

  mc::Ellipsoid Ef = Ex.get(2, F);
  std::cout << Ef << std::endl;

  require(Ef.n() == 2, "documented example projection dimension");
  require(psd(Ef.Q()), "documented example projected shape is PSD");

  std::cout << "ellimage Doxygen example test passed.\n";
  return 0;
}
