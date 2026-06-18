#define PEAK_TANH_40L4  // <-- select TorchScript here
#define SAVE_RESULTS    // <-- specify whether to save results to file
#undef ANALYSE_RATE     // <-- specify whether to analyse rate of convergence
#undef ANALYSE_TIME     // <-- specify whether to analyse computational time

#undef MC__FFMLP_CHECK
#undef MC__FFMLP_DEBUG

////////////////////////////////////////////////////////////////////////

#include <chrono>
#include <fstream>
#include <iomanip>

#if defined(MC__USE_PROFIL)
#include "mcprofil.hpp"
typedef INTERVAL I;
#elif defined(MC__USE_FILIB)
#include "mcfilib.hpp"
typedef filib::interval<double, filib::native_switched, filib::i_mode_extended>
    I;
#elif defined(MC__USE_BOOST)
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

#include "ffmlp.hpp"

typedef mc::McCormick<I> MC;
typedef mc::SupModel<mc::PWCU> PWCSM;
typedef mc::SupVar<mc::PWCU> PWCSV;
typedef mc::SupModel<mc::PWLU> PWLSM;
typedef mc::SupVar<mc::PWLU> PWLSV;

///////////////////////////////////////////////////////////////////////////////

int
main()
{
  try
  {
    // Import ANN from TorchScript file
    mc::MLP<I> ANN;
    ANN.options.EVALTORCH = false;

#if defined(PEAK_TANH_20L2)
    size_t const NX  = 2;
    double const X1L = -3., X1U = 3.;
    double const X2L = -3., X2U = 3.;
    if (!ANN.read_data("peak_Tanh_20L2.pt"))  //, true ) )
      return -1;

#elif defined(PEAK_RELU_20L2)
    size_t const NX  = 2;
    double const X1L = -3., X1U = 3.;
    double const X2L = -3., X2U = 3.;
    if (!ANN.read_data("peak_ReLU_20L2.pt"))  //, true ) )
      return -1;

#elif defined(PEAK_TANH_40L4)
    size_t const NX  = 2;
    double const X1L = -3., X1U = 3.;
    double const X2L = -3., X2U = 3.;
    std::map<unsigned, double> XREF{{0, 2.276872e-01},
                                    {1, -1.626030e+00}};  // Global min
    if (!ANN.read_data("peak_Tanh_40L4.pt"))              //, true ) )
      return -1;

#elif defined(PEAK_RELU_40L4)
    size_t const NX  = 2;
    double const X1L = -3., X1U = 3.;
    double const X2L = -3., X2U = 3.;
    std::map<unsigned, double> XREF{{0, 2.685182e-01},
                                    {1, -1.648824e+00}};  // Global min
    if (!ANN.read_data("peak_ReLU_40L4.pt"))              //, true ) )
      return -1;
#endif

    // Create DAG
    mc::FFGraph DAG;
    std::vector<mc::FFVar> X = DAG.add_vars(NX, "X");

    mc::FFMLP<I> OpNN;
    std::vector<mc::FFVar> Y{OpNN(0, X, &ANN, OpNN.COPY)};

    // Create and display subgraph
    auto SgY = DAG.subgraph(Y);
    // DAG.output( SgY, " OF Y" );
    auto StrY = mc::FFExpr::subgraph(&DAG, SgY);
    std::cout << "Y: " << StrY[0] << std::endl;

    // Evaluate DAG
    std::vector<double> DX{0, 0}, DY(1);
    ANN.eval(DX.data(), DY.data());
    std::cout << "Y(" << DX[0] << "," << DX[1] << ") = " << DY[0] << std::endl;
    DAG.eval(SgY, Y, DY, X, DX);
    std::cout << "Y(" << DX[0] << "," << DX[1] << ") = " << DY[0] << std::endl;

    // McCormick relaxation
    std::vector<I> IX{I(X1L, X1U), I(X2L, X2U)};
    std::vector<MC> MCX{MC(IX[0], DX[0]), MC(IX[1], DX[1])}, MCY;
    DAG.eval(SgY, Y, MCY, X, MCX);
    std::cout << "Y = " << MCY[0] << std::endl;

    // Piecewise-constant superposition relaxation
    PWCSM pwcmod(NX);
    mc::PWCU::options.SLOPEUSE = 0;
    pwcmod.options.PROD_METH   = PWLSM::Options::PARTIAL;  // FULL;//LOG;//NONE;
    pwcmod.options.PROD_CUT    = 0;
    pwcmod.options.USE_SHADOW  = 0;
    pwcmod.options.OPT_SHADOW  = 1;
    pwcmod.options.USE_CVXCCV  = 1;
    pwcmod.options.USE_ENDRAY  = 1;
    pwcmod.options.MAX_SUBDIV  = 32;  // 16;

    std::vector<PWCSV> PWCSVX{
        PWCSV(pwcmod, 0, I(X1L, X1U), pwcmod.options.MAX_SUBDIV),
        PWCSV(pwcmod, 1, I(X2L, X2U), pwcmod.options.MAX_SUBDIV)},
        PWCSVY;
    DAG.eval(SgY, Y, PWCSVY, X, PWCSVX);
    // std::cout << "Y = " << PWCSVY[0] << std::endl;

    // Piecewise-linear superposition relaxation
    PWLSM pwlmod(NX);
    pwlmod.options = pwcmod.options;

    // std::vector<PWLSV> PWLSVX{ PWLSV( pwlmod, 0, I(X1L,X1U),
    // pwlmod.options.MAX_SUBDIV ),
    //                            PWLSV( pwlmod, 1, I(X2L,X2U),
    //                            pwlmod.options.MAX_SUBDIV ) }, PWLSVY;
    std::vector<PWLSV> PWLSVX{PWLSV(pwlmod, 0, I(X1L, X1U), 1),
                              PWLSV(pwlmod, 1, I(X2L, X2U), 1)},
        PWLSVY;
    DAG.eval(SgY, Y, PWLSVY, X, PWLSVX);
    // std::cout << "Y = " << PWLSVY[0] << std::endl;

#ifdef SAVE_RESULTS
    // Repeated calculations at grid points
    std::ofstream resfile("test_Torch.out", std::ios_base::out);
    resfile << std::scientific << std::setprecision(5) << std::right;

    int const NPTS = 40;
    for (int iX1 = 0; iX1 < NPTS; iX1++)
    {
      for (int iX2 = 0; iX2 < NPTS; iX2++)
      {
        std::vector<double> DX{X1L + iX1 * (X1U - X1L) / (NPTS - 1.),
                               X2L + iX2 * (X2U - X2L) / (NPTS - 1.)};
        std::vector<double> DY;
        DAG.eval(SgY, Y, DY, X, DX);

        std::vector<MC> MCX{MC(I(X1L, X1U), DX[0]), MC(I(X2L, X2U), DX[1])},
            MCY;
        DAG.eval(SgY, Y, MCY, X, MCX);

        // Calculate relaxations + propagate all subgradient components
        resfile << std::setw(14) << DX[0] << std::setw(14) << DX[1]
                << std::setw(14) << DY[0] << std::setw(14) << MCY[0].l()
                << std::setw(14) << MCY[0].u() << std::setw(14) << MCY[0].cv()
                << std::setw(14) << MCY[0].cc() << std::setw(14)
                << PWLSVY[0].l() << std::setw(14) << PWLSVY[0].u()
                << std::setw(14) << PWLSVY[0].uval({{0, DX[0]}, {1, DX[1]}})
                << std::setw(14) << PWLSVY[0].oval({{0, DX[0]}, {1, DX[1]}})
                << std::endl;
      }
      resfile << std::endl;
    }
#endif

    std::vector<I> const XBND0{I(X1L, X1U), I(X2L, X2U)};
    // std::map<unsigned,double> XREF{ {0,-1.497271e-02}, {1,1.617562e+00} }; //
    // Global max std::map<unsigned,double> XREF{ {0,2.034074e-01},
    // {1,-1.637012e+00} }; // Global min

#ifdef ANALYSE_RATE
    auto const& red = [=](const I& bnd, const double& ref, const double& r)
    { return r * bnd + (1 - r) * ref; };
    auto const& min = [=](const double& x, const double& y)
    { return x < y ? x : y; };
    auto const& max = [=](const double& x, const double& y)
    { return x > y ? x : y; };

    size_t const NGRID  = 100;
    double const rate   = 0.95;
    double const rhomin = 1e-5;
    for (double rho = 1.; rho > rhomin; rho *= rate)
    {
      double min_Y = DY[0], max_Y = DY[0];
      double min_MCYcv = DY[0], max_MCYcc = DY[0];
      double distmax_MCYcv = 0., distmax_MCYcc = 0.;
      double distmax_PWC16SVYu = 0., distmax_PWC16SVYo = 0.;
      double distmax_PWC32SVYu = 0., distmax_PWC32SVYo = 0.;
      double distmax_PWC64SVYu = 0., distmax_PWC64SVYo = 0.;
      double distmax_PWC128SVYu = 0., distmax_PWC128SVYo = 0.;
      double distmax_PWL1SVYu = 0., distmax_PWL1SVYo = 0.;
      double distmax_PWL2SVYu = 0., distmax_PWL2SVYo = 0.;
      double distmax_PWL4SVYu = 0., distmax_PWL4SVYo = 0.;
      double distmax_PWL8SVYu = 0., distmax_PWL8SVYo = 0.;

      std::vector<I> XBND{red(XBND0[0], XREF[0], rho),
                          red(XBND0[1], XREF[1], rho)};

      std::vector<PWCSV> PWC16SVX(NX), PWC16SVY;
      PWC16SVX[0].set(pwcmod, 0, XBND[0], 16);
      PWC16SVX[1].set(pwcmod, 1, XBND[1], 16);
      DAG.eval(SgY, Y, PWC16SVY, X, PWC16SVX);

      std::vector<PWCSV> PWC32SVX(NX), PWC32SVY;
      PWC32SVX[0].set(pwcmod, 0, XBND[0], 32);
      PWC32SVX[1].set(pwcmod, 1, XBND[1], 32);
      DAG.eval(SgY, Y, PWC32SVY, X, PWC32SVX);

      std::vector<PWCSV> PWC64SVX(NX), PWC64SVY;
      PWC64SVX[0].set(pwcmod, 0, XBND[0], 64);
      PWC64SVX[1].set(pwcmod, 1, XBND[1], 64);
      DAG.eval(SgY, Y, PWC64SVY, X, PWC64SVX);

      std::vector<PWCSV> PWC128SVX(NX), PWC128SVY;
      PWC128SVX[0].set(pwcmod, 0, XBND[0], 128);
      PWC128SVX[1].set(pwcmod, 1, XBND[1], 128);
      DAG.eval(SgY, Y, PWC128SVY, X, PWC128SVX);

      std::vector<PWLSV> PWL1SVX(NX), PWL1SVY;
      PWL1SVX[0].set(pwlmod, 0, XBND[0], 1);
      PWL1SVX[1].set(pwlmod, 1, XBND[1], 1);
      DAG.eval(SgY, Y, PWL1SVY, X, PWL1SVX);

      //      std::vector<PWLSV> PWL2SVX(NX), PWL2SVY;
      //      PWL2SVX[0].set( pwlmod, 0, XBND[0], 2 );
      //      PWL2SVX[1].set( pwlmod, 1, XBND[1], 2 );
      //      DAG.eval( SgY, Y, PWL2SVY, X, PWL2SVX );

      //      std::vector<PWLSV> PWL4SVX(NX), PWL4SVY;
      //      PWL4SVX[0].set( pwlmod, 0, XBND[0], 4 );
      //      PWL4SVX[1].set( pwlmod, 1, XBND[1], 4 );
      //      DAG.eval( SgY, Y, PWL4SVY, X, PWL4SVX );

      //      std::vector<PWLSV> PWL8SVX(NX), PWL8SVY;
      //      PWL8SVX[0].set( pwlmod, 0, XBND[0], 8 );
      //      PWL8SVX[1].set( pwlmod, 1, XBND[1], 8 );
      //      DAG.eval( SgY, Y, PWL8SVY, X, PWL8SVX );

      for (unsigned iX1 = 0; iX1 < NGRID; iX1++)
      {
        for (unsigned iX2 = 0; iX2 < NGRID; iX2++)
        {
          std::vector<double> DX{
              mc::Op<I>::l(XBND[0]) +
                  iX1 * mc::Op<I>::diam(XBND[0]) / (NGRID - 1.),
              mc::Op<I>::l(XBND[1]) +
                  iX2 * mc::Op<I>::diam(XBND[1]) / (NGRID - 1.)};
          std::vector<double> DY;
          DAG.eval(SgY, Y, DY, X, DX);

          min_Y = min(min_Y, DY[0]);
          max_Y = max(max_Y, DY[0]);

          std::vector<MC> MCX{MC(XBND[0], DX[0]), MC(XBND[1], DX[1])}, MCY;
          DAG.eval(SgY, Y, MCY, X, MCX);

          min_MCYcv = min(min_MCYcv, MCY[0].cv());
          max_MCYcc = max(max_MCYcc, MCY[0].cc());

          distmax_MCYcv = max(distmax_MCYcv, DY[0] - MCY[0].cv());
          distmax_MCYcc = max(distmax_MCYcc, MCY[0].cc() - DY[0]);

          double const PWC16SVYu = PWC16SVY[0].uval({{0, DX[0]}, {1, DX[1]}});
          double const PWC16SVYo = PWC16SVY[0].oval({{0, DX[0]}, {1, DX[1]}});

          distmax_PWC16SVYu = max(distmax_PWC16SVYu, DY[0] - PWC16SVYu);
          distmax_PWC16SVYo = max(distmax_PWC16SVYo, PWC16SVYo - DY[0]);

          double const PWC32SVYu = PWC32SVY[0].uval({{0, DX[0]}, {1, DX[1]}});
          double const PWC32SVYo = PWC32SVY[0].oval({{0, DX[0]}, {1, DX[1]}});

          distmax_PWC32SVYu = max(distmax_PWC32SVYu, DY[0] - PWC32SVYu);
          distmax_PWC32SVYo = max(distmax_PWC32SVYo, PWC32SVYo - DY[0]);

          double const PWC64SVYu = PWC64SVY[0].uval({{0, DX[0]}, {1, DX[1]}});
          double const PWC64SVYo = PWC64SVY[0].oval({{0, DX[0]}, {1, DX[1]}});

          distmax_PWC64SVYu = max(distmax_PWC64SVYu, DY[0] - PWC64SVYu);
          distmax_PWC64SVYo = max(distmax_PWC64SVYo, PWC64SVYo - DY[0]);

          double const PWC128SVYu = PWC128SVY[0].uval({{0, DX[0]}, {1, DX[1]}});
          double const PWC128SVYo = PWC128SVY[0].oval({{0, DX[0]}, {1, DX[1]}});

          distmax_PWC128SVYu = max(distmax_PWC128SVYu, DY[0] - PWC128SVYu);
          distmax_PWC128SVYo = max(distmax_PWC128SVYo, PWC128SVYo - DY[0]);

          double const PWL1SVYu = PWL1SVY[0].uval({{0, DX[0]}, {1, DX[1]}});
          double const PWL1SVYo = PWL1SVY[0].oval({{0, DX[0]}, {1, DX[1]}});

          distmax_PWL1SVYu = max(distmax_PWL1SVYu, DY[0] - PWL1SVYu);
          distmax_PWL1SVYo = max(distmax_PWL1SVYo, PWL1SVYo - DY[0]);

          //          double const PWL2SVYu =
          //          PWL2SVY[0].uval({{0,DX[0]},{1,DX[1]}}); double const
          //          PWL2SVYo = PWL2SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          //          distmax_PWL2SVYu = max( distmax_PWL2SVYu, DY[0] - PWL2SVYu
          //          ); distmax_PWL2SVYo = max( distmax_PWL2SVYo, PWL2SVYo -
          //          DY[0] );

          //          double const PWL4SVYu =
          //          PWL4SVY[0].uval({{0,DX[0]},{1,DX[1]}}); double const
          //          PWL4SVYo = PWL4SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          //          distmax_PWL4SVYu = max( distmax_PWL4SVYu, DY[0] - PWL4SVYu
          //          ); distmax_PWL4SVYo = max( distmax_PWL4SVYo, PWL4SVYo -
          //          DY[0] );

          //          double const PWL8SVYu =
          //          PWL8SVY[0].uval({{0,DX[0]},{1,DX[1]}}); double const
          //          PWL8SVYo = PWL8SVY[0].oval({{0,DX[0]},{1,DX[1]}});

          //          distmax_PWL8SVYu = max( distmax_PWL8SVYu, DY[0] - PWL8SVYu
          //          ); distmax_PWL8SVYo = max( distmax_PWL8SVYo, PWL8SVYo -
          //          DY[0] );
        }
      }

      std::cout << std::scientific << std::setprecision(5) << std::right
                << std::setw(14) << rho << std::setw(14)
                << mc::Op<I>::l(XBND[0]) << std::setw(14)
                << mc::Op<I>::u(XBND[0]) << std::setw(14)
                << mc::Op<I>::l(XBND[1]) << std::setw(14)
                << mc::Op<I>::u(XBND[1])
                //<< std::setw(14) << min_MCYcv << std::setw(14) << max_MCYcc
                << std::setw(14) << max(distmax_MCYcv, distmax_MCYcc)
                << std::setw(14) << max(min_Y - min_MCYcv, max_MCYcc - max_Y)
                << std::setw(14) << max(distmax_PWC16SVYu, distmax_PWC16SVYo)
                << std::setw(14)
                << max(min_Y - PWC16SVY[0].l(), PWC16SVY[0].u() - max_Y)
                << std::setw(14) << max(distmax_PWC32SVYu, distmax_PWC32SVYo)
                << std::setw(14)
                << max(min_Y - PWC32SVY[0].l(), PWC32SVY[0].u() - max_Y)
                << std::setw(14) << max(distmax_PWC64SVYu, distmax_PWC64SVYo)
                << std::setw(14)
                << max(min_Y - PWC64SVY[0].l(), PWC64SVY[0].u() - max_Y)
                << std::setw(14) << max(distmax_PWC128SVYu, distmax_PWC128SVYo)
                << std::setw(14)
                << max(min_Y - PWC128SVY[0].l(), PWC128SVY[0].u() - max_Y)
                << std::setw(14) << max(distmax_PWL1SVYu, distmax_PWL1SVYo)
                << std::setw(14)
                << max(min_Y - PWL1SVY[0].l(), PWL1SVY[0].u() - max_Y)
                //                << std::setw(14) << max( distmax_PWL2SVYu,
                //                distmax_PWL2SVYo )
                //                << std::setw(14) << max( min_Y -
                //                PWL2SVY[0].l(),   PWL2SVY[0].u() - max_Y )
                //                << std::setw(14) << max( distmax_PWL4SVYu,
                //                distmax_PWL4SVYo )
                //                << std::setw(14) << max( min_Y -
                //                PWL4SVY[0].l(),   PWL4SVY[0].u() - max_Y )
                //                << std::setw(14) << max( distmax_PWL8SVYu,
                //                distmax_PWL8SVYo )
                //                << std::setw(14) << max( min_Y -
                //                PWL8SVY[0].l(),   PWL8SVY[0].u() - max_Y )
                << std::endl;
    }
#endif

#ifdef ANALYSE_TIME
    std::chrono::time_point<std::chrono::system_clock> start;
    std::chrono::microseconds walltime;
    size_t NREPEAT = 1000;

    // std::vector<I> XBND = XBND0;
    std::vector<I> XBND{XREF[0] + I(-0.1, 0.1), XREF[1] + I(-0.1, 0.1)};

    MCX[0].c(XREF[0]);
    MCX[1].c(XREF[1]);

    std::vector<MC> MCsubX{MC(XBND[0], XREF[0]), MC(XBND[1], XREF[1])},
        MCsubY(1);
    MCsubX[0].sub(2, 0);
    MCsubX[1].sub(2, 1);

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, MCsubY, X, MCsubX );
      ANN.eval(MCsubX.data(), MCsubY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "McCormick subgradient walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    start = std::chrono::system_clock::now();
    std::vector<MC> MCwk;
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, MCwk, Y, MCY, X, MCX );
      ANN.eval(MCX.data(), MCY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "McCormick walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC16SVX(NX), PWC16SVY(1), PWCwk;
    PWC16SVX[0].set(pwcmod, 0, XBND[0], 16);
    PWC16SVX[1].set(pwcmod, 1, XBND[1], 16);
    // DAG.eval( SgY, PWCwk, Y, PWC16SVY, X, PWC16SVX );
    ANN.eval(PWC16SVX.data(), PWC16SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWC16SVY, X, PWC16SVX );
      ANN.eval(PWC16SVX.data(), PWC16SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWC16 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC32SVX(NX), PWC32SVY(1);
    PWC32SVX[0].set(pwcmod, 0, XBND[0], 32);
    PWC32SVX[1].set(pwcmod, 1, XBND[1], 32);
    // DAG.eval( SgY, PWCwk, Y, PWC32SVY, X, PWC32SVX );
    ANN.eval(PWC32SVX.data(), PWC32SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWC32SVY, X, PWC32SVX );
      ANN.eval(PWC32SVX.data(), PWC32SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWC32 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC64SVX(NX), PWC64SVY(1);
    PWC64SVX[0].set(pwcmod, 0, XBND[0], 64);
    PWC64SVX[1].set(pwcmod, 1, XBND[1], 64);
    // DAG.eval( SgY, PWCwk, Y, PWC64SVY, X, PWC64SVX );
    ANN.eval(PWC64SVX.data(), PWC64SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWC64SVY, X, PWC64SVX );
      ANN.eval(PWC64SVX.data(), PWC64SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWC64 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWC128SVX(NX), PWC128SVY(1);
    PWC128SVX[0].set(pwcmod, 0, XBND[0], 128);
    PWC128SVX[1].set(pwcmod, 1, XBND[1], 128);
    // DAG.eval( SgY, PWCwk, Y, PWC128SVY, X, PWC128SVX );
    ANN.eval(PWC128SVX.data(), PWC128SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, PWCwk, Y, PWC128SVY, X, PWC128SVX );
      ANN.eval(PWC128SVX.data(), PWC128SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWC128 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    mc::PWCU::options.SLOPEUSE = 2;
    std::vector<PWCSV> PWCS4SVX(NX), PWCS4SVY(1);
    PWCS4SVX[0].set(pwcmod, 0, XBND[0], 4);
    PWCS4SVX[1].set(pwcmod, 1, XBND[1], 4);
    // DAG.eval( SgY, PWCwk, Y, PWCS4SVY, X, PWCS4SVX );
    ANN.eval(PWCS4SVX.data(), PWCS4SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWC4SVY, X, PWC4SVX );
      ANN.eval(PWCS4SVX.data(), PWCS4SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWCS4 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWCS8SVX(NX), PWCS8SVY(1);
    PWCS8SVX[0].set(pwcmod, 0, XBND[0], 8);
    PWCS8SVX[1].set(pwcmod, 1, XBND[1], 8);
    // DAG.eval( SgY, PWCwk, Y, PWCS8SVY, X, PWCS8SVX );
    ANN.eval(PWCS8SVX.data(), PWCS8SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWCS8SVY, X, PWCS8SVX );
      ANN.eval(PWCS8SVX.data(), PWCS8SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWCS8 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWCSV> PWCS16SVX(NX), PWCS16SVY(1);
    PWCS16SVX[0].set(pwcmod, 0, XBND[0], 16);
    PWCS16SVX[1].set(pwcmod, 1, XBND[1], 16);
    // DAG.eval( SgY, PWCwk, Y, PWCS16SVY, X, PWCS16SVX );
    ANN.eval(PWCS16SVX.data(), PWCS16SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, Y, PWCS16SVY, X, PWCS16SVX );
      ANN.eval(PWCS16SVX.data(), PWCS16SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWCS16 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    // NREPEAT = 200;
    std::vector<PWLSV> PWL1SVX(NX), PWL1SVY(1), PWLwk;
    PWL1SVX[0].set(pwlmod, 0, XBND[0], 1);
    PWL1SVX[1].set(pwlmod, 1, XBND[1], 1);
    DAG.eval(SgY, Y, PWL1SVY, X, PWL1SVX);
    ANN.eval(PWL1SVX.data(), PWL1SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, PWLwk, Y, PWL1SVY, X, PWL1SVX );
      ANN.eval(PWL1SVX.data(), PWL1SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWL1 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL2SVX(NX), PWL2SVY(1);
    PWL2SVX[0].set(pwlmod, 0, XBND[0], 2);
    PWL2SVX[1].set(pwlmod, 1, XBND[1], 2);
    DAG.eval(SgY, Y, PWL2SVY, X, PWL2SVX);
    ANN.eval(PWL2SVX.data(), PWL2SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, PWLwk, Y, PWL2SVY, X, PWL2SVX );
      ANN.eval(PWL2SVX.data(), PWL2SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWL2 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL4SVX(NX), PWL4SVY(1);
    PWL4SVX[0].set(pwlmod, 0, XBND[0], 4);
    PWL4SVX[1].set(pwlmod, 1, XBND[1], 4);
    DAG.eval(SgY, Y, PWL4SVY, X, PWL4SVX);
    ANN.eval(PWL4SVX.data(), PWL4SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, PWLwk, Y, PWL4SVY, X, PWL4SVX );
      ANN.eval(PWL4SVX.data(), PWL4SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWL4 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;

    std::vector<PWLSV> PWL8SVX(NX), PWL8SVY(1);
    PWL8SVX[0].set(pwlmod, 0, XBND[0], 8);
    PWL8SVX[1].set(pwlmod, 1, XBND[1], 8);
    DAG.eval(SgY, Y, PWL8SVY, X, PWL8SVX);
    ANN.eval(PWL8SVX.data(), PWL8SVY.data());

    start = std::chrono::system_clock::now();
    for (unsigned i = 0; i < NREPEAT; ++i)
      // DAG.eval( SgY, PWLwk, Y, PWL8SVY, X, PWL8SVX );
      ANN.eval(PWL8SVX.data(), PWL8SVY.data());
    walltime = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now() - start);
    std::cout << "Superposition PWL8 walltime: "
              << (walltime.count() * 1e-6) / (double)NREPEAT << std::endl;
#endif
  }

  catch (mc::FFBase::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in factorable function manipulation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && \
    !defined(MC__USE_BOOST)
  catch (I::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in natural interval extension:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif

  catch (MC::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in McCormick relaxation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  catch (PWCSM::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in superposition relaxation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  catch (PWLSM::Exceptions& eObj)
  {
    std::cerr << "Error " << eObj.ierr()
              << " in superposition relaxation:" << std::endl
              << eObj.what() << std::endl
              << "Aborts." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
