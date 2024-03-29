#define TEST4   	// <-- select test interval matrix here
////////////////////////////////////////////////////////////////////////

#include <fstream>
#include <iomanip>

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
   typedef boost::numeric::interval_lib::save_state<boost::numeric::interval_lib::rounded_transc_opp<double>> T_boost_round;
   typedef boost::numeric::interval_lib::checking_base<double> T_boost_check;
   typedef boost::numeric::interval_lib::policies<T_boost_round,T_boost_check> T_boost_policy;
   typedef boost::numeric::interval<double,T_boost_policy> I;
  #else
   #include "interval.hpp"
   typedef mc::Interval I;
  #endif
 #endif
#endif

#include "specbnd.hpp"
typedef mc::Specbnd<I> SB;

////////////////////////////////////////////////////////////////////////
// This is Example 2 in:
// M. Hladik, Applied Mathematics & Computation 219 (2013) 5584-5591
#if defined( TEST1 )
const unsigned N = 3;
const I A[N*N] = { I(0.), I(-1.),           I(-1.),
                   I(2.), I(-1.399,-0.001), I(0.),
                   I(1.), I(0.5),           I(-1.)  };

////////////////////////////////////////////////////////////////////////
// This is Example 3 in:
// M. Hladik, Applied Mathematics & Computation 219 (2013) 5584-5591
#elif defined( TEST2 )
const unsigned N = 3;
const I A[N*N] = { I(-1.),    I(0.),     I(-1.,1.),
                   I(0.),     I(-1),     I(-1.,1.),
                   I(-1.,1.), I(-1.,1.), I(0.1)     };

////////////////////////////////////////////////////////////////////////
// This is Example 4 in:
// M. Hladik, Applied Mathematics & Computation 219 (2013) 5584-5591
#elif defined( TEST3 )
const unsigned N = 4;
const I A[N*N] = { I(-3.,-2.), I(4.,5.),   I(4.,6.),   I(-1.,1.5),
                   I(-4.,-3.), I(-4.,-3.), I(-4.,-3.), I(1.,2.),
                   I(-5.,-4.), I(2.,3.),   I(-5.,-4.), I(-1.,0.),
                   I(-1.,0.1), I(0.,1.),   I(1.,2.),   I(-4.,2.5)  };

////////////////////////////////////////////////////////////////////////
// This is Example 5 in:
// M. Hladik, Applied Mathematics & Computation 219 (2013) 5584-5591
#elif defined( TEST4 )
const unsigned N = 3;
const I A[N*N] = { I(-1.8,-1.2), I(0.4,0.6),   I(0.8,1.2),
                   I(-1.2,-0.8), I(-3.6,-2.4), I(0.8,1.2),
                   I(-0.6,-0.4), I(-1.8,-1.2), I(-3.,-2.)  };
#endif

////////////////////////////////////////////////////////////////////////
int main()
////////////////////////////////////////////////////////////////////////
{  
  try{ 
    // Compute spectral bounds of interval matrix: Gerhsgorin method
    SB::options.HESSBND = SB::Options::GERSHGORIN;
    std::pair<double,double> rebndG = SB::spectral_bound_re( N, A );
    std::cout << "\nBound on real part of spectrum (Gershgorin): "
              << I(rebndG.first,rebndG.second) << std::endl;
    std::pair<double,double> imbndG = SB::spectral_bound_im( N, A );
    std::cout << "\nBound on imaginary part of spectrum (Gershgorin): "
              << I(imbndG.first,imbndG.second) << std::endl;

    // Compute spectral bounds of interval matrix: Herz & Rohn method
    SB::options.HESSBND = SB::Options::HERTZROHN;
    std::pair<double,double> rebndHR = SB::spectral_bound_re( N, A );
    std::cout << "\nBound on real part of spectrum (Hertz&Rohn): "
              << I(rebndHR.first,rebndHR.second) << std::endl;
    std::pair<double,double> imbndHR = SB::spectral_bound_im( N, A );
    std::cout << "\nBound on imaginary part of spectrum (Hertz&Rohn): "
              << I(imbndHR.first,imbndHR.second) << std::endl;
  }

#if !defined(MC__USE_PROFIL) && !defined(MC__USE_FILIB) && !defined(MC__USE_BOOST)
  catch( I::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
         << " in natural interval extension:" << std::endl
	 << eObj.what() << std::endl
         << "Aborts." << std::endl;
    return eObj.ierr();
  }
#endif
  catch( SB::Exceptions &eObj ){
    std::cerr << "Error " << eObj.ierr()
         << " in spectral bound computation:" << std::endl
	 << eObj.what() << std::endl
         << "Aborts." << std::endl;
    return eObj.ierr();
  }

  return 0;
}
