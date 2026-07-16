// Demonstrates that the reduced derivative forms (sqr(x) instead of x*x,
// a/sqr(b) instead of (a/b)/b) yield tighter McCormick relaxations -- i.e.
// they suffer less from the dependency problem -- while remaining sound.
#include <cstdio>
#include <cmath>
#include "mccormick.hpp"
#include "interval.hpp"
#include "fdiff.hpp"

typedef mc::Interval        I;
typedef mc::McCormick<I>    MC;

static double gap( const MC& x ){ return x.cc() - x.cv(); }  // relaxation width at ref point

int main()
{
  std::printf("Relaxation gap (cc-cv): smaller = tighter (less dependency blow-up)\n\n");

  // ---- self-product vs sqr ------------------------------------------------
  // x over [-1,2], reference point 0.5, seeded subgradient.
  {
    MC X( I(-1.,2.), 0.5 ); X.sub(1,0);
    MC p_xx  = X*X;                    // loose: bilinear product rule
    MC p_sqr = mc::Op<MC>::sqr(X);     // tight: uses x^2 envelope
    std::printf("[x*x vs sqr(x)] over x in [-1,2] @0.5\n");
    std::printf("  x*x   : cv=% .4f cc=% .4f  gap=%.4f\n", p_xx.cv(),  p_xx.cc(),  gap(p_xx));
    std::printf("  sqr(x): cv=% .4f cc=% .4f  gap=%.4f  -> %s\n",
                p_sqr.cv(), p_sqr.cc(), gap(p_sqr),
                gap(p_sqr) <= gap(p_xx)+1e-12 ? "tighter/equal" : "WORSE");
  }

  // ---- inv derivative coefficient: -(fv*fv) vs -sqr(fv) -------------------
  {
    MC A( I(0.5,3.), 1.5 ); A.sub(1,0);
    MC fv = mc::Op<MC>::inv(A);        // 1/a
    MC c_old = -( fv*fv );             // old coefficient
    MC c_new = -mc::Op<MC>::sqr(fv);   // new coefficient
    std::printf("[d(1/a): -(fv*fv) vs -sqr(fv)] a in [0.5,3] @1.5\n");
    std::printf("  old   : gap=%.4f\n", gap(c_old));
    std::printf("  new   : gap=%.4f  -> %s\n", gap(c_new),
                gap(c_new) <= gap(c_old)+1e-12 ? "tighter/equal" : "WORSE");
  }

  // ---- division coefficient: -(a/b)/b vs -a/sqr(b) ------------------------
  {
    const double a = 2.0;
    MC B( I(0.5,3.), 1.5 ); B.sub(1,0);
    MC c_old = -( (a/B)/B );            // old: reuses b
    MC c_new = -( a/mc::Op<MC>::sqr(B) );
    std::printf("[d(a/b): -(a/b)/b vs -a/sqr(b)] b in [0.5,3] @1.5, a=2\n");
    std::printf("  old   : cv=% .4f cc=% .4f gap=%.4f\n", c_old.cv(), c_old.cc(), gap(c_old));
    std::printf("  new   : cv=% .4f cc=% .4f gap=%.4f  -> %s\n",
                c_new.cv(), c_new.cc(), gap(c_new),
                gap(c_new) <= gap(c_old)+1e-12 ? "tighter/equal" : "WORSE");
  }

  // ---- soundness of the AD path over F<MC>: derivative relaxation must
  //      bracket the true point derivative -----------------------------------
  {
    typedef mc::F<MC> FMC;
    MC bval( I(0.5,3.), 1.5 ); bval.sub(1,0);
    FMC B = bval; B.diff(0,1);
    FMC R = 2.0 / B;                    // 2/b, exercises operator/(scalar,F)
    MC d = R.deriv(0);                  // relaxation of d(2/b)=-2/b^2 at b=1.5
    double truth = -2.0/(1.5*1.5);      // -0.8889
    bool ok = d.cv() <= truth+1e-9 && truth <= d.cc()+1e-9
              && d.l() <= truth+1e-9 && truth <= d.u()+1e-9;
    std::printf("\n[F<MC> soundness] d(2/b)@1.5 truth=% .5f  relax cv/cc=[% .4f,% .4f] I=[% .4f,% .4f] : %s\n",
                truth, d.cv(), d.cc(), d.l(), d.u(), ok?"contains":"FAIL");
  }
  return 0;
}
