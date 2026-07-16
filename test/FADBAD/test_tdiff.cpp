// Systematic test suite for mc::T Taylor-mode AD.
//   - every supported function: Taylor coefficients (orders 0..P) vs genuine
//     fadbad::T (native) and vs closed-form analytic coefficients
//   - binary/scalar operators; pow (int/real/T/scalar); xlog, arh
//   - higher-order derivative read-off (k! * coeff[k] == f^(k))
//   - interval soundness; incremental eval & reset() caching consistency
//   - deep composition / shared subexpressions
//   - multi-variable directional expansion and a small function system
// NB: raw fadbad::T's pow(T,int) is unreliable here, so integer powers are
// checked against analytic/explicit-multiplication instead of fadbad.
#include <cstdio>
#include <cmath>
#include <string>
#include "tdiff.hpp"
#include "interval.hpp"
#include "fadbad.h"
#include "tadiff.h"

typedef mc::Interval I;
static const unsigned P=9;

static int g_fail=0,g_tot=0,s_fail=0,s_tot=0; static std::string g_sec;
static void sec(const char* n){ if(!g_sec.empty()) std::printf("  -> %-40s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  g_sec=n; s_fail=s_tot=0; std::printf("[%s]\n",n); }
static void CHK(const std::string& nm,double got,double ref,double tol){ ++g_tot;++s_tot; double e=std::fabs(got-ref);
  if(e>tol*(1.+std::fabs(ref))){ ++g_fail;++s_fail; std::printf("  FAIL %-40s got=% .10e ref=% .10e |e|=%.2e\n",nm.c_str(),got,ref,e);} }
static void CHKb(const std::string& nm,bool ok){ ++g_tot;++s_tot; if(!ok){ ++g_fail;++s_fail; std::printf("  FAIL %-40s\n",nm.c_str()); } }
static void done(){ if(!g_sec.empty()) std::printf("  -> %-40s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  std::printf("\n%s  (%d/%d checks passed)\n", g_fail?"*** FAILURES ***":"ALL PASSED", g_tot-g_fail, g_tot); }

// coefficient-by-coefficient cross-check against genuine fadbad::T
template<class ADf> void cfb(const char* nm, ADf adf, double x0){
  mc::T<double> x=x0; x[1]=1.; auto f=adf(x); f.eval(P);
  fadbad::T<double> fx=x0; fx[1]=1.; auto ff=adf(fx); ff.eval(P);
  for(unsigned k=0;k<=P;k++) CHK(std::string(nm)+" c["+std::to_string(k)+"]", f[k], ff[k], 1e-11);
}

int main()
{
  // ======================= 1. functions vs genuine fadbad::T ==============
  sec("Taylor coefficients vs fadbad::T (native functions), orders 0..9");
  cfb("exp", [](auto&z){return exp(z);}, 0.3);
  cfb("log", [](auto&z){return log(z);}, 1.7);
  cfb("sqrt",[](auto&z){return sqrt(z);},1.7);
  cfb("sqr", [](auto&z){return sqr(z);}, 0.6);
  cfb("sin", [](auto&z){return sin(z);}, 0.6);
  cfb("cos", [](auto&z){return cos(z);}, 0.6);
  cfb("tan", [](auto&z){return tan(z);}, 0.6);
  cfb("asin",[](auto&z){return asin(z);},0.3);
  cfb("acos",[](auto&z){return acos(z);},0.3);
  cfb("atan",[](auto&z){return atan(z);},0.6);
  cfb("sinh",[](auto&z){return sinh(z);},0.6);
  cfb("cosh",[](auto&z){return cosh(z);},0.6);
  cfb("tanh",[](auto&z){return tanh(z);},0.6);
  cfb("1/x", [](auto&z){return 1./z;},   1.7);        // operator/(double,T) == inv path

  // ======================= 2. analytic closed forms =======================
  sec("Taylor coefficients vs analytic closed forms");
  {
    // exp(x0+t): coeff[k] = e^x0 / k!
    mc::T<double> x=0.3; x[1]=1.; mc::T<double> f=exp(x); f.eval(P);
    double fact=1.; for(unsigned k=0;k<=P;k++){ if(k)fact*=k; CHK("exp c["+std::to_string(k)+"]", f[k], std::exp(0.3)/fact, 1e-12); }
    // 1/(a-x): coeff[k] = 1/(a-x0)^{k+1}
    mc::T<double> y=0.3; y[1]=1.; mc::T<double> g=1./(2.-y); g.eval(P);
    double a=2.-0.3; for(unsigned k=0;k<=P;k++) CHK("1/(2-x) c["+std::to_string(k)+"]", g[k], std::pow(a,-(double)(k+1)), 1e-12);
    // sin(x0+t): coeff[k] = sin(x0 + k pi/2)/k!
    mc::T<double> z=0.6; z[1]=1.; mc::T<double> s=sin(z); s.eval(P);
    double f2=1.; for(unsigned k=0;k<=P;k++){ if(k)f2*=k; CHK("sin c["+std::to_string(k)+"]", s[k], std::sin(0.6+k*M_PI/2)/f2, 1e-12); }
    // inv(x) == 1/x
    mc::T<double> w=1.7; w[1]=1.; mc::T<double> iw=inv(w), rw=1./w; iw.eval(P); rw.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("inv==1/x c["+std::to_string(k)+"]", iw[k], rw[k], 1e-13);
  }

  // ======================= 3. binary & scalar operators ===================
  sec("binary / scalar operators vs fadbad::T");
  cfb("a+b*sin",[](auto&z){return z + z*sin(z);}, 0.5);
  cfb("a-b",    [](auto&z){return exp(z) - sqrt(z+2.);}, 0.5);
  cfb("a*b",    [](auto&z){return sin(z)*exp(z);}, 0.5);
  cfb("a/b",    [](auto&z){return (z+2.)/(cos(z)+2.);}, 0.5);
  cfb("scal+",  [](auto&z){return z+3.;}, 0.5);
  cfb("scal-",  [](auto&z){return 3.-z;}, 0.5);
  cfb("scal*",  [](auto&z){return 2.5*z;}, 0.5);
  cfb("scal/",  [](auto&z){return z/4.;}, 0.5);
  cfb("neg",    [](auto&z){return -exp(z);}, 0.5);

  // ======================= 4. powers ======================================
  sec("pow: integer (analytic), real, T^T, scalar^T");
  {
    // integer power (x+2)^3 vs analytic binomial and explicit multiplication
    mc::T<double> x=0.3; x[1]=1.; mc::T<double> p=pow(x+2.,3), pm=(x+2.)*(x+2.)*(x+2.);
    p.eval(P); pm.eval(P); double aa=2.3, pr[4]={aa*aa*aa,3*aa*aa,3*aa,1.};
    for(unsigned k=0;k<=P;k++){ double ref=k<4?pr[k]:0.;
      CHK("pow(x+2,3) c["+std::to_string(k)+"]", p[k], ref, 1e-11);
      CHK("pow==mult c["+std::to_string(k)+"]", p[k], pm[k], 1e-12); }
    // real power x^2.5 == exp(2.5 log x)
    mc::T<double> y=1.4; y[1]=1.; mc::T<double> pr2=pow(y,2.5), pe=exp(2.5*log(y));
    pr2.eval(P); pe.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("pow(x,2.5)==explog c["+std::to_string(k)+"]", pr2[k], pe[k], 1e-11);
    // T^T : x^y with x=x0+t, y=const-series
    mc::T<double> u=1.3; u[1]=1.; mc::T<double> vv=0.7; mc::T<double> pt=pow(u,vv), pd=exp(vv*log(u));
    pt.eval(P); pd.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("pow(x,y)==explog c["+std::to_string(k)+"]", pt[k], pd[k], 1e-11);
    // scalar^T : 2^x
    mc::T<double> s=0.4; s[1]=1.; mc::T<double> ps=pow(2.,s), pse=exp(s*std::log(2.));
    ps.eval(P); pse.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("pow(2,x)==explog c["+std::to_string(k)+"]", ps[k], pse[k], 1e-11);
  }

  // ======================= 5. xlog / arh ==================================
  sec("xlog / arh vs definitional forms");
  {
    mc::T<double> x=0.4; x[1]=1.; mc::T<double> a=xlog(x+2.), ar=(x+2.)*log(x+2.); a.eval(P); ar.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("xlog c["+std::to_string(k)+"]", a[k], ar[k], 1e-12);
    mc::T<double> y=0.4; y[1]=1.; mc::T<double> b=arh(y+2.,0.5), be=exp(-0.5/(y+2.)); b.eval(P); be.eval(P);
    for(unsigned k=0;k<=P;k++) CHK("arh c["+std::to_string(k)+"]", b[k], be[k], 1e-12);
  }

  // ======================= 6. higher-order derivative read-off ============
  sec("derivative read-off: k! * coeff[k] == f^(k)(x0)");
  {
    // f = x^2 e^x : f' = (2x+x^2)e^x, f'' = (2+4x+x^2)e^x, f''' = (6+6x+x^2)e^x
    mc::T<double> x=0.5; x[1]=1.; mc::T<double> f=sqr(x)*exp(x); f.eval(3);
    double e=std::exp(0.5);
    CHK("f",   f[0]*1.,  0.25*e, 1e-10);
    CHK("f'",  f[1]*1.,  (2*0.5+0.25)*e, 1e-10);
    CHK("f''", f[2]*2.,  (2+4*0.5+0.25)*e, 1e-10);
    CHK("f'''",f[3]*6.,  (6+6*0.5+0.25)*e, 1e-10);
  }

  // ======================= 7. interval soundness ==========================
  sec("T<Interval> soundness: coeff enclosures contain point coeffs");
  {
    mc::T<I> x=I(0.25,0.35); x[1]=1.; mc::T<I> f=exp(sin(x))+atan(x)+sqrt(x+2.); f.eval(6);
    mc::T<double> xp=0.3; xp[1]=1.; mc::T<double> fp=exp(sin(xp))+atan(xp)+sqrt(xp+2.); fp.eval(6);
    bool ok=true; for(unsigned k=0;k<=6;k++){ I c=f[k]; double p=fp[k]; ok&=(c.l()<=p+1e-9 && p<=c.u()+1e-9); }
    CHKb("all 0..6 coeff enclosures contain truth", ok);
  }

  // ======================= 8. incremental eval & reset ====================
  sec("incremental eval & reset() caching consistency");
  {
    mc::T<double> x=0.3; x[1]=1.; mc::T<double> f=exp(sin(x))+tan(x)+x*x*x;
    f.eval(3); double c[4]; for(int k=0;k<4;k++)c[k]=f[k];
    f.eval(P);                                   // extend
    for(int k=0;k<4;k++) CHK("stable c["+std::to_string(k)+"] after extend", f[k], c[k], 0);
    double cP=f[P];
    f.reset(); f.eval(P);                         // recompute from scratch
    CHK("reset reproduces c[P]", f[P], cP, 1e-14);
    for(int k=0;k<4;k++) CHK("reset reproduces c["+std::to_string(k)+"]", f[k], c[k], 1e-14);
    CHKb("length()==P+1", f.length()==P+1);
  }

  // ======================= 9. deep composition ============================
  sec("deep nested composition vs fadbad::T");
  cfb("nest6", [](auto&z){return tanh(atan(exp(sin(sqrt(z+2.)))));}, 0.4);

  // ======================= 10. directional & function system ==============
  sec("multi-variable directional expansion + function system vs fadbad::T");
  {
    const double x0=0.5,y0=0.8;
    auto F0=[](auto&x,auto&y){return exp(x*y);};
    auto F1=[](auto&x,auto&y){return sin(x)+cos(y);};
    auto F2=[](auto&x,auto&y){return sqrt(x*x+y*y+1.);};
    // expand in the x-direction: x[1]=1, y[1]=0
    auto run=[&](auto Fk, const char* nm){
      mc::T<double> x=x0,y=y0; x[1]=1.; auto f=Fk(x,y); f.eval(P);
      fadbad::T<double> fx=x0,fy=y0; fx[1]=1.; auto ff=Fk(fx,fy); ff.eval(P);
      for(unsigned k=0;k<=P;k++) CHK(std::string(nm)+" c["+std::to_string(k)+"]", f[k], ff[k], 1e-11);
    };
    run(F0,"sys f0=exp(xy)");
    run(F1,"sys f1=sin+cos");
    run(F2,"sys f2=norm");
    // expand same system in the y-direction
    auto runy=[&](auto Fk, const char* nm){
      mc::T<double> x=x0,y=y0; y[1]=1.; auto f=Fk(x,y); f.eval(P);
      fadbad::T<double> fx=x0,fy=y0; fy[1]=1.; auto ff=Fk(fx,fy); ff.eval(P);
      for(unsigned k=0;k<=P;k++) CHK(std::string(nm)+" c["+std::to_string(k)+"]", f[k], ff[k], 1e-11);
    };
    runy(F0,"sys-y f0");
  }

  done();
  return g_fail?1:0;
}
