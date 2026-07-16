// Systematic test suite for mc::F (dense) and mc::SF (sparse) forward-mode AD.
//   - every elementary function: value + derivative vs analytic, finite
//     differences, and (where native) fadbad::F, plus SF==F consistency
//   - binary/scalar operators, pow/min/max/lmtd/rlmtd/cheb/erf/xlog/arh
//   - DAG reuse, constant/mixed-operand handling
//   - VECTOR-VALUED system R^4 -> R^5: full Jacobian via F, SF, fadbad::F and
//     central differences, all cross-checked; sparse pattern verified
//   - higher-order (F<F> Hessian), interval & McCormick soundness
#include <cstdio>
#include <cmath>
#include <string>
#include "fdiff.hpp"
#include "interval.hpp"
#include "mccormick.hpp"
#include "fadbad.h"
#include "fadiff.h"

typedef mc::Interval             I;
typedef mc::McCormick<I>         MC;
static const double H = 1e-6;

static int g_fail=0,g_tot=0,s_fail=0,s_tot=0; static std::string g_sec;
static void sec(const char* n){ if(!g_sec.empty()) std::printf("  -> %-40s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  g_sec=n; s_fail=s_tot=0; std::printf("[%s]\n",n); }
static void CHK(const std::string& nm,double got,double ref,double tol){
  ++g_tot;++s_tot; double e=std::fabs(got-ref);
  if(e>tol*(1.+std::fabs(ref))){ ++g_fail;++s_fail;
    std::printf("  FAIL %-40s got=% .10e ref=% .10e |e|=%.2e\n",nm.c_str(),got,ref,e); } }
static void CHKb(const std::string& nm,bool ok){ ++g_tot;++s_tot; if(!ok){ ++g_fail;++s_fail; std::printf("  FAIL %-40s\n",nm.c_str()); } }
static std::string N2(const char* a,const char* b){ return std::string(a)+" "+b; }
static void done(){ if(!g_sec.empty()) std::printf("  -> %-40s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  std::printf("\n%s  (%d/%d checks passed)\n", g_fail?"*** FAILURES ***":"ALL PASSED", g_tot-g_fail, g_tot); }

// ---- generic unary checks (mc::F, mc::SF, analytic, finite diff) ------------
template<class ADf,class VFn,class DFn>
void uni(const char* nm, ADf adf, VFn vf, DFn df, double x0){
  mc::F<double> x=x0; x.diff(0,1); auto f=adf(x);
  mc::SF<double> sx=x0; sx.diff(0,1); auto sf=adf(sx);
  double v=vf(x0), an=df(x0), fd=(vf(x0+H)-vf(x0-H))/(2*H);
  CHK(N2(nm,"value"),        f.val(), v, 1e-12);
  CHK(N2(nm,"d/dx analytic"),f.deriv(0), an, 1e-9);
  CHK(N2(nm,"d/dx findiff"), f.deriv(0), fd, 1e-5);
  CHK(N2(nm,"SF value"),     sf.val(), v, 1e-12);
  CHK(N2(nm,"SF deriv"),     sf.deriv(0), f.deriv(0), 1e-13);
  CHKb(N2(nm,"SF nnz==1"),   sf.size()==1);
}
// ---- unary cross-check against genuine fadbad::F (native funcs only) --------
template<class ADf> void unifb(const char* nm, ADf adf, double x0){
  mc::F<double> x=x0; x.diff(0,1); auto f=adf(x);
  fadbad::F<double> fx=x0; fx.diff(0,1); auto ff=adf(fx);
  CHK(N2(nm,"F vs fadbad val"), f.val(), ff.val(), 1e-13);
  CHK(N2(nm,"F vs fadbad der"), f.deriv(0), ff.deriv(0), 1e-13);
}

// ---- vector-valued system  R^4 -> R^5 (all funcs native to fadbad & std) ----
template<class X> void SYS(const X* v, X* o){
  o[0] = v[0]*v[1] - v[2]/v[3] + sin(v[0]);
  o[1] = exp(0.5*v[1]) + sqrt(v[2]*v[2]+1.) - v[3];
  o[2] = v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3];
  o[3] = atan(v[0]*v[3]) + tanh(v[1]) - log(v[2]+2.);
  o[4] = v[0]/(v[1]+2.) + v[2]*cos(v[3]);
}
static const int NIN=4, MOUT=5;

int main()
{
  // ======================= 1. unary elementary functions ==================
  sec("unary elementary functions: value + derivative");
  uni("inv", [](auto&z){return inv(z);},  [](double z){return 1./z;},        [](double z){return -1./(z*z);},          1.7);
  uni("sqr", [](auto&z){return sqr(z);},  [](double z){return z*z;},         [](double z){return 2*z;},                0.7);
  uni("sqrt",[](auto&z){return sqrt(z);}, [](double z){return std::sqrt(z);},[](double z){return 0.5/std::sqrt(z);},   1.7);
  uni("exp", [](auto&z){return exp(z);},  [](double z){return std::exp(z);}, [](double z){return std::exp(z);},        0.6);
  uni("log", [](auto&z){return log(z);},  [](double z){return std::log(z);}, [](double z){return 1./z;},               1.7);
  uni("xlog",[](auto&z){return xlog(z);}, [](double z){return z*std::log(z);},[](double z){return std::log(z)+1.;},    1.7);
  uni("sin", [](auto&z){return sin(z);},  [](double z){return std::sin(z);}, [](double z){return std::cos(z);},        0.6);
  uni("cos", [](auto&z){return cos(z);},  [](double z){return std::cos(z);}, [](double z){return -std::sin(z);},       0.6);
  uni("tan", [](auto&z){return tan(z);},  [](double z){return std::tan(z);}, [](double z){double c=std::cos(z);return 1./(c*c);}, 0.6);
  uni("asin",[](auto&z){return asin(z);}, [](double z){return std::asin(z);},[](double z){return 1./std::sqrt(1-z*z);},0.4);
  uni("acos",[](auto&z){return acos(z);}, [](double z){return std::acos(z);},[](double z){return -1./std::sqrt(1-z*z);},0.4);
  uni("atan",[](auto&z){return atan(z);}, [](double z){return std::atan(z);},[](double z){return 1./(1+z*z);},         0.6);
  uni("sinh",[](auto&z){return sinh(z);}, [](double z){return std::sinh(z);},[](double z){return std::cosh(z);},       0.6);
  uni("cosh",[](auto&z){return cosh(z);}, [](double z){return std::cosh(z);},[](double z){return std::sinh(z);},       0.6);
  uni("tanh",[](auto&z){return tanh(z);}, [](double z){return std::tanh(z);},[](double z){double t=std::tanh(z);return 1-t*t;}, 0.6);
  uni("fabs",[](auto&z){return fabs(z);}, [](double z){return std::fabs(z);},[](double z){return z<0?-1.:1.;},         0.6);
  uni("erf", [](auto&z){return erf(z);},  [](double z){return std::erf(z);}, [](double z){return 2./std::sqrt(M_PI)*std::exp(-z*z);}, 0.6);
  uni("erfc",[](auto&z){return erfc(z);}, [](double z){return std::erfc(z);},[](double z){return -2./std::sqrt(M_PI)*std::exp(-z*z);},0.6);

  // ======================= 2. cross-check vs fadbad::F ====================
  sec("unary vs genuine fadbad::F (native functions)");
  unifb("sqr", [](auto&z){return sqr(z);}, 0.7);
  unifb("sqrt",[](auto&z){return sqrt(z);},1.7);
  unifb("exp", [](auto&z){return exp(z);}, 0.6);
  unifb("log", [](auto&z){return log(z);}, 1.7);
  unifb("sin", [](auto&z){return sin(z);}, 0.6);
  unifb("cos", [](auto&z){return cos(z);}, 0.6);
  unifb("tan", [](auto&z){return tan(z);}, 0.6);
  unifb("asin",[](auto&z){return asin(z);},0.4);
  unifb("acos",[](auto&z){return acos(z);},0.4);
  unifb("atan",[](auto&z){return atan(z);},0.6);
  unifb("sinh",[](auto&z){return sinh(z);},0.6);
  unifb("cosh",[](auto&z){return cosh(z);},0.6);
  unifb("tanh",[](auto&z){return tanh(z);},0.6);

  // ======================= 3. binary & scalar operators ===================
  sec("binary/scalar operators: value + gradient");
  {
    const double a0=1.3,b0=0.7;
    auto grad2=[&](auto fn, double da, double db, const char* nm){
      mc::F<double> a=a0,b=b0; a.diff(0,2); b.diff(1,2); auto f=fn(a,b);
      CHK(N2(nm,"d/da"), f.deriv(0), da, 1e-9);
      CHK(N2(nm,"d/db"), f.deriv(1), db, 1e-9);
      mc::SF<double> sa=a0,sb=b0; sa.diff(0,2); sb.diff(1,2); auto sf=fn(sa,sb);
      CHK(N2(nm,"SF d/da"), sf.deriv(0), f.deriv(0), 1e-13);
      CHK(N2(nm,"SF d/db"), sf.deriv(1), f.deriv(1), 1e-13);
    };
    grad2([](auto&a,auto&b){return a+b;}, 1., 1., "a+b");
    grad2([](auto&a,auto&b){return a-b;}, 1., -1., "a-b");
    grad2([](auto&a,auto&b){return a*b;}, b0, a0, "a*b");
    grad2([](auto&a,auto&b){return a/b;}, 1./b0, -a0/(b0*b0), "a/b");
    grad2([](auto&a,auto&b){return pow(a,b);}, b0*std::pow(a0,b0-1), std::pow(a0,b0)*std::log(a0), "pow(a,b)");
    // scalar variants
    mc::F<double> a=a0; a.diff(0,1);
    CHK("a+3 d/da",  (a+3.).deriv(0), 1., 1e-12);
    CHK("3-a d/da",  (3.-a).deriv(0), -1., 1e-12);
    CHK("a*3 d/da",  (a*3.).deriv(0), 3., 1e-12);
    CHK("a/4 d/da",  (a/4.).deriv(0), 0.25, 1e-12);
    CHK("2/a d/da",  (2./a).deriv(0), -2./(a0*a0), 1e-9);
    CHK("pow(a,3)",  pow(a,3).deriv(0), 3*a0*a0, 1e-9);
    CHK("pow(2,a)",  pow(2.,a).deriv(0), std::pow(2.,a0)*std::log(2.), 1e-9);
  }

  // ======================= 4. min/max/lmtd/rlmtd/cheb =====================
  sec("min / max / lmtd / rlmtd / cheb");
  {
    const double a0=1.3,b0=0.7;
    mc::F<double> a=a0,b=b0; a.diff(0,2); b.diff(1,2);
    CHK("max grad-a", max(a,b).deriv(0), 1., 1e-12);   // a>b
    CHK("max grad-b", max(a,b).deriv(1), 0., 1e-12);
    CHK("min grad-a", min(a,b).deriv(0), 0., 1e-12);
    CHK("min grad-b", min(a,b).deriv(1), 1., 1e-12);
    // lmtd non-equal: g=(a-b)/(log a-log b); check gradient vs finite diff
    auto lm=[](double p,double q){return (p-q)/(std::log(p)-std::log(q));};
    CHK("lmtd d/da FD", lmtd(a,b).deriv(0), (lm(a0+H,b0)-lm(a0-H,b0))/(2*H), 1e-6);
    CHK("lmtd d/db FD", lmtd(a,b).deriv(1), (lm(a0,b0+H)-lm(a0,b0-H))/(2*H), 1e-6);
    // lmtd equal branch: a==b -> value a, grad 1/2 each
    mc::F<double> c=a0,d=a0; c.diff(0,2); d.diff(1,2);
    CHK("lmtd eq value", lmtd(c,d).val(), a0, 1e-12);
    CHK("lmtd eq d/da",  lmtd(c,d).deriv(0), 0.5, 1e-12);
    CHK("lmtd eq d/db",  lmtd(c,d).deriv(1), 0.5, 1e-12);
    auto rl=[](double p,double q){return (std::log(p)-std::log(q))/(p-q);};
    CHK("rlmtd d/da FD", rlmtd(a,b).deriv(0), (rl(a0+H,b0)-rl(a0-H,b0))/(2*H), 1e-5);
    CHK("rlmtd eq value",rlmtd(c,d).val(), 1./a0, 1e-12);
    // cheb: T_n' via analytic (Chebyshev U)
    for(unsigned n=2;n<=8;n++){
      mc::F<double> z=0.4; z.diff(0,1);
      double dn=(std::cos(n*std::acos(0.4+H))-std::cos(n*std::acos(0.4-H)))/(2*H);
      CHK(N2("cheb'",std::to_string(n).c_str()), cheb(z,n).deriv(0), dn, 1e-4);
    }
    // arh(a,k)=exp(-k/a)
    mc::F<double> e=1.5; e.diff(0,1);
    CHK("arh value", arh(e,0.5).val(), std::exp(-0.5/1.5), 1e-12);
    CHK("arh deriv", arh(e,0.5).deriv(0), std::exp(-0.5/1.5)*0.5/(1.5*1.5), 1e-9);
  }

  // ======================= 5. DAG reuse & constants =======================
  sec("shared subexpressions & constant/mixed operands");
  {
    mc::F<double> x=0.8; x.diff(0,1);
    // f = x*x + sin(x) + x*exp(x)  (x reused 5x) : df = 2x+cos x+e^x+x e^x
    auto f=x*x+sin(x)+x*exp(x);
    CHK("reuse deriv", f.deriv(0), 2*0.8+std::cos(0.8)+std::exp(0.8)+0.8*std::exp(0.8), 1e-9);
    mc::F<double> k=3.0;                       // constant (unseeded)
    CHKb("const depend()==false", !k.depend());
    CHKb("x depend()==true", x.depend());
    auto g=x+k; CHK("x+const deriv", g.deriv(0), 1., 1e-12);
    auto h=k*x; CHK("const*x deriv", h.deriv(0), 3., 1e-12);
    CHKb("F<double> size after diff==1", x.size()==1);
  }

  // ======================= 6. VECTOR-VALUED JACOBIAN ======================
  sec("vector-valued system R^4->R^5: Jacobian F / SF / fadbad::F / finite-diff");
  {
    const double X0[NIN]={0.7,-0.4,1.3,0.9};
    double Jf[MOUT][NIN], Js[MOUT][NIN], Jb[MOUT][NIN], Jd[MOUT][NIN];
    // dense forward mc::F
    { mc::F<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++){v[i]=X0[i];v[i].diff(i,NIN);} SYS(v,o);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Jf[k][i]=o[k].deriv(i); }
    // sparse forward mc::SF
    { mc::SF<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++){v[i]=X0[i];v[i].diff(i,NIN);} SYS(v,o);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Js[k][i]=o[k].deriv(i); }
    // fadbad::F reference
    { fadbad::F<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++){v[i]=X0[i];v[i].diff(i,NIN);} SYS(v,o);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Jb[k][i]=o[k].deriv(i); }
    // central finite differences
    for(int i=0;i<NIN;i++){ double xp[NIN],xm[NIN]; for(int j=0;j<NIN;j++)xp[j]=xm[j]=X0[j];
      xp[i]+=H; xm[i]-=H; double op[MOUT],om[MOUT]; SYS(xp,op); SYS(xm,om);
      for(int k=0;k<MOUT;k++) Jd[k][i]=(op[k]-om[k])/(2*H); }
    for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++){
      std::string t="J["+std::to_string(k)+"]["+std::to_string(i)+"]";
      CHK(t+" F vs fadbad",  Jf[k][i], Jb[k][i], 1e-12);
      CHK(t+" SF vs F",      Js[k][i], Jf[k][i], 1e-12);
      CHK(t+" F vs findiff", Jf[k][i], Jd[k][i], 5e-5);
    }
    // sparse pattern: f_4 = x0/(x1+2) + x2*cos(x3) depends on {0,1,2,3} but not
    // via x-independent terms; check nnz of each output row equals #nonzero cols
    mc::SF<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++){v[i]=X0[i];v[i].diff(i,NIN);} SYS(v,o);
    int nz3=0; for(int i=0;i<NIN;i++) if(std::fabs(Jf[3][i])>1e-14) nz3++;
    CHKb("SF row3 nnz matches dense pattern", (int)o[3].size()==nz3);
    CHKb("SF dimension()==NIN", o[0].dimension()==(unsigned)NIN);
  }

  // ======================= 7. higher-order F<F> Hessian ===================
  sec("higher-order F<F>: Hessian of x^2 y + sin(x)");
  {
    typedef mc::F<double> FD; typedef mc::F<FD> FFD;
    const double x0=0.7,y0=1.3;
    FFD x; x.x()=FD(x0); x.x().diff(0,2); x.diff(0,2);
    FFD y; y.x()=FD(y0); y.x().diff(1,2); y.diff(1,2);
    FFD f=sqr(x)*y+sin(x);
    CHK("d/dx",   f.deriv(0).val(), 2*x0*y0+std::cos(x0), 1e-10);
    CHK("d/dy",   f.deriv(1).val(), x0*x0, 1e-10);
    CHK("d2/dx2", f.deriv(0).deriv(0), 2*y0-std::sin(x0), 1e-10);
    CHK("d2/dxdy",f.deriv(0).deriv(1), 2*x0, 1e-10);
    CHK("d2/dydx",f.deriv(1).deriv(0), 2*x0, 1e-10);
    CHK("symmetry",f.deriv(0).deriv(1), f.deriv(1).deriv(0), 1e-14);
  }

  // ======================= 8. interval soundness ==========================
  sec("F<Interval> soundness: derivative enclosure contains truth");
  {
    typedef mc::F<I> FI;
    const double x0=0.6;
    FI x=I(x0-0.05,x0+0.05); x.diff(0,1);
    FI f=exp(x)*sin(x)+atan(x);
    I d=f.deriv(0);
    double td=(std::exp(x0+H)*std::sin(x0+H)+std::atan(x0+H)-(std::exp(x0-H)*std::sin(x0-H)+std::atan(x0-H)))/(2*H);
    CHKb("deriv encloses point", d.l()<=td+1e-9 && td<=d.u()+1e-9);
  }

  // ======================= 9. McCormick soundness =========================
  sec("F<McCormick> soundness: relaxation & interval contain truth");
  {
    typedef mc::F<MC> FMC;
    MC xv(I(0.4,0.8),0.6); xv.sub(1,0);
    FMC x=xv; x.diff(0,1);
    FMC f=exp(x)*x;                 // d/dx = e^x (x+1)
    MC d=f.deriv(0);
    double td=std::exp(0.6)*(0.6+1.);
    CHKb("relax cv<=truth<=cc", d.cv()<=td+1e-9 && td<=d.cc()+1e-9);
    CHKb("interval contains truth", d.I().l()<=td+1e-9 && td<=d.I().u()+1e-9);
  }

  done();
  return g_fail?1:0;
}
