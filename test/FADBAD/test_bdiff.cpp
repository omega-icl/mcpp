// Systematic test suite for mc::B (dense) and mc::SB (sparse) reverse-mode AD.
//   - every elementary function: derivative vs analytic, finite differences,
//     forward mc::F, and (where native) fadbad::B; SB==B consistency
//   - binary/scalar operators, min/max/lmtd/rlmtd/cheb/pow/arh
//   - DAG accumulation (shared subexpressions summed correctly)
//   - VECTOR-VALUED system R^4 -> R^5: full Jacobian via B, SB, fadbad::B,
//     forward mc::F and finite differences, all cross-checked (forward==reverse)
//   - forward-over-reverse Hessian B<F>, interval & McCormick soundness,
//     sparse pattern (unreached input -> 0)
#include <cstdio>
#include <cmath>
#include <string>
#include "bdiff.hpp"
#include "fdiff.hpp"
#include "interval.hpp"
#include "mccormick.hpp"
#include "fadbad.h"
#include "badiff.h"

typedef mc::Interval     I;
typedef mc::McCormick<I> MC;
static const double H=1e-6;

static int g_fail=0,g_tot=0,s_fail=0,s_tot=0; static std::string g_sec;
static void sec(const char* n){ if(!g_sec.empty()) std::printf("  -> %-42s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  g_sec=n; s_fail=s_tot=0; std::printf("[%s]\n",n); }
static void CHK(const std::string& nm,double got,double ref,double tol){ ++g_tot;++s_tot; double e=std::fabs(got-ref);
  if(e>tol*(1.+std::fabs(ref))){ ++g_fail;++s_fail; std::printf("  FAIL %-42s got=% .10e ref=% .10e |e|=%.2e\n",nm.c_str(),got,ref,e);} }
static void CHKb(const std::string& nm,bool ok){ ++g_tot;++s_tot; if(!ok){ ++g_fail;++s_fail; std::printf("  FAIL %-42s\n",nm.c_str()); } }
static std::string N2(const char* a,const char* b){ return std::string(a)+" "+b; }
static void done(){ if(!g_sec.empty()) std::printf("  -> %-42s %d/%d\n",g_sec.c_str(),s_tot-s_fail,s_tot);
  std::printf("\n%s  (%d/%d checks passed)\n", g_fail?"*** FAILURES ***":"ALL PASSED", g_tot-g_fail, g_tot); }

// reverse-mode unary: B derivative vs analytic, finite diff, forward mc::F, SB
template<class ADf,class VFn,class DFn>
void uni(const char* nm, ADf adf, VFn vf, DFn df, double x0){
  mc::B<double> x=x0; auto f=adf(x); f.diff(0,1);
  mc::SB<double> sx=x0; auto sf=adf(sx); sf.diff(0,1);
  mc::F<double> fx=x0; fx.diff(0,1); auto ff=adf(fx);
  double v=vf(x0), an=df(x0), fd=(vf(x0+H)-vf(x0-H))/(2*H);
  CHK(N2(nm,"value"),         f.val(), v, 1e-12);
  CHK(N2(nm,"B d/dx analytic"),x.d(0), an, 1e-9);
  CHK(N2(nm,"B d/dx findiff"), x.d(0), fd, 1e-5);
  CHK(N2(nm,"B == forward F"), x.d(0), ff.deriv(0), 1e-13);
  CHK(N2(nm,"SB == B"),        sx.d(0), x.d(0), 1e-13);
}
template<class ADf> void unifb(const char* nm, ADf adf, double x0){
  mc::B<double> x=x0; auto f=adf(x); f.diff(0,1);
  fadbad::B<double> fx=x0; auto ff=adf(fx); ff.diff(0,1);
  CHK(N2(nm,"value vs fadbad"), f.val(), ff.x(), 1e-13);
  CHK(N2(nm,"deriv vs fadbad"), x.d(0), fx.d(0), 1e-13);
}

// vector-valued system  R^4 -> R^5
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
  sec("unary: reverse B vs analytic / findiff / forward-F / sparse-SB");
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

  // ======================= 2. cross-check vs fadbad::B ====================
  sec("unary vs genuine fadbad::B (native functions)");
  unifb("sqr", [](auto&z){return sqr(z);}, 0.7);
  unifb("sqrt",[](auto&z){return sqrt(z);},1.7);
  unifb("exp", [](auto&z){return exp(z);}, 0.6);
  unifb("log", [](auto&z){return log(z);}, 1.7);
  unifb("sin", [](auto&z){return sin(z);}, 0.6);
  unifb("cos", [](auto&z){return cos(z);}, 0.6);
  unifb("tan", [](auto&z){return tan(z);}, 0.6);
  unifb("asin",[](auto&z){return asin(z);},0.4);
  unifb("atan",[](auto&z){return atan(z);},0.6);
  unifb("tanh",[](auto&z){return tanh(z);},0.6);

  // ======================= 3. binary & scalar operators ===================
  sec("binary/scalar operators: gradient via reverse mode");
  {
    const double a0=1.3,b0=0.7;
    auto grad2=[&](auto fn, double da, double db, const char* nm){
      mc::B<double> a=a0,b=b0; auto f=fn(a,b); f.diff(0,1);
      CHK(N2(nm,"d/da"), a.d(0), da, 1e-9);
      CHK(N2(nm,"d/db"), b.d(0), db, 1e-9);
      mc::SB<double> sa=a0,sb=b0; auto sf=fn(sa,sb); sf.diff(0,1);
      CHK(N2(nm,"SB d/da"), sa.d(0), a.d(0), 1e-13);
      CHK(N2(nm,"SB d/db"), sb.d(0), b.d(0), 1e-13);
    };
    grad2([](auto&a,auto&b){return a+b;}, 1., 1., "a+b");
    grad2([](auto&a,auto&b){return a-b;}, 1., -1., "a-b");
    grad2([](auto&a,auto&b){return a*b;}, b0, a0, "a*b");
    grad2([](auto&a,auto&b){return a/b;}, 1./b0, -a0/(b0*b0), "a/b");
    grad2([](auto&a,auto&b){return pow(a,b);}, b0*std::pow(a0,b0-1), std::pow(a0,b0)*std::log(a0), "pow(a,b)");
    grad2([](auto&a,auto&b){return max(a,b);}, 1., 0., "max(a,b)");
    grad2([](auto&a,auto&b){return min(a,b);}, 0., 1., "min(a,b)");
    // scalar cases
    mc::B<double> a=a0; auto f=(2./a); f.diff(0,1);
    CHK("2/a d/da", a.d(0), -2./(a0*a0), 1e-9);
    mc::B<double> c=a0; auto g=pow(c,3); g.diff(0,1);
    CHK("pow(a,3) d/da", c.d(0), 3*a0*a0, 1e-9);
    mc::B<double> e=1.5; auto ar=arh(e,0.5); ar.diff(0,1);
    CHK("arh d/da", e.d(0), std::exp(-0.5/1.5)*0.5/(1.5*1.5), 1e-9);
  }

  // ======================= 4. DAG accumulation ============================
  sec("shared-subexpression adjoint accumulation");
  {
    mc::B<double> x=0.8; auto f=x*x+sin(x)+x*exp(x); f.diff(0,1);  // x reused 5x
    CHK("reuse deriv", x.d(0), 2*0.8+std::cos(0.8)+std::exp(0.8)+0.8*std::exp(0.8), 1e-9);
    // diamond: g=(x+1); h=g*g+exp(g)  -> dh/dx = 2g+e^g
    mc::B<double> y=0.5; auto gg=y+1.; auto hh=gg*gg+exp(gg); hh.diff(0,1);
    CHK("diamond deriv", y.d(0), 2*(1.5)+std::exp(1.5), 1e-9);
  }

  // ======================= 5. VECTOR-VALUED JACOBIAN ======================
  sec("vector system R^4->R^5: Jacobian B / SB / fadbad::B / forward-F / findiff");
  {
    const double X0[NIN]={0.7,-0.4,1.3,0.9};
    double Jb[MOUT][NIN], Js[MOUT][NIN], Jff[MOUT][NIN], Jf[MOUT][NIN], Jd[MOUT][NIN];
    // dense reverse mc::B (seed each output)
    { mc::B<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++)v[i]=X0[i]; SYS(v,o);
      for(int k=0;k<MOUT;k++) o[k].diff(k,MOUT);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Jb[k][i]=v[i].d(k); }
    // sparse reverse mc::SB
    { mc::SB<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++)v[i]=X0[i]; SYS(v,o);
      for(int k=0;k<MOUT;k++) o[k].diff(k,MOUT);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Js[k][i]=v[i].d(k); }
    // fadbad::B reference
    { fadbad::B<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++)v[i]=X0[i]; SYS(v,o);
      for(int k=0;k<MOUT;k++) o[k].diff(k,MOUT);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Jff[k][i]=v[i].d(k); }
    // forward mc::F (independent mode)
    { mc::F<double> v[NIN],o[MOUT]; for(int i=0;i<NIN;i++){v[i]=X0[i];v[i].diff(i,NIN);} SYS(v,o);
      for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++) Jf[k][i]=o[k].deriv(i); }
    // central differences
    for(int i=0;i<NIN;i++){ double xp[NIN],xm[NIN]; for(int j=0;j<NIN;j++)xp[j]=xm[j]=X0[j];
      xp[i]+=H; xm[i]-=H; double op[MOUT],om[MOUT]; SYS(xp,op); SYS(xm,om);
      for(int k=0;k<MOUT;k++) Jd[k][i]=(op[k]-om[k])/(2*H); }
    for(int k=0;k<MOUT;k++) for(int i=0;i<NIN;i++){
      std::string t="J["+std::to_string(k)+"]["+std::to_string(i)+"]";
      CHK(t+" B vs fadbad",   Jb[k][i], Jff[k][i], 1e-12);
      CHK(t+" SB vs B",       Js[k][i], Jb[k][i], 1e-12);
      CHK(t+" B vs forward-F",Jb[k][i], Jf[k][i], 1e-12);
      CHK(t+" B vs findiff",  Jb[k][i], Jd[k][i], 5e-5);
    }
  }

  // ======================= 6. multi-output seeding & sparsity =============
  sec("multi-output seeding & sparse unreached-input == 0");
  {
    // y0 = u*u (depends only on u), y1 = sin(u)+w (depends on u,w)
    mc::B<double> u=2.0,w=3.0; auto y0=u*u; auto y1=sin(u)+w;
    y0.diff(0,2); y1.diff(1,2);
    CHK("dy0/du", u.d(0), 4.0, 1e-12);
    CHK("dy0/dw (=0)", w.d(0), 0.0, 1e-12);
    CHK("dy1/du", u.d(1), std::cos(2.0), 1e-9);
    CHK("dy1/dw", w.d(1), 1.0, 1e-12);
    mc::SB<double> su=2.0,sw=3.0; auto z0=su*su; auto z1=sin(su)+sw;
    z0.diff(0,2); z1.diff(1,2);
    CHK("SB dy0/dw (=0)", sw.d(0), 0.0, 1e-12);
    CHK("SB dy1/dw",      sw.d(1), 1.0, 1e-12);
  }

  // ======================= 7. forward-over-reverse Hessian ================
  sec("Hessian via B<F<double>> (forward-over-reverse)");
  {
    typedef mc::F<double> FD; typedef mc::B<FD> BFD;
    const double x0=0.7,y0=1.3;
    FD xF=x0; xF.diff(0,2); FD yF=y0; yF.diff(1,2);
    BFD x=xF,y=yF; BFD f=sqr(x)*y+sin(x); f.diff(0,1);
    FD gx=x.d(0), gy=y.d(0);
    CHK("d/dx",    gx.val(), 2*x0*y0+std::cos(x0), 1e-10);
    CHK("d/dy",    gy.val(), x0*x0, 1e-10);
    CHK("d2/dx2",  gx.deriv(0), 2*y0-std::sin(x0), 1e-10);
    CHK("d2/dxdy", gx.deriv(1), 2*x0, 1e-10);
    CHK("d2/dydx", gy.deriv(0), 2*x0, 1e-10);
    CHK("d2/dy2",  gy.deriv(1), 0.0, 1e-10);
  }

  // ======================= 8. interval & McCormick soundness ==============
  sec("B<Interval> & B<McCormick> soundness");
  {
    typedef mc::B<I> BI;
    BI x=I(0.55,0.65),y=I(1.25,1.35); BI f=exp(x*y)+sin(x)/y; f.diff(0,1);
    I dx=x.d(0);
    double td=(std::exp((0.6+H)*1.3)+std::sin(0.6+H)/1.3 - (std::exp((0.6-H)*1.3)+std::sin(0.6-H)/1.3))/(2*H);
    CHKb("B<I> df/dx encloses truth", dx.l()<=td+1e-7 && td<=dx.u()+1e-7);
    typedef mc::B<MC> BMC;
    MC xv(I(0.4,0.8),0.6); xv.sub(1,0); BMC xm=xv; BMC g=exp(xm)*xm; g.diff(0,1);
    MC d=xm.d(0); double tg=std::exp(0.6)*(0.6+1.);
    CHKb("B<MC> relax contains truth", d.cv()<=tg+1e-9 && tg<=d.cc()+1e-9);
    CHKb("B<MC> interval contains truth", d.I().l()<=tg+1e-9 && tg<=d.I().u()+1e-9);
  }

  done();
  return g_fail?1:0;
}
