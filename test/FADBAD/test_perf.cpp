// Correctness of rewritten paths (cheb recurrence, _chain2 mixed-depend) plus
// a microbenchmark of the hot dense/sparse paths.
#include <cstdio>
#include <cmath>
#include <chrono>
#include <vector>
#include "fdiff.hpp"

using clk = std::chrono::high_resolution_clock;
static double sec( clk::time_point a, clk::time_point b )
{ return std::chrono::duration<double>(b-a).count(); }

int fails = 0;
static void check( const char* nm, double got, double ref, double tol )
{
  double e = std::fabs(got-ref);
  bool ok = e <= tol*(1.+std::fabs(ref));
  if(!ok) ++fails;
  std::printf("  %-28s got=% .10e ref=% .10e |e|=%.2e %s\n", nm, got, ref, e, ok?"ok":"FAIL");
}

int main()
{
  typedef mc::F<double>  FD;
  typedef mc::SF<double> SFD;

  // ---- cheb derivative across degrees vs central finite differences --------
  std::printf("[cheb T_n'(x) recurrence vs finite differences]\n");
  const double x0 = 0.37, h = 1e-6;
  for( unsigned n=0; n<=12; ++n ){
    FD X=x0; X.diff(0,1);
    FD Y=mc::cheb(X,n);
    double adx = Y.deriv(0);
    double fd  = ( mc::Op<double>::cheb(x0+h,n) - mc::Op<double>::cheb(x0-h,n) )/(2*h);
    char nm[32]; std::snprintf(nm,sizeof nm,"d/dx cheb(x,%u)",n);
    check(nm, adx, fd, 1e-4);
  }

  // ---- mixed-depend min/max/pow (constant F operand): must not throw -------
  std::printf("[min/max/pow with a constant-F operand]\n");
  {
    FD X=1.3; X.diff(0,2);          // depends on seed 0
    FD K=2.0;                        // constant F (no gradient)
    FD m = mc::max(X,K);            // X<K? -> value K, dX=0 ... X>K? value X, dX=1; here X<K
    check("max(x,const) value", m.val(), 2.0, 0);
    check("max(x,const) d/dx",  m.deriv(0), 0.0, 0);
    FD p = mc::pow(X,K);           // x^2 -> d/dx = 2x
    check("pow(x,const) value", p.val(), std::pow(1.3,2.0), 1e-12);
    check("pow(x,const) d/dx",  p.deriv(0), 2*1.3, 1e-12);
  }

  // ---- microbenchmark ------------------------------------------------------
  std::printf("[microbenchmark]\n");
  const int    N    = 64;     // gradient dimension
  const long   iter = 20000; // repetitions

  // dense: build N-dim seeds, run a chain of unary+binary ops
  std::vector<FD> xs(N);
  for(int i=0;i<N;++i){ xs[i]=0.5+0.001*i; xs[i].diff(i,N); }
  volatile double sink=0;
  auto t0=clk::now();
  for(long it=0; it<iter; ++it){
    FD acc = xs[it%N];
    acc = mc::exp(acc) + mc::sin(acc)*xs[(it+1)%N] - mc::sqrt(acc*acc+1.);
    sink += acc.val() + acc.deriv(0);
  }
  auto t1=clk::now();
  std::printf("  dense  F<double> N=%d : %.3f s  (%.0f ns/iter)\n",
              N, sec(t0,t1), sec(t0,t1)/iter*1e9);

  // sparse: each intermediate touches few seeds
  std::vector<SFD> ss(N);
  for(int i=0;i<N;++i){ ss[i]=0.5+0.001*i; ss[i].diff(i,N); }
  auto t2=clk::now();
  for(long it=0; it<iter; ++it){
    SFD acc = ss[it%N];
    acc = mc::exp(acc) + mc::sin(acc)*ss[(it+1)%N] - mc::sqrt(acc*acc+1.);
    sink += acc.val() + acc.deriv(it%N);
  }
  auto t3=clk::now();
  std::printf("  sparse SF<double> N=%d : %.3f s  (%.0f ns/iter, ~2 seeds/term)\n",
              N, sec(t2,t3), sec(t2,t3)/iter*1e9);

  // cheb high-degree benchmark: exercises the recurrence
  auto t4=clk::now();
  for(long it=0; it<iter; ++it){
    FD X=0.3+0.0000001*(it%7); X.diff(0,1);
    FD Y=mc::cheb(X,16);
    sink += Y.deriv(0);
  }
  auto t5=clk::now();
  std::printf("  cheb(x,16) deriv       : %.3f s  (%.0f ns/iter)\n",
              sec(t4,t5), sec(t4,t5)/iter*1e9);

  std::printf("\n%s (%d failures) sink=%.3e\n", fails? "FAILURES":"ALL CHECKS PASSED", fails, (double)sink);
  return fails?1:0;
}
