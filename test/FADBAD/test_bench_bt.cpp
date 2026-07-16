#include <cstdio>
#include <chrono>
#include BDIFF_HDR
#include TDIFF_HDR
using clk=std::chrono::high_resolution_clock;
template <typename X> X f(const X* v, int N){    // moderate expression over N inputs
  X s = v[0];
  for(int i=1;i<N;i++) s = s + v[i]*v[i-1];
  return exp(s/(double)N) + sin(s) - sqrt(s*s+1.0);
}
int main(){
  const int N=64, IT=20000;
  // ---- reverse-mode B: single-output gradient ----
  { volatile double sink=0; auto t0=clk::now();
    for(int it=0; it<IT; ++it){
      mc::B<double> v[64];
      for(int i=0;i<N;i++) v[i]=0.5+0.001*i+1e-6*it;
      mc::B<double> y=f(v,N); y.diff(0,1);
      for(int i=0;i<N;i++) sink+=v[i].d(0);
    }
    auto t1=clk::now(); double ns=std::chrono::duration<double,std::nano>(t1-t0).count()/IT;
    std::printf("B<double> N=%d grad(1 output): %.0f ns/iter  (sink=%.3e)\n",N,ns,(double)sink);
  }
  // ---- Taylor T to order P ----
  { const int P=12; volatile double sink=0; auto t0=clk::now();
    for(int it=0; it<IT; ++it){
      mc::T<double> x=0.3+1e-7*it; x[1]=1.;
      mc::T<double> y=exp(sin(x))+sqrt(x+2.)-log(x+2.)+atan(x)+tanh(x)+x*x*x;
      y.eval(P); for(int k=0;k<=P;k++) sink+=y[k];
    }
    auto t1=clk::now(); double ns=std::chrono::duration<double,std::nano>(t1-t0).count()/IT;
    std::printf("T<double> Taylor order %d: %.0f ns/iter  (sink=%.3e)\n",P,ns,(double)sink);
  }
  return 0;
}
