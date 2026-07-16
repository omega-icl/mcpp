// Which algebraic form of the derivative coefficient propagates tightest?
// Compares candidate forms for -a/b^2 (division) across several boxes/points,
// under both interval and McCormick arithmetic.
#include <cstdio>
#include <cmath>
#include "mccormick.hpp"
#include "interval.hpp"

typedef mc::Interval     I;
typedef mc::McCormick<I> MC;

int main()
{
  struct Box { double bl, bu, bp, a; };
  Box boxes[] = {
    {0.5, 3.0, 1.5, 2.0}, {0.5, 3.0, 2.5, 2.0}, {1.0, 5.0, 2.0, 1.0},
    {0.2, 1.0, 0.6, 3.0}, {2.0, 8.0, 5.0,-1.5}, {0.5, 4.0, 1.0, 2.0},
  };
  std::printf("division coeff -a/b^2 : interval width  |  McCormick gap (cc-cv)\n");
  std::printf("box[bl,bu]@bp a      (a/b)/b   a/sqr(b) | (a/b)/b   a/sqr(b)   a*inv(sqr(b))\n");
  for( auto& B : boxes ){
    // interval widths
    I bi(B.bl,B.bu);
    I i_old = -( (B.a/bi)/bi );
    I i_new = -( B.a/mc::Op<I>::sqr(bi) );
    double wI_old = mc::Op<I>::diam(i_old), wI_new = mc::Op<I>::diam(i_new);
    // McCormick gaps
    MC bm( I(B.bl,B.bu), B.bp ); bm.sub(1,0);
    MC m_old = -( (B.a/bm)/bm );
    MC m_new = -( B.a/mc::Op<MC>::sqr(bm) );
    MC m_alt = -( B.a*mc::Op<MC>::inv(mc::Op<MC>::sqr(bm)) );
    double gO = m_old.cc()-m_old.cv(), gN = m_new.cc()-m_new.cv(), gA = m_alt.cc()-m_alt.cv();
    std::printf("[%.1f,%.1f]@%.1f a=%+.1f  %7.4f  %7.4f | %7.4f  %7.4f  %7.4f  %s\n",
                B.bl,B.bu,B.bp,B.a, wI_old,wI_new, gO,gN,gA,
                gO<=gN+1e-9? "old<=new(MC)":"new<old(MC)");
  }
  return 0;
}
