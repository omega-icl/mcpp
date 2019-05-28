resfile1  = 'POLIMG-1D.out'
resfile2 = 'POLIMG-PM2-1D.out'
resfile3 = 'POLIMG-PM3-1D.out'
resfile4 = 'POLIMG-PM10-1D.out'

set term post eps enh solid color 12
set out 'POLIMG-1D.eps'

set multiplot layout 2,2 columnsfirst margins 0.1,0.9,0.1,0.9 spacing 0.1

set xlabel 'x'
set ylabel 'f(x)'
#unset key

set yrange [-2.2:2]
plot resfile1 u 1:5:6 tit '' w filledcurves lc 5 lt 1, \
  '' u 1:2 tit '' w l lt 1 lc -1 lw 3, \
  '' u 1:3 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:4 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:5 tit 'polyhedral relax.' w l lt 1 lc 7 lw 1, \
  '' u 1:6 tit '' w l lt 1 lc 7 lw 1

set yrange [-2.2:2]
plot resfile2 u 1:5:6 tit '' w filledcurves lc 5 lt 1, \
  '' u 1:2 tit '' w l lt 1 lc -1 lw 3, \
  '' u 1:3 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:4 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:5 tit '2nd-order Cheb. model ' w l lt 1 lc 7 lw 1, \
  '' u 1:6 tit '' w l lt 1 lc 7 lw 1

set yrange [-2.2:2]
plot resfile3 u 1:5:6 tit '' w filledcurves lc 5 lt 1, \
  '' u 1:2 tit '' w l lt 1 lc -1 lw 3, \
  '' u 1:3 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:4 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:5 tit '3rd-order Cheb. model ' w l lt 1 lc 7 lw 1, \
  '' u 1:6 tit '' w l lt 1 lc 7 lw 1

set yrange [-2.2:2]
plot resfile4 u 1:5:6 tit '' w filledcurves lc 5 lt 1, \
  '' u 1:2 tit '' w l lt 1 lc -1 lw 3, \
  '' u 1:3 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:4 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:5 tit '5th-order Cheb. model ' w l lt 1 lc 7 lw 1, \
  '' u 1:6 tit '' w l lt 1 lc 7 lw 1

unset multiplot
set term x11
set out
!ps2eps -B -f -l POLIMG-1D.eps
!mv POLIMG-1D.eps.eps POLIMG-1D.eps
!gv POLIMG-1D.eps
