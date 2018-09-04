resfile  = 'LINREL_1d.out'

set xlabel 'x'
set key below

plot resfile u 1:5:6 tit '' w filledcurves lc 5 lt 1, \
  '' u 1:2 tit 'function' w l lt 1 lc -1 lw 3, \
  '' u 1:3 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:4 tit '' w l lt 3 lc 7 lw 3, \
  '' u 1:5 tit 'relaxations' w l lt 1 lc 7 lw 1, \
  '' u 1:6 tit '' w l lt 1 lc 7 lw 1
 
pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 21
set out 'LINREL_1d.eps'
rep


