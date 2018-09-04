resfile  = 'LINREL_2d.out'

set xlabel 'x'
set xlabel 'y'
set key below
unset hidden3d

splot resfile u 1:2:3 tit 'function' w l lt 1 lc -1 lw 3, \
  '' u 1:2:6 tit 'relaxations' w l lt 7 lc 1 lw 1, \
  '' u 1:2:7 tit '' w l lt 1 lc 7 lw 1
#  '' u 1:2:4 tit '' w l lt 3 lc 7 lw 3, \
#  '' u 1:2:5 tit '' w l lt 3 lc 7 lw 3, \
 
pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 21
set out 'LINREL_2d.eps'
rep

