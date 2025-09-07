resfile  = 'FFDAG_2d.out'

set xlabel "x_1"
set ylabel "x_2"
set zlabel "f"
set key below
unset hidden3d

splot resfile u 1:2:3 tit 'function' w l lt 1 lc -1 lw 2, \
  '' u 1:2:6 tit 'underestimator' w l lt 7 lc 4 lw 1, \
  '' u 1:2:7 tit 'overestimator' w l lt 1 lc 7 lw 1
#  '' u 1:2:4 tit '' w l lt 3 lc 7 lw 3, \
#  '' u 1:2:5 tit '' w l lt 3 lc 7 lw 3, \
 
pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 15
set out 'FFDAG_2d.eps'
rep
set term qt

splot resfile u 1:2:($7-$3) tit 'gap' w l lt 1 lc 4 lw 1, \
  '' u 1:2:($6-$3) tit '' w l lt 7 lc 7 lw 1

pause -1 "<ENTER> TO CONTINUE"

