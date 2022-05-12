set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 65,330

splot 'INT-2D.out' u 1:2:3 tit 'function' w l lt 1, \
  '' u 1:2:4 tit 'lower bound' w l lt 3, \
  '' u 1:2:5 tit 'upper bound' w l lt 3

pause -1 "<ENTER> TO CONTINUE"

