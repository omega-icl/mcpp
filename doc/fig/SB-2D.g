set xlabel 'x'
set xrange [-0.5:0.5]
set ylabel 'y'
set yrange [-0.5:0.5]
set key below
set hidden3d
set view 65,30

splot 'SB-2D.out' u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:5 tit 'spectral bounds' w l lt 1 lc 3, \
  '' u 1:2:6 tit '' w l lt 1 lc 3

set terminal pngcairo size 350,262 enhanced font 'Verdana,7'
set out "SB-2D_spectral.png"
rep
set term wxt

pause -1 "<ENTER> TO CONTINUE"


splot 'SB-2D.out' u 1:2:7 tit 'function' w l lt 1 lc 1, \
  '' u 1:2:8 tit 'function bounds' w l lt 1 lc 3, \
  '' u 1:2:9 tit '' w l lt 1 lc 3

set terminal pngcairo size 350,262 enhanced font 'Verdana,7'
set out "SB-2D_function.png"
rep
set term wxt

pause -1 "<ENTER> TO CONTINUE"
