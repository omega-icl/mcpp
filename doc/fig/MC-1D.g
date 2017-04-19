set xlabel 'x'
set ylabel 'f(x)'
set key below

plot 'MC-1D.out' u 1:2 tit 'function' w l lt 1 lw 2, \
  '' u 1:3 tit 'interval bounds' w l lt 2 lw 2, \
  '' u 1:4 tit '' w l lt 2 lw 2, \
  '' u 1:5 tit '' w l lt 3 lw 2, \
  '' u 1:6 tit '' w l lt 3 lw 2

set terminal pngcairo size 253,189 enhanced font 'Verdana,7'
#set term png small enhanced
set out "MC-1D_relax.png"
rep

pause -1 "<ENTER> TO CONTINUE"

plot 'MC-1D.out' u 1:2 tit '' w l lt 1 lw 2, \
  '' u 1:5 tit 'convex/concave relaxations' w l lt 3 lw 2, \
  '' u 1:6 tit '' w l lt 3 lw 2, \
  '' u 1:9 tit 'affine relaxations' w l lt -1 lw 2, \
  '' u 1:10 tit '' w l lt -1 lw 2
  
set terminal pngcairo size 253,189 enhanced font 'Verdana,7'
#set term png small enhanced
set out "MC-1D_linearize.png"
rep

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 540,189 enhanced font 'Verdana,7'
set out "MC-1D.png"
set multiplot layout 1,2
#
plot 'MC-1D.out' u 1:2 tit 'function' w l lt 1 lw 2, \
  '' u 1:3 tit 'interval bounds' w l lt 2 lw 2, \
  '' u 1:4 tit '' w l lt 2 lw 2, \
  '' u 1:5 tit '' w l lt 3 lw 2, \
  '' u 1:6 tit '' w l lt 3 lw 2
#
plot 'MC-1D.out' u 1:2 tit '' w l lt 1 lw 2, \
  '' u 1:5 tit 'convex/concave relaxations' w l lt 3 lw 2, \
  '' u 1:6 tit '' w l lt 3 lw 2, \
  '' u 1:9 tit 'affine relaxations' w l lt -1 lw 2, \
  '' u 1:10 tit '' w l lt -1 lw 2
#
unset multiplot
pause -1 "<ENTER> TO CONTINUE"
