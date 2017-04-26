#set size ratio 0.7
#set key outside rm t Left reverse spacing 2
set key below

set xlabel 'x'
set xrange [-0.5:1]
set ylabel 'f(x)'
set yrange [-1.5:1.5]

### TAYLOR MODEL

plot 'TM-MC-1D.out' u 1:2 tit 'function' w l lt 1 lc 1 lw 2, \
  '' u 1:7 tit 'Taylor model' w l lt 3 lc 3 lw 2, \
  '' u 1:8 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:11 tit '' w l lt 1 lc 2 lw 2, \
  '' u 1:12 tit '' w l lt 1 lc 2 lw 2

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set out 'TM-1D.png'
rep
set term wxt

### MCCORMICK-TAYLOR MODEL

plot 'TM-MC-1D.out' u 1:2 tit '' w l lt 1 lc 1 lw 2, \
  '' u 1:9 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:10 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:13 tit 'Taylor model-based bounds' w l lt 1 lc 2 lw 2, \
  '' u 1:14 tit '' w l lt 1 lc 2 lw 2

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set out 'MCTM-1D.png'
rep
set term wxt

set terminal pngcairo size 540,189 enhanced font 'Verdana,7'
set out "TM+MCTM-1D.png"
set multiplot layout 1,2
#
plot 'TM-MC-1D.out' u 1:2 tit 'function' w l lt 1 lc 1 lw 2, \
  '' u 1:7 tit 'Taylor model' w l lt 3 lc 3 lw 2, \
  '' u 1:8 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:11 tit '' w l lt 1 lc 2 lw 2, \
  '' u 1:12 tit '' w l lt 1 lc 2 lw 2
#
plot 'TM-MC-1D.out' u 1:2 tit '' w l lt 1 lc 1 lw 2, \
  '' u 1:9 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:10 tit '' w l lt 3 lc 3 lw 2, \
  '' u 1:13 tit 'Taylor model-based bounds' w l lt 1 lc 2 lw 2, \
  '' u 1:14 tit '' w l lt 1 lc 2 lw 2
#
unset multiplot
pause -1 "<ENTER> TO CONTINUE"
