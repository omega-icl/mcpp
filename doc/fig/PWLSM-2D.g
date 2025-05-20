set xlabel 'x'
set xrange [1:2]
set ylabel 'y'
set yrange [0:1]
set key above
set hidden3d
set view 65,10,1,1.4

splot 'PWLSM-2D.out' u 1:2:3 tit 'f' w l lt 1 lc 8, \
  '' u 1:2:4 tit 'f^u' w l lt 1 lc 4, \
  '' u 1:2:5 tit 'f^o' w l lt 1 lc 7

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 490,367 enhanced font 'Verdana,9'
set out "PWLSM-2D.png"
rep

