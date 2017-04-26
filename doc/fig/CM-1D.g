#set size ratio 0.7
#set key outside rm t Left reverse spacing 2
set key below

set xlabel 'x'
set xrange [-0.5:1]
set ylabel 'f(x)'
set yrange [-0.6:0.6]

### CHEBYSHEV MODEL

plot 'CM-1D.out' u 1:2 tit 'function' w l lt 1 lc rgb 'red' lw 2, \
  '' u 1:4 tit 'Chebyshev model' w l lt 3 lc rgb 'blue' lw 2, \
  '' u 1:5 tit '' w l lt 3 lc rgb 'blue' lw 2, \
  '' u 1:8 tit 'Chebyshev bounds' w l lt 1 lc rgb 'green' lw 2, \
  '' u 1:9 tit '' w l lt 1 lc rgb 'green' lw 2

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 350,262 enhanced font 'Verdana,10'
set out 'CM-1D.png'
rep
set term wxt




