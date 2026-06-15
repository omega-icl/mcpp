resfile  = "SCM-Bernstein.out"
hullfile = "SCM-Bernstein-hull.out"

set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 55,215

stats resfile u 1:2
xL = STATS_min_x
xU = STATS_max_x
yL = STATS_min_y
yU = STATS_max_y

set key below

splot resfile u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Chebyshev approximant' w l lt 2, \
  hullfile u (xL+$1*(xU-xL)):(yL+$2*(yU-yL)):3 tit 'Bernstein Hull' w p ps 2

pause -1 "<ENTER> TO CONTINUE"  

set term post eps enh solid color 18
set out 'SCM-Bernstein.eps'
rep
set term wxt

