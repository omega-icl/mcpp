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

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 490,367 enhanced font 'Verdana,9'
set out "SB-2D_spectral.png"
rep
set term qt

set isosample 40,40
#set palette cubehelix
#set pm3d interpolate 2,2
#unset colorbox
#unset key

set view 66, 241

splot 1.+x-sin(2.*x+3.*y)-cos(3.*x-5.*y) tit 'function' w l lt 1 lc 7

pause -1 "<ENTER> TO CONTINUE"

set terminal pngcairo size 490,367 enhanced font 'Verdana,9'
set out "SB-2D_function.png"
rep
set term qt

