resfile = 'SB-2D.out'

set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 74,13

set key below

splot resfile u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:5 tit 'spectral bound' w l lt 1 lc 3, \
  '' u 1:2:6 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot resfile u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:5 tit 'spectral bound' w l lt 1 lc 3, \
  '' u 1:2:6 tit '' w l lt 1 lc 3, \
  '' u 1:2:7 tit 'Gershgorin' w l lt 1 lc 4, \
  '' u 1:2:8 tit '' w l lt 1 lc 4, \
  '' u 1:2:9 tit 'Hertz & Rohn' w l lt 1 lc 5, \
  '' u 1:2:10 tit '' w l lt 1 lc 5

pause -1 "<ENTER> TO CONTINUE"  

