set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 74,13

set key below

splot 'SB-2D.out' u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:5 tit 'spectral bound' w l lt 1 lc 3, \
  '' u 1:2:6 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:5 tit 'spectral bound' w l lt 1 lc 3, \
  '' u 1:2:6 tit '' w l lt 1 lc 3, \
  '' u 1:2:24 tit 'Gershgorin' w l lt 1 lc 4, \
  '' u 1:2:25 tit '' w l lt 1 lc 4, \
  '' u 1:2:26 tit 'Hertz & Rohn' w l lt 1 lc 5, \
  '' u 1:2:27 tit '' w l lt 1 lc 5

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:3 tit 'lambda min' w l lt 1 lc 1, \
  '' u 1:2:4 tit 'lambda max' w l lt 1 lc 2, \
  '' u 1:2:7 tit 'spectral relaxation' w l lt 1 lc 3, \
  '' u 1:2:8 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:9 tit 'function' w l lt 1 lc 1, \
  '' u 1:2:10 tit 'bounds' w l lt 1 lc 3, \
  '' u 1:2:11 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:9 tit 'function' w l lt 1 lc 1, \
  '' u 1:2:12 tit 'relaxations' w l lt 1 lc 3, \
  '' u 1:2:13 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:14 tit 'deriv w.r.t. x' w l lt 1 lc 1, \
  '' u 1:2:15 tit 'bounds' w l lt 1 lc 3, \
  '' u 1:2:16 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:14 tit 'deriv w.r.t. x' w l lt 1 lc 1, \
  '' u 1:2:17 tit 'relaxations' w l lt 1 lc 3, \
  '' u 1:2:18 tit '' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'SB-2D.out' u 1:2:19 tit 'deriv w.r.t. x' w l lt 1 lc 1, \
  '' u 1:2:20 tit 'bounds' w l lt 1 lc 3, \
  '' u 1:2:21 tit '' w l lt 1 lc 3
  
pause -1 "<ENTER> TO CONTINUE"

splot 'SB-2D.out' u 1:2:19 tit 'deriv w.r.t. x' w l lt 1 lc 1, \
  '' u 1:2:22 tit 'relaxations' w l lt 1 lc 3, \
  '' u 1:2:23 tit '' w l lt 1 lc 3
  
pause -1 "<ENTER> TO CONTINUE"
