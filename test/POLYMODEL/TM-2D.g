set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 74,13

set key below

splot 'TM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Taylor model' w l lt 2, \
  '' u 1:2:5 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

splot 'TM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:6 tit 'McCormick-Taylor model' w l lt 2, \
  '' u 1:2:7 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

splot 'TM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:8 tit 'Taylor-based interval bounds' w l lt 2, \
  '' u 1:2:9 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

splot 'TM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:10 tit 'Taylor-based McCormick relaxation' w l lt 2, \
  '' u 1:2:11 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

splot 'TM-2D.out' u 1:2:12 tit 'subgrad. along (1,0) for remainder underest.' w l lt 1, \
  '' u 1:2:13 tit 'subgrad. along (1,0) for remainder overest.' w l lt 3

pause -1 "<ENTER> TO CONTINUE"  

splot 'TM-2D.out' u 1:2:14 tit 'subgrad. along (0,1) for remainder underest.' w l lt 1, \
  '' u 1:2:15 tit 'subgrad. along (0,1) for remainder overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"
