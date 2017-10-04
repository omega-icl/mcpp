set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 65,330

splot 'MC-2D.out' u 1:2:3 tit 'function' w l lt 1, \
  '' u 1:2:4 tit 'lower bound' w l lt 3, \
  '' u 1:2:5 tit 'upper bound' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"

splot 'MC-2D.out' u 1:2:3 tit 'function' w l lt 1, \
  '' u 1:2:6 tit 'convex underst.' w l lt 3, \
  '' u 1:2:7 tit 'concave overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"

splot 'MC-2D.out' u 1:2:8 tit 'subgrad. along (1,0) for underest.' w l lt 1, \
  '' u 1:2:9 tit 'subgrad. along (1,0) for overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"

splot 'MC-2D.out' u 1:2:10 tit 'subgrad. along (0,1) for underest.' w l lt 1, \
  '' u 1:2:11 tit 'subgrad. along (0,1) for overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"

splot 'MC-2D.out' u 1:2:8 tit 'directional subgrad. for underest.' w l lt 1, \
  '' u 1:2:9 tit 'directional subgrad. for overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"
