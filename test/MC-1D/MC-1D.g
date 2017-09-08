set xlabel 'x'
set key below

plot 'MC-1D.out' u 1:2 tit 'function' w l lt 1
 
pause -1 "<ENTER> TO CONTINUE"

plot 'MC-1D.out' u 1:2 tit 'function' w l lt 1, \
  '' u 1:3 tit 'bounds' w l lt 2, \
  '' u 1:4 tit '' w l lt 2, \
  '' u 1:5 tit 'relaxations' w l lt 3, \
  '' u 1:6 tit '' w l lt 3
 
pause -1 "<ENTER> TO CONTINUE"

plot 'MC-1D.out' u 1:2 tit 'function' w l lt 1, \
  '' u 1:5 tit 'relaxations' w l lt 3, \
  '' u 1:6 tit '' w l lt 3, \
  '' u 1:9 tit 'linearizations' w l lt -1, \
  '' u 1:10 tit '' w l lt -1
  
pause -1 "<ENTER> TO CONTINUE"

plot 'MC-1D.out' u 1:7 tit 'subgrad. for underest.' w l lt 1, \
  '' u 1:8 tit 'subgrad. for overest.' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"
