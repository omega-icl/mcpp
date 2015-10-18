set xlabel 'x'
set key below

plot 'TM-1D.out' u 1:2 tit 'function' w l lc 1 lt 1, \
  '' u 1:3 tit 'Taylor model' w l lt 1 lc 3

pause -1 "<ENTER> TO CONTINUE"

plot 'TM-1D.out' u 1:2 tit 'function' w l lc 1 lt 1, \
  '' u 1:4 tit 'Taylor model' w l lt 1 lc 3, \
  '' u 1:5 tit '' w l lt 1 lc 3, \
  '' u 1:3 tit '' w l lt 0 lc 3, \
  '' u 1:8 tit 'Taylor-derived bounds' w l lt 1 lc 2, \
  '' u 1:9 tit '' w l lt 1 lc 2

set term post eps enh solid color 18
set out 'TM-1D.eps'
rep
set term x11
set out

pause -1 "<ENTER> TO CONTINUE"

plot 'TM-1D.out' u 1:2 tit 'function' w l lc 1 lt 1, \
  '' u 1:6 tit 'McCormick-Taylor model' w l lt 1 lc 3, \
  '' u 1:7 tit '' w l lt 1 lc 3, \
  '' u 1:3 tit '' w l lt 0 lc 3, \
  '' u 1:10 tit 'McCormick-Taylor-derived bounds' w l lt 1 lc 2, \
  '' u 1:11 tit '' w l lt 1 lc 2

pause -1 "<ENTER> TO CONTINUE"  

plot 'TM-1D.out' u 1:12 tit 'subgrad. for remainder underest.' w l lt 1, \
  '' u 1:13 tit 'subgrad. for remainder overest.' w l lt 3

pause -1 "<ENTER> TO CONTINUE"
