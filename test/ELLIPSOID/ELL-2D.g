set xlabel "f_1"
set ylabel "f_2"

file0 = 'Fsam.out'
file1 = 'Esam.out'
file2 = 'Isam.out'

set grid front
set key below 
set size ratio 1

set autoscale yfix

plot  file0  tit ""          w d lt 1 lc rgb "#FFA346", \
      file1  tit ""   w l lt 1 lc 3 lw 3, \
      file2  tit ""   w l lt 1 lc 1 lw 3

pause -1 "PRESS <ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'ELL-2D.eps'
rep
set term x11
set out
!ps2eps -B -l -f ELL-2D.eps
!mv ELL-2D.eps.eps ELL-2D.eps

