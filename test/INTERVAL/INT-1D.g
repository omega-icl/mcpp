set xlabel 'x'
set ylabel 'f'
set key below

plot 'INT-1D.out' u 1:2 tit 'function' w l lt 1, \
  '' u 1:3 tit 'bounds' w l lt 2, \
  '' u 1:4 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"
   
set term post eps enh solid color 18
set out 'INT-1D.eps'
rep
set term pop
set out

