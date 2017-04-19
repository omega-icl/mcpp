set xlabel "f_1"
set ylabel "f_2"
set yrange [5:35]
file0 = 'Fsam.out'
file1 = 'Esam.out'
file2 = 'Isam.out'
#file2 = 'Esam_DAG.out'

set grid front
set key below 
#set term wxt
#set term post eps enh 12
set size ratio 1
#set out 'comparison1.eps'

#set terminal png size 400,400 enhanced font "Helvetica,20"
#set output 'ELIMG1.png'
#set nokey


set autoscale yfix

plot  file0  tit ""          w filledcurves lt 1 lc rgb "#FFA346", \
      file1  tit ""   w l lt 1 lc 3 lw 3, \
      file2  tit ""   w l lt 1 lc 1 lw 3
#      file2  tit "DAG-Enclosure"   w l lt 1 lc  rgb "#00CC00", \

pause -1 "PRESS <ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'test1.eps'
rep
set term x11
set out
!ps2eps -B -l -f test1.eps
!mv test1.eps.eps test1.eps

