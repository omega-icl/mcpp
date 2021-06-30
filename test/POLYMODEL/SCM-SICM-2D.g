set xlabel 'x'
set ylabel 'y'
set zlabel 'f(x,y)' rotate parallel
#set view 59,332
set hidden3d

splot 'SICM-2D.out' u 1:2:3 tit 'Function' w l lc -1, \
      'SICM-2D.out' u 1:2:4 tit 'Interval Chebyshev predictor' w l lc 7

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'SICM-2D_pred.eps'
rep
set term qt
!ps2eps -B -f -l SICM-2D_pred.eps
!mv SICM-2D_pred.eps.eps SICM-2D_pred.eps
!gv SICM-2D_pred.eps &

splot 'SICM-2D.out' u 1:2:3 tit 'Function' w l lc -1, \
      'SICM-2D.out' u 1:2:5 tit 'Interval Chebyshev range' w l lc 7, \
      'SICM-2D.out' u 1:2:6 tit '' w l lc 7

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'SICM-2D_mod.eps'
rep
set term qt
!ps2eps -B -f -l SICM-2D_mod.eps
!mv SICM-2D_mod.eps.eps SICM-2D_mod.eps
!gv SICM-2D_mod.eps &

splot 'SCM-2D.out' u 1:2:3 tit 'Function' w l lc -1, \
      'SCM-2D.out' u 1:2:5 tit 'Interval Chebyshev range' w l lc 7, \
      'SCM-2D.out' u 1:2:6 tit '' w l lc 7

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'SCM-2D_mod.eps'
rep
set term qt
!ps2eps -B -f -l SCM-2D_mod.eps
!mv SCM-2D_mod.eps.eps SCM-2D_mod.eps
!gv SCM-2D_mod.eps &

set zlabel 'max{ P(x,y) - f(x,y) }' rotate parallel
set zrange [0:]

splot 'SCM-2D.out'  u 1:2:($6-$3) tit 'Standard Chebyshev model' w l lc -1, \
      'SICM-2D.out' u 1:2:($6-$3) tit 'Interval Chebyshev model' w l lc 7

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'SCM-SICM-2D_err.eps'
rep
set term qt
!ps2eps -B -f -l SCM-SICM-2D_err.eps
!mv SCM-SICM-2D_err.eps.eps SCM-SICM-2D_err.eps
!gv SCM-SICM-2D_err.eps &

