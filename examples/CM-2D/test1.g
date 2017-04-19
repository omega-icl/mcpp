resfile = 'test1.out'

set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 55,215

set key below

splot resfile u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Chebyshev approximant' w l lt 2

set term post eps enh solid color 18
set out 'test1_app.eps'
rep
set term wxt
set out

pause -1 "<ENTER> TO CONTINUE"  

splot resfile u 1:2:($4-$3) tit 'Chebyshev approximation error' w l lt 2

set term post eps enh solid color 18
set out 'test1_err.eps'
rep
set term wxt
set out

pause -1 "<ENTER> TO CONTINUE"  

splot resfile u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:5 tit 'Chebyshev model' w l lt 2, \
  '' u 1:2:6 tit '' w l lt 2

set term post eps enh solid color 18
set out 'test1_cmod.eps'
rep
set term wxt
set out

pause -1 "<ENTER> TO CONTINUE"  

splot resfile u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:7 tit 'Chebyshev-derived bounds' w l lt 2, \
  '' u 1:2:8 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

