set xlabel 'x'
set ylabel 'y'
set hidden3d
set view 55,215

set key below

splot 'CM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Chebyshev model' w l lt 2, \
  '' u 1:2:5 tit '' w l lt 2

set term post eps enh solid color 18
set out 'CM-2D.eps'
rep
set term wxt
set out

pause -1 "<ENTER> TO CONTINUE"  

splot 'CM-2D.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:6 tit 'Chebyshev-derived bounds' w l lt 2, \
  '' u 1:2:7 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"  

