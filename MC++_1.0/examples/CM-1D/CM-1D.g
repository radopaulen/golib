set xlabel 'x'
set key below

plot 'CM-1D.out' u 1:2 tit 'function' w l lc 1 lt 1, \
  '' u 1:4 tit 'Chebyshev model' w l lt 1 lc 3, \
  '' u 1:5 tit '' w l lt 1 lc 3, \
  '' u 1:3 tit '' w l lt 0 lc 3, \
  '' u 1:6 tit 'Chebyshev-derived bounds' w l lt 1 lc 2, \
  '' u 1:7 tit '' w l lt 1 lc 2

set term post eps enh solid color 18
set out 'CM-1D.eps'
rep
set term wxt
set out

pause -1 "<ENTER> TO CONTINUE"
