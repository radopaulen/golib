resfile = 'FPREL-1D.out'

set xlabel 'x'
set key below

plot resfile u 1:2 tit 'function' w l lt 1, \
  '' u 1:3 tit 'bounds' w l lt 2, \
  '' u 1:4 tit '' w l lt 2, \
  '' u 1:5 tit 'relaxations' w l lt 3, \
  '' u 1:6 tit '' w l lt 3
 
pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 21
set out 'FPREL-1D.eps'
rep
set term x11
