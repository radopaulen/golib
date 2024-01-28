resfile = 'FPREL-2D.out'

set xlabel 'x'
set ylabel 'y'
set zlabel 'f(x,y)'
set hidden3d

### BIVARIATE FUNCTION

splot resfile u 1:2:3 tit '' w l lt 1
  
pause -1 "<ENTER> TO CONTINUE"

### CONVEX/CONCAVE RELAXATIONS

splot resfile u 1:2:3 tit '' w l lt 1, '' u 1:2:6 tit '' w l lt 2,'' u 1:2:7 tit '' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"

### INTERVAL BOUNDS

splot resfile u 1:2:3 tit '' w l lt 1, '' u 1:2:4 tit '' w l lt 2,'' u 1:2:5 tit '' w l lt 3
  
pause -1 "<ENTER> TO CONTINUE"
