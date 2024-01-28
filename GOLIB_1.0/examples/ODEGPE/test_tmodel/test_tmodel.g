set xlabel 'x'
set ylabel 'y'
#set hidden3d
set view 74,13
set xrange [-1:1]
set yrange [-1:1]

set key below

splot 'TM-orig.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Original model' w l lt 2, \
  '' u 1:2:5 tit '' w l lt 2, \
  'TM-updt.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Reduced model' w l lt 2, \
  '' u 1:2:5 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"

#set hidden3d

splot 'TM-updt.out' u 1:2:3 tit 'Function' w l lt 1, \
  '' u 1:2:4 tit 'Reduced model' w l lt 2, \
  '' u 1:2:5 tit '' w l lt 2

pause -1 "<ENTER> TO CONTINUE"
