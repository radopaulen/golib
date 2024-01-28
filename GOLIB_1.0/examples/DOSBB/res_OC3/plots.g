set xlabel 'Nodes'
set ylabel 'Gap'
set logscale y
set yrange [7e-4:1e3]
set key spacing 1.7

plot 'NS=5_ORDER=3_TM_EIGEN.out' u 1:($4-$3) tit '3rd-order' w l lt 1 lc 2 lw 4, \
     'NS=5_ORDER=4_TM_EIGEN.out' u 1:($4-$3) tit '4th-order' w l lt 1 lc 1 lw 4, \
     'NS=5_ORDER=5_TM_EIGEN.out' u 1:($4-$3) tit '5th-order' w l lt 1 lc 3 lw 4

pause -1 "ENTER TO CONTINUE"

set term post eps enh 21
set out 'NS=5_ORDER=3-5_TM_EIGEN.eps'
rep
reset
