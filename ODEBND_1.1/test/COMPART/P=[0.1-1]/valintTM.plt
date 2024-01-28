set xlabel 'Time'
set xrange [:15]
set ylabel 'x1'
set yrange [-5:5]

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./gslint_TM6+ELL1.out' u 1:2 tit 'q=6 ell=1' w l lt 1 lc 1 lw 2, '' u 1:3 tit '' w l lt 1 lc 1 lw 2,\
'./gslint_TM4+ELL1.out' u 1:2 tit 'q=4 ell=1' w l lt 1 lc 3 lw 2, '' u 1:3 tit '' w l lt 1 lc 3 lw 2,\
'./gslint_TM2+ELL1.out' u 1:2 tit 'q=2 ell=1' w l lt 1 lc 2 lw 2, '' u 1:3 tit '' w l lt 1 lc 2 lw 2

pause -1

