set xlabel 'Time'
set xrange [:100]
set ylabel 'x1'
set yrange [-0.1:1]

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./gslint_TM6+ELL1.out' u 1:2 tit 'q=6+ell=1' w l lt 1 lc 1 lw 2, '' u 1:3 tit '' w l lt 1 lc 1 lw 2,\
'./gslint_ELL1.out' u 1:2 tit 'ell=1' w l lt 1 lc 3 lw 2, '' u 1:3 tit '' w l lt 1 lc 3 lw 2, \
'./gslint_TM6+DI.out' u 1:2 tit 'q=6+di' w l lt 1 lc 4 lw 2, '' u 1:3 tit '' w l lt 1 lc 4 lw 2,\
'./gslint_DI.out' u 1:2 tit 'di' w l lt 1 lc 0 lw 2, '' u 1:3 tit '' w l lt 1 lc 0 lw 2

pause -1

set xlabel 'Time'
set xrange [:100]
set ylabel 'x2'
#set yrange [-5:5]

plot './exact.out' u 1:4:5 tit '' w filledcurves lt 9, \
'./gslint_TM6+ELL1.out' u 1:4 tit 'q=6 ell=1' w l lt 1 lc 1 lw 2, '' u 1:5 tit '' w l lt 1 lc 1 lw 2,\
'./gslint_ELL1.out' u 1:4 tit 'ell=1' w l lt 1 lc 3 lw 2, '' u 1:5 tit '' w l lt 1 lc 3 lw 2, \
'./gslint_TM6+DI.out' u 1:4 tit 'q=6+di' w l lt 1 lc 4 lw 2, '' u 1:5 tit '' w l lt 1 lc 4 lw 2,\
'./gslint_DI.out' u 1:4 tit 'di' w l lt 1 lc 0 lw 2, '' u 1:5 tit '' w l lt 1 lc 0 lw 2

pause -1
