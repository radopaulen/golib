set xlabel 'Time'
set xrange [:50]
set ylabel 'x1'

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTM5.out' u 1:2 tit 'q=5' w l lt 1 lc 3 lw 2, '' u 1:3 tit '' w l lt 1 lc 3 lw 2, './valintTM8.out' u 1:2 tit 'q=8' w l lt 1 lc 2 lw 2, '' u 1:3 tit '' w l lt 1 lc 2 lw 2

pause -1

