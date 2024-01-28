set xlabel 'Time'
set xrange [:2.0]
set yrange [-0.4:1.2]

set grid front
set term post eps enh 12
set out 'valintTM.eps'
set multiplot layout 2,2
#
set xlabel 'Time'
set ylabel 'A'
plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTM_q=10.out' u 1:2 tit 'q=10' w l lt 1 lc 3 lw 2, '' u 1:3 tit '' w l lt 1 lc 3 lw 2, './valintTM_q=10_nbnd.out' u 1:2 tit 'q=10' w l lt 3 lc 3 lw 2, '' u 1:3 tit '' w l lt 3 lc 3 lw 2, './valintTM_q=8.out' u 1:2 tit 'q=8' w l lt 1 lc 2 lw 2, '' u 1:3 tit '' w l lt 1 lc 2 lw 2, './valintTM_q=6.out' u 1:2 tit 'q=6' w l lt 1 lc 4 lw 2, '' u 1:3 tit '' w l lt 1 lc 4 lw 2, './valintTM_q=4.out' u 1:2 tit 'q=4' w l lt 1 lc 1 lw 2, '' u 1:3 tit '' w l lt 1 lc 1 lw 2, './valintTM_q=4_nbnd.out' u 1:2 tit 'q=4' w l lt 3 lc 1 lw 2, '' u 1:3 tit '' w l lt 3 lc 1 lw 2
#
set xlabel 'Time'
set ylabel 'B'
plot './exact.out' u 1:4:5 tit '' w filledcurves lt 9, \
'./valintTM_q=10.out' u 1:4 tit 'q=10' w l lt 1 lc 3 lw 2, '' u 1:5 tit '' w l lt 1 lc 3 lw 2, './valintTM_q=10_nbnd.out' u 1:4 tit 'q=10' w l lt 3 lc 3 lw 2, '' u 1:5 tit '' w l lt 3 lc 3 lw 2, './valintTM_q=8.out' u 1:4 tit 'q=8' w l lt 1 lc 2 lw 2, '' u 1:5 tit '' w l lt 1 lc 2 lw 2, './valintTM_q=6.out' u 1:4 tit 'q=6' w l lt 1 lc 4 lw 2, '' u 1:5 tit '' w l lt 1 lc 4 lw 2, './valintTM_q=4.out' u 1:4 tit 'q=4' w l lt 1 lc 1 lw 2, '' u 1:5 tit '' w l lt 1 lc 1 lw 2, './valintTM_q=4_nbnd.out' u 1:4 tit 'q=4' w l lt 3 lc 1 lw 2, '' u 1:5 tit '' w l lt 3 lc 1 lw 2
#
set xlabel 'Time'
set ylabel 'C'
plot './exact.out' u 1:6:7 tit '' w filledcurves lt 9, \
'./valintTM_q=10.out' u 1:6 tit 'q=10' w l lt 1 lc 3 lw 2, '' u 1:7 tit '' w l lt 1 lc 3 lw 2, './valintTM_q=10_nbnd.out' u 1:6 tit 'q=10' w l lt 3 lc 3 lw 2, '' u 1:7 tit '' w l lt 3 lc 3 lw 2, './valintTM_q=8.out' u 1:6 tit 'q=8' w l lt 1 lc 2 lw 2, '' u 1:7 tit '' w l lt 1 lc 2 lw 2, './valintTM_q=6.out' u 1:6 tit 'q=6' w l lt 1 lc 4 lw 2, '' u 1:7 tit '' w l lt 1 lc 4 lw 2, './valintTM_q=4.out' u 1:6 tit 'q=4' w l lt 1 lc 1 lw 2, '' u 1:7 tit '' w l lt 1 lc 1 lw 2, './valintTM_q=4_nbnd.out' u 1:6 tit 'q=4' w l lt 3 lc 1 lw 2, '' u 1:7 tit '' w l lt 3 lc 1 lw 2
#
set xlabel 'Time'
set ylabel 'D'
plot './exact.out' u 1:8:9 tit '' w filledcurves lt 9, \
'./valintTM_q=10.out' u 1:8 tit 'q=10' w l lt 1 lc 3 lw 2, '' u 1:9 tit '' w l lt 1 lc 3 lw 2, './valintTM_q=10_nbnd.out' u 1:8 tit 'q=10' w l lt 3 lc 3 lw 2, '' u 1:9 tit '' w l lt 3 lc 3 lw 2, './valintTM_q=8.out' u 1:8 tit 'q=8' w l lt 1 lc 2 lw 2, '' u 1:9 tit '' w l lt 1 lc 2 lw 2, './valintTM_q=6.out' u 1:8 tit 'q=6' w l lt 1 lc 4 lw 2, '' u 1:9 tit '' w l lt 1 lc 4 lw 2, './valintTM_q=4.out' u 1:8 tit 'q=4' w l lt 1 lc 1 lw 2, '' u 1:9 tit '' w l lt 1 lc 1 lw 2, './valintTM_q=4_nbnd.out' u 1:8 tit 'q=4' w l lt 3 lc 1 lw 2, '' u 1:9 tit '' w l lt 3 lc 1 lw 2
#
unset multiplot
set term wxt

pause -1

