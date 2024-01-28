set term post eps enh solid color 16 size 7,4 
set out "valintTM.eps"

#set size 0.9,0.9
set origin 0,0
set multiplot layout 2,3 columnsfirst 
#scale 1.1,0.9
#
set xlabel 'Time [day]'
set xrange [:40]
set ylabel 'X1 [g/L]'
set yrange [0.1:0.6]
unset key
#set key spacing 1.5
#
plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:2 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:3 tit '' w l lt 1 lc 3 lw 2, \
'./valintTE2.out' u 1:2 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:3 tit '' w l lt 1 lc 1 lw 2
#
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'X2 [g/L]'
set yrange [0.6:1.2]
unset key
#set key spacing 1.5
#
plot './exact.out' u 1:4:5 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:4 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:5 tit '' w l lt 1 lc 3 lw 2, \
'./valintTE2.out' u 1:4 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:5 tit '' w l lt 1 lc 1 lw 2
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'S1 [g/L]'
set autoscale y
set yrange [0.5:1.5]
unset key
#set key spacing 1.5
#
plot './exact.out' u 1:6:7 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:6 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:7 tit '' w l lt 1 lc 3 lw 2, \
'./valintTE2.out' u 1:6 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:7 tit '' w l lt 1 lc 1 lw 2
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'S2 [mmol/L]'
set yrange [2:5]
#unset key
set key spacing 3 at 70,4.5
#set key outside
#
plot './exact.out' u 1:8:9 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:8 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:9 tit '' w l lt 1 lc 3 lw 2, \
'./valintTE2.out' u 1:8 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:9 tit '' w l lt 1 lc 1 lw 2
#
#set xlabel 'Time'
#set xrange [:200]
#set ylabel 'Z [mmol/L]'
#set yrange [0:100]
#set key spacing 1.5
#
#plot './exact.out' u 1:10:11 tit '' w filledcurves lt 9, \
#'./valintTE3.out' u 1:10 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:11 tit '' w l lt 1 lc 3 lw 2, \
#'./valintTE2.out' u 1:10 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:11 tit '' w l lt 1 lc 1 lw 2
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'C [mmol/L]'
set yrange [30:60]
set autoscale y
unset key
#set key spacing 1.5
#
plot './exact.out' u 1:12:13 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:12 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:13 tit '' w l lt 1 lc 3 lw 2, \
'./valintTE2.out' u 1:12 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:13 tit '' w l lt 1 lc 1 lw 2
#
#unset axis
#unset label
unset multiplot

set term wxt
set out
!ps2eps -f -B -l valintTM.eps
!mv valintTM.eps.eps valintTM.eps
