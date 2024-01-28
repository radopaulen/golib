set term post eps enh dashed color 16 size 6,2.5
set out "valintTM2.eps"

#set size 0.9,0.9
set origin 0,0
set multiplot layout 1,2 columnsfirst 
#scale 1.1,0.9
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'X2 [g/L]'
set yrange [0.6:1.1]
unset key
set key spacing 3 top right
#
plot './exact.out' u 1:4:5 tit '' w filledcurves lt 1 lc 6, \
'./valintTE3.out' u 1:4 tit 'q=3' w l lt 1 lc 3 lw 3, '' u 1:5 tit '' w l lt 1 lc 3 lw 3, \
'./valintTE2.out' u 1:4 tit 'q=2' w l lt 3 lc 1 lw 3, '' u 1:5 tit '' w l lt 3 lc 1 lw 3
#
set xlabel 'Time'
#set xrange [:200]
set ylabel 'S2 [mmol/L]'
set yrange [2:3.5]
#unset key
set key spacing 3 bottom right
#set key outside
#
plot './exact.out' u 1:8:9 tit '' w filledcurves lt 1 lc 6, \
'./valintTE3.out' u 1:8 tit 'q=3' w l lt 1 lc 3 lw 3, '' u 1:9 tit '' w l lt 1 lc 3 lw 3, \
'./valintTE2.out' u 1:8 tit 'q=2' w l lt 3 lc 1 lw 3, '' u 1:9 tit '' w l lt 3 lc 1 lw 3
#
#set xlabel 'Time'
##set xrange [:200]
#set ylabel 'C [mmol/L]'
#set yrange [30:60]
#set autoscale y
#unset key
##set key spacing 1.5
##
#plot './exact.out' u 1:12:13 tit '' w filledcurves lt 9, \
#'./valintTE3.out' u 1:12 tit 'q=3' w l lt 1 lc 3 lw 2, '' u 1:13 tit '' w l lt 1 lc 3 lw 2, \
#'./valintTE2.out' u 1:12 tit 'q=2' w l lt 1 lc 1 lw 2, '' u 1:13 tit '' w l lt 1 lc 1 lw 2
##
##unset axis
##unset label
unset multiplot

set term wxt
set out
!ps2eps -f -B -l valintTM2.eps
!mv valintTM2.eps.eps valintTM2.eps
