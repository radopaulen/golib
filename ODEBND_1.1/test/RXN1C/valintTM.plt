set term post eps enh dashed color 16 size 6,2.5
set out "valintTM.eps"

#set size 0.9,0.9
set origin 0,0
set multiplot layout 1,2 columnsfirst 
#scale 1.1,0.9
#
set xlabel 'Time'
set xrange [0:2]
set ylabel 'xB'
set yrange [0.:1.]
unset key
set key spacing 3 top left reverse Left
#
plot './exact.out' u 1:4:5 tit '' w filledcurves lt 1 lc 6, \
'./valintTE3.out' u 1:4 tit 'q=3' w l lt 3 lc 1 lw 3, '' u 1:5 tit '' w l lt 3 lc 1 lw 3, \
'./valintTE3_inv.out' u 1:4 tit 'q=3 w/ inv' w l lt 1 lc 3 lw 3, '' u 1:5 tit '' w l lt 1 lc 3 lw 3
#
set ylabel 'xC'
set yrange [0.:1.]
#unset key
set key spacing 3 top left reverse Left
#set key outside
#
plot './exact.out' u 1:6:7 tit '' w filledcurves lt 1 lc 6, \
'./valintTE3.out' u 1:6 tit 'q=3' w l lt 3 lc 1 lw 3, '' u 1:7 tit '' w l lt 3 lc 1 lw 3, \
'./valintTE3_inv.out' u 1:6 tit 'q=3 w/ inv' w l lt 1 lc 3 lw 3, '' u 1:7 tit '' w l lt 1 lc 3 lw 3
unset multiplot

set term wxt
set out
!ps2eps -f -B -l valintTM.eps
!mv valintTM.eps.eps valintTM.eps
