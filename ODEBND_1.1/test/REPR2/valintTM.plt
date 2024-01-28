set xlabel 'Time'
set xrange [:200]
set ylabel 'x1'
set key spacing 1.5

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTM1.out' u 1:2 tit 'q=1' w l lt 1 lc 3 lw 3, '' u 1:3 tit '' w l lt 1 lc 3 lw 3

set term post eps enh solid color 18
set out "valintTM1.eps"
rep
set out
set term wxt
!ps2eps -B -f -l valintTM1.eps
!mv valintTM1.eps.eps valintTM1.eps

pause -1
set xlabel 'Time'
set xrange [:200]
set ylabel 'x1'
set key spacing 1.5

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTM2.out' u 1:2 tit 'q=2' w l lt 1 lc 2 lw 3, '' u 1:3 tit '' w l lt 1 lc 2 lw 3, './valintTM1.out' u 1:2 tit 'q=1' w l lt 1 lc 3 lw 3, '' u 1:3 tit '' w l lt 1 lc 3 lw 3

set term post eps enh solid color 18
set out "valintTM2.eps"
rep
set out
set term wxt
!ps2eps -B -f -l valintTM2.eps
!mv valintTM2.eps.eps valintTM2.eps

pause -1
set xlabel 'Time'
set xrange [:200]
set ylabel 'x1'
set key spacing 1.5

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTM3.out' u 1:2 tit 'q=3' w l lt 1 lc 1 lw 3, '' u 1:3 tit '' w l lt 1 lc 1 lw 3, './valintTM2.out' u 1:2 tit 'q=2' w l lt 1 lc 2 lw 3, '' u 1:3 tit '' w l lt 1 lc 2 lw 3, './valintTM1.out' u 1:2 tit 'q=1' w l lt 1 lc 3 lw 3, '' u 1:3 tit '' w l lt 1 lc 3 lw 3

set term post eps enh solid color 18
set out "valintTM.eps"
rep
set out
set term wxt
!ps2eps -B -f -l valintTM.eps
!mv valintTM.eps.eps valintTM.eps

pause -1

set xlabel 'Time'
set xrange [150:250]
set ylabel 'x1'
unset key
set size ratio 1

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./valintTE3.out' u 1:2 tit 'q=3' w l lt 1 lc 1 lw 3, '' u 1:3 tit '' w l lt 1 lc 1 lw 3

set term post eps enh solid color 18
set out "valintTM_zoom.eps"
rep
set out
set term wxt
!ps2eps -B -f -l valintTM_zoom.eps
!mv valintTM_zoom.eps.eps valintTM_zoom.eps

pause -1

