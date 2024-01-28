set xlabel 't [day]'
set xrange [0:4]
set ylabel 'X1 [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:2 tit '' w l lt 1 lw 2, ''  u 1:3 tit '' w l lt 1 lw 2

pause -1

set ylabel 'X2 [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:4:5 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:4 tit '' w l lt 1 lw 2, ''  u 1:5 tit '' w l lt 1 lw 2

pause -1

set ylabel 'S1 [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:6:7 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:6 tit '' w l lt 1 lw 2, ''  u 1:7 tit '' w l lt 1 lw 2

pause -1

set ylabel 'S2 [mmol/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:8:9 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:8 tit '' w l lt 1 lw 2, ''  u 1:9 tit '' w l lt 1 lw 2

pause -1

set ylabel 'Z [mmol/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:10:11 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:10 tit '' w l lt 1 lw 2, ''  u 1:11 tit '' w l lt 1 lw 2

pause -1

set ylabel 'C [mmol/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:12:13 tit '' w filledcurves lt 9, \
'./gslint_TM3+ELL1.out' u 1:12 tit '' w l lt 1 lw 2, ''  u 1:13 tit '' w l lt 1 lw 2

pause -1

