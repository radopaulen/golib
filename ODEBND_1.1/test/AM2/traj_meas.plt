kVFA = 64e-3
Kb   = 6.5e-4

set xlabel 't [day]'
set xrange [0:4]
set ylabel 'X1 [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:2:3 tit '' w filledcurves lt 9, \
'./exact.out' u 1:2 tit '' w l lt 1 lw 2

pause -1

set ylabel 'X2 [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:4:5 tit '' w filledcurves lt 9, \
'./exact.out' u 1:4 tit '' w l lt 1 lw 2

pause -1

set ylabel 'CODs [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:($6+$8*kVFA):($7+$9*kVFA) tit '' w filledcurves lt 9, \
'./exact.out' u 1:($6+$8*kVFA) tit '' w l lt 1 lw 2, \
'data/inra_psfb_020404_CODs_sampled.txt' u 1:2:($2*0.1) w yerrorbars lc -1

pause -1

set ylabel 'VFA [g/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:($8*kVFA):($9*kVFA) tit '' w filledcurves lt 9, \
'./exact.out' u 1:($8*kVFA) tit '' w l lt 1 lw 2, \
'data/inra_psfb_020404_VFA_sampled.txt' u 1:2:($2*0.1) w yerrorbars lc -1

pause -1

set ylabel 'ALK [mmol/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:10:11 tit '' w filledcurves lt 9, \
'./exact.out' u 1:10 tit '' w l lt 1 lw 2, \
'data/inra_psfb_020404_Z_sampled.txt' u 1:2:($2*0.1) w yerrorbars lc -1

pause -1

set ylabel 'TIC [mmol/L]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:12:13 tit '' w filledcurves lt 9, \
'./exact.out' u 1:12 tit '' w l lt 1 lw 2, \
'data/inra_psfb_020404_TIC_sampled.txt' u 1:2:($2*0.1) w yerrorbars lc -1

pause -1

set ylabel 'pH [-]'
#set yrange [0.1:0.6]

plot './exact.out' u 1:(-log10(Kb*1e-3*($12/($10-$8)-1.))):(-log10(Kb*1e-3*($13/($11-$9)-1.))) tit '' w filledcurves lt 9, \
'./exact.out' u 1:(-log10(Kb*1e-3*($12/($10-$8)-1.))) tit '' w l lt 1 lw 2, \
'data/inra_psfb_020404_pH_sampled.txt' u 1:2:(0.1) w yerrorbars lc -1

pause -1
