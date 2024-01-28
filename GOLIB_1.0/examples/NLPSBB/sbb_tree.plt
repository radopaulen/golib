set size ratio 1
set xlabel 'p1'
set ylabel 'p2'
unset key

set parametric
p1opt = 6.000000e+00
p2opt = 0.666667e+00

plot 'sbb_tree.out' u 1:2 w l, p1opt,p2opt w p pt 7

pause -1 "<ENTER> TO CONTINUE"

set term post eps enh solid color 18
set out 'sbb_tree.eps'
rep
set term wxt
!ps2eps -B -f -l sbb_tree.eps
!mv sbb_tree.eps.eps sbb_tree.eps
