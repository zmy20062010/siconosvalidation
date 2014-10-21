result = "resultRodImpact.dat"

set terminal postscript eps enhanced color
set outpu "phase"
set multiplot layout 2,3 rowsfirst

plot result \
u 1:2 t "t--theta_1" lc rgb "red"
u 1:5 t "t--omega_1" lc rgb "blue"
u 2:5 t "theta_1--theta_1" lc rgb "magenta"

plot result \
u 1:3 t "t--x" lc rgb "red"
u 1:6 t "t--v" lc rgb "blue"
u 3:6 t "x--v" lc rgb "magenta"

unset  multiplot
