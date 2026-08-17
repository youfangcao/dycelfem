# Wound closure, adapted from the archived plot_woundclosure.gnu.
if (!exists("datafile")) datafile = 'WoundSizeHours.dat'
set term postscript color enhanced eps font "Arial,22"
set out 'WoundSizeHours.eps'
set xlabel 'Hours'
set ylabel 'Epithelial gap'
set border 4095
set grid
set yrange [0:*]
set style line 1 lt 1 lc rgb "red" lw 4
plot datafile u 1:2:3 notitle w yerrorbars ls 1, \
     datafile u 1:2   notitle w l ls 1
