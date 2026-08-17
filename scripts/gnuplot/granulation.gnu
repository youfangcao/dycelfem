# Granulation tissue, adapted from the archived plot_granulation.gnu.
if (!exists("datafile")) datafile = 'Granulation.dat'
set term postscript color enhanced eps font "Arial,22"
set out 'Granulation.eps'
set xlabel 'Hours'
set ylabel 'Area in wound bed'
set border 4095
set grid
set key inside left top spacing 1.3
set style line 1 lt 1 lc rgb "magenta" lw 4
set style line 3 lt 1 lc rgb "blue"    lw 4
set style line 5 lt 1 lc rgb "grey"    lw 4
plot datafile u 1:8:9   w yerrorlines ls 1 title "Fibroblast", \
     datafile u 1:10:11 w yerrorlines ls 3 title "Keratinocyte", \
     datafile u 1:12:13 w yerrorlines ls 5 title "ECM"
