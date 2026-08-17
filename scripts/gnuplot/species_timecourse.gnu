# Species time courses, adapted from the archived plot_speciestimedyns_4.gnu.
# Reads the column layout written by dycelfem_analyze.py:
#   1:Hour  2:S1_mean 3:S1_sd  4:S2_mean 5:S2_sd ...
# Usage:  gnuplot -e "datafile='WoundCytokinesTempoDyns.dat'" species_timecourse.gnu
if (!exists("datafile")) datafile = 'WoundCytokinesTempoDyns.dat'
if (!exists("outfile"))  outfile  = datafile.'.eps'
set term postscript color enhanced eps font "Arial,22"
set out outfile
set xlabel 'Hours'
set ylabel 'mean concentration per cell (a.u.)'
set border 4095
set key inside right top spacing 1.3
set grid
set style line 1 lt 1 lc rgb "grey"    lw 5
set style line 2 lt 1 lc rgb "red"     lw 5
set style line 3 lt 1 lc rgb "magenta" lw 5
set style line 4 lt 1 lc rgb "blue"    lw 5
set style line 5 lt 1 lc rgb "orange"  lw 5
set style line 6 lt 1 lc rgb "green"   lw 5
set style line 7 lt 1 lc rgb "purple"  lw 5
plot datafile u 1:6  w l ls 1 title "Collagen", \
     datafile u 1:8  w l ls 2 title "PDGF", \
     datafile u 1:10 w l ls 4 title "KGF", \
     datafile u 1:12 w l ls 5 title "IL1", \
     datafile u 1:14 w l ls 6 title "TGFb1", \
     datafile u 1:16 w l ls 7 title "Procollagen"
