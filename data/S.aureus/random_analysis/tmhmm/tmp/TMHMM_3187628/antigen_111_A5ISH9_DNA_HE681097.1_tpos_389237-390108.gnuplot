set arrow from 1,1.11 to 872,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_111|A5ISH9|DNA|HE681097.1|tpos:389237-390108"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:872]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_111_A5ISH9_DNA_HE681097.1_tpos_389237-390108.eps"
plot "./TMHMM_3187628/antigen_111_A5ISH9_DNA_HE681097.1_tpos_389237-390108.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
