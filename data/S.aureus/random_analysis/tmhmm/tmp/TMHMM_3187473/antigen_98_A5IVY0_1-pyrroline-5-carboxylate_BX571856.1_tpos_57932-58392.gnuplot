set arrow from 1,1.11 to 461,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_98|A5IVY0|1-pyrroline-5-carboxylate|BX571856.1|tpos:57932-58392"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:461]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_98_A5IVY0_1-pyrroline-5-carboxylate_BX571856.1_tpos_57932-58392.eps"
plot "./TMHMM_3187473/antigen_98_A5IVY0_1-pyrroline-5-carboxylate_BX571856.1_tpos_57932-58392.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
