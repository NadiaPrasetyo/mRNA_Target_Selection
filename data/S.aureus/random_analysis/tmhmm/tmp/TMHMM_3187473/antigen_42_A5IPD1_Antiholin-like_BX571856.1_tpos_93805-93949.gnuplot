set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 30,1.09 nohead lt 1 lw 40
set arrow from 31,1.11 to 39,1.11 nohead lt 4 lw 10
set arrow from 40,1.09 to 62,1.09 nohead lt 1 lw 40
set arrow from 63,1.07 to 68,1.07 nohead lt 3 lw 10
set arrow from 69,1.09 to 91,1.09 nohead lt 1 lw 40
set arrow from 92,1.11 to 94,1.11 nohead lt 4 lw 10
set arrow from 95,1.09 to 117,1.09 nohead lt 1 lw 40
set arrow from 118,1.07 to 145,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_42|A5IPD1|Antiholin-like|BX571856.1|tpos:93805-93949"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:145]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_42_A5IPD1_Antiholin-like_BX571856.1_tpos_93805-93949.eps"
plot "./TMHMM_3187473/antigen_42_A5IPD1_Antiholin-like_BX571856.1_tpos_93805-93949.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
