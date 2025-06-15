set arrow from 1,1.11 to 269,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_125|A6U0B2|NAD|CP000253.1|tpos:281693-281961"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:269]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_125_A6U0B2_NAD_CP000253.1_tpos_281693-281961.eps"
plot "./TMHMM_3187502/antigen_125_A6U0B2_NAD_CP000253.1_tpos_281693-281961.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
