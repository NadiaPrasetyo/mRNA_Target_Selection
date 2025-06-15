set arrow from 1,1.11 to 96,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_179|Q5HEW2|Serine|BX571856.1|tpos:581712-581807"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:96]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_179_Q5HEW2_Serine_BX571856.1_tpos_581712-581807.eps"
plot "./TMHMM_3187473/antigen_179_Q5HEW2_Serine_BX571856.1_tpos_581712-581807.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
