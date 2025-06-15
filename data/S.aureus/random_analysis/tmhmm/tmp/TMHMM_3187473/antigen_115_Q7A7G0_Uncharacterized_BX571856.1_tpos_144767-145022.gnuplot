set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 24,1.09 nohead lt 1 lw 40
set arrow from 25,1.11 to 256,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_115|Q7A7G0|Uncharacterized|BX571856.1|tpos:144767-145022"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:256]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_115_Q7A7G0_Uncharacterized_BX571856.1_tpos_144767-145022.eps"
plot "./TMHMM_3187473/antigen_115_Q7A7G0_Uncharacterized_BX571856.1_tpos_144767-145022.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
