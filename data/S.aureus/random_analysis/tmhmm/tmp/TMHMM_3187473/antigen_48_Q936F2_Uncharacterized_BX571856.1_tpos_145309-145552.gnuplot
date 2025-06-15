set arrow from 1,1.11 to 244,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_48|Q936F2|Uncharacterized|BX571856.1|tpos:145309-145552"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:244]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_48_Q936F2_Uncharacterized_BX571856.1_tpos_145309-145552.eps"
plot "./TMHMM_3187473/antigen_48_Q936F2_Uncharacterized_BX571856.1_tpos_145309-145552.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
