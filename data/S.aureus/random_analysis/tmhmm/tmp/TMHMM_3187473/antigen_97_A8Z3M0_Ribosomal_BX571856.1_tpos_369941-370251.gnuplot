set arrow from 1,1.11 to 311,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_97|A8Z3M0|Ribosomal|BX571856.1|tpos:369941-370251"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:311]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_97_A8Z3M0_Ribosomal_BX571856.1_tpos_369941-370251.eps"
plot "./TMHMM_3187473/antigen_97_A8Z3M0_Ribosomal_BX571856.1_tpos_369941-370251.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
