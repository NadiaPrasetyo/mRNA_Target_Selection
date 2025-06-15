set arrow from 1,1.07 to 98,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_67|A5IUL1|Serine-protein|BX571856.1|tpos:220861-220958"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:98]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_67_A5IUL1_Serine-protein_BX571856.1_tpos_220861-220958.eps"
plot "./TMHMM_3187473/antigen_67_A5IUL1_Serine-protein_BX571856.1_tpos_220861-220958.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
