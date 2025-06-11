set arrow from 1,1.11 to 590,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_49|Q2G1F5|Solute-binding|CP000253.1|tpos:56998-57587"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:590]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096803/antigen_49_Q2G1F5_Solute-binding_CP000253.1_tpos_56998-57587.eps"
plot "./TMHMM_1096803/antigen_49_Q2G1F5_Solute-binding_CP000253.1_tpos_56998-57587.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
