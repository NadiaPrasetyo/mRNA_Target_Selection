set arrow from 1,1.11 to 341,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_53|Q2FIK4|Extracellular|BX571856.1|tpos:275462-275802"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:341]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_53_Q2FIK4_Extracellular_BX571856.1_tpos_275462-275802.eps"
plot "./TMHMM_3187473/antigen_53_Q2FIK4_Extracellular_BX571856.1_tpos_275462-275802.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
