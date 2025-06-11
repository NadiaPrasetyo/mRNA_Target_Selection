set arrow from 1,1.11 to 428,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_19|Q2FWW3|Aminotransferase|BA000018.3|tpos:618899-619326"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:428]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_19_Q2FWW3_Aminotransferase_BA000018.3_tpos_618899-619326.eps"
plot "./TMHMM_1096793/antigen_19_Q2FWW3_Aminotransferase_BA000018.3_tpos_618899-619326.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
