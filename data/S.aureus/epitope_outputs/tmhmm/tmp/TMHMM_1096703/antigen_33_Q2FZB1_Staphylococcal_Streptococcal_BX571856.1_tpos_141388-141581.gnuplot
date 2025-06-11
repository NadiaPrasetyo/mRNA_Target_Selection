set arrow from 1,1.11 to 194,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_33|Q2FZB1|Staphylococcal/Streptococcal|BX571856.1|tpos:141388-141581"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:194]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096703/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_BX571856.1_tpos_141388-141581.eps"
plot "./TMHMM_1096703/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_BX571856.1_tpos_141388-141581.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
