set arrow from 1,1.07 to 89,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_33|Q2FZB1|Staphylococcal/Streptococcal|CP002114.3|tpos:132131-132219"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:89]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_CP002114.3_tpos_132131-132219.eps"
plot "./TMHMM_1211751/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_CP002114.3_tpos_132131-132219.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
