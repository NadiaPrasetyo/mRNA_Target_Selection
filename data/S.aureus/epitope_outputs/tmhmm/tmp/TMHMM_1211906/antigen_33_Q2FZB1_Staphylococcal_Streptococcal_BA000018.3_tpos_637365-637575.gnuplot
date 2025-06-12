set arrow from 1,1.11 to 211,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_33|Q2FZB1|Staphylococcal/Streptococcal|BA000018.3|tpos:637365-637575"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:211]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_BA000018.3_tpos_637365-637575.eps"
plot "./TMHMM_1211906/antigen_33_Q2FZB1_Staphylococcal_Streptococcal_BA000018.3_tpos_637365-637575.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
