set arrow from 1,1.07 to 77,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_65|Q2G2Y0|Staphylococcal|BA000018.3|tpos:128994-129070"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:77]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_65_Q2G2Y0_Staphylococcal_BA000018.3_tpos_128994-129070.eps"
plot "./TMHMM_1211906/antigen_65_Q2G2Y0_Staphylococcal_BA000018.3_tpos_128994-129070.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
