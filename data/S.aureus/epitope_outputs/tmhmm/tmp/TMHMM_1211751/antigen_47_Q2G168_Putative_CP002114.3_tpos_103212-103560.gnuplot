set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 35,1.09 nohead lt 1 lw 40
set arrow from 36,1.11 to 229,1.11 nohead lt 4 lw 10
set arrow from 230,1.09 to 252,1.09 nohead lt 1 lw 40
set arrow from 253,1.07 to 272,1.07 nohead lt 3 lw 10
set arrow from 273,1.09 to 295,1.09 nohead lt 1 lw 40
set arrow from 296,1.11 to 314,1.11 nohead lt 4 lw 10
set arrow from 315,1.09 to 337,1.09 nohead lt 1 lw 40
set arrow from 338,1.07 to 349,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_47|Q2G168|Putative|CP002114.3|tpos:103212-103560"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:349]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_47_Q2G168_Putative_CP002114.3_tpos_103212-103560.eps"
plot "./TMHMM_1211751/antigen_47_Q2G168_Putative_CP002114.3_tpos_103212-103560.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
