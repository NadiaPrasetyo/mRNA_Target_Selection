set arrow from 1,1.07 to 10,1.07 nohead lt 3 lw 10
set arrow from 11,1.09 to 33,1.09 nohead lt 1 lw 40
set arrow from 34,1.11 to 42,1.11 nohead lt 4 lw 10
set arrow from 43,1.09 to 65,1.09 nohead lt 1 lw 40
set arrow from 66,1.07 to 71,1.07 nohead lt 3 lw 10
set arrow from 72,1.09 to 94,1.09 nohead lt 1 lw 40
set arrow from 95,1.11 to 97,1.11 nohead lt 4 lw 10
set arrow from 98,1.09 to 120,1.09 nohead lt 1 lw 40
set arrow from 121,1.07 to 131,1.07 nohead lt 3 lw 10
set arrow from 132,1.09 to 154,1.09 nohead lt 1 lw 40
set arrow from 155,1.11 to 159,1.11 nohead lt 4 lw 10
set arrow from 160,1.09 to 177,1.09 nohead lt 1 lw 40
set arrow from 178,1.07 to 178,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_64|Q2G2W1|Quinolone|CP002114.3|tpos:226234-226411"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:178]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_64_Q2G2W1_Quinolone_CP002114.3_tpos_226234-226411.eps"
plot "./TMHMM_1211751/antigen_64_Q2G2W1_Quinolone_CP002114.3_tpos_226234-226411.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
