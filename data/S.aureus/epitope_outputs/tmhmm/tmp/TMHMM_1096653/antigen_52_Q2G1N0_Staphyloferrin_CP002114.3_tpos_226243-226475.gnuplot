set arrow from 1,1.07 to 1,1.07 nohead lt 3 lw 10
set arrow from 2,1.09 to 24,1.09 nohead lt 1 lw 40
set arrow from 25,1.11 to 33,1.11 nohead lt 4 lw 10
set arrow from 34,1.09 to 56,1.09 nohead lt 1 lw 40
set arrow from 57,1.07 to 62,1.07 nohead lt 3 lw 10
set arrow from 63,1.09 to 85,1.09 nohead lt 1 lw 40
set arrow from 86,1.11 to 88,1.11 nohead lt 4 lw 10
set arrow from 89,1.09 to 111,1.09 nohead lt 1 lw 40
set arrow from 112,1.07 to 122,1.07 nohead lt 3 lw 10
set arrow from 123,1.09 to 145,1.09 nohead lt 1 lw 40
set arrow from 146,1.11 to 149,1.11 nohead lt 4 lw 10
set arrow from 150,1.09 to 172,1.09 nohead lt 1 lw 40
set arrow from 173,1.07 to 192,1.07 nohead lt 3 lw 10
set arrow from 193,1.09 to 215,1.09 nohead lt 1 lw 40
set arrow from 216,1.11 to 233,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_52|Q2G1N0|Staphyloferrin|CP002114.3|tpos:226243-226475"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:233]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_52_Q2G1N0_Staphyloferrin_CP002114.3_tpos_226243-226475.eps"
plot "./TMHMM_1096653/antigen_52_Q2G1N0_Staphyloferrin_CP002114.3_tpos_226243-226475.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
