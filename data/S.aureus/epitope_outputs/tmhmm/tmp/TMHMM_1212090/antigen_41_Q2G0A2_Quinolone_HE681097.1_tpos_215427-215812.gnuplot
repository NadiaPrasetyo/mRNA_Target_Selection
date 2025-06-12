set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.11 to 36,1.11 nohead lt 4 lw 10
set arrow from 37,1.09 to 56,1.09 nohead lt 1 lw 40
set arrow from 57,1.07 to 68,1.07 nohead lt 3 lw 10
set arrow from 69,1.09 to 88,1.09 nohead lt 1 lw 40
set arrow from 89,1.11 to 91,1.11 nohead lt 4 lw 10
set arrow from 92,1.09 to 114,1.09 nohead lt 1 lw 40
set arrow from 115,1.07 to 133,1.07 nohead lt 3 lw 10
set arrow from 134,1.09 to 156,1.09 nohead lt 1 lw 40
set arrow from 157,1.11 to 159,1.11 nohead lt 4 lw 10
set arrow from 160,1.09 to 182,1.09 nohead lt 1 lw 40
set arrow from 183,1.07 to 208,1.07 nohead lt 3 lw 10
set arrow from 209,1.09 to 231,1.09 nohead lt 1 lw 40
set arrow from 232,1.11 to 245,1.11 nohead lt 4 lw 10
set arrow from 246,1.09 to 268,1.09 nohead lt 1 lw 40
set arrow from 269,1.07 to 274,1.07 nohead lt 3 lw 10
set arrow from 275,1.09 to 294,1.09 nohead lt 1 lw 40
set arrow from 295,1.11 to 299,1.11 nohead lt 4 lw 10
set arrow from 300,1.09 to 322,1.09 nohead lt 1 lw 40
set arrow from 323,1.07 to 328,1.07 nohead lt 3 lw 10
set arrow from 329,1.09 to 348,1.09 nohead lt 1 lw 40
set arrow from 349,1.11 to 357,1.11 nohead lt 4 lw 10
set arrow from 358,1.09 to 380,1.09 nohead lt 1 lw 40
set arrow from 381,1.07 to 386,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_41|Q2G0A2|Quinolone|HE681097.1|tpos:215427-215812"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:386]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_41_Q2G0A2_Quinolone_HE681097.1_tpos_215427-215812.eps"
plot "./TMHMM_1212090/antigen_41_Q2G0A2_Quinolone_HE681097.1_tpos_215427-215812.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
