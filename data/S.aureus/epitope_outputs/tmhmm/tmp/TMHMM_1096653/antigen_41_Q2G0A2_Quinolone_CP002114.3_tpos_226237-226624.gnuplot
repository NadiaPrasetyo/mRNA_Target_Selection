set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 29,1.09 nohead lt 1 lw 40
set arrow from 30,1.11 to 39,1.11 nohead lt 4 lw 10
set arrow from 40,1.09 to 62,1.09 nohead lt 1 lw 40
set arrow from 63,1.07 to 68,1.07 nohead lt 3 lw 10
set arrow from 69,1.09 to 91,1.09 nohead lt 1 lw 40
set arrow from 92,1.11 to 94,1.11 nohead lt 4 lw 10
set arrow from 95,1.09 to 117,1.09 nohead lt 1 lw 40
set arrow from 118,1.07 to 128,1.07 nohead lt 3 lw 10
set arrow from 129,1.09 to 151,1.09 nohead lt 1 lw 40
set arrow from 152,1.11 to 155,1.11 nohead lt 4 lw 10
set arrow from 156,1.09 to 178,1.09 nohead lt 1 lw 40
set arrow from 179,1.07 to 202,1.07 nohead lt 3 lw 10
set arrow from 203,1.09 to 225,1.09 nohead lt 1 lw 40
set arrow from 226,1.11 to 239,1.11 nohead lt 4 lw 10
set arrow from 240,1.09 to 259,1.09 nohead lt 1 lw 40
set arrow from 260,1.07 to 265,1.07 nohead lt 3 lw 10
set arrow from 266,1.09 to 288,1.09 nohead lt 1 lw 40
set arrow from 289,1.11 to 291,1.11 nohead lt 4 lw 10
set arrow from 292,1.09 to 314,1.09 nohead lt 1 lw 40
set arrow from 315,1.07 to 326,1.07 nohead lt 3 lw 10
set arrow from 327,1.09 to 349,1.09 nohead lt 1 lw 40
set arrow from 350,1.11 to 352,1.11 nohead lt 4 lw 10
set arrow from 353,1.09 to 375,1.09 nohead lt 1 lw 40
set arrow from 376,1.07 to 388,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_41|Q2G0A2|Quinolone|CP002114.3|tpos:226237-226624"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:388]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_41_Q2G0A2_Quinolone_CP002114.3_tpos_226237-226624.eps"
plot "./TMHMM_1096653/antigen_41_Q2G0A2_Quinolone_CP002114.3_tpos_226237-226624.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
