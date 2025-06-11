set arrow from 1,1.11 to 9,1.11 nohead lt 4 lw 10
set arrow from 10,1.09 to 32,1.09 nohead lt 1 lw 40
set arrow from 33,1.07 to 38,1.07 nohead lt 3 lw 10
set arrow from 39,1.09 to 58,1.09 nohead lt 1 lw 40
set arrow from 59,1.11 to 61,1.11 nohead lt 4 lw 10
set arrow from 62,1.09 to 84,1.09 nohead lt 1 lw 40
set arrow from 85,1.07 to 90,1.07 nohead lt 3 lw 10
set arrow from 91,1.09 to 113,1.09 nohead lt 1 lw 40
set arrow from 114,1.11 to 127,1.11 nohead lt 4 lw 10
set arrow from 128,1.09 to 150,1.09 nohead lt 1 lw 40
set arrow from 151,1.07 to 162,1.07 nohead lt 3 lw 10
set arrow from 163,1.09 to 185,1.09 nohead lt 1 lw 40
set arrow from 186,1.11 to 204,1.11 nohead lt 4 lw 10
set arrow from 205,1.09 to 223,1.09 nohead lt 1 lw 40
set arrow from 224,1.07 to 243,1.07 nohead lt 3 lw 10
set arrow from 244,1.09 to 266,1.09 nohead lt 1 lw 40
set arrow from 267,1.11 to 285,1.11 nohead lt 4 lw 10
set arrow from 286,1.09 to 308,1.09 nohead lt 1 lw 40
set arrow from 309,1.07 to 341,1.07 nohead lt 3 lw 10
set arrow from 342,1.09 to 361,1.09 nohead lt 1 lw 40
set arrow from 362,1.11 to 364,1.11 nohead lt 4 lw 10
set arrow from 365,1.09 to 382,1.09 nohead lt 1 lw 40
set arrow from 383,1.07 to 398,1.07 nohead lt 3 lw 10
set arrow from 399,1.09 to 421,1.09 nohead lt 1 lw 40
set arrow from 422,1.11 to 425,1.11 nohead lt 4 lw 10
set arrow from 426,1.09 to 444,1.09 nohead lt 1 lw 40
set arrow from 445,1.07 to 450,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_26|Q2FXN2|Amino|BA000018.3|tpos:831487-831936"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:450]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_26_Q2FXN2_Amino_BA000018.3_tpos_831487-831936.eps"
plot "./TMHMM_1096793/antigen_26_Q2FXN2_Amino_BA000018.3_tpos_831487-831936.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
