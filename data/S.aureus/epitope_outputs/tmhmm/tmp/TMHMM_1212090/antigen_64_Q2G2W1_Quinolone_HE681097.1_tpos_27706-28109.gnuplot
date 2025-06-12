set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 24,1.09 nohead lt 1 lw 40
set arrow from 25,1.11 to 38,1.11 nohead lt 4 lw 10
set arrow from 39,1.09 to 61,1.09 nohead lt 1 lw 40
set arrow from 62,1.07 to 72,1.07 nohead lt 3 lw 10
set arrow from 73,1.09 to 95,1.09 nohead lt 1 lw 40
set arrow from 96,1.11 to 98,1.11 nohead lt 4 lw 10
set arrow from 99,1.09 to 116,1.09 nohead lt 1 lw 40
set arrow from 117,1.07 to 128,1.07 nohead lt 3 lw 10
set arrow from 129,1.09 to 151,1.09 nohead lt 1 lw 40
set arrow from 152,1.11 to 154,1.11 nohead lt 4 lw 10
set arrow from 155,1.09 to 174,1.09 nohead lt 1 lw 40
set arrow from 175,1.07 to 193,1.07 nohead lt 3 lw 10
set arrow from 194,1.09 to 216,1.09 nohead lt 1 lw 40
set arrow from 217,1.11 to 219,1.11 nohead lt 4 lw 10
set arrow from 220,1.09 to 237,1.09 nohead lt 1 lw 40
set arrow from 238,1.07 to 257,1.07 nohead lt 3 lw 10
set arrow from 258,1.09 to 278,1.09 nohead lt 1 lw 40
set arrow from 279,1.11 to 287,1.11 nohead lt 4 lw 10
set arrow from 288,1.09 to 310,1.09 nohead lt 1 lw 40
set arrow from 311,1.07 to 322,1.07 nohead lt 3 lw 10
set arrow from 323,1.09 to 342,1.09 nohead lt 1 lw 40
set arrow from 343,1.11 to 345,1.11 nohead lt 4 lw 10
set arrow from 346,1.09 to 368,1.09 nohead lt 1 lw 40
set arrow from 369,1.07 to 404,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_64|Q2G2W1|Quinolone|HE681097.1|tpos:27706-28109"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:404]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_64_Q2G2W1_Quinolone_HE681097.1_tpos_27706-28109.eps"
plot "./TMHMM_1212090/antigen_64_Q2G2W1_Quinolone_HE681097.1_tpos_27706-28109.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
