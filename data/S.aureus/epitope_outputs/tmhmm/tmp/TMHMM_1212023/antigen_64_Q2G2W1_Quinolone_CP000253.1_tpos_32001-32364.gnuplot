set arrow from 1,1.11 to 14,1.11 nohead lt 4 lw 10
set arrow from 15,1.09 to 34,1.09 nohead lt 1 lw 40
set arrow from 35,1.07 to 40,1.07 nohead lt 3 lw 10
set arrow from 41,1.09 to 63,1.09 nohead lt 1 lw 40
set arrow from 64,1.11 to 66,1.11 nohead lt 4 lw 10
set arrow from 67,1.09 to 89,1.09 nohead lt 1 lw 40
set arrow from 90,1.07 to 100,1.07 nohead lt 3 lw 10
set arrow from 101,1.09 to 123,1.09 nohead lt 1 lw 40
set arrow from 124,1.11 to 127,1.11 nohead lt 4 lw 10
set arrow from 128,1.09 to 150,1.09 nohead lt 1 lw 40
set arrow from 151,1.07 to 162,1.07 nohead lt 3 lw 10
set arrow from 163,1.09 to 180,1.09 nohead lt 1 lw 40
set arrow from 181,1.11 to 184,1.11 nohead lt 4 lw 10
set arrow from 185,1.09 to 202,1.09 nohead lt 1 lw 40
set arrow from 203,1.07 to 222,1.07 nohead lt 3 lw 10
set arrow from 223,1.09 to 245,1.09 nohead lt 1 lw 40
set arrow from 246,1.11 to 254,1.11 nohead lt 4 lw 10
set arrow from 255,1.09 to 274,1.09 nohead lt 1 lw 40
set arrow from 275,1.07 to 285,1.07 nohead lt 3 lw 10
set arrow from 286,1.09 to 305,1.09 nohead lt 1 lw 40
set arrow from 306,1.11 to 308,1.11 nohead lt 4 lw 10
set arrow from 309,1.09 to 331,1.09 nohead lt 1 lw 40
set arrow from 332,1.07 to 343,1.07 nohead lt 3 lw 10
set arrow from 344,1.09 to 363,1.09 nohead lt 1 lw 40
set arrow from 364,1.11 to 364,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_64|Q2G2W1|Quinolone|CP000253.1|tpos:32001-32364"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:364]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212023/antigen_64_Q2G2W1_Quinolone_CP000253.1_tpos_32001-32364.eps"
plot "./TMHMM_1212023/antigen_64_Q2G2W1_Quinolone_CP000253.1_tpos_32001-32364.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
