set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 22,1.09 nohead lt 1 lw 40
set arrow from 23,1.07 to 28,1.07 nohead lt 3 lw 10
set arrow from 29,1.09 to 51,1.09 nohead lt 1 lw 40
set arrow from 52,1.11 to 73,1.11 nohead lt 4 lw 10
set arrow from 74,1.09 to 96,1.09 nohead lt 1 lw 40
set arrow from 97,1.07 to 104,1.07 nohead lt 3 lw 10
set arrow from 105,1.09 to 124,1.09 nohead lt 1 lw 40
set arrow from 125,1.11 to 127,1.11 nohead lt 4 lw 10
set arrow from 128,1.09 to 150,1.09 nohead lt 1 lw 40
set arrow from 151,1.07 to 156,1.07 nohead lt 3 lw 10
set arrow from 157,1.09 to 179,1.09 nohead lt 1 lw 40
set arrow from 180,1.11 to 202,1.11 nohead lt 4 lw 10
set arrow from 203,1.09 to 225,1.09 nohead lt 1 lw 40
set arrow from 226,1.07 to 231,1.07 nohead lt 3 lw 10
set arrow from 232,1.09 to 254,1.09 nohead lt 1 lw 40
set arrow from 255,1.11 to 266,1.11 nohead lt 4 lw 10
set arrow from 267,1.09 to 289,1.09 nohead lt 1 lw 40
set arrow from 290,1.07 to 295,1.07 nohead lt 3 lw 10
set arrow from 296,1.09 to 318,1.09 nohead lt 1 lw 40
set arrow from 319,1.11 to 322,1.11 nohead lt 4 lw 10
set arrow from 323,1.09 to 345,1.09 nohead lt 1 lw 40
set arrow from 346,1.07 to 364,1.07 nohead lt 3 lw 10
set arrow from 365,1.09 to 387,1.09 nohead lt 1 lw 40
set arrow from 388,1.11 to 389,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_14|A5IRD0|Na(+)/H(+)|BX571857.1|tpos:206559-206947"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:389]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_14_A5IRD0_Na(+)_H(+)_BX571857.1_tpos_206559-206947.eps"
plot "./TMHMM_3187503/antigen_14_A5IRD0_Na(+)_H(+)_BX571857.1_tpos_206559-206947.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
