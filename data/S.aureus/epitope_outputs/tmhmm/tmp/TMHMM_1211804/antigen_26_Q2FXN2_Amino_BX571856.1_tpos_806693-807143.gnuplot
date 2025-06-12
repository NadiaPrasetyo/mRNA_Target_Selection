set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 35,1.09 nohead lt 1 lw 40
set arrow from 36,1.11 to 44,1.11 nohead lt 4 lw 10
set arrow from 45,1.09 to 67,1.09 nohead lt 1 lw 40
set arrow from 68,1.07 to 79,1.07 nohead lt 3 lw 10
set arrow from 80,1.09 to 102,1.09 nohead lt 1 lw 40
set arrow from 103,1.11 to 121,1.11 nohead lt 4 lw 10
set arrow from 122,1.09 to 144,1.09 nohead lt 1 lw 40
set arrow from 145,1.07 to 148,1.07 nohead lt 3 lw 10
set arrow from 149,1.09 to 171,1.09 nohead lt 1 lw 40
set arrow from 172,1.11 to 190,1.11 nohead lt 4 lw 10
set arrow from 191,1.09 to 213,1.09 nohead lt 1 lw 40
set arrow from 214,1.07 to 233,1.07 nohead lt 3 lw 10
set arrow from 234,1.09 to 256,1.09 nohead lt 1 lw 40
set arrow from 257,1.11 to 270,1.11 nohead lt 4 lw 10
set arrow from 271,1.09 to 293,1.09 nohead lt 1 lw 40
set arrow from 294,1.07 to 327,1.07 nohead lt 3 lw 10
set arrow from 328,1.09 to 350,1.09 nohead lt 1 lw 40
set arrow from 351,1.11 to 359,1.11 nohead lt 4 lw 10
set arrow from 360,1.09 to 382,1.09 nohead lt 1 lw 40
set arrow from 383,1.07 to 402,1.07 nohead lt 3 lw 10
set arrow from 403,1.09 to 420,1.09 nohead lt 1 lw 40
set arrow from 421,1.11 to 424,1.11 nohead lt 4 lw 10
set arrow from 425,1.09 to 442,1.09 nohead lt 1 lw 40
set arrow from 443,1.07 to 451,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_26|Q2FXN2|Amino|BX571856.1|tpos:806693-807143"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:451]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211804/antigen_26_Q2FXN2_Amino_BX571856.1_tpos_806693-807143.eps"
plot "./TMHMM_1211804/antigen_26_Q2FXN2_Amino_BX571856.1_tpos_806693-807143.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
