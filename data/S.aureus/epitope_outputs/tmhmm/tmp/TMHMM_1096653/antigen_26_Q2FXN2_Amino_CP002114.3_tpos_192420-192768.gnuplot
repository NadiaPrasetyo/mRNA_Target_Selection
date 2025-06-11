set arrow from 1,1.07 to 8,1.07 nohead lt 3 lw 10
set arrow from 9,1.09 to 28,1.09 nohead lt 1 lw 40
set arrow from 29,1.11 to 37,1.11 nohead lt 4 lw 10
set arrow from 38,1.09 to 60,1.09 nohead lt 1 lw 40
set arrow from 61,1.07 to 80,1.07 nohead lt 3 lw 10
set arrow from 81,1.09 to 103,1.09 nohead lt 1 lw 40
set arrow from 104,1.11 to 125,1.11 nohead lt 4 lw 10
set arrow from 126,1.09 to 148,1.09 nohead lt 1 lw 40
set arrow from 149,1.07 to 160,1.07 nohead lt 3 lw 10
set arrow from 161,1.09 to 183,1.09 nohead lt 1 lw 40
set arrow from 184,1.11 to 197,1.11 nohead lt 4 lw 10
set arrow from 198,1.09 to 217,1.09 nohead lt 1 lw 40
set arrow from 218,1.07 to 236,1.07 nohead lt 3 lw 10
set arrow from 237,1.09 to 259,1.09 nohead lt 1 lw 40
set arrow from 260,1.11 to 278,1.11 nohead lt 4 lw 10
set arrow from 279,1.09 to 301,1.09 nohead lt 1 lw 40
set arrow from 302,1.07 to 349,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_26|Q2FXN2|Amino|CP002114.3|tpos:192420-192768"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:349]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_26_Q2FXN2_Amino_CP002114.3_tpos_192420-192768.eps"
plot "./TMHMM_1096653/antigen_26_Q2FXN2_Amino_CP002114.3_tpos_192420-192768.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
