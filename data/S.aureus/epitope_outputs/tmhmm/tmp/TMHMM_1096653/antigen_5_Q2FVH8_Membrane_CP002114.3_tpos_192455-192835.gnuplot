set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 26,1.09 nohead lt 1 lw 40
set arrow from 27,1.07 to 46,1.07 nohead lt 3 lw 10
set arrow from 47,1.09 to 69,1.09 nohead lt 1 lw 40
set arrow from 70,1.11 to 90,1.11 nohead lt 4 lw 10
set arrow from 91,1.09 to 113,1.09 nohead lt 1 lw 40
set arrow from 114,1.07 to 125,1.07 nohead lt 3 lw 10
set arrow from 126,1.09 to 148,1.09 nohead lt 1 lw 40
set arrow from 149,1.11 to 162,1.11 nohead lt 4 lw 10
set arrow from 163,1.09 to 182,1.09 nohead lt 1 lw 40
set arrow from 183,1.07 to 201,1.07 nohead lt 3 lw 10
set arrow from 202,1.09 to 224,1.09 nohead lt 1 lw 40
set arrow from 225,1.11 to 243,1.11 nohead lt 4 lw 10
set arrow from 244,1.09 to 266,1.09 nohead lt 1 lw 40
set arrow from 267,1.07 to 298,1.07 nohead lt 3 lw 10
set arrow from 299,1.09 to 318,1.09 nohead lt 1 lw 40
set arrow from 319,1.11 to 327,1.11 nohead lt 4 lw 10
set arrow from 328,1.09 to 350,1.09 nohead lt 1 lw 40
set arrow from 351,1.07 to 356,1.07 nohead lt 3 lw 10
set arrow from 357,1.09 to 379,1.09 nohead lt 1 lw 40
set arrow from 380,1.11 to 381,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_5|Q2FVH8|Membrane|CP002114.3|tpos:192455-192835"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:381]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_5_Q2FVH8_Membrane_CP002114.3_tpos_192455-192835.eps"
plot "./TMHMM_1096653/antigen_5_Q2FVH8_Membrane_CP002114.3_tpos_192455-192835.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
