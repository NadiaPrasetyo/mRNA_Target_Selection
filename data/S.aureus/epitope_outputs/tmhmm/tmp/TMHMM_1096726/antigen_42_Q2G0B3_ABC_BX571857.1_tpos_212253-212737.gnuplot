set arrow from 1,1.07 to 70,1.07 nohead lt 3 lw 10
set arrow from 71,1.09 to 90,1.09 nohead lt 1 lw 40
set arrow from 91,1.11 to 94,1.11 nohead lt 4 lw 10
set arrow from 95,1.09 to 117,1.09 nohead lt 1 lw 40
set arrow from 118,1.07 to 174,1.07 nohead lt 3 lw 10
set arrow from 175,1.09 to 197,1.09 nohead lt 1 lw 40
set arrow from 198,1.11 to 206,1.11 nohead lt 4 lw 10
set arrow from 207,1.09 to 229,1.09 nohead lt 1 lw 40
set arrow from 230,1.07 to 485,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_42|Q2G0B3|ABC|BX571857.1|tpos:212253-212737"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:485]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096726/antigen_42_Q2G0B3_ABC_BX571857.1_tpos_212253-212737.eps"
plot "./TMHMM_1096726/antigen_42_Q2G0B3_ABC_BX571857.1_tpos_212253-212737.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
