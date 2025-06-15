set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 23,1.09 nohead lt 1 lw 40
set arrow from 24,1.07 to 29,1.07 nohead lt 3 lw 10
set arrow from 30,1.09 to 48,1.09 nohead lt 1 lw 40
set arrow from 49,1.11 to 51,1.11 nohead lt 4 lw 10
set arrow from 52,1.09 to 69,1.09 nohead lt 1 lw 40
set arrow from 70,1.07 to 81,1.07 nohead lt 3 lw 10
set arrow from 82,1.09 to 104,1.09 nohead lt 1 lw 40
set arrow from 105,1.11 to 113,1.11 nohead lt 4 lw 10
set arrow from 114,1.09 to 133,1.09 nohead lt 1 lw 40
set arrow from 134,1.07 to 149,1.07 nohead lt 3 lw 10
set arrow from 150,1.09 to 172,1.09 nohead lt 1 lw 40
set arrow from 173,1.11 to 176,1.11 nohead lt 4 lw 10
set arrow from 177,1.09 to 195,1.09 nohead lt 1 lw 40
set arrow from 196,1.07 to 201,1.07 nohead lt 3 lw 10
set arrow from 202,1.09 to 224,1.09 nohead lt 1 lw 40
set arrow from 225,1.11 to 228,1.11 nohead lt 4 lw 10
set arrow from 229,1.09 to 251,1.09 nohead lt 1 lw 40
set arrow from 252,1.07 to 274,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_86|P60946|Putative|BX571857.1|tpos:713407-713680"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:274]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_86_P60946_Putative_BX571857.1_tpos_713407-713680.eps"
plot "./TMHMM_3187503/antigen_86_P60946_Putative_BX571857.1_tpos_713407-713680.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
