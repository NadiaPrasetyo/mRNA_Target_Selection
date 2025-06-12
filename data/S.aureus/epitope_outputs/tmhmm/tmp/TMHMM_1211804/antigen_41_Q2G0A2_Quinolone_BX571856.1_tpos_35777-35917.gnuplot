set arrow from 1,1.11 to 4,1.11 nohead lt 4 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.07 to 39,1.07 nohead lt 3 lw 10
set arrow from 40,1.09 to 62,1.09 nohead lt 1 lw 40
set arrow from 63,1.11 to 65,1.11 nohead lt 4 lw 10
set arrow from 66,1.09 to 83,1.09 nohead lt 1 lw 40
set arrow from 84,1.07 to 94,1.07 nohead lt 3 lw 10
set arrow from 95,1.09 to 117,1.09 nohead lt 1 lw 40
set arrow from 118,1.11 to 120,1.11 nohead lt 4 lw 10
set arrow from 121,1.09 to 140,1.09 nohead lt 1 lw 40
set arrow from 141,1.07 to 141,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_41|Q2G0A2|Quinolone|BX571856.1|tpos:35777-35917"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:141]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211804/antigen_41_Q2G0A2_Quinolone_BX571856.1_tpos_35777-35917.eps"
plot "./TMHMM_1211804/antigen_41_Q2G0A2_Quinolone_BX571856.1_tpos_35777-35917.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
