set arrow from 1,1.07 to 11,1.07 nohead lt 3 lw 10
set arrow from 12,1.09 to 34,1.09 nohead lt 1 lw 40
set arrow from 35,1.11 to 48,1.11 nohead lt 4 lw 10
set arrow from 49,1.09 to 71,1.09 nohead lt 1 lw 40
set arrow from 72,1.07 to 123,1.07 nohead lt 3 lw 10
set arrow from 124,1.09 to 141,1.09 nohead lt 1 lw 40
set arrow from 142,1.11 to 145,1.11 nohead lt 4 lw 10
set arrow from 146,1.09 to 165,1.09 nohead lt 1 lw 40
set arrow from 166,1.07 to 225,1.07 nohead lt 3 lw 10
set arrow from 226,1.09 to 248,1.09 nohead lt 1 lw 40
set arrow from 249,1.11 to 542,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_42|Q2G0B3|ABC|BX571856.1|tpos:237837-238378"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:542]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096703/antigen_42_Q2G0B3_ABC_BX571856.1_tpos_237837-238378.eps"
plot "./TMHMM_1096703/antigen_42_Q2G0B3_ABC_BX571856.1_tpos_237837-238378.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
