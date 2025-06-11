set arrow from 1,1.07 to 26,1.07 nohead lt 3 lw 10
set arrow from 27,1.09 to 49,1.09 nohead lt 1 lw 40
set arrow from 50,1.11 to 424,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_25|Q2FXJ6|Serine|BX571856.1|tpos:581539-581962"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:424]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096703/antigen_25_Q2FXJ6_Serine_BX571856.1_tpos_581539-581962.eps"
plot "./TMHMM_1096703/antigen_25_Q2FXJ6_Serine_BX571856.1_tpos_581539-581962.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
