set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 26,1.09 nohead lt 1 lw 40
set arrow from 27,1.07 to 27,1.07 nohead lt 3 lw 10
set arrow from 28,1.09 to 50,1.09 nohead lt 1 lw 40
set arrow from 51,1.11 to 69,1.11 nohead lt 4 lw 10
set arrow from 70,1.09 to 92,1.09 nohead lt 1 lw 40
set arrow from 93,1.07 to 112,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_64|Q6GID8|Na(+)/H(+)|BX571856.1|tpos:209397-209508"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:112]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_64_Q6GID8_Na(+)_H(+)_BX571856.1_tpos_209397-209508.eps"
plot "./TMHMM_3187473/antigen_64_Q6GID8_Na(+)_H(+)_BX571856.1_tpos_209397-209508.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
