set arrow from 1,1.07 to 1,1.07 nohead lt 3 lw 10
set arrow from 2,1.09 to 16,1.09 nohead lt 1 lw 40
set arrow from 17,1.11 to 25,1.11 nohead lt 4 lw 10
set arrow from 26,1.09 to 45,1.09 nohead lt 1 lw 40
set arrow from 46,1.07 to 57,1.07 nohead lt 3 lw 10
set arrow from 58,1.09 to 80,1.09 nohead lt 1 lw 40
set arrow from 81,1.11 to 99,1.11 nohead lt 4 lw 10
set arrow from 100,1.09 to 118,1.09 nohead lt 1 lw 40
set arrow from 119,1.07 to 157,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_91|Q6GIE0|Na(+)/H(+)|HE681097.1|tpos:200454-200610"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:157]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_91_Q6GIE0_Na(+)_H(+)_HE681097.1_tpos_200454-200610.eps"
plot "./TMHMM_3187628/antigen_91_Q6GIE0_Na(+)_H(+)_HE681097.1_tpos_200454-200610.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
