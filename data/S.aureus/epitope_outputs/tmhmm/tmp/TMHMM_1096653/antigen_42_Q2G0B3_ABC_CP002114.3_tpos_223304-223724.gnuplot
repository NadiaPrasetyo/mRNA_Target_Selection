set arrow from 1,1.07 to 11,1.07 nohead lt 3 lw 10
set arrow from 12,1.09 to 34,1.09 nohead lt 1 lw 40
set arrow from 35,1.11 to 37,1.11 nohead lt 4 lw 10
set arrow from 38,1.09 to 60,1.09 nohead lt 1 lw 40
set arrow from 61,1.07 to 126,1.07 nohead lt 3 lw 10
set arrow from 127,1.09 to 149,1.09 nohead lt 1 lw 40
set arrow from 150,1.11 to 421,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_42|Q2G0B3|ABC|CP002114.3|tpos:223304-223724"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:421]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_42_Q2G0B3_ABC_CP002114.3_tpos_223304-223724.eps"
plot "./TMHMM_1096653/antigen_42_Q2G0B3_ABC_CP002114.3_tpos_223304-223724.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
