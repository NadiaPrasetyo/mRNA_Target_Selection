set arrow from 1,1.11 to 382,1.11 nohead lt 4 lw 10
set arrow from 383,1.09 to 405,1.09 nohead lt 1 lw 40
set arrow from 406,1.07 to 411,1.07 nohead lt 3 lw 10
set arrow from 412,1.09 to 430,1.09 nohead lt 1 lw 40
set arrow from 431,1.11 to 433,1.11 nohead lt 4 lw 10
set arrow from 434,1.09 to 453,1.09 nohead lt 1 lw 40
set arrow from 454,1.07 to 465,1.07 nohead lt 3 lw 10
set arrow from 466,1.09 to 484,1.09 nohead lt 1 lw 40
set arrow from 485,1.11 to 885,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_33|Q2FWH7|Sensor|CP000253.1|tpos:663567-664451"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:885]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_33_Q2FWH7_Sensor_CP000253.1_tpos_663567-664451.eps"
plot "./TMHMM_3187502/antigen_33_Q2FWH7_Sensor_CP000253.1_tpos_663567-664451.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
