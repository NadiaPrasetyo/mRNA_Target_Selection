set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 21,1.09 nohead lt 1 lw 40
set arrow from 22,1.07 to 360,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_25|Q2FXJ6|Serine|BA000018.3|tpos:308669-309028"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:360]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_25_Q2FXJ6_Serine_BA000018.3_tpos_308669-309028.eps"
plot "./TMHMM_1211906/antigen_25_Q2FXJ6_Serine_BA000018.3_tpos_308669-309028.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
