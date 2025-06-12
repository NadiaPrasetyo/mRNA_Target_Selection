set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 32,1.09 nohead lt 1 lw 40
set arrow from 33,1.11 to 69,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_2|Q2FUY2|Clumping|CP000253.1|tpos:97046-97114"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:69]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212023/antigen_2_Q2FUY2_Clumping_CP000253.1_tpos_97046-97114.eps"
plot "./TMHMM_1212023/antigen_2_Q2FUY2_Clumping_CP000253.1_tpos_97046-97114.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
