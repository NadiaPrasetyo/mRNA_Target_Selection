set arrow from 1,1.11 to 220,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_37|Q2YY04|Signal|CP000253.1|tpos:8300-8519"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:220]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_37_Q2YY04_Signal_CP000253.1_tpos_8300-8519.eps"
plot "./TMHMM_3187502/antigen_37_Q2YY04_Signal_CP000253.1_tpos_8300-8519.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
