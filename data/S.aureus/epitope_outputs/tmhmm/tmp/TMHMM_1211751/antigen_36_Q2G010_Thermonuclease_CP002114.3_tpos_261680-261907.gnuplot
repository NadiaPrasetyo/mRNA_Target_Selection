set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.11 to 36,1.11 nohead lt 4 lw 10
set arrow from 37,1.09 to 59,1.09 nohead lt 1 lw 40
set arrow from 60,1.07 to 228,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_36|Q2G010|Thermonuclease|CP002114.3|tpos:261680-261907"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:228]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_36_Q2G010_Thermonuclease_CP002114.3_tpos_261680-261907.eps"
plot "./TMHMM_1211751/antigen_36_Q2G010_Thermonuclease_CP002114.3_tpos_261680-261907.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
