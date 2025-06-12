set arrow from 1,1.11 to 299,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_24|Q2FXB0|Leucotoxin|CP002114.3|tpos:772608-772906"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:299]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_24_Q2FXB0_Leucotoxin_CP002114.3_tpos_772608-772906.eps"
plot "./TMHMM_1211751/antigen_24_Q2FXB0_Leucotoxin_CP002114.3_tpos_772608-772906.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
