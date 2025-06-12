set arrow from 1,1.11 to 246,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_55|Q2G1X0|Alpha-hemolysin|CP002114.3|tpos:773449-773694"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:246]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_55_Q2G1X0_Alpha-hemolysin_CP002114.3_tpos_773449-773694.eps"
plot "./TMHMM_1211751/antigen_55_Q2G1X0_Alpha-hemolysin_CP002114.3_tpos_773449-773694.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
