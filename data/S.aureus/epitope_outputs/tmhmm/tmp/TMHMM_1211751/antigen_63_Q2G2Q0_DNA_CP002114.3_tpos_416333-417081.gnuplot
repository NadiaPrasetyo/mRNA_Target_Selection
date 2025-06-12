set arrow from 1,1.11 to 749,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_63|Q2G2Q0|DNA|CP002114.3|tpos:416333-417081"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:749]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_63_Q2G2Q0_DNA_CP002114.3_tpos_416333-417081.eps"
plot "./TMHMM_1211751/antigen_63_Q2G2Q0_DNA_CP002114.3_tpos_416333-417081.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
