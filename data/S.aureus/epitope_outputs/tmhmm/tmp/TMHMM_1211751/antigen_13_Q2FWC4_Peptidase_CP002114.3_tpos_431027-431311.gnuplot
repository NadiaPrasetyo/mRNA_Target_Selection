set arrow from 1,1.11 to 285,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_13|Q2FWC4|Peptidase|CP002114.3|tpos:431027-431311"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:285]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_13_Q2FWC4_Peptidase_CP002114.3_tpos_431027-431311.eps"
plot "./TMHMM_1211751/antigen_13_Q2FWC4_Peptidase_CP002114.3_tpos_431027-431311.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
