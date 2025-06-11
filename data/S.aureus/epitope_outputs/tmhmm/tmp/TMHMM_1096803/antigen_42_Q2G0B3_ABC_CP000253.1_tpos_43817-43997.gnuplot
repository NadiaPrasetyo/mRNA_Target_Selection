set arrow from 1,1.11 to 181,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_42|Q2G0B3|ABC|CP000253.1|tpos:43817-43997"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:181]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096803/antigen_42_Q2G0B3_ABC_CP000253.1_tpos_43817-43997.eps"
plot "./TMHMM_1096803/antigen_42_Q2G0B3_ABC_CP000253.1_tpos_43817-43997.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
