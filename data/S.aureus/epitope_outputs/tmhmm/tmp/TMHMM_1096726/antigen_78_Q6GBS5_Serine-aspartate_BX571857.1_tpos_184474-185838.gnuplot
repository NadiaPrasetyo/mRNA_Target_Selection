set arrow from 1,1.11 to 1340,1.11 nohead lt 4 lw 10
set arrow from 1341,1.09 to 1358,1.09 nohead lt 1 lw 40
set arrow from 1359,1.07 to 1365,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_78|Q6GBS5|Serine-aspartate|BX571857.1|tpos:184474-185838"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:1365]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096726/antigen_78_Q6GBS5_Serine-aspartate_BX571857.1_tpos_184474-185838.eps"
plot "./TMHMM_1096726/antigen_78_Q6GBS5_Serine-aspartate_BX571857.1_tpos_184474-185838.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
