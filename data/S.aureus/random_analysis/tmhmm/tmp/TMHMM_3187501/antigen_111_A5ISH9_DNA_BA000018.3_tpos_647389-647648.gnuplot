set arrow from 1,1.11 to 260,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_111|A5ISH9|DNA|BA000018.3|tpos:647389-647648"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:260]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187501/antigen_111_A5ISH9_DNA_BA000018.3_tpos_647389-647648.eps"
plot "./TMHMM_3187501/antigen_111_A5ISH9_DNA_BA000018.3_tpos_647389-647648.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
