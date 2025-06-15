set arrow from 1,1.11 to 187,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_184|A7X1D2|Cell|CP000253.1|tpos:340368-340554"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:187]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_184_A7X1D2_Cell_CP000253.1_tpos_340368-340554.eps"
plot "./TMHMM_3187502/antigen_184_A7X1D2_Cell_CP000253.1_tpos_340368-340554.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
