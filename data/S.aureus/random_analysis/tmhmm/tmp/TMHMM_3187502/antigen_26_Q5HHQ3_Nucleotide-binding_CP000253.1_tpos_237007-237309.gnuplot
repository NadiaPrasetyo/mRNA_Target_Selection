set arrow from 1,1.11 to 303,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_26|Q5HHQ3|Nucleotide-binding|CP000253.1|tpos:237007-237309"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:303]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_26_Q5HHQ3_Nucleotide-binding_CP000253.1_tpos_237007-237309.eps"
plot "./TMHMM_3187502/antigen_26_Q5HHQ3_Nucleotide-binding_CP000253.1_tpos_237007-237309.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
