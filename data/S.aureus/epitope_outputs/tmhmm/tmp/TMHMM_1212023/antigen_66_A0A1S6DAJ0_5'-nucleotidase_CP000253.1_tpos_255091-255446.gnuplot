set arrow from 1,1.11 to 356,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_66|A0A1S6DAJ0|5'-nucleotidase|CP000253.1|tpos:255091-255446"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:356]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212023/antigen_66_A0A1S6DAJ0_5'-nucleotidase_CP000253.1_tpos_255091-255446.eps"
plot "./TMHMM_1212023/antigen_66_A0A1S6DAJ0_5'-nucleotidase_CP000253.1_tpos_255091-255446.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
