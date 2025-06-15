set arrow from 1,1.11 to 450,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_32|Q5HFH6|Elongation|CP000253.1|tpos:286492-286941"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:450]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_32_Q5HFH6_Elongation_CP000253.1_tpos_286492-286941.eps"
plot "./TMHMM_3187502/antigen_32_Q5HFH6_Elongation_CP000253.1_tpos_286492-286941.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
