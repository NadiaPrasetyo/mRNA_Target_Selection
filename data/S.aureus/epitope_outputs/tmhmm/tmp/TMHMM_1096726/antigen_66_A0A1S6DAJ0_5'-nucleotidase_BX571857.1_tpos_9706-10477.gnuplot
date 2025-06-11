set arrow from 1,1.11 to 772,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_66|A0A1S6DAJ0|5'-nucleotidase|BX571857.1|tpos:9706-10477"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:772]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096726/antigen_66_A0A1S6DAJ0_5'-nucleotidase_BX571857.1_tpos_9706-10477.eps"
plot "./TMHMM_1096726/antigen_66_A0A1S6DAJ0_5'-nucleotidase_BX571857.1_tpos_9706-10477.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
