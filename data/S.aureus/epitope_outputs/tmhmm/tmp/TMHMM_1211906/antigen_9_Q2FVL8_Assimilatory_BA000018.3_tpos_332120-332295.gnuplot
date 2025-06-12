set arrow from 1,1.11 to 176,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_9|Q2FVL8|Assimilatory|BA000018.3|tpos:332120-332295"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:176]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_9_Q2FVL8_Assimilatory_BA000018.3_tpos_332120-332295.eps"
plot "./TMHMM_1211906/antigen_9_Q2FVL8_Assimilatory_BA000018.3_tpos_332120-332295.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
