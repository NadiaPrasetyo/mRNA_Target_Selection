set arrow from 1,1.11 to 228,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_85|A0A2S1FUJ3|Superantigen-like|CP002114.3|tpos:133080-133307"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:228]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_85_A0A2S1FUJ3_Superantigen-like_CP002114.3_tpos_133080-133307.eps"
plot "./TMHMM_1211751/antigen_85_A0A2S1FUJ3_Superantigen-like_CP002114.3_tpos_133080-133307.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
