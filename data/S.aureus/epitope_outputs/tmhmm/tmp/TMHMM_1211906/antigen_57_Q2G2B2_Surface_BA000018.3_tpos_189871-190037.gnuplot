set arrow from 1,1.07 to 19,1.07 nohead lt 3 lw 10
set arrow from 20,1.09 to 42,1.09 nohead lt 1 lw 40
set arrow from 43,1.11 to 167,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_57|Q2G2B2|Surface|BA000018.3|tpos:189871-190037"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:167]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_57_Q2G2B2_Surface_BA000018.3_tpos_189871-190037.eps"
plot "./TMHMM_1211906/antigen_57_Q2G2B2_Surface_BA000018.3_tpos_189871-190037.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
