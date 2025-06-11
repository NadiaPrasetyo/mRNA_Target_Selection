set arrow from 1,1.11 to 237,1.11 nohead lt 4 lw 10
set arrow from 238,1.09 to 255,1.09 nohead lt 1 lw 40
set arrow from 256,1.07 to 262,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_37|Q2G015|Clumping|BA000018.3|tpos:189488-189749"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:262]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_37_Q2G015_Clumping_BA000018.3_tpos_189488-189749.eps"
plot "./TMHMM_1096793/antigen_37_Q2G015_Clumping_BA000018.3_tpos_189488-189749.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
