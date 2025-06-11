set arrow from 1,1.11 to 9,1.11 nohead lt 4 lw 10
set arrow from 10,1.09 to 29,1.09 nohead lt 1 lw 40
set arrow from 30,1.07 to 41,1.07 nohead lt 3 lw 10
set arrow from 42,1.09 to 64,1.09 nohead lt 1 lw 40
set arrow from 65,1.11 to 68,1.11 nohead lt 4 lw 10
set arrow from 69,1.09 to 86,1.09 nohead lt 1 lw 40
set arrow from 87,1.07 to 101,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_10|Q2FVV4|CPBP|BA000018.3|tpos:641647-641747"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:101]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_10_Q2FVV4_CPBP_BA000018.3_tpos_641647-641747.eps"
plot "./TMHMM_1096793/antigen_10_Q2FVV4_CPBP_BA000018.3_tpos_641647-641747.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
