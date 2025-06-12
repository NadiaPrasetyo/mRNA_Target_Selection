set arrow from 1,1.07 to 26,1.07 nohead lt 3 lw 10
set arrow from 27,1.09 to 49,1.09 nohead lt 1 lw 40
set arrow from 50,1.11 to 52,1.11 nohead lt 4 lw 10
set arrow from 53,1.09 to 70,1.09 nohead lt 1 lw 40
set arrow from 71,1.07 to 76,1.07 nohead lt 3 lw 10
set arrow from 77,1.09 to 94,1.09 nohead lt 1 lw 40
set arrow from 95,1.11 to 96,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_10|Q2FVV4|CPBP|HE681097.1|tpos:789458-789553"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:96]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_10_Q2FVV4_CPBP_HE681097.1_tpos_789458-789553.eps"
plot "./TMHMM_1212090/antigen_10_Q2FVV4_CPBP_HE681097.1_tpos_789458-789553.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
