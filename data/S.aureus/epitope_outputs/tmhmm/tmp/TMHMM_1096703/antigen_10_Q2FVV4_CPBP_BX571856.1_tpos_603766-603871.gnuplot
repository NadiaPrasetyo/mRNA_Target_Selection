set arrow from 1,1.11 to 14,1.11 nohead lt 4 lw 10
set arrow from 15,1.09 to 37,1.09 nohead lt 1 lw 40
set arrow from 38,1.07 to 43,1.07 nohead lt 3 lw 10
set arrow from 44,1.09 to 61,1.09 nohead lt 1 lw 40
set arrow from 62,1.11 to 64,1.11 nohead lt 4 lw 10
set arrow from 65,1.09 to 84,1.09 nohead lt 1 lw 40
set arrow from 85,1.07 to 106,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_10|Q2FVV4|CPBP|BX571856.1|tpos:603766-603871"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:106]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096703/antigen_10_Q2FVV4_CPBP_BX571856.1_tpos_603766-603871.eps"
plot "./TMHMM_1096703/antigen_10_Q2FVV4_CPBP_BX571856.1_tpos_603766-603871.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
