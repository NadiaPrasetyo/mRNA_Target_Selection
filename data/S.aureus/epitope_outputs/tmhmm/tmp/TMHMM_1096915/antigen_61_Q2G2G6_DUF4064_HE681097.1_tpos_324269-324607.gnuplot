set arrow from 1,1.07 to 19,1.07 nohead lt 3 lw 10
set arrow from 20,1.09 to 42,1.09 nohead lt 1 lw 40
set arrow from 43,1.11 to 83,1.11 nohead lt 4 lw 10
set arrow from 84,1.09 to 106,1.09 nohead lt 1 lw 40
set arrow from 107,1.07 to 112,1.07 nohead lt 3 lw 10
set arrow from 113,1.09 to 147,1.09 nohead lt 1 lw 40
set arrow from 148,1.11 to 339,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_61|Q2G2G6|DUF4064|HE681097.1|tpos:324269-324607"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:339]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096915/antigen_61_Q2G2G6_DUF4064_HE681097.1_tpos_324269-324607.eps"
plot "./TMHMM_1096915/antigen_61_Q2G2G6_DUF4064_HE681097.1_tpos_324269-324607.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
