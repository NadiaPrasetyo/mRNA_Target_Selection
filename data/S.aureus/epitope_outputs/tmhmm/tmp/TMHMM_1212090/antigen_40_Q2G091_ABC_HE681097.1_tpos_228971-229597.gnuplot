set arrow from 1,1.11 to 627,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_40|Q2G091|ABC|HE681097.1|tpos:228971-229597"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:627]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_40_Q2G091_ABC_HE681097.1_tpos_228971-229597.eps"
plot "./TMHMM_1212090/antigen_40_Q2G091_ABC_HE681097.1_tpos_228971-229597.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
