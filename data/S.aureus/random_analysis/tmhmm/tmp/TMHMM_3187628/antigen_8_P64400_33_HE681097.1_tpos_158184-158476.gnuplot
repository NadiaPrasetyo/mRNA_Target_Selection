set arrow from 1,1.11 to 293,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_8|P64400|33|HE681097.1|tpos:158184-158476"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:293]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_8_P64400_33_HE681097.1_tpos_158184-158476.eps"
plot "./TMHMM_3187628/antigen_8_P64400_33_HE681097.1_tpos_158184-158476.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
