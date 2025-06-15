set arrow from 1,1.07 to 1,1.07 nohead lt 3 lw 10
set arrow from 2,1.09 to 20,1.09 nohead lt 1 lw 40
set arrow from 21,1.11 to 29,1.11 nohead lt 4 lw 10
set arrow from 30,1.09 to 52,1.09 nohead lt 1 lw 40
set arrow from 53,1.07 to 60,1.07 nohead lt 3 lw 10
set arrow from 61,1.09 to 83,1.09 nohead lt 1 lw 40
set arrow from 84,1.11 to 86,1.11 nohead lt 4 lw 10
set arrow from 87,1.09 to 109,1.09 nohead lt 1 lw 40
set arrow from 110,1.07 to 110,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_182|Q2FFV9|Fluoride-specific|HE681097.1|tpos:579511-579620"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:110]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_182_Q2FFV9_Fluoride-specific_HE681097.1_tpos_579511-579620.eps"
plot "./TMHMM_3187628/antigen_182_Q2FFV9_Fluoride-specific_HE681097.1_tpos_579511-579620.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
