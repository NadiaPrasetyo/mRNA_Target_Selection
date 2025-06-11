set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 29,1.09 nohead lt 1 lw 40
set arrow from 30,1.11 to 32,1.11 nohead lt 4 lw 10
set arrow from 33,1.09 to 55,1.09 nohead lt 1 lw 40
set arrow from 56,1.07 to 56,1.07 nohead lt 3 lw 10
set arrow from 57,1.09 to 76,1.09 nohead lt 1 lw 40
set arrow from 77,1.11 to 90,1.11 nohead lt 4 lw 10
set arrow from 91,1.09 to 113,1.09 nohead lt 1 lw 40
set arrow from 114,1.07 to 125,1.07 nohead lt 3 lw 10
set arrow from 126,1.09 to 148,1.09 nohead lt 1 lw 40
set arrow from 149,1.11 to 162,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_5|Q2FVH8|Membrane|BA000018.3|tpos:831524-831685"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:162]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_5_Q2FVH8_Membrane_BA000018.3_tpos_831524-831685.eps"
plot "./TMHMM_1096793/antigen_5_Q2FVH8_Membrane_BA000018.3_tpos_831524-831685.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
