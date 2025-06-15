set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 29,1.09 nohead lt 1 lw 40
set arrow from 30,1.11 to 38,1.11 nohead lt 4 lw 10
set arrow from 39,1.09 to 61,1.09 nohead lt 1 lw 40
set arrow from 62,1.07 to 73,1.07 nohead lt 3 lw 10
set arrow from 74,1.09 to 96,1.09 nohead lt 1 lw 40
set arrow from 97,1.11 to 110,1.11 nohead lt 4 lw 10
set arrow from 111,1.09 to 130,1.09 nohead lt 1 lw 40
set arrow from 131,1.07 to 136,1.07 nohead lt 3 lw 10
set arrow from 137,1.09 to 159,1.09 nohead lt 1 lw 40
set arrow from 160,1.11 to 173,1.11 nohead lt 4 lw 10
set arrow from 174,1.09 to 196,1.09 nohead lt 1 lw 40
set arrow from 197,1.07 to 211,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_14|A5IRD0|Na(+)/H(+)|CP002114.3|tpos:140117-140327"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:211]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_14_A5IRD0_Na(+)_H(+)_CP002114.3_tpos_140117-140327.eps"
plot "./TMHMM_3187290/antigen_14_A5IRD0_Na(+)_H(+)_CP002114.3_tpos_140117-140327.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
