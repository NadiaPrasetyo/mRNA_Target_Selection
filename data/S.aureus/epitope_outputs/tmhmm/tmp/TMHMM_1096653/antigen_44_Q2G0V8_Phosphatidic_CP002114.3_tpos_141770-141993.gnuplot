set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 32,1.09 nohead lt 1 lw 40
set arrow from 33,1.11 to 73,1.11 nohead lt 4 lw 10
set arrow from 74,1.09 to 93,1.09 nohead lt 1 lw 40
set arrow from 94,1.07 to 97,1.07 nohead lt 3 lw 10
set arrow from 98,1.09 to 116,1.09 nohead lt 1 lw 40
set arrow from 117,1.11 to 143,1.11 nohead lt 4 lw 10
set arrow from 144,1.09 to 161,1.09 nohead lt 1 lw 40
set arrow from 162,1.07 to 167,1.07 nohead lt 3 lw 10
set arrow from 168,1.09 to 190,1.09 nohead lt 1 lw 40
set arrow from 191,1.11 to 194,1.11 nohead lt 4 lw 10
set arrow from 195,1.09 to 217,1.09 nohead lt 1 lw 40
set arrow from 218,1.07 to 224,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_44|Q2G0V8|Phosphatidic|CP002114.3|tpos:141770-141993"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:224]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_44_Q2G0V8_Phosphatidic_CP002114.3_tpos_141770-141993.eps"
plot "./TMHMM_1096653/antigen_44_Q2G0V8_Phosphatidic_CP002114.3_tpos_141770-141993.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
