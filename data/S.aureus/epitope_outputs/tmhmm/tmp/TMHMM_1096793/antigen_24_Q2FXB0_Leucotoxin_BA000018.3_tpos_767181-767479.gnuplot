set arrow from 1,1.11 to 299,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_24|Q2FXB0|Leucotoxin|BA000018.3|tpos:767181-767479"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:299]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_24_Q2FXB0_Leucotoxin_BA000018.3_tpos_767181-767479.eps"
plot "./TMHMM_1096793/antigen_24_Q2FXB0_Leucotoxin_BA000018.3_tpos_767181-767479.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
