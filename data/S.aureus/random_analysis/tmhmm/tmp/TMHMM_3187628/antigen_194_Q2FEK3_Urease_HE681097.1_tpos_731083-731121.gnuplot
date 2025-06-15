set arrow from 1,1.07 to 39,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_194|Q2FEK3|Urease|HE681097.1|tpos:731083-731121"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:39]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_194_Q2FEK3_Urease_HE681097.1_tpos_731083-731121.eps"
plot "./TMHMM_3187628/antigen_194_Q2FEK3_Urease_HE681097.1_tpos_731083-731121.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
