set arrow from 1,1.11 to 72,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_174|Q2FHM8|DNA-directed|HE681097.1|tpos:358227-358298"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:72]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_174_Q2FHM8_DNA-directed_HE681097.1_tpos_358227-358298.eps"
plot "./TMHMM_3187628/antigen_174_Q2FHM8_DNA-directed_HE681097.1_tpos_358227-358298.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
