set arrow from 1,1.11 to 436,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_77|Q5HDD7|Immunoglobulin-binding|HE681097.1|tpos:773380-773815"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:436]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_77_Q5HDD7_Immunoglobulin-binding_HE681097.1_tpos_773380-773815.eps"
plot "./TMHMM_1212090/antigen_77_Q5HDD7_Immunoglobulin-binding_HE681097.1_tpos_773380-773815.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
