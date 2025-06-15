set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.11 to 255,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_48|Q936F2|Uncharacterized|HE681097.1|tpos:26085-26339"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:255]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_48_Q936F2_Uncharacterized_HE681097.1_tpos_26085-26339.eps"
plot "./TMHMM_3187628/antigen_48_Q936F2_Uncharacterized_HE681097.1_tpos_26085-26339.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
