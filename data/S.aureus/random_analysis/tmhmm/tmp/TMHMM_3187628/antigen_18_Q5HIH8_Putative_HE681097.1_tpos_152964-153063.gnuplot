set arrow from 1,1.07 to 100,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_18|Q5HIH8|Putative|HE681097.1|tpos:152964-153063"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:100]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_18_Q5HIH8_Putative_HE681097.1_tpos_152964-153063.eps"
plot "./TMHMM_3187628/antigen_18_Q5HIH8_Putative_HE681097.1_tpos_152964-153063.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
