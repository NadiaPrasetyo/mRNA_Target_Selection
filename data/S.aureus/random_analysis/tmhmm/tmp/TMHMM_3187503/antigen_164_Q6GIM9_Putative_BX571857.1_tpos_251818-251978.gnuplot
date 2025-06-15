set arrow from 1,1.11 to 161,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_164|Q6GIM9|Putative|BX571857.1|tpos:251818-251978"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:161]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_164_Q6GIM9_Putative_BX571857.1_tpos_251818-251978.eps"
plot "./TMHMM_3187503/antigen_164_Q6GIM9_Putative_BX571857.1_tpos_251818-251978.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
