set arrow from 1,1.11 to 364,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_36|Q6GAB5|Spermidine/putrescine|HE681097.1|tpos:323039-323402"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:364]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_36_Q6GAB5_Spermidine_putrescine_HE681097.1_tpos_323039-323402.eps"
plot "./TMHMM_3187628/antigen_36_Q6GAB5_Spermidine_putrescine_HE681097.1_tpos_323039-323402.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
