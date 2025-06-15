set arrow from 1,1.11 to 220,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_138|Q5HE63|Deoxyribose-phosphate|BA000018.3|tpos:47919-48138"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:220]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187501/antigen_138_Q5HE63_Deoxyribose-phosphate_BA000018.3_tpos_47919-48138.eps"
plot "./TMHMM_3187501/antigen_138_Q5HE63_Deoxyribose-phosphate_BA000018.3_tpos_47919-48138.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
