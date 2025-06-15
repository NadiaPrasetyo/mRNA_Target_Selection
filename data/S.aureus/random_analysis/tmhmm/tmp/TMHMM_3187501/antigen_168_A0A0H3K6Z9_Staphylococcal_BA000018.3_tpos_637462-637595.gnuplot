set arrow from 1,1.07 to 134,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_168|A0A0H3K6Z9|Staphylococcal|BA000018.3|tpos:637462-637595"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:134]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187501/antigen_168_A0A0H3K6Z9_Staphylococcal_BA000018.3_tpos_637462-637595.eps"
plot "./TMHMM_3187501/antigen_168_A0A0H3K6Z9_Staphylococcal_BA000018.3_tpos_637462-637595.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
