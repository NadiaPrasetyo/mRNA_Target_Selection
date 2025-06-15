set arrow from 1,1.11 to 311,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_110|A6U0Z9|Ribosomal|BA000018.3|tpos:357267-357577"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:311]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187501/antigen_110_A6U0Z9_Ribosomal_BA000018.3_tpos_357267-357577.eps"
plot "./TMHMM_3187501/antigen_110_A6U0Z9_Ribosomal_BA000018.3_tpos_357267-357577.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
