set arrow from 1,1.11 to 152,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_72|A0AAX2YX78|Extracellular|BX571857.1|tpos:704357-704508"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:152]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096726/antigen_72_A0AAX2YX78_Extracellular_BX571857.1_tpos_704357-704508.eps"
plot "./TMHMM_1096726/antigen_72_A0AAX2YX78_Extracellular_BX571857.1_tpos_704357-704508.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
