set arrow from 1,1.11 to 109,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_20|Q5HHK8|Glycine|BX571857.1|tpos:108044-108152"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:109]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_20_Q5HHK8_Glycine_BX571857.1_tpos_108044-108152.eps"
plot "./TMHMM_3187503/antigen_20_Q5HHK8_Glycine_BX571857.1_tpos_108044-108152.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
