set arrow from 1,1.11 to 449,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_19|P68789|Elongation|HE681097.1|tpos:296931-297379"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:449]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187628/antigen_19_P68789_Elongation_HE681097.1_tpos_296931-297379.eps"
plot "./TMHMM_3187628/antigen_19_P68789_Elongation_HE681097.1_tpos_296931-297379.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
