set arrow from 1,1.11 to 175,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_73|A8Z3N5|Bifunctional|CP002114.3|tpos:359168-359342"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:175]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_73_A8Z3N5_Bifunctional_CP002114.3_tpos_359168-359342.eps"
plot "./TMHMM_3187290/antigen_73_A8Z3N5_Bifunctional_CP002114.3_tpos_359168-359342.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
