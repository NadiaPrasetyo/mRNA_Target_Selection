set arrow from 1,1.11 to 215,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_45|Q6GJE0|Serine|CP002114.3|tpos:172783-172997"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:215]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_45_Q6GJE0_Serine_CP002114.3_tpos_172783-172997.eps"
plot "./TMHMM_3187290/antigen_45_Q6GJE0_Serine_CP002114.3_tpos_172783-172997.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
