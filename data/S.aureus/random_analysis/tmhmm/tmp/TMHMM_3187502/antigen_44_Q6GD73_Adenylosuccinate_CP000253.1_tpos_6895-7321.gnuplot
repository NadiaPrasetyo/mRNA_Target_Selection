set arrow from 1,1.11 to 427,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_44|Q6GD73|Adenylosuccinate|CP000253.1|tpos:6895-7321"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:427]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_44_Q6GD73_Adenylosuccinate_CP000253.1_tpos_6895-7321.eps"
plot "./TMHMM_3187502/antigen_44_Q6GD73_Adenylosuccinate_CP000253.1_tpos_6895-7321.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
