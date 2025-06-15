set arrow from 1,1.11 to 69,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_147|A6QHI5|UPF0735|CP000253.1|tpos:400236-400304"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:69]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187502/antigen_147_A6QHI5_UPF0735_CP000253.1_tpos_400236-400304.eps"
plot "./TMHMM_3187502/antigen_147_A6QHI5_UPF0735_CP000253.1_tpos_400236-400304.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
