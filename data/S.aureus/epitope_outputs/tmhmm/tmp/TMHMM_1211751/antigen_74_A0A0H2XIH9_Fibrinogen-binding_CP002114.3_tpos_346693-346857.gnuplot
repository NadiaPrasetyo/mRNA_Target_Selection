set arrow from 1,1.07 to 4,1.07 nohead lt 3 lw 10
set arrow from 5,1.09 to 27,1.09 nohead lt 1 lw 40
set arrow from 28,1.11 to 165,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_74|A0A0H2XIH9|Fibrinogen-binding|CP002114.3|tpos:346693-346857"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:165]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_74_A0A0H2XIH9_Fibrinogen-binding_CP002114.3_tpos_346693-346857.eps"
plot "./TMHMM_1211751/antigen_74_A0A0H2XIH9_Fibrinogen-binding_CP002114.3_tpos_346693-346857.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
