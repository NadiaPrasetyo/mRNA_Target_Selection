set arrow from 1,1.11 to 70,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_74|A0A0H2XIH9|Fibrinogen-binding|CP000253.1|tpos:65811-65880"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:70]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096803/antigen_74_A0A0H2XIH9_Fibrinogen-binding_CP000253.1_tpos_65811-65880.eps"
plot "./TMHMM_1096803/antigen_74_A0A0H2XIH9_Fibrinogen-binding_CP000253.1_tpos_65811-65880.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
