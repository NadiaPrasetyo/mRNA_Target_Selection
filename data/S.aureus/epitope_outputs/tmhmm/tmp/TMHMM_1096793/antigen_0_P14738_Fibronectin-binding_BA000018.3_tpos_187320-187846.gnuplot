set arrow from 1,1.11 to 527,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_0|P14738|Fibronectin-binding|BA000018.3|tpos:187320-187846"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:527]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096793/antigen_0_P14738_Fibronectin-binding_BA000018.3_tpos_187320-187846.eps"
plot "./TMHMM_1096793/antigen_0_P14738_Fibronectin-binding_BA000018.3_tpos_187320-187846.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
