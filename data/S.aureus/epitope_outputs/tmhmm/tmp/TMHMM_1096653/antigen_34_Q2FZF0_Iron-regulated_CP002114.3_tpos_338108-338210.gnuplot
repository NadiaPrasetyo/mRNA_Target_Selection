set arrow from 1,1.07 to 103,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_34|Q2FZF0|Iron-regulated|CP002114.3|tpos:338108-338210"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:103]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096653/antigen_34_Q2FZF0_Iron-regulated_CP002114.3_tpos_338108-338210.eps"
plot "./TMHMM_1096653/antigen_34_Q2FZF0_Iron-regulated_CP002114.3_tpos_338108-338210.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
