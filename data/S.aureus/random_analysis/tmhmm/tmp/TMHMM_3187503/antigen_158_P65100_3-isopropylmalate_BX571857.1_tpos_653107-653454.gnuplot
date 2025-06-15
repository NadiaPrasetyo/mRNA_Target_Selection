set arrow from 1,1.11 to 348,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_158|P65100|3-isopropylmalate|BX571857.1|tpos:653107-653454"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:348]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_158_P65100_3-isopropylmalate_BX571857.1_tpos_653107-653454.eps"
plot "./TMHMM_3187503/antigen_158_P65100_3-isopropylmalate_BX571857.1_tpos_653107-653454.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
