set arrow from 1,1.11 to 222,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_21|Q6GCF3|Putative|BX571857.1|tpos:105493-105714"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:222]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_21_Q6GCF3_Putative_BX571857.1_tpos_105493-105714.eps"
plot "./TMHMM_3187503/antigen_21_Q6GCF3_Putative_BX571857.1_tpos_105493-105714.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
