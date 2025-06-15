set arrow from 1,1.11 to 303,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_69|A5IQW9|Nucleotide-binding|BX571857.1|tpos:252926-253228"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:303]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_69_A5IQW9_Nucleotide-binding_BX571857.1_tpos_252926-253228.eps"
plot "./TMHMM_3187503/antigen_69_A5IQW9_Nucleotide-binding_BX571857.1_tpos_252926-253228.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
