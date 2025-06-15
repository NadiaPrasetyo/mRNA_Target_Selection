set arrow from 1,1.07 to 6,1.07 nohead lt 3 lw 10
set arrow from 7,1.09 to 26,1.09 nohead lt 1 lw 40
set arrow from 27,1.11 to 254,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_48|Q936F2|Uncharacterized|BA000018.3|tpos:33305-33558"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:254]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187501/antigen_48_Q936F2_Uncharacterized_BA000018.3_tpos_33305-33558.eps"
plot "./TMHMM_3187501/antigen_48_Q936F2_Uncharacterized_BA000018.3_tpos_33305-33558.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
