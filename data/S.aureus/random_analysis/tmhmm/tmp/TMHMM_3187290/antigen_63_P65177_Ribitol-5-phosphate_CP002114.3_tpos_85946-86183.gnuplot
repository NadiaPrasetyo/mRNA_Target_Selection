set arrow from 1,1.11 to 238,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_63|P65177|Ribitol-5-phosphate|CP002114.3|tpos:85946-86183"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:238]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_63_P65177_Ribitol-5-phosphate_CP002114.3_tpos_85946-86183.eps"
plot "./TMHMM_3187290/antigen_63_P65177_Ribitol-5-phosphate_CP002114.3_tpos_85946-86183.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
