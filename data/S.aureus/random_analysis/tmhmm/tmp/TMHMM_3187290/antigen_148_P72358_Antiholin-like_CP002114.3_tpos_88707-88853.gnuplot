set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 32,1.09 nohead lt 1 lw 40
set arrow from 33,1.11 to 41,1.11 nohead lt 4 lw 10
set arrow from 42,1.09 to 64,1.09 nohead lt 1 lw 40
set arrow from 65,1.07 to 70,1.07 nohead lt 3 lw 10
set arrow from 71,1.09 to 93,1.09 nohead lt 1 lw 40
set arrow from 94,1.11 to 96,1.11 nohead lt 4 lw 10
set arrow from 97,1.09 to 119,1.09 nohead lt 1 lw 40
set arrow from 120,1.07 to 147,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_148|P72358|Antiholin-like|CP002114.3|tpos:88707-88853"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:147]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_148_P72358_Antiholin-like_CP002114.3_tpos_88707-88853.eps"
plot "./TMHMM_3187290/antigen_148_P72358_Antiholin-like_CP002114.3_tpos_88707-88853.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
