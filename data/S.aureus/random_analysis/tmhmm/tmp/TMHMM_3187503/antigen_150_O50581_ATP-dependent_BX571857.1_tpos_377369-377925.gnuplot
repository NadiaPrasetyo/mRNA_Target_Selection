set arrow from 1,1.11 to 557,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_150|O50581|ATP-dependent|BX571857.1|tpos:377369-377925"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:557]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_150_O50581_ATP-dependent_BX571857.1_tpos_377369-377925.eps"
plot "./TMHMM_3187503/antigen_150_O50581_ATP-dependent_BX571857.1_tpos_377369-377925.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
