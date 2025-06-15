set arrow from 1,1.11 to 686,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_150|O50581|ATP-dependent|CP002114.3|tpos:369233-369918"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:686]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_150_O50581_ATP-dependent_CP002114.3_tpos_369233-369918.eps"
plot "./TMHMM_3187290/antigen_150_O50581_ATP-dependent_CP002114.3_tpos_369233-369918.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
