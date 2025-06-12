set arrow from 1,1.11 to 259,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_20|Q2FWW7|ABC|BA000018.3|tpos:406944-407202"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:259]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_20_Q2FWW7_ABC_BA000018.3_tpos_406944-407202.eps"
plot "./TMHMM_1211906/antigen_20_Q2FWW7_ABC_BA000018.3_tpos_406944-407202.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
