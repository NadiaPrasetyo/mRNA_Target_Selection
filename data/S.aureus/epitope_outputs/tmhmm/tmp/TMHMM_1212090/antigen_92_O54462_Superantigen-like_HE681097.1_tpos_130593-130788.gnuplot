set arrow from 1,1.07 to 196,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_92|O54462|Superantigen-like|HE681097.1|tpos:130593-130788"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:196]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1212090/antigen_92_O54462_Superantigen-like_HE681097.1_tpos_130593-130788.eps"
plot "./TMHMM_1212090/antigen_92_O54462_Superantigen-like_HE681097.1_tpos_130593-130788.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
