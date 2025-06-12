set arrow from 1,1.11 to 67,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_73|G0Z025|Enterotoxin-like|BX571856.1|tpos:141167-141233"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:67]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211804/antigen_73_G0Z025_Enterotoxin-like_BX571856.1_tpos_141167-141233.eps"
plot "./TMHMM_1211804/antigen_73_G0Z025_Enterotoxin-like_BX571856.1_tpos_141167-141233.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
