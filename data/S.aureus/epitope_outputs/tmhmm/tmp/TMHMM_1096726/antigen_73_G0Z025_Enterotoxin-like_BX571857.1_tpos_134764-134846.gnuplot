set arrow from 1,1.11 to 83,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_73|G0Z025|Enterotoxin-like|BX571857.1|tpos:134764-134846"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:83]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1096726/antigen_73_G0Z025_Enterotoxin-like_BX571857.1_tpos_134764-134846.eps"
plot "./TMHMM_1096726/antigen_73_G0Z025_Enterotoxin-like_BX571857.1_tpos_134764-134846.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
