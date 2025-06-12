set arrow from 1,1.11 to 309,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_8|Q2FVK3|Gamma-hemolysin|BA000018.3|tpos:767176-767484"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:309]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211906/antigen_8_Q2FVK3_Gamma-hemolysin_BA000018.3_tpos_767176-767484.eps"
plot "./TMHMM_1211906/antigen_8_Q2FVK3_Gamma-hemolysin_BA000018.3_tpos_767176-767484.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
