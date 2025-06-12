set arrow from 1,1.07 to 12,1.07 nohead lt 3 lw 10
set arrow from 13,1.09 to 35,1.09 nohead lt 1 lw 40
set arrow from 36,1.11 to 49,1.11 nohead lt 4 lw 10
set arrow from 50,1.09 to 72,1.09 nohead lt 1 lw 40
set arrow from 73,1.07 to 80,1.07 nohead lt 3 lw 10
set arrow from 81,1.09 to 100,1.09 nohead lt 1 lw 40
set arrow from 101,1.11 to 104,1.11 nohead lt 4 lw 10
set arrow from 105,1.09 to 127,1.09 nohead lt 1 lw 40
set arrow from 128,1.07 to 138,1.07 nohead lt 3 lw 10
set arrow from 139,1.09 to 161,1.09 nohead lt 1 lw 40
set arrow from 162,1.11 to 164,1.11 nohead lt 4 lw 10
set arrow from 165,1.09 to 187,1.09 nohead lt 1 lw 40
set arrow from 188,1.07 to 216,1.07 nohead lt 3 lw 10
set arrow from 217,1.09 to 239,1.09 nohead lt 1 lw 40
set arrow from 240,1.11 to 253,1.11 nohead lt 4 lw 10
set arrow from 254,1.09 to 276,1.09 nohead lt 1 lw 40
set arrow from 277,1.07 to 282,1.07 nohead lt 3 lw 10
set arrow from 283,1.09 to 305,1.09 nohead lt 1 lw 40
set arrow from 306,1.11 to 308,1.11 nohead lt 4 lw 10
set arrow from 309,1.09 to 331,1.09 nohead lt 1 lw 40
set arrow from 332,1.07 to 342,1.07 nohead lt 3 lw 10
set arrow from 343,1.09 to 365,1.09 nohead lt 1 lw 40
set arrow from 366,1.11 to 374,1.11 nohead lt 4 lw 10
set arrow from 375,1.09 to 394,1.09 nohead lt 1 lw 40
set arrow from 395,1.07 to 403,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_6|Q2FVI5|Bcr/CflA|CP002114.3|tpos:776406-776808"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:403]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_1211751/antigen_6_Q2FVI5_Bcr_CflA_CP002114.3_tpos_776406-776808.eps"
plot "./TMHMM_1211751/antigen_6_Q2FVI5_Bcr_CflA_CP002114.3_tpos_776406-776808.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
