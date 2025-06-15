set arrow from 1,1.07 to 1,1.07 nohead lt 3 lw 10
set arrow from 2,1.09 to 19,1.09 nohead lt 1 lw 40
set arrow from 20,1.11 to 22,1.11 nohead lt 4 lw 10
set arrow from 23,1.09 to 45,1.09 nohead lt 1 lw 40
set arrow from 46,1.07 to 51,1.07 nohead lt 3 lw 10
set arrow from 52,1.09 to 71,1.09 nohead lt 1 lw 40
set arrow from 72,1.11 to 96,1.11 nohead lt 4 lw 10
set arrow from 97,1.09 to 119,1.09 nohead lt 1 lw 40
set arrow from 120,1.07 to 139,1.07 nohead lt 3 lw 10
set arrow from 140,1.09 to 162,1.09 nohead lt 1 lw 40
set arrow from 163,1.11 to 176,1.11 nohead lt 4 lw 10
set arrow from 177,1.09 to 199,1.09 nohead lt 1 lw 40
set arrow from 200,1.07 to 224,1.07 nohead lt 3 lw 10
set arrow from 225,1.09 to 242,1.09 nohead lt 1 lw 40
set arrow from 243,1.11 to 246,1.11 nohead lt 4 lw 10
set arrow from 247,1.09 to 264,1.09 nohead lt 1 lw 40
set arrow from 265,1.07 to 283,1.07 nohead lt 3 lw 10
set arrow from 284,1.09 to 306,1.09 nohead lt 1 lw 40
set arrow from 307,1.11 to 315,1.11 nohead lt 4 lw 10
set arrow from 316,1.09 to 335,1.09 nohead lt 1 lw 40
set arrow from 336,1.07 to 405,1.07 nohead lt 3 lw 10
set arrow from 406,1.09 to 428,1.09 nohead lt 1 lw 40
set arrow from 429,1.11 to 429,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_192|P30329|Arsenical|CP002114.3|tpos:576369-576797"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:429]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187290/antigen_192_P30329_Arsenical_CP002114.3_tpos_576369-576797.eps"
plot "./TMHMM_3187290/antigen_192_P30329_Arsenical_CP002114.3_tpos_576369-576797.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
