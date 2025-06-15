set arrow from 1,1.07 to 19,1.07 nohead lt 3 lw 10
set arrow from 20,1.09 to 42,1.09 nohead lt 1 lw 40
set arrow from 43,1.11 to 61,1.11 nohead lt 4 lw 10
set arrow from 62,1.09 to 84,1.09 nohead lt 1 lw 40
set arrow from 85,1.07 to 90,1.07 nohead lt 3 lw 10
set arrow from 91,1.09 to 113,1.09 nohead lt 1 lw 40
set arrow from 114,1.11 to 117,1.11 nohead lt 4 lw 10
set arrow from 118,1.09 to 140,1.09 nohead lt 1 lw 40
set arrow from 141,1.07 to 159,1.07 nohead lt 3 lw 10
set arrow from 160,1.09 to 182,1.09 nohead lt 1 lw 40
set arrow from 183,1.11 to 191,1.11 nohead lt 4 lw 10
set arrow from 192,1.09 to 211,1.09 nohead lt 1 lw 40
set arrow from 212,1.07 to 241,1.07 nohead lt 3 lw 10
set arrow from 242,1.09 to 264,1.09 nohead lt 1 lw 40
set arrow from 265,1.11 to 278,1.11 nohead lt 4 lw 10
set arrow from 279,1.09 to 301,1.09 nohead lt 1 lw 40
set arrow from 302,1.07 to 312,1.07 nohead lt 3 lw 10
set arrow from 313,1.09 to 332,1.09 nohead lt 1 lw 40
set arrow from 333,1.11 to 336,1.11 nohead lt 4 lw 10
set arrow from 337,1.09 to 354,1.09 nohead lt 1 lw 40
set arrow from 355,1.07 to 374,1.07 nohead lt 3 lw 10
set arrow from 375,1.09 to 397,1.09 nohead lt 1 lw 40
set arrow from 398,1.11 to 406,1.11 nohead lt 4 lw 10
set arrow from 407,1.09 to 429,1.09 nohead lt 1 lw 40
set arrow from 430,1.07 to 466,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_146|Q8NXW9|Putative|BX571857.1|tpos:189877-190342"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:466]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_146_Q8NXW9_Putative_BX571857.1_tpos_189877-190342.eps"
plot "./TMHMM_3187503/antigen_146_Q8NXW9_Putative_BX571857.1_tpos_189877-190342.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
