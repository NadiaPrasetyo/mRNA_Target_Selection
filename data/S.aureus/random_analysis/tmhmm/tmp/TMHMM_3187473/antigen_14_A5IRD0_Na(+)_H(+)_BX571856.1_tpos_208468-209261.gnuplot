set arrow from 1,1.11 to 3,1.11 nohead lt 4 lw 10
set arrow from 4,1.09 to 21,1.09 nohead lt 1 lw 40
set arrow from 22,1.07 to 29,1.07 nohead lt 3 lw 10
set arrow from 30,1.09 to 52,1.09 nohead lt 1 lw 40
set arrow from 53,1.11 to 76,1.11 nohead lt 4 lw 10
set arrow from 77,1.09 to 99,1.09 nohead lt 1 lw 40
set arrow from 100,1.07 to 111,1.07 nohead lt 3 lw 10
set arrow from 112,1.09 to 129,1.09 nohead lt 1 lw 40
set arrow from 130,1.11 to 132,1.11 nohead lt 4 lw 10
set arrow from 133,1.09 to 155,1.09 nohead lt 1 lw 40
set arrow from 156,1.07 to 166,1.07 nohead lt 3 lw 10
set arrow from 167,1.09 to 189,1.09 nohead lt 1 lw 40
set arrow from 190,1.11 to 208,1.11 nohead lt 4 lw 10
set arrow from 209,1.09 to 231,1.09 nohead lt 1 lw 40
set arrow from 232,1.07 to 243,1.07 nohead lt 3 lw 10
set arrow from 244,1.09 to 266,1.09 nohead lt 1 lw 40
set arrow from 267,1.11 to 270,1.11 nohead lt 4 lw 10
set arrow from 271,1.09 to 293,1.09 nohead lt 1 lw 40
set arrow from 294,1.07 to 299,1.07 nohead lt 3 lw 10
set arrow from 300,1.09 to 320,1.09 nohead lt 1 lw 40
set arrow from 321,1.11 to 329,1.11 nohead lt 4 lw 10
set arrow from 330,1.09 to 352,1.09 nohead lt 1 lw 40
set arrow from 353,1.07 to 384,1.07 nohead lt 3 lw 10
set arrow from 385,1.09 to 407,1.09 nohead lt 1 lw 40
set arrow from 408,1.11 to 426,1.11 nohead lt 4 lw 10
set arrow from 427,1.09 to 449,1.09 nohead lt 1 lw 40
set arrow from 450,1.07 to 472,1.07 nohead lt 3 lw 10
set arrow from 473,1.09 to 495,1.09 nohead lt 1 lw 40
set arrow from 496,1.11 to 524,1.11 nohead lt 4 lw 10
set arrow from 525,1.09 to 547,1.09 nohead lt 1 lw 40
set arrow from 548,1.07 to 590,1.07 nohead lt 3 lw 10
set arrow from 591,1.09 to 613,1.09 nohead lt 1 lw 40
set arrow from 614,1.11 to 627,1.11 nohead lt 4 lw 10
set arrow from 628,1.09 to 647,1.09 nohead lt 1 lw 40
set arrow from 648,1.07 to 653,1.07 nohead lt 3 lw 10
set arrow from 654,1.09 to 671,1.09 nohead lt 1 lw 40
set arrow from 672,1.11 to 675,1.11 nohead lt 4 lw 10
set arrow from 676,1.09 to 698,1.09 nohead lt 1 lw 40
set arrow from 699,1.07 to 710,1.07 nohead lt 3 lw 10
set arrow from 711,1.09 to 733,1.09 nohead lt 1 lw 40
set arrow from 734,1.11 to 770,1.11 nohead lt 4 lw 10
set arrow from 771,1.09 to 790,1.09 nohead lt 1 lw 40
set arrow from 791,1.07 to 794,1.07 nohead lt 3 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_14|A5IRD0|Na(+)/H(+)|BX571856.1|tpos:208468-209261"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:794]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187473/antigen_14_A5IRD0_Na(+)_H(+)_BX571856.1_tpos_208468-209261.eps"
plot "./TMHMM_3187473/antigen_14_A5IRD0_Na(+)_H(+)_BX571856.1_tpos_208468-209261.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
