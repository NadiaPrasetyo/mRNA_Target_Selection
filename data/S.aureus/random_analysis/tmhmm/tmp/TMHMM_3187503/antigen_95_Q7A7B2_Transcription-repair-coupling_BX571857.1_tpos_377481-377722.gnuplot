set arrow from 1,1.11 to 242,1.11 nohead lt 4 lw 10
set key below
set title "TMHMM posterior probabilities for antigen_95|Q7A7B2|Transcription-repair-coupling|BX571857.1|tpos:377481-377722"
set yrange [0:1.2]
set size 2., 1.4
#set xlabel "position"
set ylabel "probability"
set xrange [1:242]
# Make the ps plot
set term postscript eps color solid "Helvetica" 30
set output "./TMHMM_3187503/antigen_95_Q7A7B2_Transcription-repair-coupling_BX571857.1_tpos_377481-377722.eps"
plot "./TMHMM_3187503/antigen_95_Q7A7B2_Transcription-repair-coupling_BX571857.1_tpos_377481-377722.plp" using 1:4 title "transmembrane" with impulses lt 1 lw 2, \
"" using 1:3 title "inside" with line lt 3 lw 2, \
"" using 1:5 title "outside" with line lt 4 lw 2
exit
