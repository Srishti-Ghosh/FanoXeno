set terminal pdfcairo size 5,4.2
set output 'qprmongraph.pdf'
set size ratio 0.80
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
#set title "Scattering phaseshift of partial vaves"
set xlabel "Scattering angle" font "Times New Roman, 28, bold" offset 0,0.25,0
set ylabel "q" font "Times New Roman, 28, bold" offset 0,0,0
set border linewidth 1.5
set xtics 30
set ytics  5  
#set xtics font "Times New Roman, 18"
#set ytics font "Times New Roman, 18"
set xrange [0:180] 
set yrange [-20:10]
set xzeroaxis lw 1.5 lt 8
#set yzeroaxis lw 1.5 lt 8
set label 1 at graph 0.93,0.85  "x(1/10)" font "Times New Roman, 18, bold" right tc lt 8 
set arrow 1 from first 133,-20 to 133,10 lt 8 lw 1.0 nohead
set multiplot
#set key at 5.9,1.8

plot 'qparameter.txt' every::1::36 w p ps 1 pt 4 lw 1.5 lc 1 notitle
plot 'qparameter.txt' every::1::37 w l lw 1.5 lc 1 notitle
plot 'qparameter.txt' every::37::64 w p ps 1 pt 8 lw 1.5 lc 2 notitle
plot 'qparameter.txt' every::36::65 w l lw 1.5 lc 2 notitle
plot 'qparameter.txt' every::65::89 w p ps 1 pt 10 lw 1.5 lc 6 notitle
plot 'qparameter.txt' every::64::90 w l lw 1.5 lc 6 notitle
plot 'qparameter.txt' every::90::132 w p ps 1 pt 12 lw 1.5 lc 7 notitle
plot 'qparameter.txt' every::89::133 w l lw 1.5 lc 7 notitle
plot 'qprm.over10.txt' every::133::180 w p ps 1 pt 6 lw 1.5 lc 8 notitle
plot 'qprm.over10.txt' every::132::180 w l lw 1.5 lc 8 notitle
#

