set terminal pdfcairo
set output 'wprmongauss.pdf'
set size ratio 0.80
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
#set title "Scattering phaseshift of partial vaves"
set xlabel "x" font "Times New Roman, 28, bold" offset 0,0.25,0
set ylabel "y" font "Times New Roman, 28, bold" offset 0,0,0
set border linewidth 1.5
set xtics 10
set ytics 10  
#set xtics font "Times New Roman, 18"
#set ytics font "Times New Roman, 18"
set xrange [-10:30] 
set yrange [-20:10]
set xzeroaxis lw 1.5 lt 8 ; set yzeroaxis lw 1.5 lt 8
set label 1 at graph 0.14,0.87  "36deg" font "Times New Roman, 18, bold" right tc lt 1 
set label 2 at graph 0.15,0.62  "64deg" font "Times New Roman, 18, bold" right tc lt 2 
set label 3 at graph 0.71,0.72  "89deg" font "Times New Roman, 18, bold" left  tc lt 6 
set label 4 at graph 0.64,0.25 "132deg" font "Times New Roman, 18, bold" left  tc lt 7 
set arrow 1 from first 10.00,-13.0 to 15.00,-17.00 lt 1 lw 2
set arrow 2 from first 18.00,5.75 to 25.00,5.00 lt 8 lw 2
set label 5 at graph 0.52,0.18   "0deg" font "Times New Roman, 18, bold" right tc lt 1
set label 6 at graph 0.85,0.90 "180deg" font "Times New Roman, 18, bold" right tc lt 1
set multiplot
#set key at 5.9,1.8
plot 'wparameter.txt' every::1::36 w p ps 1 pt 4 lw 1.5 lc 1 notitle
plot 'wparameter.txt' every::1::37 w l lw 1.5 lc 1 notitle
plot 'wparameter.txt' every::37::64 w p ps 1 pt 8 lw 1.5 lc 2 notitle
plot 'wparameter.txt' every::36::65 w l lw 1.5 lc 2 notitle
plot 'wparameter.txt' every::65::89 w p ps 1 pt 10 lw 1.5 lc 6 notitle
plot 'wparameter.txt' every::64::90 w l lw 1.5 lc 6 notitle
plot 'wparameter.txt' every::90::132 w p ps 1 pt 12 lw 1.5 lc 7 notitle
plot 'wparameter.txt' every::89::133 w l lw 1.5 lc 7 notitle
plot 'wparameter.txt' every::133::180 w p ps 1 pt 6 lw 1.5 lc 8 notitle
plot 'wparameter.txt' every::132::180 w l lw 1.5 lc 8 notitle
#
plot 'wparameter.txt' every::36::36 w p ps 1.5 pt 5 lw 1.5 lc 1 notitle
plot 'wparameter.txt' every::64::64 w p ps 1.5 pt 9 lw 1.5 lc 2 notitle
plot 'wparameter.txt' every::89::89 w p ps 1.5 pt 11 lw 1.5 lc 6 notitle
plot 'wparameter.txt' every::132::132 w p ps 1.5 pt 13 lw 1.5 lc 7 notitle
#
