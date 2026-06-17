set terminal pdfcairo size 5,4
set output 'wprmwithprofile.pdf'
#set size ratio 0.80
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
set xrange [-15:30] 
set yrange [-20:15]
set xzeroaxis lw 1.5 lt 8 ; set yzeroaxis lw 1.5 lt 8
set label 1 at graph 0.18,0.95  "36deg" font "Times New Roman, 18, bold" right tc lt 1 
set label 2 at graph 0.20,0.51  "64deg" font "Times New Roman, 18, bold" right tc lt 2 
set label 3 at graph 0.59,0.90  "89deg" font "Times New Roman, 18, bold" right tc lt 6 
set label 4 at graph 0.72,0.51 "132deg" font "Times New Roman, 18, bold" left  tc lt 7 
set arrow 1 from first 10.00,-13.0 to 15.00,-17.00 lt 1 lw 2
set arrow 2 from first 21.50,5.75 to 28.50,5.00 lt 8 lw 2
set label 5 at graph 0.58,0.15   "0deg" font "Times New Roman, 18, bold" right tc lt 1
set label 6 at graph 0.99,0.78 "180deg" font "Times New Roman, 18, bold" right tc lt 8
set multiplot
set size 1,1; set origin 0,0
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
set size ratio 1.5
set size 0.25,0.35; set origin 0.17,0.57
unset label 
unset arrow
unset xlabel; unset ylabel
unset tics
set xrange [3:5]
set yrange [10:110]
#clear
plot '36deg.txt' w l lw 3 lc 7 notitle, '36asy.txt' w l dt 3 lw 3 lc 8 notitle
set size 0.25,0.35; set origin 0.18,0.20
set xrange [3:5]
set yrange [ 0: 60]
plot '64deg.txt' w l lw 3 lc 7 notitle, '64asy.txt' w l dt 3 lw 3 lc 8 notitle
set size 0.22,0.32; set origin 0.54,0.67
set xrange [3:5]
set yrange [ 0:300] 
plot '89deg.txt' w l lw 3 lc 7 notitle, '89asy.txt' w l dt 3 lw 3 lc 8 notitle
set size 0.25,0.35; set origin 0.60,0.20
set xrange [3:5]
set yrange [ 0:250]
plot '132deg.txt' w l lw 3 lc 7 notitle, '132asy.txt' w l dt 3 lw 3 lc 8 notitle

