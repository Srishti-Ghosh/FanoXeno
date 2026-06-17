set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key at 5.9,2.0
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
set title "Phaseshift analysis of l=4 resonance"
set xlabel "Kinetic Energy (1/cm)"
set ylabel "Scattering phaseshift (rad.)
set border linewidth 2
#set xtics font "Times New Roman, 18"
#set ytics font "Times New Roman, 18"
set xrange [2:6] 
set yrange [0:3]
plot 'l=4.ps.txt' w l lt 1 lw 2 lc 7 title "VPA Calculation", atan2(1,-(x-4.13)/(0.56/2))-0.128 w p ps 0.5 lt 5 lc 6 title "Resonance fit"
pause -100   
set terminal pdfcairo
set output 'psanalysis.pdf' 
replot
