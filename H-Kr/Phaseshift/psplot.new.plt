set terminal pdfcairo
set output 'psplot.new.pdf'
set size ratio 0.55
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
#set title "Scattering phaseshift of partial vaves"
set xlabel "Center of mass kinetic energy (1/cm)"
set ylabel "Scattering phaseshift (rad.)
set border linewidth 2
set xtics 1
set ytics 1  
#set xtics font "Times New Roman, 18"
#set ytics font "Times New Roman, 18"
set xrange [2:6] 
set yrange [0.0:3.2]
set multiplot
set key at 5.9,1.8
set samples 50
plot atan2(1,-(x-4.13)/(0.56/2))-0.128 w p ps 1.2 lw 3 lt 6 lc 'red' title "Resonance fit:"
set key at 2.9,1.5
plot 'l=0.ps.txt' w l lt 2 dt 2 lw 4 lc 'black' title "l = 0:"
set key at 2.9,2.5
plot 'l=1.ps.txt' w l lt 3 dt 3 lw 6 lc 'sea-green' title "l = 1:"
set key at 3.7,2.9
plot 'l=2.ps.txt' w l lt 4 dt 4 lw 4 lc 'dark-orange' title "l = 2:"
set key at 4.5,2.7
plot 'l=3.ps.txt' w l lt 6 dt 5 lw 4 lc 'brown' title "l = 3:"
set key at 5.9,1.5
plot 'l=4.ps.txt' w l lt 7      lw 4 lc 'red' title "Res. wave l=4:"
set key at 5.3,0.7 
plot 'l=5.ps.txt' w l lt 1 dt 2 lw 4 lc 'blue' title "l = 5:"
set key at 5.3,0.4
plot 'l=6.ps.txt' w l lt 8 dt 3 lw 6 lc 'purple' title "l = 6:"
