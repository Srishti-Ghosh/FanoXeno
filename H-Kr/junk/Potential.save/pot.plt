set terminal pdfcairo
set output 'pot.pdf'
set size ratio 0.80
#set colorsequence classic
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
#set title "Scattering phaseshift of partial vaves"
set xlabel "Internuclear distance (au)"
set ylabel "Potential energy (eV)"
set border linewidth 2
set xtics 2
set ytics 0.002 
#set xtics font "Times New Roman, 18"
#set ytics font "Times New Roman, 18"
set xrange [4:20] 
set yrange [-0.007:0.005]
set multiplot
#set key at 5.9,1.8
#plot 'potential.txt' w l lt 7 dt 6 lw 3 lc 7 notitle
#set style line 1 lc rgb '#ff0000' lt 1 lw 2 pi -1 ps 1.0
#set style line 2 lc rgb '#ffffff' lt 1 lw 2 pi -1 ps 1.0
#plot 'l=0.pot.txt' w l ls 1 title "l=0:"
#plot 'l=1.pot.txt' w l ls 2 title "l=1:"
plot 'l=0.pot.txt' w l lt 8 dt 1 lw 3 notitle
plot 'l=1.pot.txt' w l lt 3 dt 3 lw 3 notitle
plot 'l=2.pot.txt' w l lt 4 dt 4 lw 3 notitle
plot 'l=3.pot.txt' w l lt 6 dt 5 lw 3 notitle
plot 'l=4.pot.txt' w l lt 1 dt 1 lw 3 notitle
plot 'l=5.pot.txt' w l lt 7 dt 7 lw 3 notitle
plot 'l=6.pot.txt' w l lt 8 dt 8 lw 3 notitle
plot 0 w l lt 8 lw 2 notitle
plot 0.00012398*4.13 w l lt 7 lw 2 notitle
