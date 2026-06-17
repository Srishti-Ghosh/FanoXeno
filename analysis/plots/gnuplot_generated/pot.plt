set terminal pdfcairo size 8,6
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/pot.pdf'
set size ratio 0.80
set tics font "Times New Roman, 18"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key font "Times New Roman,16"
set title font "Times New Roman, 24, bold"
set title "H+Xe Potential"
set xlabel "Internuclear distance (a₀)"
set ylabel "Potential energy (meV)"
set border linewidth 2
set xtics 2
set ytics 2
set xrange [4:15]      # Focus on the well and immediate long-range tail
set yrange [-10:5]     # Clip the repulsive wall to see the -7.08 meV well
plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2)*27211.386) w l lt 8 lw 3 title "V(r) H+Xe"
