set terminal pdfcairo size 8,6
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/psplot.pdf'
set tics font "Times New Roman, 16"
set xlabel font "Times New Roman, 18"
set ylabel font "Times New Roman, 18"
set key outside right top font "Times New Roman,14"
set title font "Times New Roman, 24, bold"
set title "H+Xe Phase Shifts"
set xlabel "Energy (cm⁻¹)"
set ylabel "Phase Shift η_l (radians)"
set border linewidth 2
set xtics 5
set xrange [0:50]
set yrange [-5:10]      # Adjust based on your data's specific phase range
set ytics 3.14159       # Set tick marks to multiples of Pi (π)
set format y "%.1f"
plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/phshift.txt' using 1:2 w l lw 2.2 title 'η_0(E)', \
     '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/phshift.txt' using 1:3 w l lw 2.2 title 'η_1(E)'
