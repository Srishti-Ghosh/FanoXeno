
    set terminal pdfcairo size 8,5; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/angle_vs_w_magnitude.pdf'
    set title "H+Xe: Magnitude |w(θ)| (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "|w(θ)|" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/wparameter.txt' using 3:(sqrt($1**2 + $2**2)) w lines linecolor rgb "navy" lw 2.5 notitle
    