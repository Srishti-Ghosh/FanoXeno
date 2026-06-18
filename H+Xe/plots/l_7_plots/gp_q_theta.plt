
    set terminal pdfcairo size 8,5; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/angle_vs_q_parameter.pdf'
    set title "H+Xe: Angle-Dependent Fano q-parameter (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "q(θ)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/qparameter.txt' using 1:2 w lines linecolor rgb "dark-red" lw 2.5 notitle
    