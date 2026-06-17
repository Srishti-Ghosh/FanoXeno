
    set terminal pdfcairo size 8,5
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/q_vs_theta.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Angle-Dependent Fano q-parameter (l=7)"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "q(θ)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/qparameter.txt' using 1:2 w lines linecolor rgb "dark-red" lw 2.5 notitle
    