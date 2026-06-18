
    set terminal pdfcairo size 8,5; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/angle_vs_z_background.pdf'
    set title "H+Xe: Background Factor z(θ) (l=6)" font "Times New Roman, 20, bold"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"; set ylabel "z(θ) (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [0:180]; set xtics 30; set logscale y 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/qparameter.txt' using 1:3 w lines linecolor rgb "forest-green" lw 2.5 notitle
    