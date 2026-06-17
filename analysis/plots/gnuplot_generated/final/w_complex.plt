
    set terminal pdfcairo size 6,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/real_w_vs_imag_w.pdf'
    set title "H+Xe: Complex w(θ) Parameter Trajectory (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"
    set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set size square
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3 w lines palette lw 3 title "Angle θ Evolution"
    