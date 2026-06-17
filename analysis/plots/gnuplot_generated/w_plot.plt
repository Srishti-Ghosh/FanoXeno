
    set terminal pdfcairo size 6,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/w_complex_plane.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Complex w(θ) Parameter Trajectory (l=7)"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"
    set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set size square
    # Plotting Column 1 (Real) vs Column 2 (Imaginary), color-mapped to Column 3 (Angle)
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3 w lines palette lw 3 title "Angle θ Evolution"
    