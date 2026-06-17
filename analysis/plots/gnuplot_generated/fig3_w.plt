
    set terminal pdfcairo size 7,7
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/fig3_w_complex_plane.pdf'
    set title font "Times New Roman, 20, bold"
    set title "Figure 3: Complex w(θ) Parameter Trajectory (l=7)"
    set xlabel "Re[w(θ)]" font "Times New Roman, 18"
    set ylabel "Im[w(θ)]" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    
    # Force a square aspect ratio so circles don't look like ovals
    set size square 
    
    # Plotting Column 1 (Real) vs Column 2 (Imaginary)
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2 w lines linecolor rgb "navy" lw 2.5 notitle
    