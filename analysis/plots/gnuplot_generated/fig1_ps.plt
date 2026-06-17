
    set terminal pdfcairo size 8,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/fig1_phase_shifts.pdf'
    set title font "Times New Roman, 20, bold"
    set title "Figure 1: H+Xe Phase Shifts η_l(E)"
    set xlabel "Collision Energy E (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Phase Shift η_l (radians)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    
    # Position the legend outside the graph
    set key outside right top font "Times New Roman,14"
    
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/phshift_matrix.txt' using 1:2 w lines lw 1.5 title 'l=0', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/phshift_matrix.txt' using 1:3 w lines lw 1.5 title 'l=1', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/phshift_matrix.txt' using 1:8 w lines lw 1.5 title 'l=6', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/phshift_matrix.txt' using 1:9 w lines lw 3.0 title 'l=7', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/phshift_matrix.txt' using 1:10 w lines lw 1.5 title 'l=8'
    