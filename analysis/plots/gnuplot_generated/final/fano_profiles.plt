
    set terminal pdfcairo size 8,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Beutler-Fano Resonance Profiles (l=7, Er=25.89 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [24.5:27.5] 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/slice_30deg.txt' using 1:2 w lines lw 2.5 title 'θ = 30.0°', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/slice_60deg.txt' using 1:2 w lines lw 2.5 title 'θ = 60.0°', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/slice_120deg.txt' using 1:2 w lines lw 2.5 title 'θ = 120.0°'
    