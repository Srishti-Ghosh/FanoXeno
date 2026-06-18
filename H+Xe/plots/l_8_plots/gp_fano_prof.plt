
    set terminal pdfcairo size 8,6; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Archetypal Fano Resonance Profiles (l=8, Er=39.95 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [38.45:41.45] 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/slice_0.txt' using 1:2 w lines lw 2.5 title 'θ = 0° (q ≈ -0.0)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/slice_46.txt' using 1:2 w lines lw 2.5 title 'θ = 46° (q ≈ 1.5)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/slice_79.txt' using 1:2 w lines lw 2.5 title 'θ = 79° (q ≈ -1.6)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/slice_126.txt' using 1:2 w lines lw 2.5 title 'θ = 126° (q ≈ -81.9)'
    