
    set terminal pdfcairo size 8,6; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_5_plots/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Archetypal Fano Resonance Profiles (l=5, Er=5.41 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [3.91:6.91] 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_5_plots/slice_42.txt' using 1:2 w lines lw 2.5 title 'θ = 42° (q ≈ 0.0)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_5_plots/slice_90.txt' using 1:2 w lines lw 2.5 title 'θ = 90° (q ≈ 0.9)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_5_plots/slice_93.txt' using 1:2 w lines lw 2.5 title 'θ = 93° (q ≈ -69.4)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_5_plots/slice_158.txt' using 1:2 w lines lw 2.5 title 'θ = 158° (q ≈ -1.5)'
    