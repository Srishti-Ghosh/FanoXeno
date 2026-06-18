
    set terminal pdfcairo size 8,6; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Archetypal Fano Resonance Profiles (l=7, Er=25.98 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [24.48:27.48] 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/slice_78.txt' using 1:2 w lines lw 2.5 title 'θ = 78° (q ≈ -1.5)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/slice_86.txt' using 1:2 w lines lw 2.5 title 'θ = 86° (q ≈ 0.0)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/slice_109.txt' using 1:2 w lines lw 2.5 title 'θ = 109° (q ≈ -26.0)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/slice_127.txt' using 1:2 w lines lw 2.5 title 'θ = 127° (q ≈ 1.5)'
    