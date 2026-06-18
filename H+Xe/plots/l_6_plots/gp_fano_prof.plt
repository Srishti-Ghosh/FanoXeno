
    set terminal pdfcairo size 8,6; set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/energy_vs_dcs_profiles.pdf'
    set title "H+Xe: Archetypal Fano Resonance Profiles (l=6, Er=14.475 cm⁻¹)" font "Times New Roman, 20, bold"
    set xlabel "Collision Energy (cm⁻¹)" font "Times New Roman, 18"
    set ylabel "Differential Cross Section (a.u.)" font "Times New Roman, 18"
    set grid; set border linewidth 2; set xrange [12.975:15.975] 
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/slice_71.txt' using 1:2 w lines lw 2.5 title 'θ = 71° (q ≈ -24.3)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/slice_104.txt' using 1:2 w lines lw 2.5 title 'θ = 104° (q ≈ -1.5)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/slice_134.txt' using 1:2 w lines lw 2.5 title 'θ = 134° (q ≈ 1.5)', '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/slice_162.txt' using 1:2 w lines lw 2.5 title 'θ = 162° (q ≈ -0.0)'
    