
    set terminal pdfcairo size 8,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/koike_graph_1_effective_potentials.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Effective Interaction Potentials (l=0 to 8)"
    set xlabel "Internuclear distance R (a₀)" font "Times New Roman, 18"
    set ylabel "Effective Potential Energy (meV)" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    
    # Zoom in exactly on the potential well region
    set xrange [4:20]
    set yrange [-200:20]
    
    # Draw a line at E=0 to clearly show where the well disappears
    set xzeroaxis
    
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 0*(0+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=0', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 1*(1+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=1', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 2*(2+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=2', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 3*(3+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=3', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 4*(4+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=4', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 5*(5+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=5', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 6*(6+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=6', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 7*(7+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=7', \
         '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/potential.txt' using 1:(($2 + 8*(8+1)/(2.0*1823.14*$1**2)) * 27211.386) w lines lw 2 title 'l=8'
    