
    set terminal pdfcairo size 8,5
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/w_magnitude.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: Magnitude of w-parameter |w(θ)| (l=7)"
    set xlabel "Scattering Angle θ (degrees)" font "Times New Roman, 18"
    set ylabel "|w(θ)|" font "Times New Roman, 18"
    set grid
    set border linewidth 2
    set xrange [0:180]
    set xtics 30
    plot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 3:(sqrt($1**2 + $2**2)) w lines linecolor rgb "navy" lw 2.5 notitle
    