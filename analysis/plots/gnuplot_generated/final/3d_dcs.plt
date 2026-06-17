
    set terminal pdfcairo size 5,6
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/energy_angle_vs_dcs.pdf'
    set tics font "Times New Roman, 20"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 24"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 24"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 24"
    set title "H+Xe Angular DCS" font "Times New Roman, 24, bold"
    set border linewidth 1.8
    set logscale z
    set zrange [0.1:1000]
    set xrange [1:50]
    set yrange [0:180]
    set xyplane at 0.1
    set hidden3d
    set view 50, 45, 1, 1 
    splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/angdcs.2pi.txt' w l palette lw 1.2 notitle
    