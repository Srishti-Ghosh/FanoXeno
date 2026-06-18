set terminal pdfcairo size 6,6
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/energy_angle_vs_asydcs.pdf'

    set tics font "Times New Roman, 16"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 20"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 20"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 20"
    set border linewidth 1.5; set logscale z; set zrange [0.1:1000]
    set xrange [37.95:41.95]; set yrange [0:180]; set xyplane at 0.1
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    
    set title "H+Xe Synthesized Fano DCS (l=8)" font "Times New Roman, 20, bold"
    set hidden3d front; set view 50, 45, 1, 1 
    splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_8_plots/asydcs.2pi.txt' w l palette lw 1.2 notitle
    