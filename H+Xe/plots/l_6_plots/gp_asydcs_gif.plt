set terminal gif animate delay 7 size 800,700 optimize
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/energy_angle_vs_asydcs_rotating.gif'

    set tics font "Times New Roman, 16"
    set xlabel "E (cm⁻¹)" offset 0,-1 font "Times New Roman, 20"
    set ylabel "θ (degrees)" offset 1,-1 font "Times New Roman, 20"
    set zlabel "dσ/dΩ (a.u.)" offset 3,0 font "Times New Roman, 20"
    set border linewidth 1.5; set logscale z; set zrange [0.1:1000]
    set xrange [12.475:16.475]; set yrange [0:180]; set xyplane at 0.1
    set palette defined (0 "#000040", 1 "#0000ff", 2 "#00ffff", 3 "#00ff00", 4 "#ffff00", 5 "#ff0000")
    
    set title "H+Xe Synthesized Fano DCS (l=6)" font "Times New Roman, 20, bold"
    unset pm3d; set hidden3d front; set logscale cb; set cbrange [0.1:100]
    do for [ang=0:356:4] { set view 60, ang, 1, 1.2; splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_6_plots/asydcs.2pi.txt' every 3:3 w l palette lw 1.2 notitle }
    