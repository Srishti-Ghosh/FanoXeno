set terminal pdfcairo size 8,8
set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/real_w_imag_w_vs_angle.pdf'

    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)" font "Times New Roman, 20, bold"    
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-2,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-2,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    set grid; set border linewidth 1.5; set xyplane at 0; set zrange [0:180]; set ztics 30    
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    
    set view 60, 55, 1, 1.2
    splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/wparameter.txt' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/wparameter.txt' using 1:2:3:3 w impulses palette lw 0.5 notitle, \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/H+Xe/plots/l_7_plots/wparameter.txt' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    