
    set terminal gif animate delay 5 size 800,800 optimize
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/final/real_w_imag_w_vs_angle_rotating.gif'
    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)" font "Times New Roman, 20, bold"
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    set grid
    set border linewidth 1.5
    set ticslevel 0
    set xyplane at 0.1
    set zrange [0:180]
    set ztics 30
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    do for [ang=0:356:4] {
        set view 60, ang, 1, 1.2
        splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3:3 w impulses palette lw 0.5 notitle, \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    }
    set output
    