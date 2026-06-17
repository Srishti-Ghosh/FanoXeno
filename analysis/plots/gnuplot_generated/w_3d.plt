
    set terminal pdfcairo size 8,8
    set output '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/analysis/plots/gnuplot_generated/koike_graph_3_w_trajectory_3d.pdf'
    set title font "Times New Roman, 20, bold"
    set title "H+Xe: 3D Angular Evolution of w(θ) (l=7)"
    
    set xlabel "Re[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set ylabel "Im[w(θ)]" font "Times New Roman, 16" offset 0,-1,0
    set zlabel "Scattering Angle θ (°)" font "Times New Roman, 16" offset 3,0,0
    
    set grid
    set border linewidth 1.5
    
    # Make the Z-axis (Angle) start exactly at the floor
    set xyplane at 0
    set zrange [0:180]
    set ztics 30
    
    # Adjust viewing angle to match the Koike 3D perspective
    set view 60, 55, 1, 1.2
    
    # High contrast color map
    set palette defined (0 "navy", 1 "blue", 2 "cyan", 3 "green", 4 "yellow", 5 "red")
    
    # LAYER 1: The X-Y Plane Trace (Shadow) - Forcing Z to 0
    # LAYER 2: The Drop Lines (Impulses)
    # LAYER 3: The Thick 3D Trajectory Tube
    
    splot '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:(0):3 w lines black lw 2.5 title "X-Y Trace", \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3:3 w impulses palette lw 0.5 notitle, \
          '/Users/srishtighosh/Downloads/CodesForColdCollisionProject/vpa/wparameter.txt' using 1:2:3:3 w lines palette lw 5 title "3D Trajectory"
    